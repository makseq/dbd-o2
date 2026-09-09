// loader.js — загрузка прогона ДБР (manifest.json + frames.bin + series.bin) в браузере.
//
// Это КАНОНИЧЕСКИЙ ДЕКОДЕР формата: recorder.js кодирует, этот файл декодирует,
// тест recorder.test.js гоняет round-trip именно через него.
//
// Обратные формулы квантования (парные к recorder.js):
//   log   : v = (q === 0) ? 0 : 10^(log10(min) + q/scale)      // min = max*10^-decades
//   lin   : v = q / scale                                       // scale = qmax/max
//   asinh : v = v0 * sinh(q / scale)                            // q = 0 -> ТОЧНЫЙ ноль
//
// Единственная зависимость — словарь i18n.js (тексты ошибок); в node работает
// при передаче своего fetch:
//   Playback.load(dirUrl, { fetch: makeFileFetch() })

import { t as tr } from './i18n.js';

export const SUPPORTED_FORMAT_VERSION = 1;

const TYPED = {
  u8: Uint8Array, u16: Uint16Array, i8: Int8Array, i16: Int16Array,
  f32: Float32Array, f64: Float64Array,
};

/** Декод одного квантованного поля кадра в Float32Array. */
export function decodeField(q, codec, spec, out) {
  const n = q.length;
  const dst = out && out.length === n ? out : new Float32Array(n);
  const map = spec.map;
  if (!codec || !(codec.scale > 0)) { dst.fill(0); return dst; }
  if (map === 'log') {
    const lmin = Math.log10(codec.min);
    const inv = 1 / codec.scale;
    for (let k = 0; k < n; k++) {
      const v = q[k];
      dst[k] = v === 0 ? 0 : Math.pow(10, lmin + v * inv);
    }
  } else if (map === 'lin') {
    const inv = 1 / codec.scale;
    for (let k = 0; k < n; k++) dst[k] = q[k] * inv;
  } else if (map === 'asinh') {
    const v0 = codec.v0, inv = 1 / codec.scale;
    for (let k = 0; k < n; k++) {
      const v = q[k];
      dst[k] = v === 0 ? 0 : v0 * Math.sinh(v * inv);
    }
  } else {
    throw new Error(tr('loader.badMap', { map }));
  }
  return dst;
}

function readHeader(buf, magic) {
  const dv = new DataView(buf);
  let s = '';
  for (let i = 0; i < 8; i++) s += String.fromCharCode(dv.getUint8(i));
  if (s !== magic) throw new Error(tr('loader.badMagic', { got: s.replace(/\0/g, ''), want: magic.replace(/\0/g, '') }));
  return {
    version: dv.getUint32(8, true),
    headerBytes: dv.getUint32(12, true),
    stride: dv.getUint32(16, true),
    count: dv.getUint32(20, true),
    a: dv.getUint32(24, true),
    b: dv.getUint32(28, true),
    c: dv.getUint32(32, true),
  };
}

async function fetchBuffer(fetchFn, url, onProgress, total) {
  const res = await fetchFn(url);
  if (!res.ok) throw new Error(`loader: ${url} -> HTTP ${res.status}`);
  const len = Number(res.headers && res.headers.get ? res.headers.get('content-length') : 0) || total || 0;
  // Потоковая загрузка с прогрессом (браузер); в node-шиме body отсутствует -> обычный путь.
  if (onProgress && res.body && typeof res.body.getReader === 'function') {
    const reader = res.body.getReader();
    const chunks = [];
    let got = 0;
    for (;;) {
      const { done, value } = await reader.read();
      if (done) break;
      chunks.push(value); got += value.byteLength;
      onProgress(got, len);
    }
    const out = new Uint8Array(got);
    let o = 0;
    for (const c of chunks) { out.set(c, o); o += c.byteLength; }
    return out.buffer;
  }
  const ab = await res.arrayBuffer();
  if (onProgress) onProgress(ab.byteLength, ab.byteLength);
  return ab;
}

/** LRU-кэш декодированных кадров. */
class FrameCache {
  constructor(limit) { this.limit = limit; this.map = new Map(); }
  get(key) {
    if (!this.map.has(key)) return undefined;
    const v = this.map.get(key);
    this.map.delete(key); this.map.set(key, v); // освежить
    return v;
  }
  set(key, v) {
    if (this.map.has(key)) this.map.delete(key);
    this.map.set(key, v);
    while (this.map.size > this.limit) this.map.delete(this.map.keys().next().value);
  }
  clear() { this.map.clear(); }
  get size() { return this.map.size; }
}

export class Playback {
  constructor(manifest, opts) {
    this.manifest = manifest;
    this.baseUrl = opts.baseUrl;
    this.fetchFn = opts.fetch;
    this.cache = new FrameCache(opts.cacheFrames ?? 64);
    this.fieldByName = new Map(manifest.fields.map((f) => [f.name, f]));
    this.frameCount = manifest.frames.length;
    this.frameTimes = Float64Array.from(manifest.frames.map((f) => f.t));
    this.framesBuf = null;   // ArrayBuffer целиком (режим 'eager')
    this.mode = opts.mode || 'eager';
    this.series = null;
  }

  /**
   * @param {string} url — URL каталога прогона (или самого manifest.json)
   * @param {object} opts — { fetch, mode: 'eager'|'range', cacheFrames, onProgress, withSeries }
   */
  static async load(url, opts = {}) {
    const fetchFn = opts.fetch || (typeof fetch !== 'undefined' ? fetch.bind(globalThis) : null);
    if (!fetchFn) throw new Error(tr('loader.noFetch'));
    const base = url.endsWith('manifest.json') ? url.slice(0, -'manifest.json'.length) : (url.endsWith('/') ? url : url + '/');
    const res = await fetchFn(base + 'manifest.json');
    if (!res.ok) throw new Error(tr('loader.manifestHttp', { status: res.status }));
    const manifest = await res.json();
    if (manifest.formatVersion !== SUPPORTED_FORMAT_VERSION) {
      throw new Error(tr('loader.badVersion', { got: manifest.formatVersion, want: SUPPORTED_FORMAT_VERSION }));
    }
    const pb = new Playback(manifest, { ...opts, baseUrl: base, fetch: fetchFn });
    if (pb.mode === 'eager') {
      pb.framesBuf = await fetchBuffer(fetchFn, base + 'frames.bin', opts.onProgress);
      pb.header = readHeader(pb.framesBuf, 'DBD2FRM\0');
      const expect = manifest.frameLayout.headerBytes + manifest.frameLayout.frameStride * manifest.frameLayout.frameCount;
      if (pb.framesBuf.byteLength !== expect) {
        throw new Error(tr('loader.framesSize', { got: pb.framesBuf.byteLength, want: expect }));
      }
      if (pb.header.count !== manifest.frames.length) {
        throw new Error(tr('loader.framesCount', { got: pb.header.count, want: manifest.frames.length }));
      }
    }
    if (opts.withSeries !== false) await pb.loadSeries();
    return pb;
  }

  /** Скалярные ряды -> колоночные типизированные массивы. */
  async loadSeries() {
    const s = this.manifest.series;
    const buf = await fetchBuffer(this.fetchFn, this.baseUrl + s.file);
    const h = readHeader(buf, 'DBD2SER\0');
    if (h.count !== s.count || h.stride !== s.recordStride) {
      throw new Error(tr('loader.seriesMismatch', { got: `${h.count}/${h.stride}`, want: `${s.count}/${s.recordStride}` }));
    }
    const dv = new DataView(buf);
    const n = s.count, base = s.headerBytes, st = s.recordStride;
    const out = {};
    for (let k = 0; k < s.names.length; k++) {
      out[s.names[k]] = s.dtypes[k] === 'f64' ? new Float64Array(n) : new Float32Array(n);
    }
    for (let i = 0; i < n; i++) {
      const o = base + i * st;
      for (let k = 0; k < s.names.length; k++) {
        const off = o + s.offsets[k];
        out[s.names[k]][i] = s.dtypes[k] === 'f64' ? dv.getFloat64(off, true) : dv.getFloat32(off, true);
      }
    }
    this.series = out;
    return out;
  }

  _frameBytes(i) {
    const fl = this.manifest.frameLayout;
    if (this.mode === 'eager') {
      const off = this.manifest.frames[i].byteOffset;
      if (off + fl.frameStride > this.framesBuf.byteLength) throw new Error(tr('loader.frameOutOfFile', { i }));
      return { buf: this.framesBuf, off };
    }
    const blk = this.blocks && this.blocks.get(i);
    if (!blk) throw new Error(tr('loader.frameNotLoaded', { i }));
    return { buf: blk, off: 0 };
  }

  /** Режим 'range': подкачать блок кадра одним Range-запросом (кэш блоков LRU). */
  async fetchFrameBlock(i) {
    if (this.mode === 'eager') return;
    if (!this.blocks) this.blocks = new FrameCache(this.cache.limit);
    if (this.blocks.get(i)) return;
    const fl = this.manifest.frameLayout;
    const off = this.manifest.frames[i].byteOffset;
    const res = await this.fetchFn(this.baseUrl + 'frames.bin', {
      headers: { range: `bytes=${off}-${off + fl.frameStride - 1}` },
    });
    if (!res.ok) throw new Error(tr('loader.rangeHttp', { i, status: res.status }));
    const ab = await res.arrayBuffer();
    if (ab.byteLength !== fl.frameStride) {
      throw new Error(tr('loader.rangeShort', { got: ab.byteLength, want: fl.frameStride }));
    }
    this.blocks.set(i, ab);
  }

  /** Универсальный доступ: работает в обоих режимах. */
  async getFrameAsync(i, field) {
    if (this.mode !== 'eager') await this.fetchFrameBlock(i);
    return this.getFrame(i, field);
  }

  /** Подкачать диапазон кадров вперёд (плавное проигрывание в режиме range). */
  async prefetch(i0, i1) {
    if (this.mode === 'eager') return;
    for (let i = Math.max(0, i0); i <= Math.min(this.frameCount - 1, i1); i++) await this.fetchFrameBlock(i);
  }

  /** Декодированное поле кадра (Float32Array длиной nrOut*nzOut, раскладка i*nzOut+j). */
  getFrame(i, field) {
    if (i < 0 || i >= this.frameCount) throw new Error(tr('loader.noFrame', { i, n: this.frameCount }));
    const key = `${i}|${field}`;
    const hit = this.cache.get(key);
    if (hit) return hit;
    const spec = this.fieldByName.get(field);
    if (!spec) throw new Error(tr('loader.noField', { field }));
    const { buf, off } = this._frameBytes(i);
    const q = new TYPED[spec.dtype](buf, off + spec.offsetInFrame, spec.shape[0] * spec.shape[1]);
    const codec = this.manifest.frames[i].fields[field];
    const dec = decodeField(q, codec, spec);
    this.cache.set(key, dec);
    return dec;
  }

  /** Сырые квантованные отсчёты кадра (без копии) — для быстрой отрисовки палитрой. */
  getFrameRaw(i, field) {
    const spec = this.fieldByName.get(field);
    if (!spec) throw new Error(tr('loader.noField', { field }));
    const { buf, off } = this._frameBytes(i);
    return new TYPED[spec.dtype](buf, off + spec.offsetInFrame, spec.shape[0] * spec.shape[1]);
  }

  /** Кодек кадра (min/max/scale/log) — плеер обязан показывать покадровую шкалу. */
  getCodec(i, field) { return this.manifest.frames[i].fields[field]; }

  /** Поверхностный заряд на кадре: { sigmaL, sigmaR } (Float32Array по r, Кл/м^2). */
  getSurface(i) {
    const { buf, off } = this._frameBytes(i);
    const s = this.manifest.frameLayout.surfaces;
    return {
      sigmaL: new Float32Array(buf.slice(off + s.sigmaL.offsetInFrame, off + s.sigmaL.offsetInFrame + s.sigmaL.length * 4)),
      sigmaR: new Float32Array(buf.slice(off + s.sigmaR.offsetInFrame, off + s.sigmaR.offsetInFrame + s.sigmaR.length * 4)),
    };
  }

  /** Индекс кадра, ближайшего к времени t (бинарный поиск). */
  findFrame(t) {
    const a = this.frameTimes;
    let lo = 0, hi = a.length - 1;
    if (t <= a[0]) return 0;
    if (t >= a[hi]) return hi;
    while (hi - lo > 1) { const m = (lo + hi) >> 1; if (a[m] <= t) lo = m; else hi = m; }
    return (t - a[lo] <= a[hi] - t) ? lo : hi;
  }

  get fields() { return this.manifest.fields.map((f) => f.name); }
  get grid() { return this.manifest.grid; }
  clearCache() { this.cache.clear(); }
}

export default Playback;
