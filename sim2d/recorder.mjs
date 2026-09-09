// recorder.mjs — запись результатов 2D-прогона ДБР на диск.
//
// Контейнер прогона: <dir>/manifest.json + frames.bin + series.bin
//   frames.bin  — квантованные 2D-поля по ГАЗОВОЙ части сетки + sigmaL(r), sigmaR(r) на кадр
//   series.bin  — скалярные временные ряды (в полном уровне детализации — КАЖДЫЙ шаг)
//   manifest.json — пишется ПОСЛЕДНИМ, служит маркером «прогон завершён»
//
// Раскладка полей — как в солвере: field[i*nz + j], z — быстрый индекс.
// Все квантователи и их обратные формулы задокументированы в loader.mjs (декодер),
// здесь — только кодер. Тест recorder.test.mjs проверяет round-trip через ЛОАДЕР.

import fs from 'node:fs';
import path from 'node:path';

export const FORMAT_VERSION = 1;

export const FRAMES_MAGIC = 'DBD2FRM\0';
export const SERIES_MAGIC = 'DBD2SER\0';
export const FRAMES_HEADER_BYTES = 64;
export const SERIES_HEADER_BYTES = 64;

// Скалярные ряды: t хранится Float64 (t~3e-4 с при dt~3e-12 с — Float32 не хватает),
// остальные — Float32. Запись фиксированной длины, выровнена на 8 байт.
export const SERIES_NAMES = [
  't', 'Uapp', 'Ugap', 'Icond', 'Idisp', 'Itot', 'Q',
  'maxEN', 'o3ppm', 'sigmaMax', 'photoEmitTotalL', 'photoEmitTotalR',
];
export const SERIES_STRIDE = 56; // 8 (f64 t) + 11*4 (f32) = 52 -> паддинг до 56

export const DTYPES = {
  u8:  { bytes: 1, signed: false, qmax: 255,   Arr: Uint8Array },
  u16: { bytes: 2, signed: false, qmax: 65535, Arr: Uint16Array },
  i8:  { bytes: 1, signed: true,  qmax: 127,   Arr: Int8Array },
  i16: { bytes: 2, signed: true,  qmax: 32767, Arr: Int16Array },
};

const align = (x, a) => Math.ceil(x / a) * a;

// ---------------------------------------------------------------------------
// Описание полей кадра
// ---------------------------------------------------------------------------
// map: 'log'   — знакоположительное, огромный динамический диапазон (плотности,
//                свечение, фотоскорости). Покадровая нормировка на vmax, DEC декад вниз.
//      'lin'   — знакоположительное, умеренный диапазон (|E|, E/N).
//      'asinh' — ЗНАКОПЕРЕМЕННОЕ (объёмный заряд): линейно около нуля, логарифмически
//                на хвостах; q=0 декодируется в ТОЧНЫЙ ноль.
//
// ВАЖНО про разрядность (честно, без замалчивания):
//   uint8 + 6 декад даёт шаг 6/255 декады => предельная относительная ошибка
//   10^(6/510) - 1 = 2.75 %. Требование ТЗ «< 1 % для лог-полей» с uint8 НЕДОСТИЖИМО
//   в принципе. Поэтому:
//     level='full'    — лог-поля в uint16 (ошибка 1.05e-5), это научный продукт;
//     level='compact' — лог-поля в uint8 (ошибка <= 2.75 %), это картинка для плеера.
//   Тест проверяет строгие пороги на 'full' и теоретическую границу на 'compact'.
function defaultFields(level) {
  const q = level === 'compact' ? 'u8' : 'u16';
  return [
    { name: 'n_e',             src: 'n.e',             map: 'log',   decades: 6, dtype: q,
      unit: 'м^-3',      label: 'Плотность электронов' },
    { name: 'ionizRate',       src: 'ionizRate',       map: 'log',   decades: 5, dtype: q,
      unit: 'м^-3 с^-1', label: 'Скорость ионизации (свечение)' },
    { name: 'rho',             src: 'rho',             map: 'asinh', v0rel: 1e-4, dtype: 'i16',
      unit: 'Кл/м^3',    label: 'Объёмный заряд' },
    { name: 'Emag',            src: 'E',               map: 'lin',   dtype: 'u16',
      unit: 'В/м',       label: '|E|' },
    { name: 'EN',              src: 'EN',              map: 'lin',   dtype: level === 'compact' ? 'u8' : 'u16',
      unit: 'Тд',        label: 'Приведённое поле E/N' },
    { name: 'photoIonRate',    src: 'photoIonRate',    map: 'log',   decades: 5, dtype: q,
      unit: 'м^-3 с^-1', label: 'Скорость фотоионизации' },
    { name: 'photoDetachRate', src: 'photoDetachRate', map: 'log',   decades: 5, dtype: q,
      unit: 'м^-3 с^-1', label: 'Скорость фотоотлипания' },
    { name: 'n_O3m',           src: 'n.O3m',           map: 'log',   decades: 6, dtype: q,
      unit: 'м^-3',      label: 'Плотность O3-' },
    { name: 'n_O3',            src: 'n.O3',            map: 'log',   decades: 4, dtype: q,
      unit: 'м^-3',      label: 'Плотность O3' },
  ];
}

// Бюджеты подобраны так, чтобы КАТАЛОГ ПРОГОНА ЦЕЛИКОМ (frames+series+manifest) влезал:
//   full    <= ~730 МБ  (frames 700 + series ~25 + manifest ~1)  -> 4 прогона < 3 ГБ
//   compact <= ~24 МБ   (frames 20 + series 3.4 + manifest ~0.5) -> < 30 МБ, цель плеера
export const LEVEL_DEFAULTS = {
  full:    { dsR: 1, dsZ: 1, maxBytes: 700e6,  seriesMaxRecords: Infinity, interestScale: 1 },
  compact: { dsR: 2, dsZ: 2, maxBytes: 20e6,   seriesMaxRecords: 60000,    interestScale: 2 },
};

// ---------------------------------------------------------------------------
// Кодеры (обратные формулы — в loader.mjs, менять только парой)
// ---------------------------------------------------------------------------

// Лог: q = 0 <=> v <= vmin (или v <= 0). Иначе v = 10^(lmin + q/scale).
function encodeLog(src, out, spec) {
  const n = src.length;
  let vmax = 0;
  for (let k = 0; k < n; k++) { const v = src[k]; if (v > vmax && Number.isFinite(v)) vmax = v; }
  const qmax = DTYPES[spec.dtype].qmax;
  if (!(vmax > 0)) { out.fill(0); return { min: 0, max: 0, scale: 0, log: true, decades: spec.decades }; }
  const lmax = Math.log10(vmax);
  const lmin = lmax - spec.decades;
  const scale = qmax / spec.decades;
  for (let k = 0; k < n; k++) {
    const v = src[k];
    if (!(v > 0) || !Number.isFinite(v)) { out[k] = 0; continue; }
    let q = Math.round((Math.log10(v) - lmin) * scale);
    if (q < 0) q = 0; else if (q > qmax) q = qmax;
    out[k] = q;
  }
  return { min: Math.pow(10, lmin), max: vmax, scale, log: true, decades: spec.decades };
}

// Линейный: v = q / scale, scale = qmax / vmax.
function encodeLin(src, out, spec) {
  const n = src.length;
  let vmax = 0;
  for (let k = 0; k < n; k++) { const v = src[k]; if (Number.isFinite(v) && v > vmax) vmax = v; }
  const qmax = DTYPES[spec.dtype].qmax;
  if (!(vmax > 0)) { out.fill(0); return { min: 0, max: 0, scale: 0, log: false }; }
  const scale = qmax / vmax;
  for (let k = 0; k < n; k++) {
    const v = src[k];
    if (!Number.isFinite(v) || v <= 0) { out[k] = 0; continue; }
    let q = Math.round(v * scale);
    if (q > qmax) q = qmax;
    out[k] = q;
  }
  return { min: 0, max: vmax, scale, log: false };
}

// asinh (ЗНАКОПЕРЕМЕННОЕ): v = v0 * sinh(q / scale). q = 0 -> ровно 0.
function encodeAsinh(src, out, spec) {
  const n = src.length;
  let vmax = 0;
  for (let k = 0; k < n; k++) { const a = Math.abs(src[k]); if (Number.isFinite(a) && a > vmax) vmax = a; }
  const qmax = DTYPES[spec.dtype].qmax; // положительный предел знакового типа
  if (!(vmax > 0)) { out.fill(0); return { min: 0, max: 0, scale: 0, v0: 0, log: 'asinh' }; }
  const v0 = vmax * spec.v0rel;
  const scale = qmax / Math.asinh(vmax / v0);
  for (let k = 0; k < n; k++) {
    const v = src[k];
    if (!Number.isFinite(v) || v === 0) { out[k] = 0; continue; }
    let q = Math.round(scale * Math.asinh(v / v0));
    if (q > qmax) q = qmax; else if (q < -qmax) q = -qmax;
    out[k] = q;
  }
  return { min: -vmax, max: vmax, scale, v0, log: 'asinh' };
}

function encodeField(src, out, spec) {
  if (spec.map === 'log') return encodeLog(src, out, spec);
  if (spec.map === 'lin') return encodeLin(src, out, spec);
  if (spec.map === 'asinh') return encodeAsinh(src, out, spec);
  throw new Error(`recorder: неизвестный map '${spec.map}'`);
}

// ---------------------------------------------------------------------------
// Вспомогательное: разбор сетки из состояния солвера
// ---------------------------------------------------------------------------
function resolveGrid(state) {
  const probe = state.E || state.rho || (state.n && state.n.e);
  if (!probe) throw new Error('recorder: в state нет ни E, ни rho, ни n.e — не могу определить сетку');
  const total = probe.length;
  const r = state.r, z = state.z;
  if (!r || !z) throw new Error('recorder: в state нет r/z');
  const cands = [[r.length, z.length], [r.length - 1, z.length], [r.length, z.length - 1], [r.length - 1, z.length - 1]];
  for (const [nr, nz] of cands) if (nr > 0 && nz > 0 && nr * nz === total) return { nr, nz };
  throw new Error(`recorder: не сходится сетка: |r|=${r.length}, |z|=${z.length}, cells=${total}`);
}

// Грани из центров (если солвер не отдал rf/zf явно).
function facesFrom(centers, n, first) {
  if (centers.length === n + 1) return Float64Array.from(centers);
  const f = new Float64Array(n + 1);
  f[0] = first !== undefined ? first : centers[0] - 0.5 * (centers[1] - centers[0]);
  for (let i = 1; i < n; i++) f[i] = 0.5 * (centers[i - 1] + centers[i]);
  f[n] = centers[n - 1] + (centers[n - 1] - f[n - 1]);
  return f;
}

function gasRows(state, nr, nz) {
  const m = state.gasMask;
  let j0 = -1, j1 = -1;
  if (m && m.length === nz) {
    for (let j = 0; j < nz; j++) if (m[j]) { if (j0 < 0) j0 = j; j1 = j; }
  } else if (m && m.length === nr * nz) {
    for (let j = 0; j < nz; j++) if (m[j]) { if (j0 < 0) j0 = j; j1 = j; } // столбец i=0
  }
  if (j0 < 0) { j0 = 0; j1 = nz - 1; }
  return { j0, j1 };
}

function pick(state, srcPath) {
  if (srcPath.startsWith('n.')) return state.n ? state.n[srcPath.slice(2)] : undefined;
  return state[srcPath];
}

// ---------------------------------------------------------------------------
// Recorder
// ---------------------------------------------------------------------------
export class Recorder {
  /**
   * @param {object} opts
   *   dir              — каталог прогона (создаётся)
   *   runId            — имя прогона
   *   level            — 'full' | 'compact'
   *   params           — параметры прогона (пишутся в манифест целиком)
   *   meta             — произвольные метаданные (версии, git-sha, комментарий)
   *   tEnd             — ожидаемое время конца прогона (для бюджетного темпа выборки)
   *   dtFrameMin/Max   — 1 нс / 1 мкс по умолчанию
   *   maxBytes         — бюджет frames.bin
   *   seriesMaxRecords — бюджет записей series.bin (Infinity в 'full')
   *   dsR, dsZ         — прореживание сетки (площадное усреднение) для компактного уровня
   */
  constructor(opts = {}) {
    const level = opts.level || 'full';
    const L = LEVEL_DEFAULTS[level];
    if (!L) throw new Error(`recorder: неизвестный level '${level}'`);
    this.level = level;
    this.dir = opts.dir || path.join('data', opts.runId || 'run');
    this.runId = opts.runId || path.basename(this.dir);
    this.params = opts.params || {};
    this.meta = opts.meta || {};
    this.tEnd = opts.tEnd || 0;
    this.dtFrameMin = opts.dtFrameMin ?? 1e-9;
    this.dtFrameMax = opts.dtFrameMax ?? 1e-6;
    this.maxBytes = opts.maxBytes ?? L.maxBytes;
    this.seriesMaxRecords = opts.seriesMaxRecords ?? L.seriesMaxRecords;
    this.dsR = opts.dsR ?? L.dsR;
    this.dsZ = opts.dsZ ?? L.dsZ;
    this.interestScale = opts.interestScale ?? L.interestScale;
    this.interestScale0 = this.interestScale;      // базовое значение, ниже регулятор не опускает
    this.dtFrameMinEff = this.dtFrameMin;          // текущий минимальный интервал (растёт под бюджет)
    this._lastCtrlN = 0;
    this.fieldSpecs = opts.fields || defaultFields(level);

    this.inited = false;
    this.closed = false;
    this.frames = [];
    this.frameCount = 0;
    this.seriesCount = 0;
    this.stepCount = 0;
    this.nonFiniteHits = 0;

    // состояние адаптивного эмиттера
    this.lastFrameT = -Infinity;
    this.lastLogNe = -Infinity;
    this.lastI = 0;
    this.lastUgap = 0;
    this.IpeakRun = 0;
  }

  // ---- инициализация по первому состоянию --------------------------------
  _init(state) {
    const { nr, nz } = resolveGrid(state);
    this.nr = nr; this.nz = nz;
    const { j0, j1 } = gasRows(state, nr, nz);
    this.jGas0 = j0; this.jGas1 = j1;
    this.nzGas = j1 - j0 + 1;

    this.rf = state.rf ? Float64Array.from(state.rf) : facesFrom(state.r, nr, 0);
    this.zf = state.zf ? Float64Array.from(state.zf) : facesFrom(state.z, nz, 0);

    // веса площадного усреднения: Acell[i] = pi*(rf[i+1]^2 - rf[i]^2), dz[j]
    this.wR = new Float64Array(nr);
    for (let i = 0; i < nr; i++) this.wR[i] = Math.PI * (this.rf[i + 1] ** 2 - this.rf[i] ** 2);
    this.wZ = new Float64Array(nz);
    for (let j = 0; j < nz; j++) this.wZ[j] = this.zf[j + 1] - this.zf[j];

    this.nrOut = Math.ceil(nr / this.dsR);
    this.nzOut = Math.ceil(this.nzGas / this.dsZ);
    this.nOut = this.nrOut * this.nzOut;

    // раскладка кадра: поля по порядку, каждый блок выровнен на 4 байта,
    // затем sigmaL[nrOut], sigmaR[nrOut] (Float32), stride выровнен на 8
    let off = 0;
    this.fields = [];
    for (const spec of this.fieldSpecs) {
      const dt = DTYPES[spec.dtype];
      if (!dt) throw new Error(`recorder: неизвестный dtype '${spec.dtype}' у поля ${spec.name}`);
      const byteLength = this.nOut * dt.bytes;
      const present = pick(state, spec.src) != null;
      this.fields.push({
        ...spec, offsetInFrame: off, byteLength, present,
        shape: [this.nrOut, this.nzOut], layout: 'i-major (z contiguous)',
      });
      off = align(off + byteLength, 4);
    }
    this.sigLoOffset = off; off += this.nrOut * 4;
    this.sigHiOffset = off; off += this.nrOut * 4;
    this.frameStride = align(off, 8);

    // буферы
    this.frameBuf = Buffer.alloc(this.frameStride);
    this.workF64 = new Float64Array(this.nOut);       // поле после прореживания
    this.accW = new Float64Array(this.nOut);          // сумма весов для усреднения
    this.qViews = this.fields.map((f) => {
      const dt = DTYPES[f.dtype];
      return new dt.Arr(this.frameBuf.buffer, this.frameBuf.byteOffset + f.offsetInFrame, this.nOut);
    });
    this.sigLoView = new Float32Array(this.frameBuf.buffer, this.frameBuf.byteOffset + this.sigLoOffset, this.nrOut);
    this.sigHiView = new Float32Array(this.frameBuf.buffer, this.frameBuf.byteOffset + this.sigHiOffset, this.nrOut);

    this.maxFrames = Math.max(1, Math.floor(this.maxBytes / this.frameStride));

    fs.mkdirSync(this.dir, { recursive: true });
    // manifest.json пишется последним; если остался от прошлого прогона — убрать
    const mpath = path.join(this.dir, 'manifest.json');
    if (fs.existsSync(mpath)) fs.unlinkSync(mpath);

    this.framesFd = fs.openSync(path.join(this.dir, 'frames.bin'), 'w');
    this.seriesFd = fs.openSync(path.join(this.dir, 'series.bin'), 'w');
    fs.writeSync(this.framesFd, this._framesHeader(), 0, FRAMES_HEADER_BYTES, 0);
    fs.writeSync(this.seriesFd, this._seriesHeader(), 0, SERIES_HEADER_BYTES, 0);
    this.framesBytes = FRAMES_HEADER_BYTES;
    this.seriesBytes = SERIES_HEADER_BYTES;

    this.seriesChunk = Buffer.alloc(SERIES_STRIDE * 1024);
    this.seriesChunkN = 0;

    this.tFirst = state.t; this.tLast = state.t;
    this.inited = true;
  }

  _framesHeader() {
    const b = Buffer.alloc(FRAMES_HEADER_BYTES);
    b.write(FRAMES_MAGIC, 0, 'latin1');
    b.writeUInt32LE(FORMAT_VERSION, 8);
    b.writeUInt32LE(FRAMES_HEADER_BYTES, 12);
    b.writeUInt32LE(this.frameStride, 16);
    b.writeUInt32LE(this.frameCount, 20);
    b.writeUInt32LE(this.nrOut, 24);
    b.writeUInt32LE(this.nzOut, 28);
    b.writeUInt32LE(this.fields.length, 32);
    b.writeUInt32LE(0, 36);
    b.writeDoubleLE(this.tFirst ?? 0, 40);
    b.writeDoubleLE(this.tLast ?? 0, 48);
    return b;
  }

  _seriesHeader() {
    const b = Buffer.alloc(SERIES_HEADER_BYTES);
    b.write(SERIES_MAGIC, 0, 'latin1');
    b.writeUInt32LE(FORMAT_VERSION, 8);
    b.writeUInt32LE(SERIES_HEADER_BYTES, 12);
    b.writeUInt32LE(SERIES_STRIDE, 16);
    b.writeUInt32LE(this.seriesCount, 20);
    b.writeUInt32LE(SERIES_NAMES.length, 24);
    b.writeUInt32LE(0, 28);
    return b;
  }

  // ---- прореживание + вырезка газовой части ------------------------------
  // Площадное усреднение с весом Acell[i]*dz[j] ДО квантования (см. NUMERICS_2D §9.5).
  _reduce(src) {
    const { nr, nz, jGas0, jGas1, dsR, dsZ, nrOut, nzOut } = this;
    const out = this.workF64, w = this.accW;
    if (dsR === 1 && dsZ === 1) {
      for (let i = 0; i < nr; i++) {
        const base = i * nz + jGas0, ob = i * nzOut;
        for (let j = 0; j < nzOut; j++) {
          const v = src[base + j];
          out[ob + j] = Number.isFinite(v) ? v : (this.nonFiniteHits++, 0);
        }
      }
      return out;
    }
    out.fill(0); w.fill(0);
    for (let i = 0; i < nr; i++) {
      const io = (i / dsR) | 0;
      const wr = this.wR[i];
      for (let j = jGas0; j <= jGas1; j++) {
        const jo = ((j - jGas0) / dsZ) | 0;
        let v = src[i * nz + j];
        if (!Number.isFinite(v)) { v = 0; this.nonFiniteHits++; }
        const ww = wr * this.wZ[j];
        out[io * nzOut + jo] += v * ww;
        w[io * nzOut + jo] += ww;
      }
    }
    for (let k = 0; k < out.length; k++) if (w[k] > 0) out[k] /= w[k];
    return out;
  }

  _reduceSurf(src, dst) {
    const { nr, dsR, nrOut } = this;
    if (dsR === 1) { for (let i = 0; i < nrOut; i++) dst[i] = src[i]; return; }
    const acc = new Float64Array(nrOut), w = new Float64Array(nrOut);
    for (let i = 0; i < nr; i++) {
      const io = (i / dsR) | 0, ww = this.wR[i];
      const v = Number.isFinite(src[i]) ? src[i] : 0;
      acc[io] += v * ww; w[io] += ww;
    }
    for (let i = 0; i < nrOut; i++) dst[i] = w[i] > 0 ? acc[i] / w[i] : 0;
  }

  // ---- критерий «интересности» (NUMERICS_2D §9.3) -------------------------
  _wantFrame(state) {
    const t = state.t;
    const dtSince = t - this.lastFrameT;
    if (dtSince < this.dtFrameMinEff) return false;
    if (dtSince >= this.dtFrameMax) return true;
    let neMax = 0;
    if (Number.isFinite(state.neMax)) {
      neMax = state.neMax;                       // солвер уже посчитал — не сканируем массив
    } else {
      const ne = state.n && state.n.e;
      if (ne) for (let k = 0; k < ne.length; k++) if (ne[k] > neMax) neMax = ne[k];
    }
    const s = this.interestScale;
    const dLogNe = neMax > 0 ? Math.abs(Math.log10(neMax) - this.lastLogNe) : 0;
    const dI = Math.abs((state.Icond || 0) - this.lastI) / Math.max(this.IpeakRun, 1e-12);
    const U0 = Math.abs(this.params.U0kV ? this.params.U0kV * 1e3 : (state.Uapp || 1));
    const dU = Math.abs((state.Ugap || 0) - this.lastUgap) / Math.max(U0, 1);
    return (dLogNe > 0.04 * s) || (dI > 0.02 * s) || (dU > 0.005 * s);
  }

  // Бюджет НЕ распределяется равномерно по времени: это убило бы всю адаптивность
  // (кадры обязаны идти пачкой во время импульса и редко в темноте). Вместо квоты —
  // регулятор ПОРОГА: если по текущему темпу прогон не влезает в maxFrames, критерий
  // «интересности» и минимальный интервал ужесточаются, сохраняя пачечную структуру.
  _updateBudgetControl(t) {
    if (!this.tEnd || this.frameCount < 32) return;
    if (this.frameCount - this._lastCtrlN < 32) return;
    this._lastCtrlN = this.frameCount;
    const frac = (t - this.tFirst) / this.tEnd;
    if (frac < 0.02) return;
    const projected = this.frameCount / frac;
    if (projected > 0.90 * this.maxFrames) {
      this.interestScale *= 1.3;
      this.dtFrameMinEff = Math.min(this.dtFrameMax / 10, this.dtFrameMinEff * 1.3);
    } else if (projected < 0.45 * this.maxFrames) {
      this.interestScale = Math.max(this.interestScale0, this.interestScale / 1.15);
      this.dtFrameMinEff = Math.max(this.dtFrameMin, this.dtFrameMinEff / 1.15);
    }
  }

  _budgetOk(t) {
    this._updateBudgetControl(t);
    return this.frameCount < this.maxFrames;
  }

  // Решение об эмиссии кэшируется по t: shouldRecord() и record() на одном шаге
  // обязаны дать ОДИН И ТОТ ЖЕ ответ, иначе продюсер запишет незаполненные поля.
  _decide(state) {
    if (this._pendT === state.t) return this._pendWant;
    const w = this._wantFrame(state) && this._budgetOk(state.t);
    this._pendT = state.t; this._pendWant = w;
    return w;
  }

  /**
   * Будет ли этот шаг записан кадром. Позволяет продюсеру (солверу) заполнять
   * тяжёлые диагностические 2D-поля ТОЛЬКО перед реальной записью кадра.
   * Вызывать до record() на том же state (тот же state.t).
   */
  shouldRecord(state) {
    if (!this.inited) return true;
    const Ic = Math.abs(state.Icond || 0);
    if (Ic > this.IpeakRun) this.IpeakRun = Ic;
    return this._decide(state);
  }

  /**
   * Вызывать ПОСЛЕ каждого шага солвера.
   * @param {object} state — s.state
   * @param {boolean} force — принудительно записать кадр (первый/последний/событие)
   * @returns {boolean} записан ли кадр
   */
  record(state, force = false) {
    if (this.closed) throw new Error('recorder: запись после close()');
    if (!this.inited) { this._init(state); force = true; }
    this.tLast = state.t;
    this.stepCount++;
    this._appendSeries(state);
    const Ic = Math.abs(state.Icond || 0);
    if (Ic > this.IpeakRun) this.IpeakRun = Ic;
    if (!force && !this._decide(state)) return false;
    if (this.frameCount >= this.maxFrames) return false;
    this._writeFrame(state);
    return true;
  }

  _writeFrame(state) {
    const rec = { index: this.frameCount, t: state.t, step: this.stepCount,
                  byteOffset: this.framesBytes, fields: {} };
    this.frameBuf.fill(0);
    for (let fi = 0; fi < this.fields.length; fi++) {
      const f = this.fields[fi];
      const src = pick(state, f.src);
      if (!src) { rec.fields[f.name] = { min: 0, max: 0, scale: 0, log: f.map === 'log' }; continue; }
      const red = this._reduce(src);
      rec.fields[f.name] = encodeField(red, this.qViews[fi], f);
    }
    if (state.sigmaL) this._reduceSurf(state.sigmaL, this.sigLoView); else this.sigLoView.fill(0);
    if (state.sigmaR) this._reduceSurf(state.sigmaR, this.sigHiView); else this.sigHiView.fill(0);

    fs.writeSync(this.framesFd, this.frameBuf, 0, this.frameStride, this.framesBytes);
    this.framesBytes += this.frameStride;
    this.frames.push(rec);
    this.frameCount++;

    // обновление опорных значений эмиттера
    this.lastFrameT = state.t;
    let neMax = 0; const ne = state.n && state.n.e;
    if (ne) for (let k = 0; k < ne.length; k++) if (ne[k] > neMax) neMax = ne[k];
    this.lastLogNe = neMax > 0 ? Math.log10(neMax) : -Infinity;
    this.lastI = state.Icond || 0;
    this.lastUgap = state.Ugap || 0;
  }

  _seriesBudgetOk() {
    if (!Number.isFinite(this.seriesMaxRecords)) return true;
    if (this.seriesCount >= this.seriesMaxRecords) return false;
    if (!this.tEnd) return true;
    const frac = Math.min(1, Math.max(0, (this.tLast - this.tFirst) / this.tEnd));
    return this.seriesCount < Math.max(1, Math.ceil(this.seriesMaxRecords * frac) + 1);
  }

  _appendSeries(state) {
    if (!this._seriesBudgetOk()) return;
    const b = this.seriesChunk;
    const o = this.seriesChunkN * SERIES_STRIDE;
    b.writeDoubleLE(state.t ?? 0, o);
    for (let k = 1; k < SERIES_NAMES.length; k++) {
      const v = state[SERIES_NAMES[k]];
      b.writeFloatLE(Number.isFinite(v) ? v : 0, o + 8 + (k - 1) * 4);
    }
    b.fill(0, o + 52, o + SERIES_STRIDE);
    this.seriesChunkN++;
    this.seriesCount++;
    if (this.seriesChunkN * SERIES_STRIDE >= b.length) this._flushSeries();
  }

  _flushSeries() {
    if (this.seriesChunkN === 0) return;
    const n = this.seriesChunkN * SERIES_STRIDE;
    fs.writeSync(this.seriesFd, this.seriesChunk, 0, n, this.seriesBytes);
    this.seriesBytes += n;
    this.seriesChunkN = 0;
  }

  /** Закрыть прогон: дописать заголовки и manifest.json (ПОСЛЕДНИМ). */
  /**
   * Промежуточная фиксация: заголовки + manifest.json на диск БЕЗ закрытия файлов.
   * Нужна для длинных прогонов: если процесс убьют, индекс кадров (он живёт в памяти)
   * иначе потеряется и все записанные байты станут мусором. Вызывать редко —
   * manifest с индексом кадров может весить мегабайты.
   */
  checkpoint(extra = {}) {
    if (this.closed || !this.inited) return null;
    this._flushSeries();
    fs.writeSync(this.framesFd, this._framesHeader(), 0, FRAMES_HEADER_BYTES, 0);
    fs.writeSync(this.seriesFd, this._seriesHeader(), 0, SERIES_HEADER_BYTES, 0);
    const tmp = path.join(this.dir, 'manifest.json.tmp');
    fs.writeFileSync(tmp, JSON.stringify(this._manifest({ ...extra, checkpoint: true })));
    fs.renameSync(tmp, path.join(this.dir, 'manifest.json'));   // атомарная замена
    return path.join(this.dir, 'manifest.json');
  }

  close(extra = {}) {
    if (this.closed) return this.manifestPath;
    if (!this.inited) throw new Error('recorder: close() без единого record()');
    this._flushSeries();
    fs.writeSync(this.framesFd, this._framesHeader(), 0, FRAMES_HEADER_BYTES, 0);
    fs.writeSync(this.seriesFd, this._seriesHeader(), 0, SERIES_HEADER_BYTES, 0);
    fs.closeSync(this.framesFd);
    fs.closeSync(this.seriesFd);

    const manifest = this._manifest(extra);
    this.manifestPath = path.join(this.dir, 'manifest.json');
    fs.writeFileSync(this.manifestPath, JSON.stringify(manifest));
    this.manifest = manifest;
    this.closed = true;
    return this.manifestPath;
  }

  _manifest(extra) {
    const MM = 1e3; // м -> мм
    const rfOut = new Float64Array(this.nrOut + 1);
    for (let i = 0; i <= this.nrOut; i++) rfOut[i] = this.rf[Math.min(this.nr, i * this.dsR)] * MM;
    const zfOut = new Float64Array(this.nzOut + 1);
    for (let j = 0; j <= this.nzOut; j++) zfOut[j] = this.zf[Math.min(this.nz, this.jGas0 + j * this.dsZ)] * MM;
    const rcOut = [], zcOut = [];
    for (let i = 0; i < this.nrOut; i++) rcOut.push(0.5 * (rfOut[i] + rfOut[i + 1]));
    for (let j = 0; j < this.nzOut; j++) zcOut.push(0.5 * (zfOut[j] + zfOut[j + 1]));

    const gasMaskZ = [];
    for (let j = 0; j < this.nz; j++) gasMaskZ.push(j >= this.jGas0 && j <= this.jGas1 ? 1 : 0);

    return {
      formatVersion: FORMAT_VERSION,
      runId: this.runId,
      level: this.level,
      createdAt: new Date().toISOString(),
      params: this.params,
      meta: this.meta,
      limits: {
        // честное ограничение модели — дублируется в README
        note: 'Осесимметричная модель описывает только ЦЕНТРАЛЬНЫЙ канал; ' +
              'внеосевой канал становится кольцом, азимутальные филаментационные моды отсутствуют принципиально.',
      },
      grid: {
        units: 'mm',
        nr: this.nr, nz: this.nz,
        jGas0: this.jGas0, jGas1: this.jGas1, nzGas: this.nzGas,
        dsR: this.dsR, dsZ: this.dsZ,
        nrOut: this.nrOut, nzOut: this.nzOut,
        rFaces: Array.from(rfOut), zFaces: Array.from(zfOut),
        rCenters: rcOut, zCenters: zcOut,
        gasMaskZ,
        dielectrics: [
          { name: 'diel_lo', z0: this.zf[0] * MM, z1: this.zf[this.jGas0] * MM },
          { name: 'gap',     z0: this.zf[this.jGas0] * MM, z1: this.zf[this.jGas1 + 1] * MM, gas: true },
          { name: 'diel_hi', z0: this.zf[this.jGas1 + 1] * MM, z1: this.zf[this.nz] * MM },
        ],
        electrodes: { powered_z: this.zf[0] * MM, ground_z: this.zf[this.nz] * MM },
      },
      frameLayout: {
        file: 'frames.bin',
        headerBytes: FRAMES_HEADER_BYTES,
        frameStride: this.frameStride,
        frameCount: this.frameCount,
        surfaces: {
          sigmaL: { offsetInFrame: this.sigLoOffset, dtype: 'f32', length: this.nrOut },
          sigmaR: { offsetInFrame: this.sigHiOffset, dtype: 'f32', length: this.nrOut },
        },
      },
      fields: this.fields.map((f) => ({
        name: f.name, dtype: f.dtype, map: f.map, decades: f.decades, v0rel: f.v0rel,
        offsetInFrame: f.offsetInFrame, byteLength: f.byteLength,
        shape: f.shape, layout: f.layout, unit: f.unit, label: f.label, present: f.present,
      })),
      frames: this.frames,
      series: {
        file: 'series.bin',
        headerBytes: SERIES_HEADER_BYTES,
        recordStride: SERIES_STRIDE,
        count: this.seriesCount,
        decimated: this.seriesCount < this.stepCount,
        steps: this.stepCount,
        names: SERIES_NAMES,
        dtypes: SERIES_NAMES.map((n) => (n === 't' ? 'f64' : 'f32')),
        offsets: SERIES_NAMES.map((n, k) => (k === 0 ? 0 : 8 + (k - 1) * 4)),
      },
      sampling: {
        dtFrameMin: this.dtFrameMin, dtFrameMax: this.dtFrameMax,
        dtFrameMinFinal: this.dtFrameMinEff,
        maxFrames: this.maxFrames, maxBytes: this.maxBytes,
        interestScale0: this.interestScale0, interestScaleFinal: this.interestScale,
        seriesMaxRecords: Number.isFinite(this.seriesMaxRecords) ? this.seriesMaxRecords : null,
      },
      stats: {
        steps: this.stepCount, frames: this.frameCount, seriesRecords: this.seriesCount,
        framesBytes: this.framesBytes, seriesBytes: this.seriesBytes,
        tFirst: this.tFirst, tLast: this.tLast,
        nonFiniteHits: this.nonFiniteHits,
      },
      ...extra,
    };
  }
}

/**
 * Пара рекордеров: полный (научный) + компактный (для плеера).
 * Полный  -> <dataDir>/<runId>/          (ряды БЕЗ прореживания)
 * Компакт -> <dataDir>/<runId>-compact/  (2x2 даунсэмплинг, uint8 лог-поля, ряды прорежены)
 */
export function createRecorders(opts = {}) {
  const dataDir = opts.dataDir || 'data';
  const runId = opts.runId || 'run';
  const base = { ...opts };
  delete base.dataDir; delete base.levels; delete base.perLevel;
  const levels = opts.levels || ['full', 'compact'];
  const perLevel = opts.perLevel || {};      // { full: {maxBytes}, compact: {...} }
  const recs = levels.map((level) => new Recorder({
    ...base, ...(perLevel[level] || {}), level,
    runId: level === 'full' ? runId : `${runId}-${level}`,
    dir: path.join(dataDir, level === 'full' ? runId : `${runId}-${level}`),
  }));
  return {
    recorders: recs,
    // true, если ХОТЯ БЫ один уровень запишет кадр -> продюсеру пора заполнять поля
    shouldRecord(state) { return recs.some((r) => r.shouldRecord(state)); },
    record(state, force = false) { return recs.map((r) => r.record(state, force)); },
    checkpoint(extra) { return recs.map((r) => r.checkpoint(extra)); },
    close(extra) { return recs.map((r) => r.close(extra)); },
  };
}

/**
 * fetch-подобная функция поверх файловой системы — чтобы loader.mjs (браузерный)
 * можно было прогнать в node (тесты, проверка прогона).
 */
export function makeFileFetch() {
  return async function fileFetch(url, init = {}) {
    let p = url;
    if (p.startsWith('file://')) p = decodeURIComponent(new URL(p).pathname);
    const range = init.headers && (init.headers.range || init.headers.Range);
    let buf;
    if (range) {
      const m = /bytes=(\d+)-(\d+)?/.exec(range);
      const start = Number(m[1]);
      const end = m[2] !== undefined ? Number(m[2]) : fs.statSync(p).size - 1;
      const len = end - start + 1;
      buf = Buffer.alloc(len);
      const fd = fs.openSync(p, 'r');
      fs.readSync(fd, buf, 0, len, start);
      fs.closeSync(fd);
    } else {
      buf = fs.readFileSync(p);
    }
    const ab = buf.buffer.slice(buf.byteOffset, buf.byteOffset + buf.byteLength);
    return {
      ok: true, status: range ? 206 : 200,
      headers: { get: (h) => (h.toLowerCase() === 'content-length' ? String(buf.byteLength) : null) },
      async arrayBuffer() { return ab; },
      async json() { return JSON.parse(buf.toString('utf8')); },
      async text() { return buf.toString('utf8'); },
      body: null,
    };
  };
}

export default Recorder;
