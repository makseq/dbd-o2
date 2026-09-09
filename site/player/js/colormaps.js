// colormaps.js — палитры для сцены ДБР. Зависимость одна: словарь i18n.js.

import { t as tr } from './i18n.js';
//
// Каждая палитра — список опорных точек [pos, r, g, b]; из него строится
// LUT на 256 уровней (Uint8ClampedArray, 4 байта на уровень, alpha = 255).
// LUT кэшируется по имени: сцена дёргает их каждый кадр.

const STOPS = {
  // ── последовательные ───────────────────────────────────────────────────
  inferno: [
    [0.000, 0, 0, 4], [0.125, 20, 11, 52], [0.250, 58, 9, 99], [0.375, 96, 19, 110],
    [0.500, 133, 33, 107], [0.625, 172, 50, 88], [0.750, 206, 75, 54],
    [0.875, 231, 113, 18], [0.940, 246, 165, 26], [0.975, 250, 205, 66],
    [1.000, 252, 255, 164],
  ],
  magma: [
    [0.000, 0, 0, 4], [0.125, 19, 11, 52], [0.250, 59, 15, 112], [0.375, 101, 21, 110],
    [0.500, 140, 41, 129], [0.625, 183, 55, 121], [0.750, 222, 73, 104],
    [0.875, 247, 112, 92], [0.940, 253, 152, 110], [1.000, 252, 253, 191],
  ],
  plasma: [
    [0.000, 13, 8, 135], [0.125, 75, 3, 161], [0.250, 125, 3, 168], [0.375, 168, 34, 150],
    [0.500, 203, 70, 121], [0.625, 229, 107, 93], [0.750, 248, 148, 65],
    [0.875, 253, 195, 40], [1.000, 240, 249, 33],
  ],
  viridis: [
    [0.000, 68, 1, 84], [0.125, 72, 36, 117], [0.250, 65, 68, 135], [0.375, 53, 95, 141],
    [0.500, 42, 120, 142], [0.625, 33, 145, 140], [0.750, 39, 173, 129],
    [0.875, 92, 200, 99], [0.940, 170, 220, 50], [1.000, 253, 231, 37],
  ],
  turbo: [
    [0.000, 48, 18, 59], [0.100, 65, 88, 201], [0.200, 35, 150, 243], [0.300, 25, 199, 220],
    [0.400, 48, 229, 166], [0.500, 114, 246, 105], [0.600, 182, 242, 61],
    [0.700, 232, 209, 45], [0.800, 253, 155, 45], [0.900, 232, 86, 17], [1.000, 122, 4, 3],
  ],
  // холодная «фотонная» — для фотоионизации/фотоотлипания
  ice: [
    [0.000, 2, 5, 16], [0.180, 8, 40, 76], [0.360, 10, 84, 128], [0.540, 18, 140, 168],
    [0.720, 70, 198, 202], [0.880, 150, 235, 226], [1.000, 236, 255, 252],
  ],
  // «свечение» — почти чёрный ноль, затем фиолетово-оранжевый разряд
  glow: [
    [0.000, 2, 2, 8], [0.120, 24, 8, 48], [0.280, 62, 14, 96], [0.440, 118, 26, 110],
    [0.600, 176, 48, 90], [0.760, 224, 96, 50], [0.880, 248, 162, 46],
    [0.950, 253, 214, 110], [1.000, 255, 250, 214],
  ],

  // ── расходящиеся (0.5 — нейтральный «ноль») ────────────────────────────
  diverging: [
    [0.000, 96, 176, 255], [0.180, 46, 116, 226], [0.360, 26, 62, 140],
    [0.500, 22, 26, 34],
    [0.640, 140, 42, 46], [0.820, 214, 66, 54], [1.000, 255, 148, 118],
  ],
  // расходящаяся с яркими краями (для sigma поверх диэлектрика)
  sigma: [
    [0.000, 120, 196, 255], [0.250, 56, 128, 236], [0.460, 30, 48, 78],
    [0.500, 26, 32, 42],
    [0.540, 84, 36, 44], [0.750, 226, 78, 62], [1.000, 255, 176, 140],
  ],
};

export const COLORMAP_NAMES = Object.keys(STOPS);

export const COLORMAP_LABELS = {
  inferno: 'inferno', magma: 'magma', plasma: 'plasma', viridis: 'viridis',
  turbo: 'turbo', ice: 'cmap.ice', glow: 'cmap.glow',
  diverging: 'diverging (±)', sigma: 'sigma (±)',
};

/** Подпись палитры на активном языке (ключ словаря или готовое латинское имя). */
export function colormapLabel(name) {
  const v = COLORMAP_LABELS[name];
  if (!v) return name;
  return v.startsWith('cmap.') ? tr(v) : v;
}

const CACHE = new Map();

/** LUT палитры: Uint8ClampedArray(256*4), RGBA. */
export function getLUT(name) {
  const key = STOPS[name] ? name : 'inferno';
  let lut = CACHE.get(key);
  if (lut) return lut;
  const stops = STOPS[key];
  lut = new Uint8ClampedArray(256 * 4);
  let s = 0;
  for (let i = 0; i < 256; i++) {
    const t = i / 255;
    while (s < stops.length - 2 && t > stops[s + 1][0]) s++;
    const a = stops[s], b = stops[s + 1];
    const span = b[0] - a[0];
    const u = span > 0 ? Math.min(1, Math.max(0, (t - a[0]) / span)) : 0;
    lut[i * 4 + 0] = a[1] + (b[1] - a[1]) * u;
    lut[i * 4 + 1] = a[2] + (b[2] - a[2]) * u;
    lut[i * 4 + 2] = a[3] + (b[3] - a[3]) * u;
    lut[i * 4 + 3] = 255;
  }
  CACHE.set(key, lut);
  return lut;
}

/** Цвет в виде css rgb() для отметки t∈[0,1]. */
export function cssColor(name, t) {
  const lut = getLUT(name);
  const i = Math.max(0, Math.min(255, Math.round(t * 255))) * 4;
  return `rgb(${lut[i]},${lut[i + 1]},${lut[i + 2]})`;
}

/** Палитра по умолчанию для поля. */
export function defaultColormap(field) {
  switch (field) {
    case 'ionizRate': return 'glow';
    case 'rho': return 'diverging';
    case 'Emag': return 'viridis';
    case 'EN': return 'turbo';
    case 'photoIonRate': return 'ice';
    case 'photoDetachRate': return 'plasma';
    case 'n_O3m': return 'magma';
    case 'n_O3': return 'viridis';
    default: return 'inferno';
  }
}

/** Расходящаяся ли палитра (нужен симметричный диапазон и отметка 0). */
export function isDiverging(name) { return name === 'diverging' || name === 'sigma'; }

export default { getLUT, cssColor, defaultColormap, isDiverging, COLORMAP_NAMES, COLORMAP_LABELS };
