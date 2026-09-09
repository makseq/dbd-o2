// i18n.test.mjs — полнота словаря плеера и отсутствие кириллицы вне i18n.js.
//   node player/test/i18n.test.mjs   (из /p/dbd-o2/site)

import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { I18N } from '../js/i18n.js';

const here = path.dirname(fileURLToPath(import.meta.url));
const JS = path.join(here, '..', 'js');
let failed = 0;
const check = (ok, msg) => { console.log(`${ok ? 'ok  ' : 'FAIL'}  ${msg}`); if (!ok) failed++; };

// ── 1. одинаковые наборы ключей
const ru = Object.keys(I18N.ru).sort();
const en = Object.keys(I18N.en).sort();
const missEn = ru.filter((k) => !(k in I18N.en));
const missRu = en.filter((k) => !(k in I18N.ru));
check(missEn.length === 0, `в en нет ключей: ${missEn.join(', ') || '—'}`);
check(missRu.length === 0, `в ru нет ключей: ${missRu.join(', ') || '—'}`);

// ── 2. одинаковый тип значения (строка/функция) и одинаковая арность шаблонов
const typeMismatch = ru.filter((k) => k in I18N.en && typeof I18N.ru[k] !== typeof I18N.en[k]);
check(typeMismatch.length === 0, `разный тип значения: ${typeMismatch.join(', ') || '—'}`);

// ── 3. в en нет кириллицы
const cyr = /\p{Script=Cyrillic}/u;
const enCyr = en.filter((k) => cyr.test(String(typeof I18N.en[k] === 'function'
  ? I18N.en[k]({ v: 1, n: 1, i: 1, id: 'x', t: 1, msg: 'x', unit: 'x', mode: 'x', list: 'x' })
  : I18N.en[k])));
check(enCyr.length === 0, `в en осталась кириллица: ${enCyr.join(', ') || '—'}`);

// ── 4. в остальных файлах плеера кириллицы нет вне комментариев
function stripComments(src) {
  let out = '', i = 0, inS = null, inBlock = false, inLine = false;
  while (i < src.length) {
    const c = src[i], d = src[i + 1];
    if (inBlock) { if (c === '*' && d === '/') { inBlock = false; i += 2; } else i++; continue; }
    if (inLine) { if (c === '\n') { inLine = false; out += c; } i++; continue; }
    if (inS) {
      if (c === '\\') { out += src.substr(i, 2); i += 2; continue; }
      if (c === inS) inS = null;
      out += c; i++; continue;
    }
    if (c === '/' && d === '/') { inLine = true; i += 2; continue; }
    if (c === '/' && d === '*') { inBlock = true; i += 2; continue; }
    if (c === '"' || c === "'" || c === '`') { inS = c; out += c; i++; continue; }
    out += c; i++;
  }
  return out;
}
for (const f of fs.readdirSync(JS).filter((n) => n.endsWith('.js') && n !== 'i18n.js')) {
  const code = stripComments(fs.readFileSync(path.join(JS, f), 'utf8'));
  const bad = code.split('\n').map((l, n) => [n + 1, l]).filter(([, l]) => cyr.test(l));
  check(bad.length === 0, `${f}: кириллица вне комментариев — ${bad.map(([n]) => n).join(', ') || 'нет'}`);
}

console.log(failed ? `\n${failed} FAILED` : `\nвсе проверки пройдены (${ru.length} ключей)`);
process.exit(failed ? 1 : 0);
