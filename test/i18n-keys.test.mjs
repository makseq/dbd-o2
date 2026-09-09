// test/i18n-keys.test.mjs — полнота и согласованность словарей ru/en.
// Запуск: node test/i18n-keys.test.mjs   (exit code 0 = все проверки прошли)
//
// Требование конвенции I18N.md §4: «Оба языка — полный набор ключей
// (проверка: одинаковые наборы ключей в ru и en, тест в test/)», плюс
// §4: «в I18N.en нет кириллицы».
//
// Тест не знает про конкретные единицы сайта: он сам находит все словари
// site/**/i18n.js и проверяет каждый. Новый словарь подхватывается сам.

import { readFileSync, readdirSync, statSync } from 'node:fs';
import { fileURLToPath, pathToFileURL } from 'node:url';
import { dirname, join, relative } from 'node:path';

const ROOT = join(dirname(fileURLToPath(import.meta.url)), '..');
const SITE = join(ROOT, 'site');

// Словари читают window.LANG на верхнем уровне — в node окна нет.
globalThis.window = globalThis.window || {};

let passed = 0, failed = 0;
const ok = (cond, name, detail) => {
  if (cond) { passed++; return; }
  failed++;
  console.error(`FAIL  ${name}${detail ? '\n      ' + detail : ''}`);
};

function findDicts(dir, out = []) {
  for (const e of readdirSync(dir)) {
    if (e === 'node_modules' || e.startsWith('.')) continue;
    const p = join(dir, e);
    if (statSync(p).isDirectory()) findDicts(p, out);
    else if (e === 'i18n.js') out.push(p);
  }
  return out;
}

const CYR = /\p{Script=Cyrillic}/u;

// Ключи, объявленные в исходнике: строки вида `  'some.key': …` на своей строке.
// Каждый ключ обязан встретиться РОВНО дважды — по разу в ru и в en.
// Один раз = ключа нет в одном из языков; три и больше = дубликат внутри языка
// (JS молча оставит последний, и объект Object.keys такую ошибку не покажет).
function sourceKeyCounts(src) {
  // Ключи бывают по нескольку в строке ('unit.kV': 'кВ', 'unit.V': 'В',),
  // поэтому ищем по всему тексту: начало строки/`{`/запятая, затем 'ключ':
  const counts = new Map();
  const re = /(?:^|[{,])[ \t\n]*'([^'\n]+)'[ \t]*:/gm;
  let m;
  while ((m = re.exec(src)) !== null) {
    re.lastIndex = m.index + m[0].length - 1;   // запятая может начинать следующий ключ
    counts.set(m[1], (counts.get(m[1]) || 0) + 1);
  }
  return counts;
}

const dicts = findDicts(SITE);
ok(dicts.length > 0, 'словари найдены', `в ${SITE} нет ни одного i18n.js`);

for (const file of dicts) {
  const name = relative(ROOT, file);
  const mod = await import(pathToFileURL(file).href);
  const I18N = mod.I18N;
  if (!I18N || !I18N.ru || !I18N.en) {
    ok(false, `${name}: экспортирует I18N { ru, en }`,
      'тест полноты ключей не может работать без экспорта I18N');
    continue;
  }

  const ru = Object.keys(I18N.ru), en = Object.keys(I18N.en);
  const missingEn = ru.filter((k) => !(k in I18N.en));
  const missingRu = en.filter((k) => !(k in I18N.ru));
  ok(missingEn.length === 0, `${name}: все ключи ru есть в en`,
    missingEn.length ? `нет в en: ${missingEn.join(', ')}` : '');
  ok(missingRu.length === 0, `${name}: все ключи en есть в ru`,
    missingRu.length ? `нет в ru: ${missingRu.join(', ')}` : '');
  ok(ru.length === en.length, `${name}: одинаковое число ключей`,
    `ru ${ru.length}, en ${en.length}`);

  // Тип значения обязан совпадать: если ru — шаблон-функция, а en — строка,
  // вызов t(key, params) в английской версии молча потеряет подстановки.
  const typeMismatch = ru.filter((k) => k in I18N.en && typeof I18N.ru[k] !== typeof I18N.en[k]);
  ok(typeMismatch.length === 0, `${name}: типы значений совпадают`,
    typeMismatch.map((k) => `${k}: ru ${typeof I18N.ru[k]} / en ${typeof I18N.en[k]}`).join('; '));

  // Дубликаты и односторонние ключи — по исходному тексту файла.
  // Вложенные объекты словаря (карты единиц и т.п.) в счёт не идут: берём
  // только те ключи, что реально попали в верхний уровень ru/en.
  const counts = new Map([...sourceKeyCounts(readFileSync(file, 'utf8'))]
    .filter(([k]) => k in I18N.ru || k in I18N.en));
  const dup = [...counts].filter(([, n]) => n > 2).map(([k, n]) => `${k} ×${n}`);
  ok(dup.length === 0, `${name}: нет дублирующихся ключей`, dup.join(', '));
  const lone = [...counts].filter(([k, n]) => n === 1 && (k in I18N.ru || k in I18N.en))
    .map(([k]) => k);
  ok(lone.length === 0, `${name}: каждый ключ объявлен в обоих языках`, lone.join(', '));
  // Сумма объявлений в исходнике обязана сойтись с числом ключей в объектах:
  // расхождение = либо дубликат (JS оставил последний), либо ключ, который
  // разбор исходника не увидел, — и тогда две проверки выше ничего не значат.
  const declared = [...counts.values()].reduce((a, b) => a + b, 0);
  ok(declared === ru.length + en.length, `${name}: разбор исходника сходится с объектом`,
    `в тексте ${declared} объявлений, в объектах ${ru.length + en.length}`);

  // Английский словарь не должен содержать кириллицы.
  const cyrEn = en.filter((k) => typeof I18N.en[k] === 'string' && CYR.test(I18N.en[k]));
  ok(cyrEn.length === 0, `${name}: в en нет кириллицы`,
    cyrEn.map((k) => `${k}: ${I18N.en[k]}`).join(' | '));

  // Пустые значения — почти всегда забытый перевод.
  const empty = [...ru, ...en].filter((k) => I18N.ru[k] === '' || I18N.en[k] === '');
  ok(empty.length === 0, `${name}: нет пустых строк`, [...new Set(empty)].join(', '));

  console.log(`  ${name}: ${ru.length} ключей, ru/en совпадают`);
}

console.log(`\ni18n-keys: ${passed} проверок прошло, ${failed} провалено`);
process.exit(failed ? 1 : 0);
