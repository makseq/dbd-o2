// i18n.js — словарь плеера 2D-ДБР (ru/en). Подключается ПЕРВЫМ: все остальные
// модули импортируют отсюда t() и не содержат ни одной пользовательской строки.
//
// Язык берётся из window.LANG, который выставляет инлайн-скрипт в <head>
// (?lang=ru|en -> localStorage.lang -> 'ru'). В node window нет -> 'ru'.
//
// Значения — строки или функции (шаблоны с параметрами).
// Ключи одинаковы в обоих языках; полнота проверяется тестом.

const RU = {
  // ── шапка / общее
  'doc.title': '2D-плеер осесимметричной модели | ДБР в кислороде',
  'doc.desc': 'Барьерный разряд в кислороде кадр за кадром: плотность электронов, поле, свечение и фотопроцессы в двумерной осесимметричной модели.',
  'brand.name': 'ДБР · O₂',
  'brand.sub': '2D (r,z)',
  'lang.aria': 'Язык / Language',
  'nav.home': '← Барьерный разряд в кислороде',

  'ctl.run': 'прогон',
  'ctl.data': 'данные',
  'ctl.data.compact': 'компакт',
  'ctl.data.full': 'полные',
  'ctl.speed': 'скорость',
  'ctl.colormap': 'палитра',

  'btn.compare': 'сравнить фото вкл/выкл',
  'btn.compare.title': 'Сравнение run-default и run-nophoto бок о бок',
  'btn.plots': ' графики',
  'btn.plots.title': 'Показать/скрыть колонку графиков',
  'btn.png.title': 'Сохранить текущий кадр как PNG',
  'btn.webm': '⏺ WebM',
  'btn.webm.title': 'Записать анимацию в WebM (MediaRecorder)',
  'btn.webm.stop': '⏹ стоп',
  'btn.csv.title': 'Выгрузить временные ряды в CSV',

  // ── панели графиков
  'panel.metrics': 'Метрики',
  'panel.wave': 'Осциллограммы U/I',
  'panel.wave.hint': 'клик = перемотка',
  'panel.liss': 'Фигура Лиссажу Q–V',
  'panel.liss.hint': 'фит vs аналитика',
  'panel.axial': 'Плотности вдоль оси r=0',
  'panel.axial.hint': 'клик — solo, alt — скрыть',
  'panel.radial': 'σ(r) на обеих поверхностях',
  'panel.photo': 'Фотопроцессы',
  'btn.logI': 'лог I',
  'btn.logI.title': 'логарифмическая ось тока',
  'btn.winzoom': 'окно=лупа',
  'btn.winzoom.title': 'окно графика = окно лупы таймлайна',
  'btn.sigasinh.title': 'знаковая лог-подобная ось',

  // ── боковая колонка
  'panel.field': 'Поле',
  'panel.scale': 'Шкала',
  'panel.readout': 'Мгновенные величины',
  'panel.runinfo': 'Прогон',
  'panel.limit': 'Ограничение модели',
  'scale.log': 'лог',
  'scale.lin': 'лин',
  'scale.frame': 'кадр',
  'scale.global': 'весь прогон',
  'opt.afterglow': 'послесвечение',
  'opt.afterglow.title': 'скорость затухания следа',
  'opt.sigma': 'накладки σ(r)',
  'opt.probe': 'зонд под курсором',

  // ── транспорт
  'tr.first': 'в начало (Home)',
  'tr.prev': 'кадр назад (←)',
  'tr.play': 'пуск/пауза (пробел)',
  'tr.next': 'кадр вперёд (→)',
  'tr.last': 'в конец (End)',
  'tr.a': 'начало участка (клавиша [)',
  'tr.b': 'конец участка (клавиша ])',
  'tr.loopclr': '✕ цикл',
  'tr.loopclr.title': 'сбросить участок (\\)',
  'tr.hint': 'верх — весь прогон (мкс), полоса под ним — плотность записи кадров · низ — лупа на импульс (нс, колесо = масштаб) · '
    + 'пробел пуск/пауза · ←/→ кадр · Shift+←/→ ×10 · [ ] цикл',

  // ── прогоны
  'run.default': 'база · 10 кВ · 10 кГц · фото ВКЛ',
  'run.nophoto': 'без фотопроцессов · 10 кВ · 10 кГц',
  'run.low': 'низкое U · 6 кВ (диффузный)',
  'run.high': 'высокое U · 14 кВ',
  'run.fast': 'высокая f · 10 кВ · 30 кГц',
  'run.synth': 'синтетика (отладка формата)',

  // ── единицы времени
  'unit.ms': 'мс',
  'unit.us': 'мкс',
  'unit.ns': 'нс',
  'unit.ps': 'пс',
  'unit.mm': 'мм',
  'unit.um': 'мкм',
  'unit.V': 'В',
  'unit.A': 'А',
  'unit.MB': 'МБ',

  // ── статус и загрузка
  'status.loading': 'загрузка…',
  'status.loadingRun': (p) => `${p.id}: ${p.mb} МБ`,
  'status.loadingRunOf': (p) => `${p.id}: ${p.mb} МБ / ${p.total} МБ`,
  'status.loaded': (p) => `${p.runId} · ${p.frames} кадров · сетка ${p.nr}×${p.nz} · ${p.tLast} мкс`,
  'status.pngSaved': (p) => `PNG сохранён · ${p.field} · ${p.t}`,
  'status.csv': (p) => `CSV: ${p.rows} строк × ${p.cols} колонок`,
  'status.webmUnsupported': 'MediaRecorder не поддерживается браузером',
  'status.webmRecording': 'запись WebM…',
  'status.webmDone': (p) => `WebM записан (${p.mb} МБ)`,
  'err.zeroFrames': (p) => `прогон ${p.id}: в manifest.json ноль кадров`,
  'err.noData': 'данные не загружены',
  'overlay.noData': 'Нет данных',
  'overlay.loadFailed': 'Данные не загрузились',
  'overlay.loadFailedText': (p) => `${p.err}\n\nПроверьте, что каталог ${p.dir} существует и содержит `
    + 'manifest.json / frames.bin / series.bin,\nи что страница открыта через HTTP-сервер '
    + '(file:// не работает из-за CORS).',

  // ── сцена
  'scene.photoOn': 'ФОТО ВКЛ',
  'scene.photoOff': 'ФОТО ВЫКЛ',
  'scene.noFrameData': 'нет данных',
  'scene.metalHV': 'МЕТАЛЛ · U(t)',
  'scene.metalGnd': 'МЕТАЛЛ · ⏚ 0 В',
  'scene.gas': 'газ O₂ · 1 атм · зазор 1.0 мм',
  'scene.axis0': 'ось r = 0',
  'scene.axisZ': 'z, мм',
  'scene.axisR': 'r, мм',
  'scene.profileOnAxis': 'профиль на оси r = 0',
  'scene.maxAt': (p) => `max = ${p.v} при z = ${p.z} мм`,
  'scene.rHalf': (p) => `r½ = ${p.v} мкм по r`,
  'scene.dielThick': (p) => `${p.v} мм`,
  'scene.scaleSync': ' · шкала синхр.',
  'scene.hudTime': (p) => `t = ${p.t} мкс · кадр ${p.i}/${p.n}`,
  'scene.afterglow': 'послесвечение',
  'scene.probePos': (p) => `r = ${p.r}  z = ${p.z} мм`,
  'scene.probeSigma': (p) => `   σ = ${p.v} Кл/м²`,

  // ── таймлайн
  'tl.whole': 'ВЕСЬ ПРОГОН',
  'tl.zoom': 'ЛУПА НА ИМПУЛЬС (масштаб окна)',
  'tl.uAmp': (p) => `U ±${p.v} В`,
  'tl.iAmp': (p) => `I ±${p.v} А (asinh)`,
  'tl.localNorm': 'нормировка по окну',
  'tl.frames': 'кадры',
  'tl.clock': (p) => `${p.t}  ·  кадр ${p.i}/${p.n}`,
  'tl.loop': (p) => `цикл ${p.a} … ${p.b}`,

  // ── чипы полей
  'field.chip.ionizRate': 'свечение',
  'field.chip.rho': 'заряд ρ',
  'field.chip.photoIonRate': 'фотоионизация',
  'field.chip.photoDetachRate': 'фотоотлипание',

  // ── названия полей (дублируют manifest.json, чтобы не зависеть от данных)
  'field.label.n_e': 'Плотность электронов',
  'field.label.ionizRate': 'Скорость ионизации (свечение)',
  'field.label.rho': 'Объёмный заряд',
  'field.label.Emag': '|E|',
  'field.label.EN': 'Приведённое поле E/N',
  'field.label.photoIonRate': 'Скорость фотоионизации',
  'field.label.photoDetachRate': 'Скорость фотоотлипания',
  'field.label.n_O3m': 'Плотность O3-',
  'field.label.n_O3': 'Плотность O3',
  'unit.m-3': 'м^-3',
  'unit.m-3s-1': 'м^-3 с^-1',
  'unit.C/m3': 'Кл/м^3',
  'unit.V/m': 'В/м',
  'unit.Td': 'Тд',

  // ── информация о прогоне
  'info.params': (p) => `U₀ = <b>${p.U0} кВ</b>, f = <b>${p.f} кГц</b>`,
  'info.photo': (p) => `фотопроцессы: <b>${p.on ? 'ВКЛ' : 'ВЫКЛ'}</b>`,
  'info.grid': (p) => `сетка ${p.nr}×${p.nz} (запись ${p.nrOut}×${p.nzOut})`,
  'info.counts': (p) => `кадров ${p.frames}, шагов ${p.steps}`,
  'info.tspan': (p) => `t ∈ [0, ${p.tLast} мкс]`,
  'info.crash': (p) => `⚠ прогон оборван на пробое (t = ${p.t} мкс):`,
  // Полный текст аварии из manifest.json (в ru-версии на экран идёт сам исходник, ключ — для парности наборов).
  'crash.stepRejected': (p) => `DBD2D: шаг ${p.step} (t = ${p.t} с) не проходит приёмку даже при `
    + `dt = ${p.dt} с. Причина: E/N = ${p.en} Тд. Типовой источник — разгон объёмного заряда `
    + 'в ПРИСТЕНОЧНОЙ ячейке. ВНИМАНИЕ: прежняя формулировка «лечится сгущением z-сетки» ОПРОВЕРГНУТА '
    + 'измерением (docs/DIVERGENCE_ANALYSIS §8.3: dz_wall от 15.6 до 0.184 мкм — срыв всегда, момент '
    + 'сдвигается на 6 %). Действующий диагноз (§8.4): тепловой поток электронов ¼v_th·n_e на КАТОДНЫЙ '
    + 'барьер не подавлялся отталкивающим полем, ток на барьер менял знак, и петля σ↓→E↑→ионизация↑ '
    + "уходила в разнос. Лечится ГУ Хагелаара (wallBC: 'hagelaar', дефолт) — дрейфовый член входит "
    + "со знаком (2a−1). Если авария всё же случилась при wallBC = 'hagelaar', это НЕ известный дефект: "
    + 'сгущение сетки и уменьшение dt здесь не помогают, нужна диагностика. Состояние откачено к началу '
    + 'шага и пригодно для чекпойнта.',
  'limit.default': 'Осесимметричная модель описывает только центральный канал.',
  'limit.axisymmetric': 'Осесимметричная модель описывает только ЦЕНТРАЛЬНЫЙ канал; внеосевой канал '
    + 'становится кольцом, азимутальные филаментационные моды отсутствуют принципиально.',
  'meta.syntheticWarning': 'СИНТЕТИЧЕСКИЙ ПРОГОН: правдоподобная форма сигнала, НЕ результат решения '
    + 'уравнений. Пригоден только для разработки плеера.',

  // ── строка мгновенных величин
  'ro.t': 't',
  'ro.Uapp': 'U(t), В',
  'ro.Ugap': 'U на зазоре, В',
  'ro.Icond': 'I пров., А',
  'ro.Idisp': 'I смещ., А',
  'ro.Itot': 'I полн., А',
  'ro.maxEN': 'max E/N, Тд',
  'ro.sigmaMax': 'σ max, Кл/м²',
  'ro.o3ppm': 'O₃, ppm',
  'ro.photoL': 'фотоэмиссия L',
  'ro.photoR': 'фотоэмиссия R',

  // ── метрики / статус
  'metrics.srcLeft': '← run-default (левая сцена)',
  'metrics.photoIntegrating': 'интегрирование по объёму…',
  'metrics.photoMissing': (p) => `нет каналов: ${p.list}`,
  'metrics.photoRate': 'темп по объёму газа, эв./с',
  'metrics.error': (p) => `ошибка: ${p.msg}`,
  'metrics.note.cumulative': 'накопл. / V_gas',
  'metrics.note.medianDQ': 'медиана dQ/dU',
  'metrics.note.loopOpen': '⚠ петля не замкнута',
  'metrics.note.perPeriod': 'за период',
  'metrics.limiterLiss': 'Лиссажу ⚠',

  // ── режимы кнопок
  'liss.period': 'период',
  'liss.last8': 'последние 8',
  'liss.all': 'все',
  'ax.axis': 'ось',
  'ax.mean': 'средн. по r',
  'export.caption': (p) => `${p.tag} · ${p.field} · t = ${p.t} · ${p.scaleKind} · шкала: ${p.scope}`,

  // ── графики (plots.js)
  'plots.noSeries': 'нет series.bin',
  'plots.noQU': 'в series нет Q/Uapp',
  'plots.noFrame': 'нет кадра',
  'plots.noPhotoSeries': 'фоторяды не построены',
  'plots.logNote': 'log: |I|, синяя полоса = I < 0',
  'plots.vsAnalytic': '  vs аналит.',
  'plots.r2sides': 'R² сторон',
  'plots.closure': 'замыкание',
  'plots.sigmaLo': (p) => `σ_lo (z=0.5 мм)  ось: ${p.v} nC/cm²`,
  'plots.sigmaHi': (p) => `σ_hi (z=1.5 мм)  ось: ${p.v} nC/cm²`,
  'plots.asinhNote': 'asinh-ось: линейна у нуля',
  'plots.clipped': (p) => `${p.n} точек ≤ 0 или ниже ${p.lo} — не показаны`,
  'plots.axialHint': 'клик — solo · alt+клик — скрыть · × — сорт не записан в прогоне',
  'plots.photoHint': '× — канал не записан в контейнере · клик — скрыть/показать',
  'plots.photoAxis': 'темп, событий/с (по объёму газа)',
  'plots.photoIonRate': (p) => `фотоионизация ${p.v}`,
  'plots.zeros': (p) => `${p.n} нулей — на лог-оси разрыв, а не ноль`,
  'plots.axisT': (p) => `t, ${p.unit}`,
  'plots.axisU': 'U, kV',
  'plots.axisI': 'I, mA',
  'plots.axisIlog': '|I|, A (log)',
  'plots.axisUapp': 'U_app, kV',
  'plots.axisQ': (p) => `Q, ${p.unit}`,
  'plots.axisR': 'r, mm',
  'plots.axisSigma': 'σ, nC/cm²',
  'plots.axisZ': 'z, mm',
  'plots.axisN': 'n, m⁻³',

  // ── фит Лиссажу: предупреждения
  'fit.shortPeriod': 'период не набран: точек < 16',
  'fit.noSegments': 'нет пригодных сегментов dQ/dU',
  'fit.noFourSides': 'нет четырёх сторон: петля не параллелограмм',
  'fit.loopOpen': 'петля не замкнута',
  'fit.lowR2': (p) => `R² разрядных сторон = ${p.v} < 0.98`,
  'fit.fragmented': (p) => `ёмкостная ветвь раздроблена микроразрядами (разброс b ${p.v} % от ΔQ)`,
  'fit.notSeparated': (p) => `ветви не разделены (×${p.v})`,
  'fit.fewSides': 'найдено < 4 сторон',
  'fit.deviation': (p) => `расхождение с аналитикой ${p.v} %`,
  'fit.unreliable': ' — наклоны Мэнли недостоверны',

  // ── фотоканалы
  'chan.photoIon': 'фотоионизация',
  'chan.photoEmit': 'фотоэмиссия (пов.)',
  'chan.photoDetach': 'фотоотлипание',
  'chan.impactIon': 'ударная ионизация',
  'chan.collDetach': 'столкн. отлипание',
  'chan.collDetach.note': 'не записано в контейнере (recorder не пишет detachRate)',

  // ── палитры
  'cmap.ice': 'ice (холодная)',
  'cmap.glow': 'glow (свечение)',

  // ── ошибки загрузчика (диагностика разработчика)
  'loader.badMap': (p) => `loader: неизвестный map '${p.map}'`,
  'loader.badMagic': (p) => `loader: неверная сигнатура '${p.got}', ожидалась '${p.want}'`,
  'loader.noFetch': 'loader: нет fetch — передайте opts.fetch',
  'loader.manifestHttp': (p) => `loader: manifest.json -> HTTP ${p.status} (прогон не завершён?)`,
  'loader.badVersion': (p) => `loader: версия формата ${p.got}, поддерживается ${p.want}`,
  'loader.framesSize': (p) => `loader: frames.bin ${p.got} Б, по манифесту ${p.want} Б`,
  'loader.framesCount': (p) => `loader: заголовок frames.bin: ${p.got} кадров, манифест: ${p.want}`,
  'loader.seriesMismatch': (p) => `loader: series.bin рассинхронизирован с манифестом (${p.got} против ${p.want})`,
  'loader.frameOutOfFile': (p) => `loader: кадр ${p.i} выходит за пределы frames.bin`,
  'loader.frameNotLoaded': (p) => `loader: кадр ${p.i} ещё не загружен — вызовите await getFrameAsync(${p.i}, ...) или prefetch()`,
  'loader.rangeHttp': (p) => `loader: Range-запрос кадра ${p.i} -> HTTP ${p.status}`,
  'loader.rangeShort': (p) => `loader: сервер вернул ${p.got} Б вместо ${p.want} (Range не поддержан?)`,
  'loader.noFrame': (p) => `loader: кадра ${p.i} нет (всего ${p.n})`,
  'loader.noField': (p) => `loader: поля '${p.field}' нет в прогоне`,
  'plots.noCanvas': 'plots: не передан canvas',
  'metrics.noContainer': 'MetricsPanel: не передан контейнер',

  // ── plots-demo.html
  'demo.title': 'plots-demo · графики плеера 2D-ДБР',
  'demo.sub': '2D ДБР · O₂ · осесимметрия',
  'demo.run.synth': 'run-synth-compact (2 периода, синтетика)',
  'demo.run.default': 'run-default-compact (реальный, обрыв на пробое)',
  'demo.run.nophoto': 'run-nophoto-compact',
  'demo.run.low': 'run-low-compact (6 кВ)',
  'demo.run.high': 'run-high-compact (14 кВ)',
  'demo.run.fast': 'run-fast-compact (30 кГц)',
  'demo.pause': '⏸ пауза',
  'demo.play': '▶ пуск',
  'demo.liss': (p) => `Лиссажу: ${p.mode}`,
  'demo.hWave': 'D · Waveforms',
  'demo.hWave.hint': 'min/max-декимация · клик = перемотка',
  'demo.hLiss': 'E · Q–V Лиссажу',
  'demo.hPhoto': 'Фотопроцессы',
  'demo.hPhoto.hint': 'объёмный темп, лог-ось',
  'demo.hAxial': 'F · Плотности вдоль оси',
  'demo.hAxial.hint': 'клик — solo · alt+клик — скрыть',
  'demo.hRadial': 'σ(r) на обеих поверхностях',
  'demo.loading': 'загрузка…',
  'demo.loadingRun': (p) => `загрузка ${p.id}…`,
  'demo.progress': (p) => `${p.id}: ${p.pct} %`,
  'demo.synthWarn': '⚠ СИНТЕТИЧЕСКИЙ ПРОГОН — не результат решения уравнений',
  'demo.photoRows': (p) => `фоторяды: ${p.n} кадров`,
  'demo.clock': (p) => `t = ${p.t}   кадр ${p.i}/${p.n}   период ${p.period}`,
  'demo.frozen': (p) => `   ⚠ кадры кончились на ${p.t} — поля заморожены, ряды идут дальше`,
  'demo.error': (p) => `ОШИБКА: ${p.msg}`,
  'demo.loadFailed': (p) => `не удалось загрузить прогон: ${p.msg}`,
  'demo.serveHint': 'запускать через HTTP-сервер из корня /p/dbd-o2',
};

const EN = {
  'doc.title': '2D axisymmetric player, frame by frame | DBD in oxygen',
  'doc.desc': 'Watch a barrier discharge in oxygen frame by frame: electron density, electric field, glow and photoprocesses in a 2D axisymmetric model.',
  'brand.name': 'DBD · O₂',
  'brand.sub': '2D (r,z)',
  'lang.aria': 'Language',
  'nav.home': '← Barrier discharge in oxygen',

  'ctl.run': 'run',
  'ctl.data': 'data',
  'ctl.data.compact': 'compact',
  'ctl.data.full': 'full',
  'ctl.speed': 'speed',
  'ctl.colormap': 'colormap',

  'btn.compare': 'compare photo on/off',
  'btn.compare.title': 'Show run-default and run-nophoto side by side',
  'btn.plots': ' plots',
  'btn.plots.title': 'Show or hide the plots column',
  'btn.png.title': 'Save the current frame as PNG',
  'btn.webm': '⏺ WebM',
  'btn.webm.title': 'Record the animation to WebM (MediaRecorder)',
  'btn.webm.stop': '⏹ stop',
  'btn.csv.title': 'Export the time series to CSV',

  'panel.metrics': 'Metrics',
  'panel.wave': 'U/I waveforms',
  'panel.wave.hint': 'click to seek',
  'panel.liss': 'Q–V Lissajous figure',
  'panel.liss.hint': 'fit vs analytic',
  'panel.axial': 'Densities along the r = 0 axis',
  'panel.axial.hint': 'click to solo, alt to hide',
  'panel.radial': 'σ(r) on both surfaces',
  'panel.photo': 'Photoprocesses',
  'btn.logI': 'log I',
  'btn.logI.title': 'Logarithmic current axis',
  'btn.winzoom': 'window = zoom',
  'btn.winzoom.title': 'Plot window follows the timeline zoom window',
  'btn.sigasinh.title': 'Signed, log-like axis',

  'panel.field': 'Field',
  'panel.scale': 'Scale',
  'panel.readout': 'Instantaneous values',
  'panel.runinfo': 'Run',
  'panel.limit': 'Model limitation',
  'scale.log': 'log',
  'scale.lin': 'lin',
  'scale.frame': 'frame',
  'scale.global': 'whole run',
  'opt.afterglow': 'afterglow',
  'opt.afterglow.title': 'Trail decay rate',
  'opt.sigma': 'σ(r) overlays',
  'opt.probe': 'probe under cursor',

  'tr.first': 'To start (Home)',
  'tr.prev': 'Previous frame (←)',
  'tr.play': 'Play / pause (space)',
  'tr.next': 'Next frame (→)',
  'tr.last': 'To end (End)',
  'tr.a': 'Range start (key [)',
  'tr.b': 'Range end (key ])',
  'tr.loopclr': '✕ loop',
  'tr.loopclr.title': 'Clear the range (\\)',
  'tr.hint': 'top — whole run (µs), the strip below it — frame recording density · bottom — pulse zoom (ns, wheel = zoom) · '
    + 'space play/pause · ←/→ frame · Shift+←/→ ×10 · [ ] loop',

  'run.default': 'baseline · 10 kV · 10 kHz · photo ON',
  'run.nophoto': 'no photoprocesses · 10 kV · 10 kHz',
  'run.low': 'low U · 6 kV (diffuse)',
  'run.high': 'high U · 14 kV',
  'run.fast': 'high f · 10 kV · 30 kHz',
  'run.synth': 'synthetic (format debugging)',

  'unit.ms': 'ms',
  'unit.us': 'µs',
  'unit.ns': 'ns',
  'unit.ps': 'ps',
  'unit.mm': 'mm',
  'unit.um': 'µm',
  'unit.V': 'V',
  'unit.A': 'A',
  'unit.MB': 'MB',

  'status.loading': 'loading…',
  'status.loadingRun': (p) => `${p.id}: ${p.mb} MB`,
  'status.loadingRunOf': (p) => `${p.id}: ${p.mb} MB / ${p.total} MB`,
  'status.loaded': (p) => `${p.runId} · ${p.frames} ${plural(p.frames, ['frame', 'frames'])} · `
    + `grid ${p.nr}×${p.nz} · ${p.tLast} µs`,
  'status.pngSaved': (p) => `PNG saved · ${p.field} · ${p.t}`,
  'status.csv': (p) => `CSV: ${p.rows} ${plural(p.rows, ['row', 'rows'])} × `
    + `${p.cols} ${plural(p.cols, ['column', 'columns'])}`,
  'status.webmUnsupported': 'MediaRecorder is not supported by this browser',
  'status.webmRecording': 'recording WebM…',
  'status.webmDone': (p) => `WebM recorded (${p.mb} MB)`,
  'err.zeroFrames': (p) => `run ${p.id}: manifest.json has zero frames`,
  'err.noData': 'data not loaded',
  'overlay.noData': 'No data',
  'overlay.loadFailed': 'Could not load the data',
  'overlay.loadFailedText': (p) => `${p.err}\n\nCheck that the folder ${p.dir} exists and contains `
    + 'manifest.json / frames.bin / series.bin,\nand that the page is served over HTTP '
    + '(file:// fails because of CORS).',

  'scene.photoOn': 'PHOTO ON',
  'scene.photoOff': 'PHOTO OFF',
  'scene.noFrameData': 'no data',
  'scene.metalHV': 'METAL · U(t)',
  'scene.metalGnd': 'METAL · ⏚ 0 V',
  'scene.gas': 'O₂ gas · 1 atm · 1.0 mm gap',
  'scene.axis0': 'r = 0 axis',
  'scene.axisZ': 'z, mm',
  'scene.axisR': 'r, mm',
  'scene.profileOnAxis': 'profile on the r = 0 axis',
  'scene.maxAt': (p) => `max = ${p.v} at z = ${p.z} mm`,
  'scene.rHalf': (p) => `r½ = ${p.v} µm in r`,
  'scene.dielThick': (p) => `${p.v} mm`,
  'scene.scaleSync': ' · scale synced',
  'scene.hudTime': (p) => `t = ${p.t} µs · frame ${p.i}/${p.n}`,
  'scene.afterglow': 'afterglow',
  'scene.probePos': (p) => `r = ${p.r}  z = ${p.z} mm`,
  'scene.probeSigma': (p) => `   σ = ${p.v} C/m²`,

  'tl.whole': 'WHOLE RUN',
  'tl.zoom': 'PULSE ZOOM (window scale)',
  'tl.uAmp': (p) => `U ±${p.v} V`,
  'tl.iAmp': (p) => `I ±${p.v} A (asinh)`,
  'tl.localNorm': 'normalized to window',
  'tl.frames': 'frames',
  'tl.clock': (p) => `${p.t}  ·  frame ${p.i}/${p.n}`,
  'tl.loop': (p) => `loop ${p.a} … ${p.b}`,

  'field.chip.ionizRate': 'glow',
  'field.chip.rho': 'charge ρ',
  'field.chip.photoIonRate': 'photoionization',
  'field.chip.photoDetachRate': 'photodetachment',

  'field.label.n_e': 'Electron density',
  'field.label.ionizRate': 'Ionization rate (glow)',
  'field.label.rho': 'Space charge',
  'field.label.Emag': '|E|',
  'field.label.EN': 'Reduced field E/N',
  'field.label.photoIonRate': 'Photoionization rate',
  'field.label.photoDetachRate': 'Photodetachment rate',
  'field.label.n_O3m': 'O3- density',
  'field.label.n_O3': 'O3 density',
  'unit.m-3': 'm^-3',
  'unit.m-3s-1': 'm^-3 s^-1',
  'unit.C/m3': 'C/m^3',
  'unit.V/m': 'V/m',
  'unit.Td': 'Td',

  'info.params': (p) => `U₀ = <b>${p.U0} kV</b>, f = <b>${p.f} kHz</b>`,
  'info.photo': (p) => `photoprocesses: <b>${p.on ? 'ON' : 'OFF'}</b>`,
  'info.grid': (p) => `grid ${p.nr}×${p.nz} (recorded ${p.nrOut}×${p.nzOut})`,
  'info.counts': (p) => `${p.frames} ${plural(p.frames, ['frame', 'frames'])}, `
    + `${p.steps} ${plural(p.steps, ['step', 'steps'])}`,
  'info.tspan': (p) => `t ∈ [0, ${p.tLast} µs]`,
  'info.crash': (p) => `⚠ run aborted at breakdown (t = ${p.t} µs):`,
  'crash.stepRejected': (p) => `DBD2D: step ${p.step} (t = ${p.t} s) fails acceptance even at `
    + `dt = ${p.dt} s. Cause: E/N = ${p.en} Td. Typical source — space-charge runaway in the `
    + 'NEAR-WALL cell. WARNING: the earlier claim that it is “cured by refining the z-grid” is DISPROVED '
    + 'by measurement (docs/DIVERGENCE_ANALYSIS §8.3: dz_wall from 15.6 to 0.184 µm — the blow-up always '
    + 'occurs, the timing shifts by 6 %). Current diagnosis (§8.4): the electron '
    + 'thermal flux ¼v_th·n_e onto the CATHODE barrier was not suppressed by the repelling field, the '
    + 'current into the barrier changed sign, and the loop σ↓→E↑→ionization↑ ran away. Cured by the '
    + "Hagelaar boundary condition (wallBC: 'hagelaar', the default) — the drift term enters with sign "
    + "(2a−1). If the crash still happens with wallBC = 'hagelaar', it is NOT a known defect: refining "
    + 'the grid and reducing dt do not help here, diagnostics are needed. The state was rolled back to '
    + 'the start of the step and is usable as a checkpoint.',
  'limit.default': 'The axisymmetric model describes the central channel only.',
  'limit.axisymmetric': 'The axisymmetric model describes the CENTRAL channel only; an off-axis channel '
    + 'becomes a ring, and azimuthal filamentation modes are absent by construction.',
  'meta.syntheticWarning': 'SYNTHETIC RUN: a plausible waveform, NOT a solution of the equations. '
    + 'Useful only for developing the player.',

  'ro.t': 't',
  'ro.Uapp': 'U(t), V',
  'ro.Ugap': 'U across gap, V',
  'ro.Icond': 'I cond., A',
  'ro.Idisp': 'I displ., A',
  'ro.Itot': 'I total, A',
  'ro.maxEN': 'max E/N, Td',
  'ro.sigmaMax': 'σ max, C/m²',
  'ro.o3ppm': 'O₃, ppm',
  'ro.photoL': 'photoemission L',
  'ro.photoR': 'photoemission R',

  'metrics.srcLeft': '← run-default (left scene)',
  'metrics.photoIntegrating': 'integrating over volume…',
  'metrics.photoMissing': (p) => `channels missing: ${p.list}`,
  'metrics.photoRate': 'rate over the gas volume, events/s',
  'metrics.error': (p) => `error: ${p.msg}`,
  'metrics.note.cumulative': 'cumulative / V_gas',
  'metrics.note.medianDQ': 'median dQ/dU',
  'metrics.note.loopOpen': '⚠ loop not closed',
  'metrics.note.perPeriod': 'per period',
  'metrics.limiterLiss': 'Lissajous ⚠',

  'liss.period': 'period',
  'liss.last8': 'last 8',
  'liss.all': 'all',
  'ax.axis': 'axis',
  'ax.mean': 'mean over r',
  'export.caption': (p) => `${p.tag} · ${p.field} · t = ${p.t} · ${p.scaleKind} · scale: ${p.scope}`,

  'plots.noSeries': 'no series.bin',
  'plots.noQU': 'series has no Q/Uapp',
  'plots.noFrame': 'no frame',
  'plots.noPhotoSeries': 'photo series not built',
  'plots.logNote': 'log: |I|, blue band = I < 0',
  'plots.vsAnalytic': '  vs analytic',
  'plots.r2sides': 'R² of sides',
  'plots.closure': 'closure',
  'plots.sigmaLo': (p) => `σ_lo (z = 0.5 mm)  axis: ${p.v} nC/cm²`,
  'plots.sigmaHi': (p) => `σ_hi (z = 1.5 mm)  axis: ${p.v} nC/cm²`,
  'plots.asinhNote': 'asinh axis: linear near zero',
  'plots.clipped': (p) => `${p.n} ${plural(p.n, ['point', 'points'])} ≤ 0 or below ${p.lo} — not shown`,
  'plots.axialHint': 'click to solo · alt+click to hide · × — species not recorded in the run',
  'plots.photoHint': '× — channel not recorded in the container · click to hide or show',
  'plots.photoAxis': 'rate, events/s (over the gas volume)',
  'plots.photoIonRate': (p) => `photoionization ${p.v}`,
  'plots.zeros': (p) => `${p.n} ${plural(p.n, ['zero', 'zeros'])} — a gap on the log axis, not a zero`,
  'plots.axisT': (p) => `t, ${p.unit}`,
  'plots.axisU': 'U, kV',
  'plots.axisI': 'I, mA',
  'plots.axisIlog': '|I|, A (log)',
  'plots.axisUapp': 'U_app, kV',
  'plots.axisQ': (p) => `Q, ${p.unit}`,
  'plots.axisR': 'r, mm',
  'plots.axisSigma': 'σ, nC/cm²',
  'plots.axisZ': 'z, mm',
  'plots.axisN': 'n, m⁻³',

  'fit.shortPeriod': 'period incomplete: fewer than 16 points',
  'fit.noSegments': 'no usable dQ/dU segments',
  'fit.noFourSides': 'four sides not found: the loop is not a parallelogram',
  'fit.loopOpen': 'loop not closed',
  'fit.lowR2': (p) => `R² of the discharge sides = ${p.v} < 0.98`,
  'fit.fragmented': (p) => `capacitive branch broken up by microdischarges (b spread ${p.v} % of ΔQ)`,
  'fit.notSeparated': (p) => `branches not separated (×${p.v})`,
  'fit.fewSides': 'fewer than 4 sides found',
  'fit.deviation': (p) => `deviation from the analytic value ${p.v} %`,
  'fit.unreliable': ' — Manley slopes unreliable',

  'chan.photoIon': 'photoionization',
  'chan.photoEmit': 'photoemission (surface)',
  'chan.photoDetach': 'photodetachment',
  'chan.impactIon': 'impact ionization',
  'chan.collDetach': 'collisional detachment',
  'chan.collDetach.note': 'not recorded in the container (recorder writes no detachRate)',

  'cmap.ice': 'ice (cool)',
  'cmap.glow': 'glow',

  'loader.badMap': (p) => `loader: unknown map '${p.map}'`,
  'loader.badMagic': (p) => `loader: bad signature '${p.got}', expected '${p.want}'`,
  'loader.noFetch': 'loader: no fetch — pass opts.fetch',
  'loader.manifestHttp': (p) => `loader: manifest.json -> HTTP ${p.status} (run not finished?)`,
  'loader.badVersion': (p) => `loader: format version ${p.got}, supported ${p.want}`,
  'loader.framesSize': (p) => `loader: frames.bin is ${p.got} B, manifest says ${p.want} B`,
  'loader.framesCount': (p) => `loader: frames.bin header says ${p.got} `
    + `${plural(p.got, ['frame', 'frames'])}, manifest says ${p.want}`,
  'loader.seriesMismatch': (p) => `loader: series.bin out of sync with the manifest (${p.got} vs ${p.want})`,
  'loader.frameOutOfFile': (p) => `loader: frame ${p.i} lies past the end of frames.bin`,
  'loader.frameNotLoaded': (p) => `loader: frame ${p.i} is not loaded yet — call await getFrameAsync(${p.i}, ...) or prefetch()`,
  'loader.rangeHttp': (p) => `loader: range request for frame ${p.i} -> HTTP ${p.status}`,
  'loader.rangeShort': (p) => `loader: server returned ${p.got} B instead of ${p.want} (Range unsupported?)`,
  'loader.noFrame': (p) => `loader: there is no frame ${p.i} (${p.n} in total)`,
  'loader.noField': (p) => `loader: field '${p.field}' is not in this run`,
  'plots.noCanvas': 'plots: no canvas given',
  'metrics.noContainer': 'MetricsPanel: no container given',

  'demo.title': 'plots-demo · plots of the 2D DBD player',
  'demo.sub': '2D DBD · O₂ · axisymmetric',
  'demo.run.synth': 'run-synth-compact (2 periods, synthetic)',
  'demo.run.default': 'run-default-compact (real, aborted at breakdown)',
  'demo.run.nophoto': 'run-nophoto-compact',
  'demo.run.low': 'run-low-compact (6 kV)',
  'demo.run.high': 'run-high-compact (14 kV)',
  'demo.run.fast': 'run-fast-compact (30 kHz)',
  'demo.pause': '⏸ pause',
  'demo.play': '▶ play',
  'demo.liss': (p) => `Lissajous: ${p.mode}`,
  'demo.hWave': 'D · Waveforms',
  'demo.hWave.hint': 'min/max decimation · click to seek',
  'demo.hLiss': 'E · Q–V Lissajous',
  'demo.hPhoto': 'Photoprocesses',
  'demo.hPhoto.hint': 'volume rate, log axis',
  'demo.hAxial': 'F · Densities along the axis',
  'demo.hAxial.hint': 'click to solo · alt+click to hide',
  'demo.hRadial': 'σ(r) on both surfaces',
  'demo.loading': 'loading…',
  'demo.loadingRun': (p) => `loading ${p.id}…`,
  'demo.progress': (p) => `${p.id}: ${p.pct} %`,
  'demo.synthWarn': '⚠ SYNTHETIC RUN — not a solution of the equations',
  'demo.photoRows': (p) => `photo series: ${p.n} ${plural(p.n, ['frame', 'frames'])}`,
  'demo.clock': (p) => `t = ${p.t}   frame ${p.i}/${p.n}   period ${p.period}`,
  'demo.frozen': (p) => `   ⚠ frames ran out at ${p.t} — fields are frozen, series continue`,
  'demo.error': (p) => `ERROR: ${p.msg}`,
  'demo.loadFailed': (p) => `could not load the run: ${p.msg}`,
  'demo.serveHint': 'serve it over HTTP from the /p/dbd-o2 root',
};

export const I18N = { ru: RU, en: EN };

function detect() {
  if (typeof window === 'undefined') return 'en';   // конвенция: английский по умолчанию
  const l = window.LANG;
  return (l === 'en' || l === 'ru') ? l : 'en';
}

export const LANG = detect();
export const L = I18N[LANG];

/** Строка словаря. Отсутствующий ключ падает на ru, затем на сам ключ. */
export function t(key, params) {
  let v = L[key];
  if (v === undefined) v = RU[key];
  if (v === undefined) return key;
  return typeof v === 'function' ? v(params || {}) : v;
}

/** Множественное число: ru — три формы, en — две. */
export function plural(n, forms) {
  if (LANG === 'en') return Math.abs(n) === 1 ? forms[0] : forms[1];
  const a = Math.abs(n) % 100, b = a % 10;
  if (a > 10 && a < 20) return forms[2];
  if (b > 1 && b < 5) return forms[1];
  if (b === 1) return forms[0];
  return forms[2];
}

/** Локализованные подпись и единица поля (данные лежат в manifest.json по-русски). */
export function fieldLabel(spec) {
  if (!spec) return '';
  const k = `field.label.${spec.name}`;
  return (L[k] !== undefined || RU[k] !== undefined) ? t(k) : (spec.label || spec.name);
}
const UNIT_KEY = {
  'м^-3': 'unit.m-3', 'м^-3 с^-1': 'unit.m-3s-1', 'Кл/м^3': 'unit.C/m3',
  'В/м': 'unit.V/m', 'Тд': 'unit.Td',
};
export function fieldUnit(spec) {
  if (!spec || !spec.unit) return '';
  const k = UNIT_KEY[spec.unit];
  return k ? t(k) : spec.unit;
}
/** Текст, пришедший из manifest.json: переводится, если он нам знаком. */
const DATA_KEY = {};
for (const key of ['limit.axisymmetric', 'meta.syntheticWarning']) DATA_KEY[RU[key]] = key;
// Авария солвера: одинаковый шаблон во всех прогонах, числа разные.
const CRASH_RE = /^DBD2D: шаг (\d+) \(t = ([^)]*?) с\) не проходит приёмку даже при dt = (\S+) с\. Причина: E\/N = (\S+) Тд\./;
export function dataText(s) {
  if (!s) return '';
  const str = String(s);
  if (LANG === 'ru') return str;              // русские тексты данных не трогаем
  const k = DATA_KEY[str.trim()];
  if (k) return t(k);
  const m = CRASH_RE.exec(str);
  if (m) return t('crash.stepRejected', { step: m[1], t: m[2], dt: m[3], en: m[4] });
  return str;
}

/** Статические строки шаблона: data-i18n="key", data-i18n-attr="title:key;aria-label:key". */
export function applyStatic(root, titleKey) {
  const scope = root || (typeof document !== 'undefined' ? document : null);
  if (!scope) return;
  scope.querySelectorAll('[data-i18n]').forEach((n) => { n.textContent = t(n.dataset.i18n); });
  scope.querySelectorAll('[data-i18n-attr]').forEach((n) => {
    for (const pair of n.dataset.i18nAttr.split(';')) {
      const i = pair.indexOf(':');
      if (i < 0) continue;
      n.setAttribute(pair.slice(0, i).trim(), t(pair.slice(i + 1).trim()));
    }
  });
  const tk = titleKey || 'doc.title';
  if (typeof document !== 'undefined') {
    document.title = t(tk);
    applyMeta(tk, titleKey ? null : 'doc.desc');
  }
}

/** SEO-мета под активный язык: description, og:*, twitter:* и og:locale. */
function applyMeta(titleKey, descKey) {
  const set = (sel, val) => {
    const el = document.querySelector(sel);
    if (el && val) el.setAttribute('content', val);
  };
  const title = t(titleKey || 'doc.title');
  const desc = descKey ? t(descKey) : null;
  set('meta[property="og:title"]', title);
  set('meta[name="twitter:title"]', title);
  if (desc) {
    set('meta[name="description"]', desc);
    set('meta[property="og:description"]', desc);
    set('meta[name="twitter:description"]', desc);
  }
  set('meta[property="og:locale"]', LANG === 'ru' ? 'ru_RU' : 'en_US');
}

/** Переключатель RU/EN: href той же страницы с ?lang=xx, активный помечен. */
export function initLangSwitch(root) {
  const scope = root || (typeof document !== 'undefined' ? document : null);
  if (!scope) return;
  scope.querySelectorAll('[data-lang-link]').forEach((a) => {
    const lang = a.getAttribute('data-lang-link');
    const u = new URL(window.location.href);
    u.searchParams.set('lang', lang);
    a.setAttribute('href', u.pathname + u.search + u.hash);
    if (lang === LANG) a.setAttribute('aria-current', 'true');
    else a.removeAttribute('aria-current');
  });
  const nav = scope.querySelector('.lang');
  if (nav) nav.setAttribute('aria-label', t('lang.aria'));
}

export default { I18N, LANG, L, t, plural, fieldLabel, fieldUnit, dataText, applyStatic, initLangSwitch };
