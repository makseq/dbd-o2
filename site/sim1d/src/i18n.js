// i18n.js — словарь ru/en для 1D-симулятора ДБР. Единственный файл проекта,
// в котором разрешена кириллица (объект I18N.ru). Подключается ПЕРВЫМ:
// ui.js импортирует его до всего остального, window.LANG выставляется инлайном
// в <head> (см. index.html).
//
// Значение ключа — либо строка, либо функция (шаблон) от объекта параметров.
// Символы физических величин (U_app, C_diel, E/N, ε_r, ρ) НЕ переводятся:
// это формульная нотация, одинаковая в обоих языках. Переводятся слова и
// единицы измерения (kV → кВ, ms → мс, W → Вт …).

export const LANG = (window.LANG === 'en' || window.LANG === 'ru') ? window.LANG : 'en';

/** en: plural(n, 'cycle', 'cycles') */
const plEn = (n, one, many) => (Math.abs(n) === 1 ? one : many);
/** ru: три формы — 1 пробой, 2 пробоя, 5 пробоев */
const plRu = (n, one, few, many) => {
  const a = Math.abs(n) % 100, b = a % 10;
  if (a > 10 && a < 20) return many;
  if (b > 1 && b < 5) return few;
  if (b === 1) return one;
  return many;
};

const I18N = {
  // ══════════════════════════════════════════════════════════════ РУССКИЙ ══
  ru: {
    'doc.title': '1D-симулятор дрейф-диффузии | ДБР в кислороде',
    'doc.desc': 'Барьерный разряд в кислороде вживую в браузере: дрейф-диффузионный перенос, уравнение Пуассона с поверхностным зарядом, девять компонент плазмы.',

    // ── единицы измерения ──
    'unit.kV': 'кВ', 'unit.V': 'В', 'unit.kHz': 'кГц', 'unit.mm': 'мм',
    'unit.um': 'мкм', 'unit.ms': 'мс', 'unit.us': 'мкс', 'unit.ns': 'нс',
    'unit.s': 'с', 'unit.W': 'Вт', 'unit.pF': 'пФ', 'unit.mJ': 'мДж',
    'unit.JL': 'Дж/л', 'unit.ppm': 'ppm', 'unit.Td': 'Тд', 'unit.mA': 'мА',
    'unit.A': 'А', 'unit.ohm': 'Ом', 'unit.K': 'К', 'unit.Torr': 'Торр',
    'unit.cm2': 'см²', 'unit.m3': 'м⁻³', 'unit.m3s': 'м⁻³с⁻¹',
    'unit.Cm3': 'Кл·м⁻³', 'unit.Vm': 'В/м', 'unit.C': 'Кл', 'unit.nC': 'нКл',
    'unit.nCcm2': 'нКл/см²', 'unit.ms_': 'м/с', 'unit.pct': '%',

    // ── шапка ──
    'lang.aria': 'Язык / Language',
    'hdr.home': '← Барьерный разряд в кислороде',
    'hdr.preset': 'пресет',
    'hdr.presetAria': 'Пресет параметров',
    'preset.filamentary': 'Филаментарный',
    'preset.townsend': 'Таунсендовский / однородный',
    'preset.ozonizer': 'Озонатор',
    'hdr.summaryInit': 'U₀ 10,0 кВ · f 10,0 кГц · зазор 1,00 мм',
    'hdr.summary': (p) => `U₀ ${p.u} кВ · f ${p.f} кГц · зазор ${p.gap} мм · ε_r ${p.eps}`,
    'run.running': 'СЧЁТ',
    'run.paused': 'ПАУЗА',
    'hdr.fpsInit': '– fps',
    'hdr.fps': (p) => `${p.n} fps`,
    'hdr.speedInit': '⟳ 1×',
    'hdr.helpTitle': 'Горячие клавиши (?)',
    'hdr.helpAria': 'Горячие клавиши',
    'hdr.drawerTitle': 'Управление',
    'hdr.drawerAria': 'Открыть панель управления',

    // ── B. разрядный промежуток ──
    'panel.gap': 'Разрядный промежуток',
    'gap.fieldAria': 'Отображаемое поле',
    'gap.log': 'лог',
    'gap.auto': 'авто',
    'gap.modeAria': 'Режим показа промежутка',
    'gap.snapshot': 'снимок',
    'gap.streak': 'развёртка x–t',
    'gap.canvasAria': 'Карта поля в разрядном промежутке',
    'gap.glowBadge': 'послесвечение вкл.',
    'gap.foot': '1D-модель · однородность в поперечном направлении',
    'gap.metal': 'металл',
    'gap.dielectric': 'диэлектрик',
    'gap.gasGap': 'газовый зазор',
    'gap.axisX': 'x, мм',
    'gap.now': 'сейчас',
    'gap.streakInfo': (p) => `↓ время · 1 строка = ${p.row} · v_фронта ≈ ${p.v}`,
    'gap.frontSpeed': (p) => `${p.v} м/с`,
    'gap.frontNA': 'н/д (замедлите воспроизведение)',
    'gap.sigma': (p) => `${p.v} нКл/см²`,
    'gap.colorbar': (p) => `${p.name}, ${p.unit}`,
    'gap.readout': (p) => `${p.name}(x=${p.x} мм) = ${p.v} ${p.unit}`,

    // ── поля цветовой карты ──
    'field.ne': 'n_e', 'field.ne.unit': 'м⁻³',
    'field.rho': 'ρ', 'field.rho.unit': 'Кл·м⁻³',
    'field.sion': 'S_ion (послесвечение)', 'field.sion.unit': 'м⁻³с⁻¹',
    'field.Eabs': '|E|', 'field.Eabs.unit': 'В/м',
    'field.EN': 'E/N', 'field.EN.unit': 'Тд',

    // ── C. метрики ──
    'panel.metrics': 'Метрики',
    'metrics.periodInit': 'период #0 · t = 0,000 мс',
    'metrics.period': (p) => `период #${p.i} · t = ${p.t}`,
    'metrics.periodAria': (p) => `период ${p.i}: ${p.bd} ${plRu(p.bd, 'пробой', 'пробоя', 'пробоев')}, `
      + `пиковый ток ${p.ipk} мА, мощность ${p.pw} Вт`,
    'metric.power': 'Мощность',
    'metric.spec': 'Уд. энергия',
    'metric.o3': 'Озон',
    'metric.en': 'макс. E/N',
    'metric.ipk': 'Пиковый ток',
    'metric.bd': 'Пробоев за период',
    'metric.cdiel': 'C_diel',
    'metric.ccell': 'C_cell',
    'metric.eper': 'Энергия за период',
    'solver.statusInit': 'dt – · шагов/кадр – · стенное→сим –',
    'solver.status': (p) => `dt ${p.dt} с · шагов/кадр ${p.steps} · стенное ${p.wall} мс → сим ${p.sim}`
      + ` · отбраковка ${p.rej}%${p.resid}${p.clip}`,
    'solver.resid': (p) => ` · невязка ${p.v}`,
    'solver.clip': (p) => ` · срез ${p.v} Кл`,

    // ── D. осциллограммы ──
    'panel.wave': 'Осциллограммы',
    'wave.scaleAria': 'Шкала тока',
    'wave.lin': 'лин',
    'wave.log': 'лог',
    'wave.canvasAria': 'Осциллограммы напряжения и тока',
    'wave.cursorInit': 't = –',
    'wave.cursor': (p) => `t = ${p.t}`,
    'wave.trackU': 'кВ',
    'wave.trackILog': 'А (лог)',
    'wave.trackI': 'мА',
    'wave.iBoth': 'I_disch (жирный) · I_total',
    'wave.iOne': 'I_disch',
    'wave.noSamplesWindow': (p) => `в этом окне нет отсчётов — история охватывает последние ${p.span}`,
    'wave.noSamplesYet': 'отсчётов пока нет — нажмите ⏵ пуск',
    'wave.decimated': (p) => `прореживание ×${p.n} (огибающая min–max)`,

    // ── E. фигура Лиссажу ──
    'panel.qv': 'Фигура Лиссажу Q–V',
    'panel.qv.sub': '(Мэнли)',
    'qv.windowAria': 'Окно фигуры Лиссажу',
    'qv.win1': '1 период',
    'qv.win8': 'последние 8',
    'qv.canvasAria': 'Фигура Лиссажу: заряд против напряжения',
    'qv.canvasWarnInit': '⚠ неидеальная фигура Лиссажу',
    'qv.footInit': 'C_diel – · C_cell – · E – · P –',
    'qv.axisU': 'U_app, кВ',
    'qv.axisQ': 'Q, нКл',
    'qv.foot': (p) => `C<sub>diel</sub> ${p.cd}${p.nb}пФ <span class="muted">(геом. ${p.cdg})</span> · `
      + `C<sub>cell</sub> ${p.cc}${p.nb}пФ · E ${p.e}${p.nb}мДж · P ${p.pw}${p.nb}Вт · U_burn ${p.ub}${p.nb}кВ${p.ghosts}`,
    'qv.ghosts': (p) => ` · призраков ${p.n}`,
    'qv.warnSides': (p) => `⚠ стороны не выделены (${p.on}/${p.off}) · показана геометрическая C`,
    'qv.warnShape': (p) => `⚠ это не параллелограмм · R² ${p.r2on}/${p.r2off} · дуги ${p.on}/${p.off}`,
    'qv.warnTitle': 'Наклоны Мэнли берутся из МНК-фитов горящей и тёмной дуг. '
      + 'R² < 0,98 означает, что сторона не прямая, то есть фигура не параллелограмм, '
      + 'а C_diel — наклон фита, а не наклон стороны параллелограмма.',

    // ── F. профили плотностей ──
    'panel.prof': 'Профили плотностей',
    'panel.prof.sub': 'log₁₀ м⁻³',
    'prof.modeAria': 'Режим профилей',
    'prof.x': 'профиль по x',
    'prof.t': 'временной ряд',
    'prof.cvdSafe': 'CVD-палитра',
    'prof.canvasAria': 'Профили плотностей компонентов',
    'prof.cursorInit': 'x = –',
    'prof.cursor': (p) => `x = ${p.x} мм`,
    'prof.collecting': 'идёт накопление данных…',
    'prof.axisX': 'x, мм',
    'prof.axisT': 't',
    'prof.tip': (p) => `x = ${p.x} мм\n${p.rows}\nΣq/e = ${p.q}`,
    'legend.btnTitle': 'клик = соло · alt+клик = скрыть',
    'legend.hint': 'соло: клик · скрыть: alt-клик',

    // ── G. управление ──
    'panel.ctl': 'Управление',
    'ctl.resetSoft': '⟲ сброс',
    'ctl.resetSoftTitle': 'Сбросить состояние, параметры оставить',
    'ctl.resetPreset': '⟲⟲ пресет',
    'ctl.resetPresetTitle': 'Вернуться к пресету',
    'ctl.gridInit': 'Δx – · λ_D – ',
    'grid.info': (p) => `${p.mark}Δx = ${p.dx}${p.nb}мкм · λ_D(сейчас) = ${p.ld}${p.nb}мкм · `
      + `λ_D/Δx сейчас ${p.now} / худшее ${p.worst}`,
    'ctl.U0kV': 'U₀ амплитуда',
    'ctl.freqKHz': 'f частота',
    'ctl.gapMM': 'd_зазор',
    'ctl.dielMM1': 'd_диэл ×2',
    'ctl.epsR': 'ε_r',
    'ctl.gamma': 'γ втор. эмиссия',
    'ctl.seedDensity': 'n₀ затравка',
    'ctl.pressureTorr': 'p давление',
    'ctl.tempK': 'T газа',
    'ctl.areaCM2': 'S площадь',
    'ctl.ballastOhm': 'R_посл',
    'ctl.nCells': 'N_x ячеек',
    'ctl.valueAria': (p) => `${p.label}, значение${p.unit}`,
    'ctl.rangeTitle': (p) => `${p.label}: ${p.min}…${p.max}${p.unit}`,
    'ctl.secAccuracy': 'Точность / скорость',
    'ctl.mode': 'режим',
    'ctl.modeAria': 'Пресет точности солвера',
    'ctl.mode.demo': 'demo — самый быстрый, CFL 0,8',
    'ctl.mode.default': 'default — CFL 0,4',
    'ctl.mode.accurate': 'accurate — CFL 0,2, медленно',

    // ── H. транспорт и таймлайны ──
    'tr.resetTitle': 'Вернуться к t = 0 (R)',
    'tr.play': '⏵ пуск',
    'tr.pause': '⏸ пауза',
    'tr.playTitle': 'Пуск / пауза (пробел)',
    'tr.step': 'шаг dt',
    'tr.stepTitle': 'Один шаг солвера',
    'tr.period': '1 период',
    'tr.periodTitle': 'Просчитать ровно один период',
    'tr.steady': 'до стационара',
    'tr.steadyTitle': 'Считать до выхода на стационар',
    'tr.prevEv': '⏮ событие',
    'tr.prevEvTitle': 'Предыдущий пробой (Shift+←)',
    'tr.nextEv': 'событие ⏭',
    'tr.nextEvTitle': 'Следующий пробой (Shift+→)',
    'tr.speed': 'скорость',
    'tr.speedInit': '–',
    'tr.speedBadge': (p) => `⟳ ${p.v} с/с`,
    'tr.speedLabel': (p) => `${p.v} с/с  (1 с → ${p.t})`,
    'tr.macro': 'МАКРО',
    'tr.micro': 'МИКРО',
    'tr.macroAria': 'Общий таймлайн с событиями пробоя',
    'tr.microAria': 'Лупа детального таймлайна',
    'tr.lens': 'лупа',
    'tr.lens1T': '1 T',
    'tr.lens1TTitle': 'один полный период',
    'tr.lens8T': '8 T',
    'tr.lens8TTitle': 'восемь периодов',
    'tr.lensFil': 'фил',
    'tr.lensFilTitle': 'масштаб филамента, 50 нс',
    'tr.lensLabel': (p) => `лупа ${p.w} (${p.n} T)`,
    'macro.noHistory': 'нет истории',
    'macro.total': (p) => `${p.t} (${p.n} T)`,

    // ── сообщения (live-region) ──
    'msg.geometryChanged': 'геометрия изменена — симуляция перезапущена с t = 0',
    'msg.breakdownAt': (p) => `пробой в момент ${p.t}`,
    'msg.steady': (p) => `достигнут стационар (Δ ${p.d} %)`,
    'msg.noMoreEvents': 'в буфере больше нет событий пробоя',

    // ── модальные окна: таблицы, отчёт, справка ──
    'modal.titleInit': 'Таблица',
    'modal.closeAria': 'Закрыть',
    'tool.table': 'Табличный вид',
    'table.profTitle': (p) => `Профили при t = ${p.t}`,
    'table.th.xmm': 'x, мм',
    'table.th.E': 'E, В/м',
    'table.th.EN': 'E/N, Тд',
    'table.th.rho': 'ρ',
    'table.waveTitle': 'Отсчёты осциллограмм',
    'table.th.t': 't, с',
    'table.th.Uapp': 'U_app, В',
    'table.th.Ugap': 'U_gap, В',
    'table.th.I': 'I, А',
    'table.th.Q': 'Q, Кл',
    'table.manleyTitle': 'Отчёт Мэнли',
    'table.manleyBody': (p) => `# dbd-o2, отчёт по фигуре Лиссажу Q–V (Мэнли)
C_diel     = ${p.cd} пФ
C_cell     = ${p.cc} пФ
U_burn     = ${p.ub} кВ
E/период   = ${p.e} мДж
P          = ${p.pw} Вт
параметры  = ${p.params}`,
    'help.title': 'Горячие клавиши',
    'help.body': `Пробел     пуск / пауза
←  →       шаг кадра / перемотка
Shift+← →  предыдущий / следующий пробой
+  −       скорость
1 … 5      поле цветовой карты (n_e, ρ, S_ion, |E|, E/N)
L          лин / лог по току
T          таблица профилей
E          выгрузить CSV
?          эта справка`,

    // ── экспорт ──
    'export.mock': '(ФИКТИВНЫЙ СОЛВЕР — не физика)',
    'export.header': (p) => `# ВНИМАНИЕ: выгружается только кольцевой буфер истории — последние `
      + `${p.span} мкс (~${p.periods} периодов) из ${p.total} мкс прогона\n`,
    'export.params': '# параметры',
    'png.footer': (p) => `dbd-o2${p.mock} · U₀ ${p.u} кВ · f ${p.f} кГц · зазор ${p.gap} мм · `
      + `ε_r ${p.eps} · γ ${p.gamma} · t = ${p.t} · ${p.iso}`,
  },

  // ══════════════════════════════════════════════════════════════ ENGLISH ══
  en: {
    'doc.title': '1D drift-diffusion simulator, nine species | DBD in oxygen',
    'doc.desc': 'A barrier discharge in oxygen simulated live in the browser: drift-diffusion transport, Poisson with surface charge and nine plasma species.',

    'unit.kV': 'kV', 'unit.V': 'V', 'unit.kHz': 'kHz', 'unit.mm': 'mm',
    'unit.um': 'µm', 'unit.ms': 'ms', 'unit.us': 'µs', 'unit.ns': 'ns',
    'unit.s': 's', 'unit.W': 'W', 'unit.pF': 'pF', 'unit.mJ': 'mJ',
    'unit.JL': 'J/L', 'unit.ppm': 'ppm', 'unit.Td': 'Td', 'unit.mA': 'mA',
    'unit.A': 'A', 'unit.ohm': 'Ω', 'unit.K': 'K', 'unit.Torr': 'Torr',
    'unit.cm2': 'cm²', 'unit.m3': 'm⁻³', 'unit.m3s': 'm⁻³s⁻¹',
    'unit.Cm3': 'C·m⁻³', 'unit.Vm': 'V/m', 'unit.C': 'C', 'unit.nC': 'nC',
    'unit.nCcm2': 'nC/cm²', 'unit.ms_': 'm/s', 'unit.pct': '%',

    'lang.aria': 'Language',
    'hdr.home': '← Barrier discharge in oxygen',
    'hdr.preset': 'preset',
    'hdr.presetAria': 'Parameter preset',
    'preset.filamentary': 'Filamentary',
    'preset.townsend': 'Townsend / homogeneous',
    'preset.ozonizer': 'Ozonizer',
    'hdr.summaryInit': 'U₀ 10.0 kV · f 10.0 kHz · gap 1.00 mm',
    'hdr.summary': (p) => `U₀ ${p.u} kV · f ${p.f} kHz · gap ${p.gap} mm · ε_r ${p.eps}`,
    'run.running': 'RUNNING',
    'run.paused': 'PAUSED',
    'hdr.fpsInit': '– fps',
    'hdr.fps': (p) => `${p.n} fps`,
    'hdr.speedInit': '⟳ 1×',
    'hdr.helpTitle': 'Keyboard shortcuts (?)',
    'hdr.helpAria': 'Keyboard shortcuts',
    'hdr.drawerTitle': 'Controls',
    'hdr.drawerAria': 'Open controls',

    'panel.gap': 'Discharge gap',
    'gap.fieldAria': 'Displayed field',
    'gap.log': 'log',
    'gap.auto': 'auto',
    'gap.modeAria': 'Gap view mode',
    'gap.snapshot': 'snapshot',
    'gap.streak': 'x–t streak',
    'gap.canvasAria': 'Discharge gap field map',
    'gap.glowBadge': 'afterglow ON',
    'gap.foot': '1D model · uniform in transverse direction',
    'gap.metal': 'metal',
    'gap.dielectric': 'dielectric',
    'gap.gasGap': 'gas gap',
    'gap.axisX': 'x, mm',
    'gap.now': 'now',
    'gap.streakInfo': (p) => `↓ time · 1 row = ${p.row} · v_front ≈ ${p.v}`,
    'gap.frontSpeed': (p) => `${p.v} m/s`,
    'gap.frontNA': 'n/a (slow down playback)',
    'gap.sigma': (p) => `${p.v} nC/cm²`,
    'gap.colorbar': (p) => `${p.name}, ${p.unit}`,
    'gap.readout': (p) => `${p.name}(x=${p.x} mm) = ${p.v} ${p.unit}`,

    'field.ne': 'n_e', 'field.ne.unit': 'm⁻³',
    'field.rho': 'ρ', 'field.rho.unit': 'C·m⁻³',
    'field.sion': 'S_ion (afterglow)', 'field.sion.unit': 'm⁻³s⁻¹',
    'field.Eabs': '|E|', 'field.Eabs.unit': 'V/m',
    'field.EN': 'E/N', 'field.EN.unit': 'Td',

    'panel.metrics': 'Metrics',
    'metrics.periodInit': 'period #0 · t = 0.000 ms',
    'metrics.period': (p) => `period #${p.i} · t = ${p.t}`,
    'metrics.periodAria': (p) => `period ${p.i}: ${p.bd} ${plEn(p.bd, 'breakdown', 'breakdowns')}, `
      + `peak current ${p.ipk} mA, power ${p.pw} W`,
    'metric.power': 'Power',
    'metric.spec': 'Spec. energy',
    'metric.o3': 'Ozone',
    'metric.en': 'max E/N',
    'metric.ipk': 'Peak current',
    'metric.bd': 'Breakdowns/T',
    'metric.cdiel': 'C_diel',
    'metric.ccell': 'C_cell',
    'metric.eper': 'Energy/period',
    'solver.statusInit': 'dt – · steps/frame – · wall→sim –',
    'solver.status': (p) => `dt ${p.dt} s · steps/frame ${p.steps} · wall ${p.wall} ms → sim ${p.sim}`
      + ` · rejects ${p.rej}%${p.resid}${p.clip}`,
    'solver.resid': (p) => ` · resid ${p.v}`,
    'solver.clip': (p) => ` · clip ${p.v} C`,

    'panel.wave': 'Waveforms',
    'wave.scaleAria': 'Current scale',
    'wave.lin': 'lin',
    'wave.log': 'log',
    'wave.canvasAria': 'Voltage and current waveforms',
    'wave.cursorInit': 't = –',
    'wave.cursor': (p) => `t = ${p.t}`,
    'wave.trackU': 'kV',
    'wave.trackILog': 'A (log)',
    'wave.trackI': 'mA',
    'wave.iBoth': 'I_disch (bold) · I_total',
    'wave.iOne': 'I_disch',
    'wave.noSamplesWindow': (p) => `no samples in this window — history covers the last ${p.span}`,
    'wave.noSamplesYet': 'no samples yet — press ⏵ play',
    'wave.decimated': (p) => `decimated ×${p.n} (min–max envelope)`,

    'panel.qv': 'Q–V Lissajous',
    'panel.qv.sub': '(Manley)',
    'qv.windowAria': 'Lissajous window',
    'qv.win1': '1 period',
    'qv.win8': 'last 8',
    'qv.canvasAria': 'Charge versus voltage Lissajous figure',
    'qv.canvasWarnInit': '⚠ non-ideal Lissajous',
    'qv.footInit': 'C_diel – · C_cell – · E – · P –',
    'qv.axisU': 'U_app, kV',
    'qv.axisQ': 'Q, nC',
    'qv.foot': (p) => `C<sub>diel</sub> ${p.cd}${p.nb}pF <span class="muted">(geom ${p.cdg})</span> · `
      + `C<sub>cell</sub> ${p.cc}${p.nb}pF · E ${p.e}${p.nb}mJ · P ${p.pw}${p.nb}W · U_burn ${p.ub}${p.nb}kV${p.ghosts}`,
    'qv.ghosts': (p) => ` · ghosts ${p.n}`,
    'qv.warnSides': (p) => `⚠ sides not resolved (${p.on}/${p.off}) · geometric C shown`,
    'qv.warnShape': (p) => `⚠ not a parallelogram · R² ${p.r2on}/${p.r2off} · arcs ${p.on}/${p.off}`,
    'qv.warnTitle': 'Manley slopes are taken from least-squares fits of the burning / dark arcs. '
      + 'R² < 0.98 means the side is not straight, i.e. the figure is not a parallelogram '
      + 'and C_diel is a fitted slope, not a parallelogram slope.',

    'panel.prof': 'Density profiles',
    'panel.prof.sub': 'log₁₀ m⁻³',
    'prof.modeAria': 'Profile mode',
    'prof.x': 'x-profile',
    'prof.t': 'time-series',
    'prof.cvdSafe': 'CVD-safe',
    'prof.canvasAria': 'Species density profiles',
    'prof.cursorInit': 'x = –',
    'prof.cursor': (p) => `x = ${p.x} mm`,
    'prof.collecting': 'collecting data…',
    'prof.axisX': 'x, mm',
    'prof.axisT': 't',
    'prof.tip': (p) => `x = ${p.x} mm\n${p.rows}\nΣq/e = ${p.q}`,
    'legend.btnTitle': 'click = solo · alt+click = hide',
    'legend.hint': 'solo: click · hide: alt-click',

    'panel.ctl': 'Controls',
    'ctl.resetSoft': '⟲ reset',
    'ctl.resetSoftTitle': 'Reset state, keep parameters',
    'ctl.resetPreset': '⟲⟲ preset',
    'ctl.resetPresetTitle': 'Reset to preset',
    'ctl.gridInit': 'Δx – · λ_D – ',
    'grid.info': (p) => `${p.mark}Δx = ${p.dx}${p.nb}µm · λ_D(now) = ${p.ld}${p.nb}µm · `
      + `λ_D/Δx now ${p.now} / worst ${p.worst}`,
    'ctl.U0kV': 'U₀ amplitude',
    'ctl.freqKHz': 'f frequency',
    'ctl.gapMM': 'd_gap',
    'ctl.dielMM1': 'd_diel ×2',
    'ctl.epsR': 'ε_r',
    'ctl.gamma': 'γ sec. emis.',
    'ctl.seedDensity': 'n₀ seed',
    'ctl.pressureTorr': 'p pressure',
    'ctl.tempK': 'T gas',
    'ctl.areaCM2': 'S area',
    'ctl.ballastOhm': 'R_series',
    'ctl.nCells': 'N_x cells',
    'ctl.valueAria': (p) => `${p.label} value${p.unit}`,
    'ctl.rangeTitle': (p) => `${p.label}: ${p.min}…${p.max}${p.unit}`,
    'ctl.secAccuracy': 'Accuracy / speed',
    'ctl.mode': 'mode',
    'ctl.modeAria': 'Solver accuracy preset',
    'ctl.mode.demo': 'demo — fastest, CFL 0.8',
    'ctl.mode.default': 'default — CFL 0.4',
    'ctl.mode.accurate': 'accurate — CFL 0.2, slow',

    'tr.resetTitle': 'Reset to t = 0 (R)',
    'tr.play': '⏵ play',
    'tr.pause': '⏸ pause',
    'tr.playTitle': 'Play / pause (Space)',
    'tr.step': 'step dt',
    'tr.stepTitle': 'One solver step',
    'tr.period': '1 period',
    'tr.periodTitle': 'Run exactly one period',
    'tr.steady': 'run to steady',
    'tr.steadyTitle': 'Run to steady state',
    'tr.prevEv': '⏮ event',
    'tr.prevEvTitle': 'Previous breakdown (Shift+←)',
    'tr.nextEv': 'event ⏭',
    'tr.nextEvTitle': 'Next breakdown (Shift+→)',
    'tr.speed': 'speed',
    'tr.speedInit': '–',
    'tr.speedBadge': (p) => `⟳ ${p.v} s/s`,
    'tr.speedLabel': (p) => `${p.v} s/s  (1 s → ${p.t})`,
    'tr.macro': 'MACRO',
    'tr.micro': 'MICRO',
    'tr.macroAria': 'Macro timeline with breakdown events',
    'tr.microAria': 'Micro timeline lens',
    'tr.lens': 'lens',
    'tr.lens1T': '1 T',
    'tr.lens1TTitle': 'one full period',
    'tr.lens8T': '8 T',
    'tr.lens8TTitle': 'eight periods',
    'tr.lensFil': 'fil',
    'tr.lensFilTitle': 'filament scale, 50 ns',
    'tr.lensLabel': (p) => `lens ${p.w} (${p.n} T)`,
    'macro.noHistory': 'no history',
    'macro.total': (p) => `${p.t} (${p.n} T)`,

    'msg.geometryChanged': 'geometry changed — simulation restarted from t = 0',
    'msg.breakdownAt': (p) => `breakdown at ${p.t}`,
    'msg.steady': (p) => `steady state reached (Δ ${p.d} %)`,
    'msg.noMoreEvents': 'no further breakdown event in buffer',

    'modal.titleInit': 'Table',
    'modal.closeAria': 'Close',
    'tool.table': 'Table view',
    'table.profTitle': (p) => `Profiles @ t = ${p.t}`,
    'table.th.xmm': 'x, mm',
    'table.th.E': 'E, V/m',
    'table.th.EN': 'E/N, Td',
    'table.th.rho': 'ρ',
    'table.waveTitle': 'Waveform samples',
    'table.th.t': 't, s',
    'table.th.Uapp': 'U_app, V',
    'table.th.Ugap': 'U_gap, V',
    'table.th.I': 'I, A',
    'table.th.Q': 'Q, C',
    'table.manleyTitle': 'Manley report',
    'table.manleyBody': (p) => `# dbd-o2 Q–V (Manley) report
C_diel     = ${p.cd} pF
C_cell     = ${p.cc} pF
U_burn     = ${p.ub} kV
E/period   = ${p.e} mJ
P          = ${p.pw} W
params     = ${p.params}`,
    'help.title': 'Keyboard shortcuts',
    'help.body': `Space      play / pause
←  →       step frame / seek
Shift+← →  previous / next breakdown event
+  −       speed
1 … 5      colormap field (n_e, ρ, S_ion, |E|, E/N)
L          lin / log current
T          table view of profiles
E          export CSV
?          this help`,

    'export.mock': '(MOCK SOLVER — not physical)',
    'export.header': (p) => `# NOTE: only the history ring buffer is exported — the last `
      + `${p.span} µs (~${p.periods} periods) out of ${p.total} µs of the run\n`,
    'export.params': '# params',
    'png.footer': (p) => `dbd-o2${p.mock} · U₀ ${p.u} kV · f ${p.f} kHz · gap ${p.gap} mm · `
      + `ε_r ${p.eps} · γ ${p.gamma} · t = ${p.t} · ${p.iso}`,
  },
};

export const L = I18N[LANG];

/** t('key') или t('key', {…}) для шаблонов-функций. */
export function t(key, params) {
  const v = L[key] !== undefined ? L[key] : I18N.ru[key];
  if (v === undefined) return key;
  return typeof v === 'function' ? v(params || {}) : v;
}

/** Множественное число для активного языка: n + существительное. */
export function plural(n, one, few, many) {
  return LANG === 'ru' ? plRu(n, one, few, many === undefined ? few : many) : plEn(n, one, few);
}

export const KEYS = { ru: Object.keys(I18N.ru), en: Object.keys(I18N.en) };
/** Полный словарь — нужен тесту полноты ключей (test/i18n-keys.test.mjs). */
export { I18N };

/**
 * Заполняет статическую разметку: data-i18n="key" -> textContent,
 * data-i18n-attr="aria-label:key;title:key" -> атрибуты.
 * Плюс <title> документа и подсветка активной ссылки переключателя.
 */
/** SEO-мета под активный язык: description, og:*, twitter:* и og:locale. */
function applyMeta() {
  const set = (sel, val) => {
    const el = document.querySelector(sel);
    if (el && val) el.setAttribute('content', val);
  };
  const title = t('doc.title');
  const desc = t('doc.desc');
  set('meta[name="description"]', desc);
  set('meta[property="og:title"]', title);
  set('meta[property="og:description"]', desc);
  set('meta[property="og:locale"]', LANG === 'ru' ? 'ru_RU' : 'en_US');
  set('meta[name="twitter:title"]', title);
  set('meta[name="twitter:description"]', desc);
}

export function applyStatic(root = document) {
  document.title = t('doc.title');
  applyMeta();
  for (const el of root.querySelectorAll('[data-i18n]')) el.textContent = t(el.dataset.i18n);
  for (const el of root.querySelectorAll('[data-i18n-attr]')) {
    for (const pair of el.dataset.i18nAttr.split(';')) {
      const i = pair.indexOf(':');
      if (i < 0) continue;
      el.setAttribute(pair.slice(0, i).trim(), t(pair.slice(i + 1).trim()));
    }
  }
  // ссылки переключателя: сохраняем прочие query-параметры страницы
  for (const a of root.querySelectorAll('[data-lang-link]')) {
    const lang = a.dataset.langLink;
    const q = new URLSearchParams(location.search);
    q.set('lang', lang);
    a.href = location.pathname + '?' + q.toString();
    if (lang === LANG) a.setAttribute('aria-current', 'true');
    else a.removeAttribute('aria-current');
  }
}
