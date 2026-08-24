/* =============================================================================
 * FastBN Network Viewer — cosmos.gl (window.Cosmos) を使ったネットワーク表示。
 *
 * データは serve.py の JSON API から取る:
 *   GET /api/networks           … 見つかったネットワークの一覧
 *   GET /api/network?id=<id>    … ノード・エッジ・重要度・確率・判定
 *
 * 表示の約束:
 *   * 矢印は「親 → 子」。DAG の向きはマルコフ同値類の中で決まらないことがある。
 *   * エッジの色と太さは |metric| (visualize.py と同じく絶対値) か
 *     ブートストラップ確率、または正解構造との判定。
 *   * ノードの色は BN 構造上の役割 (親のみ / 中間 / 子のみ) と注目遺伝子。
 * ========================================================================== */
(() => {
  'use strict';

  const Cosmos = window.Cosmos;
  if (!Cosmos || !Cosmos.Graph) {
    document.getElementById('loading').textContent =
      'cosmos.gl を読み込めませんでした。viewer/vendor/ のビルドを確認してください ' +
      '(python3 viewer/serve.py --fetch-vendor)。';
    return;
  }

  const el = (id) => document.getElementById(id);
  const SPACE = 4096;

  /* 失敗をユーザに見える形で出す (WebGL が無い環境など)。
     ブラウザが出す無害な警告 (ResizeObserver ループ) は無視し、
     1 つ目のネットワークを描けたあとはオーバーレイを出さず HUD に出す。 */
  const BENIGN = /ResizeObserver loop/i;
  function fail(msg) {
    if (BENIGN.test(msg)) { console.warn('[viewer] (無視)', msg); return; }
    console.error('[viewer]', msg);
    if (S.net) { el('hud').textContent = '⚠ ' + msg.split('\n')[0]; return; }
    const box = el('loading');
    box.hidden = false;
    box.textContent = msg;
    box.style.whiteSpace = 'pre-wrap';
    box.style.padding = '0 24px';
    box.style.textAlign = 'center';
  }
  window.addEventListener('error', (e) => fail('JavaScript エラー: ' + (e.message || e.error)));
  window.addEventListener('unhandledrejection', (e) =>
    fail('初期化に失敗しました: ' + (e.reason && e.reason.message ? e.reason.message : e.reason) +
         '\nWebGL2 が使えないブラウザ・環境では表示できません。'));

  const METRIC_LABEL = {
    dlogL: '|ΔlogL|', dBIC: '|ΔBIC|', dK2: '|ΔK2|', dBDeu: '|ΔBDeu|',
    mean_dlogL_per_sample: '|平均ΔlogL/サンプル|',
    std_dlogL_per_sample: '|標準偏差ΔlogL/サンプル|',
  };
  const EVAL_COLOR = {
    TP: [0.35, 0.84, 0.55], FP: [1.0, 0.36, 0.48],
    FP_reversed: [1.0, 0.66, 0.30], FN: [0.55, 0.58, 0.65],
  };
  const ROLE_COLOR = {
    src: [0.35, 0.84, 0.55], mid: [0.42, 0.66, 1.0],
    sink: [0.73, 0.55, 1.0], tgt: [1.0, 0.36, 0.48],
  };
  /* 低 → 高 の 3 色ランプ (暗いスレート → ティール → 琥珀) */
  const RAMP = [[0.24, 0.31, 0.42], [0.30, 0.79, 0.94], [1.0, 0.82, 0.40]];

  const S = {
    catalog: [], net: null, graph: null,
    positions: null, linkArr: null,
    deg: [], role: [], parents: [], children: [],
    attr: 'none', attrValues: null, attrRank: null, isNumeric: false,
    attrData: {},                       /* key -> {values, rank, min, max, n} */
    topPct: 100, selection: null, hover: null,
    filterMode: 'simple',               /* simple: 上位X% / advanced: 条件の組み合わせ */
    combine: 'and', conds: [], visibleCount: 0, histLog: false,
    labels: true, labelCount: 25, arrows: true, curved: true, targets: true,
    paused: false, labelIdx: [], pool: [], raf: 0,
    fitTimers: [], userZoomed: false, neighborSet: new Set(),
  };

  /* ---------- 色ユーティリティ ---------- */
  function ramp(t) {
    t = Math.max(0, Math.min(1, t));
    const seg = t < 0.5 ? 0 : 1;
    const u = seg === 0 ? t * 2 : (t - 0.5) * 2;
    const a = RAMP[seg], b = RAMP[seg + 1];
    return [a[0] + (b[0] - a[0]) * u, a[1] + (b[1] - a[1]) * u, a[2] + (b[2] - a[2]) * u];
  }
  const rgbaCss = (c, a) =>
    `rgba(${Math.round(c[0] * 255)},${Math.round(c[1] * 255)},${Math.round(c[2] * 255)},${a})`;

  /* 再現性のある擬似乱数 (同じネットワークなら毎回同じ初期配置) */
  function mulberry32(seed) {
    return function () {
      seed |= 0; seed = (seed + 0x6D2B79F5) | 0;
      let t = Math.imul(seed ^ (seed >>> 15), 1 | seed);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }

  /* ---------- 一覧の取得 ---------- */
  async function boot() {
    let cat;
    try {
      cat = await fetch('api/networks').then((r) => r.json());
    } catch (e) {
      el('loading').textContent = 'serve.py に接続できません: ' + e;
      return;
    }
    S.catalog = cat.networks;
    const sel = el('network');
    const groups = new Map();
    for (const n of S.catalog) {
      if (!groups.has(n.group)) groups.set(n.group, []);
      groups.get(n.group).push(n);
    }
    for (const [g, list] of groups) {
      const og = document.createElement('optgroup');
      og.label = `${g} (${list.length})`;
      for (const n of list) {
        const o = document.createElement('option');
        o.value = n.id;
        const flags = [n.has_importance ? '重要度' : null, n.has_prob ? '確率' : null,
                       n.has_eval ? '判定' : null].filter(Boolean).join('/');
        o.textContent = `${n.dir.replace(/\/out$/, '')} — ${n.kind === 'consensus' ? 'コンセンサス' : '学習'}` +
          ` · ${n.n_edges} 辺${flags ? ' · ' + flags : ''}`;
        og.appendChild(o);
      }
      sel.appendChild(og);
    }
    sel.value = cat.default;
    sel.addEventListener('change', () => load(sel.value));
    wireUI();
    await load(cat.default);
  }

  /* ---------- ネットワークの読み込み ---------- */
  async function load(id) {
    el('loading').hidden = false;
    el('loading').textContent = '読み込み中…';
    let d;
    try {
      d = await fetch('api/network?id=' + encodeURIComponent(id)).then((r) => r.json());
    } catch (e) {
      el('loading').textContent = '読み込みに失敗しました: ' + e;
      return;
    }
    if (d.error) { el('loading').textContent = d.error; return; }
    S.net = d;
    S.selection = null;
    S.hover = null;
    S.userZoomed = false;

    const n = d.nodes.length;
    S.deg = d.nodes.map((x) => x.in + x.out);
    S.role = d.nodes.map((x) => (x.t ? 'tgt' : x.in === 0 ? 'src' : x.out === 0 ? 'sink' : 'mid'));
    S.parents = Array.from({ length: n }, () => []);
    S.children = Array.from({ length: n }, () => []);
    for (let k = 0; k < d.links.length; k += 2) {
      const s = d.links[k], t = d.links[k + 1];
      S.children[s].push([t, k / 2]);
      S.parents[t].push([s, k / 2]);
    }

    /* 初期配置: 中心付近の円板にばらまく (次数が高いほど内側) */
    const rnd = mulberry32(n * 7919 + d.links.length);
    const pos = new Float32Array(n * 2);
    const maxDeg = Math.max(1, ...S.deg);
    for (let i = 0; i < n; i++) {
      const inner = 0.35 + 0.65 * (1 - S.deg[i] / maxDeg);
      const r = SPACE * 0.26 * inner * Math.sqrt(rnd());
      const a = rnd() * Math.PI * 2;
      pos[i * 2] = SPACE / 2 + r * Math.cos(a);
      pos[i * 2 + 1] = SPACE / 2 + r * Math.sin(a);
    }
    S.positions = pos;
    S.linkArr = new Float32Array(d.links);

    prepareAttrs();
    /* 前の網の条件は、同じ属性がこの網にもある分だけ引き継ぐ */
    const usable = new Set(condKeys());
    S.conds = S.conds.filter((c) => usable.has(c.key));
    fillAttrSelect();
    ensureGraph();
    S.graph.setPointPositions(S.positions);
    S.graph.setLinks(S.linkArr);
    S.graph.setPointSizes(pointSizes());
    S.graph.setPointColors(pointColors());
    applyEdgeStyles();
    S.graph.setConfigPartial({
      curvedLinks: S.curved,
      linkDefaultArrows: S.arrows,
      simulationRepulsion: parseFloat(el('repulsion').value),
      simulationLinkDistance: parseFloat(el('linkDist').value),
    });
    S.graph.render();
    S.graph.start();
    S.paused = false;
    el('btnPause').textContent = '一時停止';
    scheduleFits();

    /* 強調表示の config はネットワークを切り替えても残るので明示的に消す
       (前の網のノード番号で greyout が効いたままになる) */
    S.graph.setConfigPartial({
      focusedPointIndex: undefined,
      highlightedPointIndices: undefined,
    });

    renderInfo();
    setFilterMode(S.filterMode);
    renderSelection();
    doSearch();
    updateLabelSet();
    el('loading').hidden = true;
    console.log(`[viewer] loaded ${d.id}: ${d.n_nodes} nodes, ${d.n_edges} edges, ` +
                `attr=${S.attr}, metrics=[${d.metrics.join(',')}]`);
  }

  /* 力学レイアウトは数秒かけて広がるので、収束の途中で何度か全体表示に合わせる
     (1 回だけだと初期配置に合わせた倍率のまま網が縮んで小さく見える)。 */
  function scheduleFits() {
    S.fitTimers.forEach(clearTimeout);
    S.fitTimers = [700, 1800, 3600, 6500].map((ms) =>
      setTimeout(() => { if (!S.userZoomed) S.graph.fitView(500, 0.22); }, ms));
  }

  function ensureGraph() {
    if (S.graph) return;
    S.graph = new Cosmos.Graph(el('graph'), {
      spaceSize: SPACE,
      backgroundColor: '#0c0f14',
      pointDefaultSize: 4,
      pointOpacity: 0.95,
      pointGreyoutOpacity: 0.12,
      linkGreyoutOpacity: 0.08,
      linkWidthScale: 1,
      linkOpacity: 1,
      linkArrowsSizeScale: 0.9,
      linkDefaultArrows: true,
      curvedLinks: true,
      curvedLinkWeight: 0.55,
      scalePointsOnZoom: true,
      enableDrag: true,
      hoveredPointCursor: 'pointer',
      renderHoveredPointRing: true,
      hoveredPointRingColor: '#ffd166',
      focusedPointRingColor: '#4cc9f0',
      simulationFriction: 0.88,
      simulationGravity: 0.06,
      simulationRepulsion: 1.0,
      simulationLinkSpring: 1.0,
      simulationLinkDistance: 20,
      simulationDecay: 5000,
      transitionDuration: 0,
      rescalePositions: false,
      fitViewOnInit: false,
      onPointMouseOver: (i, p) => { S.hover = i; showTooltip(i, p); },
      onPointMouseOut: () => { S.hover = null; el('tooltip').hidden = true; },
      onPointClick: (i) => selectNode(i),
      onBackgroundClick: () => selectNode(null),
      onZoom: (e, userDriven) => { if (userDriven) S.userZoomed = true; drawLabels(); },
      onSimulationEnd: () => { if (!S.userZoomed) S.graph.fitView(500, 0.22); },
    });
    if (S.graph.ready && typeof S.graph.ready.then === 'function') {
      S.graph.ready
        .then(() => console.log('[viewer] cosmos.gl ready (WebGL2 device initialized)'))
        .catch((e) => fail('cosmos.gl の初期化に失敗しました: ' + (e && e.message ? e.message : e) +
                          '\nWebGL2 対応のブラウザで開いてください。'));
    }
    loopLabels();
  }

  /* ---------- ノードの見た目 ---------- */
  function pointSizes() {
    const n = S.net.nodes.length;
    const out = new Float32Array(n);
    for (let i = 0; i < n; i++) out[i] = 3.2 + 2.1 * Math.sqrt(S.deg[i]);
    return out;
  }

  function pointColors() {
    const n = S.net.nodes.length;
    const out = new Float32Array(n * 4);
    for (let i = 0; i < n; i++) {
      const role = S.role[i] === 'tgt' && !S.targets ? 'mid' : S.role[i];
      const c = ROLE_COLOR[role];
      out[i * 4] = c[0]; out[i * 4 + 1] = c[1]; out[i * 4 + 2] = c[2];
      out[i * 4 + 3] = 1;
    }
    return out;
  }

  /* ---------- エッジの属性 ---------- */

  /* 数値属性ごとに 値・順位 (0..1)・範囲 を用意する。
     詳細フィルタは色に使っている属性とは無関係に、どの属性でも条件に使えるので
     ネットワーク読み込み時にまとめて計算しておく。 */
  function prepareAttrs() {
    const m = S.net.n_edges;
    S.attrData = {};
    const add = (key, vals) => {
      if (!vals) return;
      const idx = [];
      for (let i = 0; i < m; i++) {
        if (vals[i] != null && Number.isFinite(vals[i])) idx.push(i);
      }
      if (!idx.length) return;
      const sorted = idx.slice().sort((a, b) => vals[a] - vals[b]);
      const rank = new Array(m).fill(null);
      sorted.forEach((i, k) => { rank[i] = sorted.length > 1 ? k / (sorted.length - 1) : 1; });
      S.attrData[key] = {
        values: vals, rank, n: idx.length,
        min: vals[sorted[0]], max: vals[sorted[sorted.length - 1]],
        /* しきい値ごとの件数を二分探索で数えるための昇順の値 */
        sortedVals: sorted.map((i) => vals[i]),
        hist: null,   /* 初回描画時に作る */
      };
    };
    for (const key of S.net.metrics) add(key, S.net.link_metrics[key]);
    if (S.net.has_prob) add('prob', S.net.link_prob);
  }

  /* ---------- ヒストグラム (しきい値を決めるための分布表示) ---------- */

  /* 昇順配列で v 以上 / 以下の件数 (二分探索) */
  function countGe(sorted, v) {
    let lo = 0, hi = sorted.length;
    while (lo < hi) { const mid = (lo + hi) >> 1; if (sorted[mid] < v) lo = mid + 1; else hi = mid; }
    return sorted.length - lo;
  }
  function countLe(sorted, v) {
    let lo = 0, hi = sorted.length;
    while (lo < hi) { const mid = (lo + hi) >> 1; if (sorted[mid] <= v) lo = mid + 1; else hi = mid; }
    return lo;
  }

  const HIST_BINS = 32;
  /* ヒストグラムの読み方 + 縦軸の目盛り切替。
     重要度の分布は強く歪む (大半が小さく、裾が長い) ので、線形目盛りだと裾のビンが
     1px になって見えない。対数にすると裾を見ながらしきい値を決められる。 */
  const histLegend = () =>
    '<div class="histlegend">' +
    '<span><i style="background:rgba(76,201,240,.85)"></i>通る</span>' +
    '<span><i style="background:rgba(147,161,181,.35)"></i>全体</span>' +
    '<span><i style="background:rgba(255,209,102,.85)"></i>残る本数 (累積)</span>' +
    '<span><i style="background:#ff5c7a"></i>しきい値 · クリックで指定</span>' +
    `<label class="hlog"><input type="checkbox" class="histlog"${S.histLog ? ' checked' : ''}>` +
    '縦軸を対数に</label>' +
    '</div>';

  function buildHist(d) {
    if (d.hist) return d.hist;
    const lo = d.min, hi = d.max;
    const span = hi - lo;
    const w = span > 0 ? span / HIST_BINS : 1;
    const counts = new Array(HIST_BINS).fill(0);
    for (const v of d.sortedVals) {
      let b = span > 0 ? Math.floor((v - lo) / w) : 0;
      if (b >= HIST_BINS) b = HIST_BINS - 1;     /* 最大値は最後のビンに入れる */
      if (b < 0) b = 0;
      counts[b]++;
    }
    d.hist = { lo, hi, w: span > 0 ? w : 0, counts, maxCount: Math.max(...counts) };
    return d.hist;
  }

  /* 条件のしきい値を「値」で返す。順位モードならその順位に対応する値に直す。
     判定は常にこの値で行うので、同じ値のエッジは必ずまとめて通る / 通らない
     (順位で切ると同値の辺が恣意的に分断され、ヒストグラムとも食い違う)。
     quantileValue は毎エッジ呼ぶには重いので条件ごとにキャッシュする。 */
  function condThresholdValue(c, d) {
    if (c.mode !== 'pct') return c.value;
    const key = `${c.op}:${c.pct}`;
    if (c._tk !== key || c._td !== d) {
      c._tk = key; c._td = d;
      c._tv = quantileValue(d, c.op === 'ge' ? 1 - c.pct / 100 : c.pct / 100);
    }
    return c._tv;
  }

  /* ヒストグラムを描く。
       bars      : 各ビンの本数 (通る側は明るく)
       cum       : 「ここにしきい値を置いたら何辺残るか」の曲線 (単調なので歪んだ分布でも読める)
       marker    : いまのしきい値
     クリック / ドラッグで値を指定できる。onPick(value) が呼ばれる。 */
  function drawHist(host, d, c, onPick) {
    const H = buildHist(d);
    const W = Math.max(160, host.clientWidth || 272);
    const h = 62, padB = 9, padT = 3;
    const plotH = h - padB - padT;
    const x = (v) => (H.w > 0 ? ((v - H.lo) / (H.hi - H.lo)) * W : W / 2);
    const invX = (px) => H.lo + (Math.max(0, Math.min(W, px)) / W) * (H.hi - H.lo);
    const thr = condThresholdValue(c, d);

    /* ビンごとの「条件を通る本数」。値は昇順・ビンも昇順なので、通る側 (op が ge なら
       上位 totalPass 本) を端から順に配れば正確に分かる。ビンの中心で判定すると、
       値が離散的なとき (ブートストラップ確率など) に通る辺を「通らない」色で塗って
       しまうため、この数え方にしている。 */
    const totalPass = c.op === 'ge' ? countGe(d.sortedVals, thr) : countLe(d.sortedVals, thr);
    const passCount = new Array(HIST_BINS).fill(0);
    let remain = totalPass;
    if (c.op === 'ge') {
      for (let b = HIST_BINS - 1; b >= 0 && remain > 0; b--) {
        passCount[b] = Math.min(H.counts[b], remain); remain -= passCount[b];
      }
    } else {
      for (let b = 0; b < HIST_BINS && remain > 0; b++) {
        passCount[b] = Math.min(H.counts[b], remain); remain -= passCount[b];
      }
    }

    /* 棒: 全体を薄い灰色で描き、そのうち条件を通る分を明るく重ねる
       (境界のビンは一部だけが通るので、色を塗り分けるより積み上げが正確) */
    const barH = (n) => (n <= 0 ? 0 : Math.max(1, (S.histLog
      ? Math.log1p(n) / Math.log1p(H.maxCount)
      : n / H.maxCount) * plotH));
    let bars = '';
    for (let b = 0; b < HIST_BINS; b++) {
      const n = H.counts[b];
      if (!n) continue;
      const x0 = (b / HIST_BINS) * W, x1 = ((b + 1) / HIST_BINS) * W;
      const bw = Math.max(1, x1 - x0 - 0.6).toFixed(1);
      const bh = barH(n);
      bars += `<rect x="${x0.toFixed(1)}" y="${(padT + plotH - bh).toFixed(1)}" ` +
        `width="${bw}" height="${bh.toFixed(1)}" fill="rgba(147,161,181,.28)"></rect>`;
      const ph = barH(passCount[b]);
      if (ph > 0) {
        bars += `<rect x="${x0.toFixed(1)}" y="${(padT + plotH - ph).toFixed(1)}" ` +
          `width="${bw}" height="${ph.toFixed(1)}" fill="rgba(76,201,240,.85)"></rect>`;
      }
    }

    /* 累積: 各ビン境界にしきい値を置いたときの残り本数 */
    const pts = [];
    for (let b = 0; b <= HIST_BINS; b++) {
      const v = H.lo + b * H.w;
      const keep = c.op === 'ge' ? countGe(d.sortedVals, v) : countLe(d.sortedVals, v);
      pts.push(`${((b / HIST_BINS) * W).toFixed(1)},` +
               `${(padT + plotH - (keep / d.n) * plotH).toFixed(1)}`);
    }
    const cum = `<polyline points="${pts.join(' ')}" fill="none" ` +
      `stroke="rgba(255,209,102,.85)" stroke-width="1.4"></polyline>`;

    const mx = x(thr);
    const marker =
      `<line x1="${mx.toFixed(1)}" y1="${padT}" x2="${mx.toFixed(1)}" y2="${padT + plotH}" ` +
      `stroke="#ff5c7a" stroke-width="1.6"></line>` +
      `<polygon points="${(mx - 3.5).toFixed(1)},${padT + plotH} ${(mx + 3.5).toFixed(1)},` +
      `${padT + plotH} ${mx.toFixed(1)},${(padT + plotH - 4.5).toFixed(1)}" fill="#ff5c7a"></polygon>`;

    const axis = `<line x1="0" y1="${padT + plotH}" x2="${W}" y2="${padT + plotH}" ` +
      `stroke="rgba(147,161,181,.35)" stroke-width="1"></line>` +
      `<text x="0" y="${h - 1}" class="ht">${esc(fmtNum(H.lo))}</text>` +
      `<text x="${W}" y="${h - 1}" class="ht" text-anchor="end">${esc(fmtNum(H.hi))}</text>`;

    host.innerHTML =
      `<svg class="hist" width="${W}" height="${h}" viewBox="0 0 ${W} ${h}">` +
      bars + cum + marker + axis + '</svg>';

    /* --- クリック / ドラッグでしきい値を指定 --- */
    const svg = host.firstChild;
    const pick = (ev) => {
      const r = svg.getBoundingClientRect();
      onPick(invX(ev.clientX - r.left));
    };
    let dragging = false;
    svg.addEventListener('pointerdown', (ev) => {
      dragging = true; svg.setPointerCapture(ev.pointerId); pick(ev); ev.preventDefault();
    });
    svg.addEventListener('pointermove', (ev) => {
      if (dragging) { pick(ev); return; }
      const r = svg.getBoundingClientRect();
      const v = invX(ev.clientX - r.left);
      const keep = c.op === 'ge' ? countGe(d.sortedVals, v) : countLe(d.sortedVals, v);
      host.title = `ここを境にすると ${keep} 辺 (${(100 * keep / d.n).toFixed(0)}%) ` +
        `· 値 ${fmtNum(v)}`;
    });
    const stop = (ev) => {
      if (!dragging) return;
      dragging = false;
      try { svg.releasePointerCapture(ev.pointerId); } catch (e) { /* noop */ }
      refreshFilter(true);        /* 離したときだけ UI を作り直す */
    };
    svg.addEventListener('pointerup', stop);
    svg.addEventListener('pointercancel', stop);
  }

  const attrLabel = (key) => (key === 'prob' ? 'ブートストラップ確率'
    : key === 'eval' ? '正解構造との判定' : (METRIC_LABEL[key] || key));

  /* 条件に使える属性 (数値属性 + 判定) */
  function condKeys() {
    const keys = Object.keys(S.attrData);
    if (S.net.has_eval) keys.push('eval');
    return keys;
  }

  function defaultCond(key) {
    if (key === 'eval') {
      return { key, on: true, statuses: ['TP'] };
    }
    const d = S.attrData[key];
    /* 既定は「上位 30%」。しきい値の絶対値は分布を見てから決める方が多いので
       まず順位で入れておく。 */
    return { key, on: true, mode: 'pct', op: 'ge', pct: 30,
             value: d ? d.min + (d.max - d.min) * 0.7 : 0 };
  }

  function condPass(c, k) {
    if (c.key === 'eval') {
      const st = S.net.link_eval ? S.net.link_eval[k] : null;
      return st != null && c.statuses.includes(st);
    }
    const d = S.attrData[c.key];
    if (!d) return false;
    const v = d.values[k];
    if (v == null) return false;                 /* 値が無いエッジは通さない */
    const thr = condThresholdValue(c, d);
    return c.op === 'ge' ? v >= thr : v <= thr;
  }

  function activeConds() { return S.conds.filter((c) => c.on !== false); }

  /* 簡易モードの「上位 X%」を条件 1 本として表す (詳細モードと同じ判定を使う) */
  function simpleCond() {
    if (!S._simple) S._simple = { key: null, mode: 'pct', op: 'ge', pct: 100 };
    S._simple.key = S.attr;
    S._simple.pct = S.topPct;
    return S._simple;
  }

  function visibleAt(k) {
    if (S.filterMode === 'simple') {
      if (!S.isNumeric || S.topPct >= 100) return true;
      return condPass(simpleCond(), k);
    }
    const cs = activeConds();
    if (!cs.length) return true;
    return S.combine === 'and' ? cs.every((c) => condPass(c, k))
                               : cs.some((c) => condPass(c, k));
  }

  /* ---------- 詳細フィルタの UI ---------- */
  function renderConds() {
    const box = el('conds');
    box.innerHTML = '';
    const keys = condKeys();
    if (!keys.length) {
      box.innerHTML = '<div class="note">この網には数値属性 (重要度・確率) が無いため、' +
        '条件を作れません。</div>';
      el('filterStat').textContent = '';
      return;
    }
    S.conds.forEach((c, ci) => {
      const d = S.attrData[c.key];
      const wrap = document.createElement('div');
      wrap.className = 'cond' + (c.on === false ? ' off' : '');

      const r1 = document.createElement('div');
      r1.className = 'r';
      r1.innerHTML =
        `<label class="tg" title="この条件を一時的に無効にする">` +
        `<input type="checkbox" class="en" ${c.on === false ? '' : 'checked'}></label>` +
        `<select class="k">${keys.map((k) =>
          `<option value="${k}" ${k === c.key ? 'selected' : ''}>${esc(attrLabel(k))}</option>`
        ).join('')}</select>` +
        '<button class="del" type="button" title="この条件を削除">✕</button>';
      wrap.appendChild(r1);

      if (c.key === 'eval') {
        const counts = {};
        for (let k = 0; k < S.net.n_edges; k++) {
          const st = S.net.link_eval[k];
          if (st != null) counts[st] = (counts[st] || 0) + 1;
        }
        const chips = document.createElement('div');
        chips.className = 'chips';
        chips.innerHTML = ['TP', 'FP', 'FP_reversed', 'FN'].map((st) =>
          `<label><input type="checkbox" class="st" value="${st}" ` +
          `${c.statuses.includes(st) ? 'checked' : ''}>` +
          `<i class="bar" style="background:${rgbaCss(EVAL_COLOR[st], 1)}"></i>` +
          `${st} <span class="num">${counts[st] || 0}</span></label>`
        ).join('');
        wrap.appendChild(chips);
        const info = document.createElement('div');
        info.className = 'cinfo';
        info.textContent = `${countPass(c)} 辺が該当`;
        wrap.appendChild(info);
      } else {
        const r2 = document.createElement('div');
        r2.className = 'r';
        const step = d ? Math.max((d.max - d.min) / 200, 1e-6) : 0.01;
        r2.innerHTML =
          `<select class="mode">
             <option value="pct" ${c.mode === 'pct' ? 'selected' : ''}>順位</option>
             <option value="value" ${c.mode === 'value' ? 'selected' : ''}>値</option>
           </select>` +
          (c.mode === 'pct'
            ? `<select class="op">
                 <option value="ge" ${c.op === 'ge' ? 'selected' : ''}>上位</option>
                 <option value="le" ${c.op === 'le' ? 'selected' : ''}>下位</option>
               </select>
               <input type="number" class="val" min="1" max="100" step="1" value="${c.pct}">
               <span class="tg">%</span>`
            : `<select class="op">
                 <option value="ge" ${c.op === 'ge' ? 'selected' : ''}>≥</option>
                 <option value="le" ${c.op === 'le' ? 'selected' : ''}>≤</option>
               </select>
               <input type="number" class="val" step="${step}" value="${fmtInput(c.value)}">`);
        wrap.appendChild(r2);

        const sl = document.createElement('input');
        sl.type = 'range';
        sl.className = 'sl';
        if (c.mode === 'pct') {
          sl.min = 1; sl.max = 100; sl.step = 1; sl.value = c.pct;
        } else if (d) {
          sl.min = d.min; sl.max = d.max; sl.step = step; sl.value = c.value;
        } else { sl.disabled = true; }
        wrap.appendChild(sl);

        const hist = document.createElement('div');
        hist.className = 'histbox';
        wrap.appendChild(hist);

        const info = document.createElement('div');
        info.className = 'cinfo';
        wrap.appendChild(info);
        condInfo(info, c, d);

        if (d) {
          /* 描画は詳細パネルが表示されてから (幅が 0 だと潰れる) */
          requestAnimationFrame(() => {
            if (!document.body.contains(hist)) return;
            drawHist(hist, d, c, (v) => setFromHist(c, d, v, wrap));
          });
        }
      }

      /* --- 行ごとの操作 --- */
      wrap.querySelector('.en').addEventListener('change', (e) => {
        c.on = e.target.checked; refreshFilter();
      });
      wrap.querySelector('.k').addEventListener('change', (e) => {
        S.conds[ci] = defaultCond(e.target.value); refreshFilter(true);
      });
      wrap.querySelector('.del').addEventListener('click', () => {
        S.conds.splice(ci, 1); refreshFilter(true);
      });
      const modeSel = wrap.querySelector('.mode');
      if (modeSel) {
        modeSel.addEventListener('change', (e) => {
          c.mode = e.target.value;
          if (c.mode === 'value' && d) {
            /* 順位 → 値: いま選ばれている順位の位置にある値を初期値にする */
            c.value = quantileValue(d, c.op === 'ge' ? 1 - c.pct / 100 : c.pct / 100);
          }
          refreshFilter(true);
        });
        wrap.querySelector('.op').addEventListener('change', (e) => {
          c.op = e.target.value; refreshFilter(true);
        });
        const num = wrap.querySelector('.val');
        const sl2 = wrap.querySelector('.sl');
        const infoNode = wrap.querySelector('.cinfo');
        const histNode = wrap.querySelector('.histbox');
        const setVal = (v, redraw) => {
          if (c.mode === 'pct') c.pct = Math.max(1, Math.min(100, Math.round(v)));
          else c.value = v;
          if (infoNode) condInfo(infoNode, c, d);
          if (histNode && d) drawHist(histNode, d, c, (nv) => setFromHist(c, d, nv, wrap));
          refreshFilter(redraw);
        };
        num.addEventListener('input', () => {
          const v = parseFloat(num.value);
          if (Number.isFinite(v)) { if (sl2) sl2.value = v; setVal(v, false); }
        });
        num.addEventListener('change', () => refreshFilter(true));
        if (sl2) {
          sl2.addEventListener('input', () => {
            const v = parseFloat(sl2.value);
            num.value = c.mode === 'pct' ? Math.round(v) : fmtInput(v);
            setVal(v, false);
          });
          sl2.addEventListener('change', () => refreshFilter(true));
        }
      }
      wrap.querySelectorAll('.st').forEach((cb) => cb.addEventListener('change', () => {
        c.statuses = [...wrap.querySelectorAll('.st')].filter((x) => x.checked)
          .map((x) => x.value);
        refreshFilter(true);
      }));

      box.appendChild(wrap);
    });

    const cs = activeConds();
    const hasHist = S.conds.some((c) => c.key !== 'eval' && S.attrData[c.key]);
    el('filterStat').innerHTML = (cs.length
      ? `結果: <b>${S.visibleCount}</b> / ${S.net.n_edges} 辺 ` +
        `(${cs.length} 条件を ${S.combine === 'and' ? 'すべて満たす' : 'いずれか満たす'})`
      : '条件が無いので全エッジを表示しています') + (hasHist ? histLegend() : '');
  }

  /* ヒストグラムで選ばれた値を条件に反映する (順位モードなら % に直す)。
     ドラッグ中は行を作り直さず、数値入力・スライダ・棒の色だけ更新する。 */
  function setFromHist(c, d, v, wrap) {
    if (c.mode === 'pct') {
      const keep = c.op === 'ge' ? countGe(d.sortedVals, v) : countLe(d.sortedVals, v);
      c.pct = Math.max(1, Math.min(100, Math.round((100 * keep) / d.n)));
    } else {
      c.value = v;
    }
    if (wrap) {
      const num = wrap.querySelector('.val');
      const sl = wrap.querySelector('.sl');
      const shown = c.mode === 'pct' ? c.pct : c.value;
      if (num) num.value = c.mode === 'pct' ? c.pct : fmtInput(c.value);
      if (sl) sl.value = shown;
      const info = wrap.querySelector('.cinfo');
      if (info) condInfo(info, c, d);
      const hist = wrap.querySelector('.histbox');
      if (hist) drawHist(hist, d, c, (nv) => setFromHist(c, d, nv, wrap));
    }
    refreshFilter(false);
  }

  /* 条件 1 本あたりの該当本数 (スライダ操作中もここだけ更新する) */
  function condInfo(node, c, d) {
    node.textContent = d
      ? `${countPass(c)} 辺が該当 · 範囲 ${fmtNum(d.min)} 〜 ${fmtNum(d.max)}` +
        (d.n < S.net.n_edges ? ` · 値なし ${S.net.n_edges - d.n} 辺` : '')
      : 'この属性はこの網にありません';
  }

  function countPass(c) {
    let n = 0;
    for (let k = 0; k < S.net.n_edges; k++) if (condPass(c, k)) n++;
    return n;
  }

  /* 順位 t (0..1) の位置にある実際の値 */
  function quantileValue(d, t) {
    const i = Math.round(Math.max(0, Math.min(1, t)) * (d.sortedVals.length - 1));
    return d.sortedVals[i];
  }

  /* フィルタを反映する。redraw=true なら条件 UI も作り直す
     (スライダ操作中に作り直すとつまみを失うので false で呼ぶ)。 */
  function refreshFilter(redraw) {
    applyEdgeStyles();
    hud();
    if (S.filterMode === 'simple') drawSimpleHist();
    else if (redraw) renderConds();
    else updateFilterStat();
  }

  function updateFilterStat() {
    const cs = activeConds();
    if (!cs.length) return;
    const hasHist = S.conds.some((c) => c.key !== 'eval' && S.attrData[c.key]);
    el('filterStat').innerHTML =
      `結果: <b>${S.visibleCount}</b> / ${S.net.n_edges} 辺 ` +
      `(${cs.length} 条件を ${S.combine === 'and' ? 'すべて満たす' : 'いずれか満たす'})` +
      (hasHist ? histLegend() : '');
  }

  /* 簡易モード: 色に使っている属性の分布と「上位 X%」の位置を出す。
     ここでも棒をクリックすれば X% を決められる。 */
  function drawSimpleHist() {
    const host = el('simpleHist');
    const d = S.attrData[S.attr];
    if (!d || S.filterMode !== 'simple') {
      host.innerHTML = '';
      el('simpleHistInfo').textContent = S.filterMode !== 'simple' ? ''
        : S.attr === 'eval' ? '判定には分布がありません (詳細モードで状態を選べます)'
        : S.attr === 'none' ? '属性を選ぶと分布が出ます'
        : '';
      return;
    }
    const c = simpleCond();
    drawHist(host, d, c, (v) => {
      const keep = countGe(d.sortedVals, v);
      S.topPct = Math.max(5, Math.min(100, Math.round((100 * keep) / d.n / 5) * 5));
      el('topPct').value = S.topPct;
      el('topPctVal').textContent = S.topPct + '%';
      refreshFilter(false);      /* simple モードではこの中で描き直される */
    });
    const thr = condThresholdValue(c, d);
    const keep = S.topPct >= 100 ? d.n : countGe(d.sortedVals, thr);
    el('simpleHistInfo').innerHTML =
      `しきい値 <b>${esc(fmtNum(thr))}</b> 付近 · 該当 <b>${keep}</b> / ` +
      `${d.n} 辺 · 範囲 ${esc(fmtNum(d.min))} 〜 ${esc(fmtNum(d.max))}` + histLegend();
  }

  function setFilterMode(mode) {
    S.filterMode = mode;
    el('fmSimple').classList.toggle('on', mode === 'simple');
    el('fmAdv').classList.toggle('on', mode === 'advanced');
    el('filterSimple').hidden = mode !== 'simple';
    el('filterAdv').hidden = mode === 'simple';
    if (mode === 'advanced' && !S.conds.length) {
      const first = condKeys()[0];
      if (first) S.conds = [defaultCond(S.attrData[S.attr] ? S.attr : first)];
    }
    refreshFilter(true);
  }

  function fillAttrSelect() {
    const sel = el('edgeAttr');
    sel.innerHTML = '';
    const opts = [];
    for (const m of S.net.metrics) opts.push([m, METRIC_LABEL[m] || m]);
    if (S.net.has_prob) opts.push(['prob', 'ブートストラップ確率']);
    if (S.net.has_eval) opts.push(['eval', '正解構造との判定']);
    opts.push(['none', '一律 (属性なし)']);
    for (const [v, t] of opts) {
      const o = document.createElement('option');
      o.value = v; o.textContent = t; sel.appendChild(o);
    }
    S.attr = opts[0][0];
    sel.value = S.attr;
    /* 数値でない属性 (判定 / 一律) では上位 X% の絞り込みができない */
    el('topPct').disabled = (S.attr === 'eval' || S.attr === 'none');
  }

  /* 色・太さに使う値 (0..1)。重要度は桁が違うので順位、確率は値そのまま。 */
  function computeAttr() {
    if (S.attr === 'eval') {
      S.attrValues = S.net.link_eval; S.isNumeric = false; S.attrRank = null; return;
    }
    if (S.attr === 'none') {
      S.attrValues = null; S.isNumeric = false; S.attrRank = null; return;
    }
    const d = S.attrData[S.attr];
    S.isNumeric = !!d;
    S.attrValues = d ? d.values : null;
    S.attrRank = !d ? null
      : (S.attr === 'prob'
        ? d.values.map((v) => (v == null ? null : Math.max(0, Math.min(1, v))))
        : d.rank);
  }

  function applyEdgeStyles() {
    computeAttr();
    const m = S.net.n_edges;
    const colors = new Float32Array(m * 4);
    const widths = new Float32Array(m);
    const arrows = new Array(m);
    let shown = 0;
    const selSet = S.selection == null ? null :
      new Set([...S.parents[S.selection].map((x) => x[1]),
               ...S.children[S.selection].map((x) => x[1])]);

    for (let k = 0; k < m; k++) {
      let c, a, w;
      if (S.attr === 'eval') {
        c = EVAL_COLOR[S.attrValues && S.attrValues[k]] || [0.5, 0.55, 0.62];
        a = 0.75; w = 1.5;
      } else if (S.isNumeric) {
        const t = S.attrRank[k];
        if (t == null) { c = [0.42, 0.46, 0.53]; a = 0.25; w = 0.8; }
        else { c = ramp(t); a = 0.22 + 0.68 * t; w = 0.6 + 3.6 * t; }
      } else {
        c = [0.45, 0.66, 0.90]; a = 0.65; w = 1.6;
      }
      /* 簡易モード: 色の属性の上位 X% / 詳細モード: 条件の組み合わせ */
      const visible = visibleAt(k);
      if (visible) shown++;
      /* 選択中はその周辺だけ強調 */
      if (selSet) {
        if (selSet.has(k)) { c = [1.0, 0.82, 0.40]; a = 0.95; w = Math.max(w, 1.8); }
        else a *= 0.10;
      }
      if (!visible) a = 0;
      colors[k * 4] = c[0]; colors[k * 4 + 1] = c[1]; colors[k * 4 + 2] = c[2];
      colors[k * 4 + 3] = a;
      widths[k] = visible ? w : 0;
      arrows[k] = S.arrows && visible;
    }
    S.visibleCount = shown;
    S.graph.setLinkColors(colors);
    S.graph.setLinkWidths(widths);
    S.graph.setLinkArrows(arrows);
    /* setter は「保留中の変更」を積むだけなので、render で GPU に反映させる
       (第 2 引数 0 = アニメーションなしで即反映、シミュレーションには触らない) */
    S.graph.render(undefined, 0);
    renderEdgeLegend();
  }

  function renderEdgeLegend() {
    const box = el('edgeLegend');
    if (S.attr === 'eval') {
      box.innerHTML = ['TP', 'FP', 'FP_reversed', 'FN'].map((k) =>
        `<span><i class="bar" style="background:${rgbaCss(EVAL_COLOR[k], 1)}"></i>${k}</span>`).join('');
      return;
    }
    if (!S.isNumeric) { box.innerHTML = '<span>属性なし (一律の色)</span>'; return; }
    const vals = S.attrValues.filter((v) => v != null);
    const lo = Math.min(...vals), hi = Math.max(...vals);
    const stops = [0, .25, .5, .75, 1].map((t) => `${rgbaCss(ramp(t), 1)} ${t * 100}%`).join(',');
    const fmt = (v) => (S.attr === 'prob' ? v.toFixed(2) :
      Math.abs(v) >= 100 ? v.toFixed(0) : v.toPrecision(3));
    box.innerHTML =
      `<div class="ramp" style="background:linear-gradient(90deg,${stops})"></div>` +
      `<span>低 ${fmt(lo)}</span><span>高 ${fmt(hi)}</span>` +
      (S.attr === 'prob' ? '' : '<span class="note">色・太さは順位で正規化</span>');
  }

  /* ---------- パネル ---------- */
  function renderInfo() {
    const d = S.net;
    const rows = [
      ['ディレクトリ', `<b>${esc(d.dir)}</b>`],
      ['種類', esc(d.kind_label)],
      ['ノード', `<b>${d.n_nodes}</b>${d.n_vars ? ` / 変数 ${d.n_vars}` : ''}`],
      ['エッジ', `<b>${d.n_edges}</b>`],
      ['注目遺伝子', d.n_targets ? `<b>${d.n_targets}</b>` : 'なし'],
    ];
    el('netinfo').innerHTML = rows.map(([k, v]) =>
      `<div><span class="k">${k}</span>${v}</div>`).join('');
    const w = el('warnings');
    if (d.warnings && d.warnings.length) {
      w.hidden = false;
      w.innerHTML = d.warnings.map(esc).join('<br>');
    } else w.hidden = true;

    const f = d.files || {};
    const items = [['エッジ', f.edges], ['名前', f.named], ['重要度', f.importance],
                   ['確率', f.prob], ['判定', f.eval], ['変数表', f.var_map],
                   ['注目遺伝子', (f.targets || []).join(', ')]];
    el('files').innerHTML = '<div class="lbl">読み込んだファイル</div>' +
      items.filter(([, v]) => v && v.length)
        .map(([k, v]) => `<div>${k}: <code>${esc(String(v))}</code></div>`).join('');
    hud();
  }

  function hud() {
    const d = S.net;
    const shown = S.visibleCount;
    el('hud').textContent =
      `${d.n_nodes} ノード / ${shown}${shown !== d.n_edges ? ` (全 ${d.n_edges})` : ''} エッジ` +
      (S.selection != null ? ` · 選択: ${d.nodes[S.selection].n}` : '');
  }

  function selectNode(i) {
    S.selection = i;
    S.neighborSet = i == null ? new Set() :
      new Set([i, ...S.parents[i].map((x) => x[0]), ...S.children[i].map((x) => x[0])]);
    S.graph.setConfigPartial({
      focusedPointIndex: i == null ? undefined : i,
      highlightedPointIndices: i == null ? undefined :
        [i, ...S.parents[i].map((x) => x[0]), ...S.children[i].map((x) => x[0])],
    });
    applyEdgeStyles();
    renderSelection();
    updateLabelSet();
    hud();
  }

  function renderSelection() {
    const box = el('selectionBox');
    if (S.selection == null) { box.hidden = true; return; }
    box.hidden = false;
    const i = S.selection, nd = S.net.nodes[i];
    const list = (arr, dir) => {
      if (!arr.length) return '<ul><li class="empty">なし</li></ul>';
      const rows = arr.map(([j, k]) => {
        const v = S.isNumeric && S.attrValues ? S.attrValues[k] : null;
        const txt = v == null ? '' : (S.attr === 'prob' ? v.toFixed(2) : fmtNum(v));
        return `<li data-node="${j}"><span>${esc(S.net.nodes[j].n)}</span><span class="v">${txt}</span></li>`;
      });
      return `<ul>${rows.join('')}</ul>`;
    };
    el('selection').innerHTML =
      `<div class="name">${esc(nd.n)}</div>` +
      `<div class="meta">列インデックス ${nd.i} · 親 ${nd.in} / 子 ${nd.out}` +
      `${nd.t ? ' · 注目遺伝子' : ''}</div>` +
      `<h4>親 (→ このノード)</h4>${list(S.parents[i])}` +
      `<h4>子 (このノード →)</h4>${list(S.children[i])}`;
    el('selection').querySelectorAll('li[data-node]').forEach((li) => {
      li.addEventListener('click', () => {
        const j = +li.dataset.node;
        selectNode(j);
        S.graph.zoomToPointByIndex(j, 350, 4, true);
      });
    });
  }

  function doSearch() {
    const q = el('search').value.trim().toLowerCase();
    const ul = el('searchResults');
    if (!q) { ul.innerHTML = ''; return; }
    const hits = [];
    S.net.nodes.forEach((nd, i) => {
      if (nd.n.toLowerCase().includes(q)) hits.push([i, nd]);
    });
    hits.sort((a, b) => S.deg[b[0]] - S.deg[a[0]]);
    if (!hits.length) { ul.innerHTML = '<li class="empty">該当なし</li>'; return; }
    ul.innerHTML = hits.slice(0, 40).map(([i, nd]) =>
      `<li data-node="${i}"><span>${esc(nd.n)}</span>` +
      `<span class="deg">親${nd.in}/子${nd.out}</span></li>`).join('');
    ul.querySelectorAll('li[data-node]').forEach((li) => {
      li.addEventListener('click', () => {
        const j = +li.dataset.node;
        selectNode(j);
        S.graph.zoomToPointByIndex(j, 400, 4, true);
      });
    });
  }

  /* ---------- ツールチップ ---------- */
  function showTooltip(i, spacePos) {
    const nd = S.net.nodes[i];
    const t = el('tooltip');
    t.innerHTML = `<div>${esc(nd.n)}</div>` +
      `<div class="t2">親 ${nd.in} / 子 ${nd.out} · 列 ${nd.i}${nd.t ? ' · 注目' : ''}</div>`;
    t.hidden = false;
    try {
      const [x, y] = S.graph.spaceToScreenPosition(spacePos);
      t.style.left = (x + 12) + 'px';
      t.style.top = (y + 10) + 'px';
    } catch (e) { /* 初期化直後は座標変換できないことがある */ }
  }

  /* ---------- ラベル ---------- */
  function updateLabelSet() {
    const set = new Set();
    if (S.labelCount > 0) {
      const byDeg = S.net.nodes.map((_, i) => i).sort((a, b) => S.deg[b] - S.deg[a]);
      byDeg.slice(0, S.labelCount).forEach((i) => set.add(i));
    }
    if (S.targets) S.net.nodes.forEach((nd, i) => { if (nd.t) set.add(i); });
    if (S.selection != null) {
      set.add(S.selection);
      S.parents[S.selection].forEach(([j]) => set.add(j));
      S.children[S.selection].forEach(([j]) => set.add(j));
    }
    S.labelIdx = [...set].slice(0, 250);
    S.graph.trackPointPositionsByIndices(S.labelIdx);
    drawLabels();
  }

  function drawLabels() {
    const box = el('labels');
    if (!S.labels || !S.net || !S.graph) { box.innerHTML = ''; S.pool = []; return; }
    let map;
    try { map = S.graph.getTrackedPointPositionsMap(); } catch (e) { return; }
    const w = box.clientWidth, h = box.clientHeight;
    let k = 0;
    for (const i of S.labelIdx) {
      const p = map.get(i);
      if (!p || Number.isNaN(p[0])) continue;
      let sp;
      try { sp = S.graph.spaceToScreenPosition(p); } catch (e) { return; }
      if (sp[0] < -60 || sp[1] < -20 || sp[0] > w + 60 || sp[1] > h + 20) continue;
      let span = S.pool[k];
      if (!span) {
        span = document.createElement('span');
        box.appendChild(span);
        S.pool[k] = span;
      }
      const nd = S.net.nodes[i];
      span.textContent = nd.n;
      span.className = i === S.selection ? 'sel' : (nd.t && S.targets ? 'tgt' : '');
      span.style.opacity = (S.selection == null || S.neighborSet.has(i)) ? '1' : '0.3';
      span.style.left = sp[0] + 'px';
      span.style.top = sp[1] + 'px';
      span.hidden = false;
      k++;
    }
    for (let j = k; j < S.pool.length; j++) S.pool[j].hidden = true;
  }

  function loopLabels() {
    const tick = () => { drawLabels(); S.raf = requestAnimationFrame(tick); };
    if (!S.raf) S.raf = requestAnimationFrame(tick);
  }

  /* ---------- UI 配線 ---------- */
  function wireUI() {
    el('edgeAttr').addEventListener('change', (e) => {
      S.attr = e.target.value;
      el('topPct').disabled = (S.attr === 'eval' || S.attr === 'none');
      refreshFilter(false);
    });
    el('topPct').addEventListener('input', (e) => {
      S.topPct = +e.target.value;
      el('topPctVal').textContent = S.topPct + '%';
      refreshFilter(false);
    });
    el('fmSimple').addEventListener('click', () => setFilterMode('simple'));
    el('fmAdv').addEventListener('click', () => setFilterMode('advanced'));
    el('cbAnd').addEventListener('click', () => {
      S.combine = 'and';
      el('cbAnd').classList.add('on'); el('cbOr').classList.remove('on');
      refreshFilter(true);
    });
    el('cbOr').addEventListener('click', () => {
      S.combine = 'or';
      el('cbOr').classList.add('on'); el('cbAnd').classList.remove('on');
      refreshFilter(true);
    });
    el('addCond').addEventListener('click', () => {
      const keys = condKeys();
      if (!keys.length) return;
      /* まだ使っていない属性を優先して追加する */
      const used = new Set(S.conds.map((c) => c.key));
      const key = keys.find((k) => !used.has(k)) || keys[0];
      S.conds.push(defaultCond(key));
      refreshFilter(true);
    });
    el('clearCond').addEventListener('click', () => {
      S.conds = []; refreshFilter(true);
    });

    el('search').addEventListener('input', doSearch);
    el('clearSel').addEventListener('click', () => selectNode(null));

    el('tgLabels').addEventListener('change', (e) => {
      S.labels = e.target.checked; el('labels').innerHTML = ''; S.pool = []; drawLabels();
    });
    el('tgArrows').addEventListener('change', (e) => {
      S.arrows = e.target.checked;
      S.graph.setConfigPartial({ linkDefaultArrows: S.arrows });
      applyEdgeStyles();
    });
    el('tgCurved').addEventListener('change', (e) => {
      S.curved = e.target.checked;
      S.graph.setConfigPartial({ curvedLinks: S.curved });
    });
    el('tgTargets').addEventListener('change', (e) => {
      S.targets = e.target.checked;
      S.graph.setPointColors(pointColors());
      S.graph.render(undefined, 0);
      updateLabelSet();
    });
    el('labelCount').addEventListener('input', (e) => {
      S.labelCount = +e.target.value;
      el('labelCountVal').textContent = S.labelCount;
      updateLabelSet();
    });

    el('btnPause').addEventListener('click', () => {
      S.paused = !S.paused;
      if (S.paused) S.graph.pause(); else S.graph.unpause();
      el('btnPause').textContent = S.paused ? '再開' : '一時停止';
    });
    el('btnFit').addEventListener('click', () => {
      S.userZoomed = false; S.graph.fitView(500, 0.22);
    });
    el('btnRestart').addEventListener('click', () => {
      S.graph.setPointPositions(S.positions);
      S.graph.render();
      S.graph.start(1);
      S.paused = false;
      S.userZoomed = false;
      el('btnPause').textContent = '一時停止';
      scheduleFits();
    });
    el('repulsion').addEventListener('input', (e) => {
      el('repulsionVal').textContent = e.target.value;
      S.graph.setConfigPartial({ simulationRepulsion: +e.target.value });
      S.graph.start(0.2);
    });
    el('linkDist').addEventListener('input', (e) => {
      el('linkDistVal').textContent = e.target.value;
      S.graph.setConfigPartial({ simulationLinkDistance: +e.target.value });
      S.graph.start(0.2);
    });

    document.addEventListener('keydown', (e) => {
      if (e.target.tagName === 'INPUT' || e.target.tagName === 'SELECT') return;
      if (e.key === 'f') S.graph.fitView(500, 0.25);
      if (e.key === 'Escape') selectNode(null);
      if (e.key === ' ') { e.preventDefault(); el('btnPause').click(); }
    });
    document.addEventListener('change', (e) => {
      if (!e.target.classList || !e.target.classList.contains('histlog')) return;
      S.histLog = e.target.checked;
      refreshFilter(S.filterMode !== 'simple');
    });

    window.addEventListener('resize', () => {
      drawLabels();
      if (S.net) refreshFilter(S.filterMode !== 'simple');
    });
  }

  /* ---------- 小物 ---------- */
  const esc = (s) => String(s).replace(/[&<>"']/g,
    (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
  const fmtInput = (v) => (Math.abs(v) >= 100 ? v.toFixed(1)
    : Math.abs(v) >= 1 ? v.toFixed(3) : v.toPrecision(3));
  const fmtNum = (v) => (Math.abs(v) >= 1000 ? v.toFixed(0)
    : Math.abs(v) >= 1 ? v.toFixed(2) : v.toPrecision(2));

  boot();
})();
