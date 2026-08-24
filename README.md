# Bayesian Network Structure Learning Tools

A collection of C++ and Python tools for scalable **Bayesian Network structure learning**, evaluation, and postprocessing.

Includes:

* High-performance **Hill-Climbing + Tabu Search** structure learning (`fast_bn`)
* Support for **BIC**, **K2**, and **BDeu** scoring
* Incremental updates, candidate parent pre-selection, and bitset optimizations for large-scale networks (D = 10³–10⁴)
* **Bootstrap structure sampling** for edge stability analysis
* **Edge importance** analysis by score perturbation
* **Python CPT generator** from learned count tables

---

## 📦 Contents

| Component                           | Language | Description                                                          |
| ----------------------------------- | -------- | -------------------------------------------------------------------- |
| `fast_bn`                           | C++17    | Main executable for structure learning, bootstrap, and edge analysis |
| `compute_cpt_from_counts.py`        | Python 3 | Postprocessor: compute CPTs from `--save-counts` output              |
| `compute_bs_prob.py`                | Python 3 | Postprocessor: integrate bootstrap runs into a consensus structure   |
| `init_edges.tsv` / `all_counts.tsv` | TSV      | Intermediate structures saved/loaded between runs                    |
| `script/`                           | Bash+Py  | Reusable pipeline for **bulk RNA expression data**: preprocessing (normalize → log → filter → discretize), learning, edge importance, bootstrap stability, per-group comparison, plots, HTML report — see `script/README.md` |
| `example_sc/`                       | Bash+Py  | End-to-end pipeline example on **pre-discretized single-cell** mouse expression data (switchable dataset / discretization / #genes; groups = tissues) — see `example_sc/README.md` |
| `example_bnlearn/`                  | Bash+Py  | **Benchmark** against the bnlearn repository (`asia`, `cancer`, `earthquake`, `sachs`, `survey`): sample from the true network, learn, and score with SHD / directed & skeleton P-R-F1 / SID / KL — see `example_bnlearn/README.md` |
| `example_sachs/`                    | Bash+Py  | **Benchmark** on the Sachs single-cell protein-signaling data (Zenodo) against the experimentally validated pathway — see `example_sachs/README.md` |
| `example_dream/`                    | Bash+Py  | **Benchmark** on the DREAM challenges (DREAM4 in silico, DREAM5; HPN-DREAM needs manual Synapse download) — see `example_dream/README.md` |
| `example_bulk/`                     | Bash+Py  | End-to-end pipeline example on generated dummy bulk RNA counts (true network known, so accuracy is measurable) — see `example_bulk/README.md` |
| `viewer/`                           | Py+JS    | Interactive network viewer built on [cosmos.gl](https://github.com/cosmosgl/graph) (GPU force graph). `python3 viewer/serve.py` finds every learned network under the repository and opens it in the browser — see `viewer/README.md` |
| `skills/`, `commands/`, `.claude-plugin/`, `.codex-plugin/` | Markdown+Bash+Py | Agent Skills (shared by Claude Code and OpenAI Codex), Claude slash commands, and both plugin/marketplace manifests — see [Use from coding agents](#-use-from-coding-agents-claude-code--claude-science--codex) |

## 🔌 Use from coding agents (Claude Code / Claude Science / Codex)

The analysis workflow is packaged as **Agent Skills** in `skills/`, which is the single source of
truth for both ecosystems (the `SKILL.md` format is common to Claude Code and OpenAI Codex):

| Skill | Role |
| --- | --- |
| `expression-network-inference` | An expression matrix arrives → inspect it, derive parameters, run preprocess → learn → importance → bootstrap → per-group → figures → `report.html`, then interpret the result |
| `network-structure-evaluation` | Score a learned DAG against a known pathway, a gold standard, a BIF network or a simulation's true DAG (SHD, directed/skeleton P-R-F1, SID, KL) |

They carry the operational knowledge that is easy to get wrong: node index = column position, the
`--iters 0` idiom, `--seed` (not `--bootstrap-seed`), how `ITERS` / `MAX_PARENTS` / `N_BINS` /
bootstrap counts follow from D, N and the number of bins, and why edge directions must not be
reported as causal claims.

### Inside this repository — no configuration needed

`.claude/skills/` and `.codex/skills/` hold relative symlinks to `skills/`, so both agents discover
the skills as project skills as soon as they are started in this directory. `AGENTS.md` (Codex) and
`CLAUDE.md` (Claude Code) carry the repository-level invariants.

> Git checkouts on filesystems without symlink support (Windows without `core.symlinks`) get plain
> text files instead; copy `skills/<name>` into `.codex/skills/` / `.claude/skills/` there.

### Claude Code / Claude Science

This repository is also a Claude Code **plugin marketplace**: `.claude-plugin/marketplace.json`
lists one plugin, `bn-analysis`, whose root is the repository itself — so installing it brings
`fast_bn.cpp`, `script/` and the runnable examples along with the skills, and
`${CLAUDE_PLUGIN_ROOT}` points at a complete FastBN checkout.

```bash
/plugin marketplace add kojima-r/FastBN
/plugin install bn-analysis@fastbn
```

The plugin adds three slash commands: `/bn-analysis:analyze <matrix> [meta] [dir]` (end-to-end
inference), `/bn-analysis:setup [--smoke-test]` (build the binary, check deps) and
`/bn-analysis:evaluate <true-edges|bif>` (compare against a reference structure).

Local checkout (development / testing):

```bash
claude --plugin-dir /path/to/FastBN     # load the plugin for one session
claude plugin validate /path/to/FastBN  # check the manifests
```

### OpenAI Codex

The same skills work in Codex CLI. Personal install:

```bash
cp -r skills/expression-network-inference "${CODEX_HOME:-$HOME/.codex}/skills/"
cp -r skills/network-structure-evaluation "${CODEX_HOME:-$HOME/.codex}/skills/"
```

Or install the whole thing as a Codex plugin — `.codex-plugin/plugin.json` is the manifest and
`.agents/plugins/marketplace.json` is the repo marketplace:

```bash
codex plugin marketplace add kojima-r/FastBN   # or a local path to this checkout
codex plugin add bn-analysis@fastbn
codex plugin list
```

Each skill ships `agents/openai.yaml` (display name, short description, `$skill` starter prompt)
for Codex's skill list. Invoke a skill explicitly with `$expression-network-inference` /
`$network-structure-evaluation`, or let Codex trigger it from the description.

> A **local** plugin install copies the working tree as it is, including downloaded datasets and
> generated run directories (this checkout is several hundred MB). For distribution prefer the Git
> source, which carries only tracked files.

Validate a change to the Codex manifests with the CLI's own system skills:

```bash
python3 "${CODEX_HOME:-$HOME/.codex}/skills/.system/plugin-creator/scripts/validate_plugin.py" .
python3 "${CODEX_HOME:-$HOME/.codex}/skills/.system/skill-creator/scripts/quick_validate.py" \
        skills/expression-network-inference
```

### Helper scripts (usable directly, no agent involved)

```bash
S=skills/expression-network-inference/scripts

bash $S/fastbn_env.sh --check                 # locate FastBN, build the binary, check deps
python3 $S/inspect_matrix.py --input counts.tsv --meta sample_meta.tsv
                                              # orientation / value type / N / D / variance
                                              # + recommended parameters with reasoning
bash $S/new_analysis.sh my_analysis --expr counts.tsv --meta sample_meta.tsv
                                              # scaffold an analysis dir with a filled config.sh
```

`inspect_matrix.py` decides raw counts vs. normalized vs. already-log-transformed vs. already
discretized, separates annotation columns (including numeric ones such as `gene_length`), and
derives `NORMALIZE` / `N_BINS` / `MAX_PARENTS` / `TOP_VAR_GENES` / `ITERS` / bootstrap counts
from N, D and the variance distribution.

---

## 🔭 Interactive viewer (cosmos.gl)

`figures/*.png` are static; `viewer/` is for exploring a network interactively. It renders the
learned DAG with [cosmos.gl](https://github.com/cosmosgl/graph) — a WebGL force graph that runs
the layout on the GPU.

```bash
python3 viewer/serve.py        # scans the repo, opens the browser, no build step, no network access
```

It starts on the learned network of `example_bulk/out`, and the dropdown lists **every** network
found under the repository (all `example_*` benchmark runs included — 350+ in a full checkout), so
the already-computed examples can be inspected right away. If nothing has been learned yet, run
`example_bulk/run_all.sh` (1–2 min) first.

* Edge color and width: `|ΔlogL|` / `|ΔBIC|` / `|ΔK2|` / `|ΔBDeu|` (same absolute values as
  `visualize.py`), bootstrap probability, or the TP/FP/FP_reversed/FN verdict from
  `eval_*_edges.tsv` — rank-normalized.
* Edge filtering has two modes: **simple** (top X% by the attribute used for color) and
  **advanced**, which combines thresholds on *several* metrics at once with AND/OR — e.g.
  "`|ΔlogL|` in the top 30% **and** bootstrap probability ≥ 0.8". Each condition takes a
  percentile or an absolute threshold (`≥`/`≤`), verdict conditions take a set of
  TP/FP/FP_reversed/FN, and every condition reports how many edges it passes on its own.
* Both modes draw a **histogram of the metric** so a threshold can be chosen from the actual
  distribution: grey bars are the distribution, the blue overlay is the part that passes, an
  amber curve shows how many edges would survive at each cut point, and a red marker sits at the
  current threshold. Click or drag on the histogram to set the threshold; importance
  distributions are heavily skewed, so a log count axis is one checkbox away.
* Node color: source (no parents) / internal / sink (no children) / target gene; size = degree.
* Click or search a gene to see its parents and children with their importance values, and to
  focus its neighborhood; arrows point parent → child.
* `python3 viewer/serve.py --root my_analysis` limits the scan to one analysis directory;
  `--list` prints what was found; `--select <id>` picks the network to open with.

Everything is Python standard library plus a vendored cosmos.gl UMD build (MIT) in
`viewer/vendor/`, so there is no npm/bundler step and the viewer works offline. Details, including
which files each network is assembled from, are in [`viewer/README.md`](viewer/README.md).

---

### Analyzing your own bulk RNA expression data

```bash
mkdir my_analysis && cd my_analysis
cp ../example_bulk/config.sh .        # edit EXPR_INPUT / SAMPLE_META / ID_COL ...
source ./config.sh
../script/run_pipeline.sh             # preprocess → learn → importance → bootstrap → plots → report
```

Each stage is also a standalone script (`../script/preprocess.sh`,
`learn_structure.sh`, `edge_importance.sh`, `bootstrap_stability.sh`,
`importance_groups.sh`, `viz*.sh`, `make_report.sh`) configured through
environment variables. See `script/README.md` for the full list.

---

## 🚀 Build

Requires **g++ 13+** and standard libraries only (no external deps).

```bash
g++ -O3 -march=native -std=c++17 fast_bn.cpp -o fast_bn
```
(see `compile.sh`)

If you want to compile with the gprof profiler, see `compile_prof.sh`

---

## 🧩 Usage Overview

### Structure Learning

```bash
./fast_bn --input data.tsv --score bic \
  --iters 3000 --tabu 20 \
  --max-parents 3 --max-children 4 \
  --topk 50 --mi-sample 8000 --mi-budget 1500 \
  --reach lazy --jindex-cache 32 \
  --verbose
```

**Supported scores:** `bic`, `k2`, `bdeu`

Output files:

* `init_edges.tsv` – learned DAG edges
* `all_counts.tsv` – node-level counts used for CPT estimation

---

### Bootstrap Mode

Run multiple structure learning rounds on bootstrap-resampled datasets to estimate edge stability.

```bash
./fast_bn --input data.tsv --score bic \
  --bootstrap 100 \
  --bootstrap-seed 42 \
  --save-bootstrap-counts results/boot_edges.tsv
```

Output:
`results/boot_edges_seed0042.tsv` — columns `u, v, count, prob` (edge bootstrap frequency).

---

### Edge Importance Mode

Evaluates each edge’s contribution by removing it and recomputing scores.

```bash
./fast_bn --input test.tsv --score bic \
  --edge-importance \
  --init init_edges.tsv \
  --counts all_counts.tsv \
  --save-edge-importance edge_importance.tsv
```

Output columns:

```
u  v  ΔlogL  ΔBIC  ΔK2  ΔBDeu  meanΔlogL_per_sample  stdΔlogL_per_sample
```

---

### Candidate Parent Selection

You can pre-select parent candidates using either **Mutual Information (MI)** or **Chi-square p-value** criteria:

```bash
--cand-metric mi|chi2
--topk <K>                    # keep top K candidates
--mi-threshold <t>            # MI ≥ t (nats)
--chi2-p-threshold <p>        # p ≤ p_th (default=0.05)
--mi-sample <n>               # subsample rows for MI/Chi2 computation
--mi-budget <n>               # max #candidates tested
```

---

### Runtime / Logging

```bash
--verbose                # print detailed logs
--reach dense|lazy       # reachability check mode
--jindex-cache <cap>     # cache size for parent configurations
--alpha <a>              # smoothing coefficient
--ess <val>              # equivalent sample size (BDeu only)
```

---

### Help and Language Mode

```bash
./fast_bn --help
```

Language of help is auto-selected based on environment:

| Variable               | Language |
| ---------------------- | -------- |
| `BN_LANG=ja` (default) | Japanese |
| `BN_LANG=en`           | English  |

---

## 🧠 File Formats

### Input Dataset (`--input`)

* CSV or TSV with **header row**
* Each column = discrete variable (already integer-encoded)

Example:

```
A   B   C
0   1   2
1   0   1
...
```

### Learned Structure (`init_edges.tsv`)

```
u   v
0   1
1   2
```

### Counts (`all_counts.tsv`)

```
v   j   k|*   n
0   0   0     120
0   0   1     45
0   0   *     165
```

---

## 🧮 Postprocessing: Compute CPTs (Python)

`compute_cpt_from_counts.py` converts `--save-counts` output to normalized conditional probabilities.

### Example

```bash
python3 compute_cpt_from_counts.py \
  --counts all_counts.tsv \
  --out cpts.tsv \
  --alpha 1.0
```

Or per-node:

```bash
python3 compute_cpt_from_counts.py \
  --counts all_counts.tsv \
  --out-dir cpts/ \
  --alpha 0.0 --skip-nan
```

### Options

| Option        | Description                               | Default      |
| ------------- | ----------------------------------------- | ------------ |
| `--counts`    | Input `all_counts.tsv`                    | *(required)* |
| `--out`       | Single output TSV (`v j k prob`)          | stdout       |
| `--out-dir`   | Output per-node CPT files (`cpt_<v>.tsv`) | —            |
| `--alpha`     | Dirichlet smoothing                       | `0.0`        |
| `--precision` | Decimal precision                         | `12`         |
| `--skip-nan`  | Skip undefined (nij=0) rows               | off          |
| `--quiet`     | Suppress logs                             | off          |

---

## ⚙️ Implementation Details

* **Language:** C++17 (gcc 13+)
* **Design:** Incremental scoring with full REMOVE updates
* **Optimization:**

  * Candidate-parent precomputation (MI/Chi2)
  * Bitset closure for reachability
  * Cached parent-configuration indexing
  * Optional tabu search and bootstrap
* **Scalability:** Designed for D up to 10⁴–3×10⁴ nodes


