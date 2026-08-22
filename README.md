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
| `.claude-plugin/`, `skills/`, `commands/` | Markdown+Bash+Py | Claude Code plugin / marketplace manifests, Agent Skills and slash commands that drive this pipeline from an agent — see [Use from Claude Code / Claude Science](#-use-from-claude-code--claude-science-plugin--skills) |

## 🔌 Use from Claude Code / Claude Science (plugin & skills)

This repository doubles as a Claude Code **plugin marketplace**. `.claude-plugin/marketplace.json`
lists one plugin, `bn-analysis`, whose root is the repository itself — so installing the plugin
brings `fast_bn.cpp`, `script/` and the runnable examples along with the skills, and
`${CLAUDE_PLUGIN_ROOT}` points at a complete FastBN checkout.

### Install as a plugin

```bash
/plugin marketplace add kojima-r/FastBN
/plugin install bn-analysis@fastbn
```

Working from a local checkout (development / testing):

```bash
claude --plugin-dir /path/to/FastBN     # load for one session
claude plugin validate /path/to/FastBN  # check the manifests
```

### What the plugin provides

| Component | Name | Role |
| --- | --- | --- |
| Skill | `expression-network-inference` | An expression matrix arrives → inspect it, derive parameters, run preprocess → learn → importance → bootstrap → per-group → figures → `report.html`, then interpret the result |
| Skill | `network-structure-evaluation` | Score a learned DAG against a known pathway, a gold standard, a BIF network or a simulation's true DAG (SHD, directed/skeleton P-R-F1, SID, KL) |
| Command | `/bn-analysis:analyze <matrix> [meta] [dir]` | End-to-end network inference from a data file |
| Command | `/bn-analysis:setup [--smoke-test]` | Resolve/build `fast_bn`, check Python deps, optional end-to-end smoke test |
| Command | `/bn-analysis:evaluate <true-edges\|bif>` | Compare a learned network with a reference structure |

The skills carry the operational knowledge that is easy to get wrong: node index = column
position, the `--iters 0` idiom, `--seed` (not `--bootstrap-seed`), how `ITERS` /
`MAX_PARENTS` / `N_BINS` / bootstrap counts follow from D, N and the number of bins, and why
edge directions must not be reported as causal claims.

### Use as a plain skills directory (Claude Science, claude.ai, other agents)

`skills/` follows the standard Agent Skills layout (`<name>/SKILL.md` + `references/` +
`scripts/`), so it can be used without the plugin machinery: copy a skill directory into
`~/.claude/skills/` (or `.claude/skills/` for one project), or upload it as a custom skill.
In that mode point the helper scripts at this repository:

```bash
export FASTBN_HOME=/path/to/FastBN
```

### Helper scripts (usable directly)

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


