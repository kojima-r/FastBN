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
| `example/`                          | Bash     | Minimal `fast_bn` walkthrough on pre-discretized mouse expression data |
| `example_bulk/`                     | Bash+Py  | End-to-end pipeline example on generated dummy bulk RNA counts (true network known, so accuracy is measurable) — see `example_bulk/README.md` |

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


