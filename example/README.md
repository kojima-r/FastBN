# Bayesian Network Structure Learning Example (Gene Expression Data)

This directory provides a sample environment for **Bayesian network structure learning**
using **mouse gene expression data**.

## Prerequisites

* This directory is assumed to be the **current working directory**
* The project must be **built in advance** by following `../README.md`
* The `fast_bn` binary must be available and executable at:

```bash
../fast_bn
```

## Directory Structure (Overview)

The directory includes shell scripts for downloading data and running structure learning, discretized datasets, and an output directory for results.

```text
.
├── 00download.sh
├── 01run.sh
├── 01run_tri.sh
├── 02run_bs.sh
├── 02run_bs_tri.sh
├── 03integrate_bs.sh
├── 04score_check.sh
├── 05importance.sh
├── data_bin/ #generated
└── data_tri/ #generated
```

## Step 00: Data Preparation

### 00download.sh
The script `00download.sh` downloads the required tabular data files.
```bash
./00download.sh
```

* Downloads the required **table-format data files**
* The data are based on **mouse gene expression datasets**

### Data Files Used for Structure Learning

The following files are mainly used for structure learning in this example:

* `data_bin/all_disc100.tsv`
* `data_tri/all_disc_tri100.tsv`

These datasets are **restricted to 100 variables**.
Other variants are also included:

| File name      | Description    |
| -------------- | -------------- |
| `all_disc10`   | 10 variables   |
| `all_disc100`  | 100 variables  |
| `all_disc1000` | 1000 variables |
| `all_disc`     | All variables  |

You can inspect the contents using tools such as `less`.
They are simple TSV (tab-separated values) files.

```bash
less data_bin/all_disc100.tsv
```

### Discretization

* The original gene expression data are **continuous**
* They are **discretized** to be used with this Bayesian network tool

| Directory  | Description                                     |
| ---------- | ----------------------------------------------- |
| `data_bin` | Discretized into **2 values (binary)**          |
| `data_tri` | Discretized into **3 values (ternary)** (`tri`) |

## Step 01: Initial Structure Learning

Run the following scripts to perform the initial structure learning:

```bash
./01run.sh
./01run_tri.sh
```

* `01run.sh`: Binary (2-state) data
* `01run_tri.sh`: Ternary (3-state) data

The first script uses binary discretized data, while the second uses ternary discretized data.
These scripts estimate an initial network structure from the data without using any prior structural information.

## Step 02: Bootstrap-Based Continuous Learning

After initial structure learning, the learned network structure can be reused as an initial structure for bootstrap-based continuous learning.

```bash
./02run_bs.sh
./02run_bs_tri.sh
```

### Notes

In this stage, the structure learned in 01run.sh is typically used as the initial network.
However, this initial learning step is not mandatory.
If you want to start bootstrap learning directly, simply comment out the line specifying the initial structure (for example, `--init ./out/edges.tsv`) in the corresponding script.

### Output

* All results are written to the `out/` directory
* Learned network structures (e.g., edge lists) are saved in TSV format

## Post-processing Workflow

The steps described below correspond to **post-processing** performed after bootstrap-based structure learning.
These procedures integrate bootstrap results, recompute model scores, and evaluate the importance of individual edges in the learned Bayesian network.

### Step 03: Integrating Bootstrap Results

Bootstrap results obtained from multiple runs are first aggregated using the following script:

```bash
./03integrate_bs.sh
```

This script performs the following tasks:

* Collects the results from all bootstrap runs
* Integrates them into consensus structures
* Outputs integrated result files with the prefix:

```text
out/integ_*
```

These files summarize the structural information obtained across bootstrap samples.

### Step 04: Recomputing Model Scores

Using the integrated bootstrap results, model scores are recomputed by running:

```bash
./04score_check.sh
```

In the current example:

* Scores are calculated using the **same dataset used for training**
* This allows direct comparison of scores before and after bootstrap integration

More generally, this step can also be applied to:

* Validation datasets
* Train/validation splits
* Other evaluation schemes commonly used in machine learning

### Step 05: Estimating Edge Importance

Next, the importance of each edge in the Bayesian network is calculated:

```bash
./05importance.sh
```

This step evaluates how frequently each edge appears across bootstrap samples and provides a quantitative measure of edge stability.

## Output Files

The final edge importance scores are written to the following file:

```text
./out/edge_importance.tsv
```

This file contains edge-level importance values and can be used for:

* Filtering unreliable edges
* Visualizing network confidence
* Downstream statistical or biological analysis


