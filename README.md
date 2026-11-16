# An Empirical Benchmarking Analysis of Generative AI Models vs. Traditional Heuristics for Codon Optimization, Structural Stability, and Cost Efficiency in mRNA Vaccine Engineering

[![Empirical mRNA Benchmark Pipeline](https://github.com/hasanmd1/mrna-ai-benchmark/actions/workflows/benchmark.yml/badge.svg)](https://github.com/hasanmd1/mrna-ai-benchmark/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Abstract & Motivation
Traditional mRNA vaccine and therapeutic development require navigating an astronomical sequence design space (e.g., $2.4 \times 10^{632}$ possible candidate sequences for a spike protein). Conventional codon-optimization software relies on rule-based heuristics that are slow, expensive to iterate through laboratory trials, and frequently fail to simultaneously optimize mRNA folding stability ($\Delta G$), GC-content balance, and translational efficiency.

This repository provides the complete, reproducible computational benchmarking framework accompanying our preprint study. It evaluates traditional rule-based heuristics against thermodynamic folding-aware tools and state-of-the-art deep generative AI models across key performance indicators (KPIs).

---

## Research Questions Addressed
- **RQ1**: How much do AI-driven generative sequence engines reduce computational wall-clock time ($T_{\text{design}}$) and wet-lab screening iterations compared to traditional heuristics?
- **RQ2**: What are the comparative trade-offs between AI-optimized sequences and standard baselines regarding minimum free energy ($\Delta G$ folding stability) and GC-content balance across varying protein lengths?
- **RQ3**: How do predicted translation efficiencies (Codon Adaptation Index, CAI) scale when utilizing deep knowledge-infused machine learning models vs. conventional codon tables?
- **RQ4**: What is the estimated economic impact ($\Delta \$$ research cost and timeline reduction) of integrating automated AI sequence optimization upfront in preclinical pipelines?

---

## Repository Architecture

```text
mrna-ai-benchmark/
├── .github/
│   └── workflows/
│       └── benchmark.yml
├── data/
│   └── target_orfs.fasta
├── src/
│   └── empirical_benchmark.py
├── results/
│   ├── empirical_benchmark_metrics.csv
│   └── summary_table.md
├── requirements.txt
└── README.md

```

---

## Experimental Methodology & Execution Matrix

The benchmarking pipeline evaluates target Open Reading Frames (ORFs) across three distinct computational tiers:

1. **Baseline Heuristic**: Conventional single-codon frequency tables.
2. **Intermediate Thermodynamic**: Folding-aware optimization considering local melting temperature windows.
3. **State-of-the-Art (SOTA) Generative AI**: Deep generative models and grammar-based lattice parsers maximizing CAI and structural resistance.

---

## Quickstart & Reproducibility Guide

### 1. Clone the Repository

```bash
git clone [https://github.com/hasanmd1/mrna-ai-benchmark.git](https://github.com/hasanmd1/mrna-ai-benchmark.git)
cd mrna-ai-benchmark

```

### 2. Install Dependencies

```bash
pip install -r requirements.txt

```

### 3. Execute the Benchmarking Script

```bash
python src/empirical_benchmark.py --fasta data/target_orfs.fasta --output_dir results

```

### 4. Inspect Results

Generated raw datasets (`empirical_benchmark_metrics.csv`) and manuscript-ready summary tables (`summary_table.md`) will be available in the `results/` directory.

---

## Citation & License

This project is open-source under the MIT License. If you utilize this benchmarking framework in your preprint or publication, please cite our repository and accompanying manuscript.

```
