# 🌽 K-mers GWAS Pipeline (based on Voichek et al. 2020)

This repository provides a SLURM-ready pipeline for running **k-mer-based genome-wide association studies (GWAS)**, using the [`kmersGWAS`](https://github.com/voichek/kmersGWAS) library developed by [Voichek and Weigel (2020)](https://www.nature.com/articles/s41587-019-0380-0).

The pipeline is optimized for **high-performance computing (HPC)** environments and was tested on whole-genome-sequenced maize sweet corn and DArT-sequenced potato accessions using UF’s SLURM cluster.

---

## 📦 Based on

- **Original library:** [`voichek/kmersGWAS`](https://github.com/voichek/kmersGWAS)  
- **Authors:** Yoav Voichek and Detlef Weigel  
- **License:** GPL-3.0  
- **Citation:**  
  > Voichek Y, Weigel D. (2020). Identifying genetic variants underlying phenotypic variation in plants without complete genomes. *Nature Biotechnology*.

---

## 🚀 Overview of Pipeline Steps

| Step | Description |
|------|-------------|
| **1** | Prepare sample input files (FASTQ or preprocessed) |
| **2** | Run `kmc` to count k-mers (canonical and non-canonical) |
| **3** | Add strand information |
| **4** | Filter k-mers found in multiple samples |
| **5** | Build presence/absence matrix |
| **6** | Compute kinship matrix |
| **7** | Convert matrix to PLINK format (optional) |
| **8** | Run GWAS using `pipeline.py` with permutation thresholds |

---

## ⚙️ Requirements

### 🧠 Software

- Linux 64-bit system
- [KMC 3](https://github.com/refresh-bio/KMC)
- GEMMA
- Python 2.7
- R ≥ 4.0
- `kmersGWAS` compiled binaries (placed in `bin/`)

### 📚 R Packages

The following R packages are auto-installed by the pipeline if not already available:

- `MASS`
- `Mvnpermute`
- `matrixcalc`

---

## 📁 Repository Structure

kmers-gwas-pipeline/
├── kmers_01.sh # SLURM job script (main pipeline)
├── bin/ # Compiled kmersGWAS binaries
├── example_input/ # Example phenotype/sample files
├── external_programs/ # KMC and GEMMA binaries
├── README.md
└── docs/ # Optional documentation or diagrams

---

## 📂 Input Requirements

- `Sample-genotypes.txt` with one sample name per line
- Paired-end FASTQ files named as:
  - `sample_R1.fastq.gz`
  - `sample_R2.fastq.gz`
- Phenotype files (`*.pheno`) with numeric values in the same sample order

---

## 📊 Output

- Filtered k-mers list (`kmers_to_use`)
- K-mers presence/absence matrix (`kmers_table`)
- Kinship matrix (`kmers_table.kinship`)
- (Optional) PLINK `.bed/.bim/.fam` files
- GWAS summary results with permutation-adjusted thresholds

---

## 💡 Examples

This pipeline was tested on:

- **WGS Maize sweet corn**
- **DArT-sequenced potato accessions**

Traits analyzed include:

`AC` (Anthers color), `COV` (Convexity), `CUR` (Ear curvature), `DTP` (Days to pollination), `DTS` (Days to silk), `EH` (Ear height), `EL` (Ear length), `EW` (Ear width), `FLH` (Flag leaf height), `GER` (Germination), `LA` (Leaf angle), `PH` (Plant height), `SC` (Silk color), `SL` (Shank length), `SOL` (Ear solidity), `TBN` (Tassel branch number), `TE` (Tassel extension), `TIL` (Tiller number), `TP` (Taper), `TPF` (Tip fill), among others.

The pipeline is adaptable to any species with FASTQ data and matching phenotype data.

---

## 🙏 Acknowledgments

This work uses the **kmersGWAS** library developed by **Yoav Voichek** and **Detlef Weigel** (2020):  
🔗 https://github.com/voichek/kmersGWAS

> We thank the authors for making the tools publicly available under the GPL-3.0 license.
![K-mers GWAS pipeline diagram](docs/pipeline_diagram.png)
*K-mer-based GWAS workflow. Adapted from Corut & Wallace, 2024*
---

## 📜 License

This repository is distributed under the same [GPL-3.0 license](https://www.gnu.org/licenses/gpl-3.0.en.html) as the original `kmersGWAS` library.

If you use this pipeline in published work, please **cite** the original authors.
