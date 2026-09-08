<p align="center"><img src="https://github.com/user-attachments/assets/8c997c9e-d74f-4521-8bf1-8173b777bb86" style="width: 50%; height: auto;"></p>

<div align="center">

**An efficient assembly toolkit for organellar genomes**

</div>

***

<div align="center">

[![Release Version](https://img.shields.io/github/v/release/aiPGAB/PMAT2?style=flat-square)](https://github.com/aiPGAB/PMAT2/releases)
[![License](https://img.shields.io/github/license/aiPGAB/PMAT2?style=flat-square)](LICENSE)
[![Last Commit](https://img.shields.io/github/last-commit/aiPGAB/PMAT2?style=flat-square)](https://github.com/aiPGAB/PMAT2/commits)
[![Wiki Documentation](https://img.shields.io/badge/docs-Wiki-blue?style=flat-square)](https://github.com/aiPGAB/PMAT2/wiki)
</div>

<a href="https://github.com/user-attachments/assets/1e4e48f9-7428-4b2f-a037-5e1f65da5b4e">
  <img src="https://github.com/user-attachments/assets/1e4e48f9-7428-4b2f-a037-5e1f65da5b4e" align="right" width="500" alt="Graphical Abstract">
</a>

PMAT2 is a specialized toolkit for the **de novo assembly of mitochondrial and chloroplast genomes** from PacBio HiFi and ONT/CLR sequencing data. It supports animal, plant, and fungal organellar genome assembly and integrates optimized graphical algorithms for resolving complex structural repeats.

If you encounter any problems using PMAT2, please contact the authors via email to join the user group (include your *name + organization + PMAT2* in the message):

- **Changwei Bi**: bichwei@njfu.edu.cn  
- **Fuchuan Han**: hanfc@caf.ac.cn

## <a name="C9">Citation</a>

Fuchuan Han, Changwei Bi, Yicun Chen, Xiaogang Dai, Zefu Wang, Huaitong Wu, Ning Sun, et al. 2025. PMAT2: An Efficient Graphical Assembly Toolkit for Comprehensive Organellar Genomes. iMeta 4: e70064. https://doi.org/10.1002/imt2.70064 (if you use PMAT2 for organellar genomes)</br>
Bi C, Shen F, Han F, Qu Y, et al. PMAT: an efficient plant mitogenome assembly toolkit using ultra-low coverage HiFi sequencing data. Horticulture Research. (2024). uhae023, https://doi.org/10.1093/hr/uhae023 (if you use PMAT for plant genomes)

## <a name="C1">Installation</a>

PMAT2 version 2.2.0 provides native binary execution on Linux x86_64 without mandatory container dependencies, as well as a prebuilt standalone container for containerized environments.

### Method 1: Download Release Archive

```sh
wget https://github.com/aiPGAB/PMAT2/archive/refs/tags/v2.2.0.tar.gz
tar -zxvf v2.2.0.tar.gz
cd PMAT2-2.2.0
./PMAT --help
```

### Method 2: Clone from GitHub

```sh
git clone --depth 1 https://github.com/aiPGAB/PMAT2.git
cd PMAT2
./PMAT --help
```
*(Optional: If you modify the C source code, simply run `make clean && make` to recompile.)*

### Method 3: Standalone SIF Container (Apptainer / Singularity)

For users who prefer containerized execution, a standalone all-in-one image (`PMAT_v2.2.0.sif`) containing PMAT2, NCBI BLAST+, and Newbler tools is available on GitHub Releases:

```sh
# Download the SIF image from GitHub Releases
wget https://github.com/aiPGAB/PMAT2/releases/download/v2.2.0/PMAT_v2.2.0.sif
chmod +x PMAT_v2.2.0.sif

# Execute directly or via apptainer
./PMAT_v2.2.0.sif autoMito --help
# Or:
apptainer exec PMAT_v2.2.0.sif PMAT autoMito --help
```

> **Note on Container Storage Binding**: If your input or output files are located on external storage mounts (such as `/pub`, `/data`, or `/scratch`), please add `-B /path` when running Apptainer, for example: `apptainer exec -B /pub PMAT_v2.2.0.sif PMAT autoMito ...`

## <a name="C2">Requirement</a>

- **Linux x86_64 Operating System**
- [**BLASTn > 2.10.0**](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download): Required in system `PATH` (already included if using `PMAT_v2.2.0.sif`).
- [**Canu > v2.0**](https://github.com/marbl/canu) or [**NextDenovo**](https://github.com/Nextomics/NextDenovo): Required only for ONT or CLR raw sequencing data error correction.
- [**zlib**](https://www.zlib.net/): Required in system `PATH` (only if compiling from source code).

## <a name="C3">Command Overview</a>

### autoMito

One-step de novo assembly of organellar genomes:

```sh
# HiFi sequencing data
PMAT autoMito -i hifi.fastq.gz -o hifi_assembly -t hifi -m -T 8

# ONT sequencing data
PMAT autoMito -i ont.fastq.gz -o ont_assembly -t ont -S nextdenovo -C canu -N nextdenovo

# CLR sequencing data
PMAT autoMito -i clr.fastq.gz -o clr_assembly -t clr -S canu -C canu
```

### graphBuild

Manual seed extension and assembly graph reconstruction when autoMito requires customized parameters:

```sh
PMAT graphBuild -i assembly_out/subsample -a assembly_out/assembly_result -o graphBuild_result -s 1 312 356 -T 8
```

## <a name="C4">Documentation and Detailed Tutorials</a>

Comprehensive documentation, step-by-step usage workflows, real-world demos, and parameter guidelines have been migrated to the **[PMAT2 GitHub Wiki](https://github.com/aiPGAB/PMAT2/wiki)**:

- [**Quick Start Guide**](https://github.com/aiPGAB/PMAT2/wiki/Quick-Start): Complete assembly workflow and best-practice recommendations.
- [**Real-World Demos and Benchmarks**](https://github.com/aiPGAB/PMAT2/wiki/Tutorials-and-Demos): Detailed walkthroughs and run times for *Arabidopsis thaliana*, *Malus domestica*, *Phaseolus vulgaris*, and *Populus deltoides*.
- [**Command-Line Parameter Reference**](https://github.com/aiPGAB/PMAT2/wiki/Command-Reference): Detailed option descriptions for `autoMito` and `graphBuild`.
- [**Output Files Explanation**](https://github.com/aiPGAB/PMAT2/wiki/Output-Files): Format specifications for assembly contigs, GFA graphs, FASTA results, and assembly assessments.
- [**Troubleshooting and FAQ**](https://github.com/aiPGAB/PMAT2/wiki/Troubleshooting-and-FAQ): Solutions for memory management, storage binding, and assembly unlooping.

## <a name="C5">Output Files</a>

```plaintext
output_dir/
├── assembly_result/
│   ├── PMATAllContigs.fna       # Assembly contigs
│   └── PMATContigGraph.txt      # Contig relationships
├── gfa_result/
│   ├── PMAT_mt_raw.gfa          # Initial mitogenome graph
│   ├── PMAT_mt_main.gfa         # Optimized mitogenome graph
│   ├── PMAT_mt.fa               # Final mitogenome assembly
│   ├── PMAT_pt_raw.gfa          # Initial chloroplast graph
│   ├── PMAT_pt_main.gfa         # Optimized chloroplast graph
│   └── PMAT_pt.fa               # Final chloroplast assembly
├── gkmer_result/
│   ├── gkmer_histo.txt          # Kmer frequency
│   └── summary.txt              # Genome size estimation
├── subsample/
│   └── PMAT_cut_seq.fa          # Subsampled reads for assembly
└── PMAT_orgAss.txt              # Organellar assembly assessment
```

## <a name="C6">Version History</a>

PMAT version 2.0.1 (24/11/21)</br>
Updates:

- Optimized assembly strategy for organellar genomes, enabling faster and more accurate sequence capture.
- Implemented support for animal, plant, and fungal organellar genomes.
- Enhanced genome graph untangling algorithms to resolve complex repeat structures.
- Parallelized key steps in the workflow to improve runtime efficiency.

PMAT version 2.1.0 (25/2/1)</br>
Updates:

- Added `orgAss` module to evaluate the completeness of assembly results.

PMAT version 2.2.0 (26/09/08)</br>
Updates:

- Simplified execution: Removed internal container routing for cleaner, faster execution and easier HPC cluster scheduler integration.
- Standalone container release: Standalone all-in-one `.sif` image provided on GitHub Releases for container-based execution.
