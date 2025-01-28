<p  align="center"><img src="https://github.com/user-attachments/assets/8c997c9e-d74f-4521-8bf1-8173b777bb86" style="width: 50%; height: auto;">

<div align="center">

**An efficient assembly toolkit for organellar genomes**

</div>

***

<div align="center">

[![Release Version](https://img.shields.io/github/v/release/aiPGAB/PMAT2?style=flat-square)](https://github.com/aiPGAB/PMAT2/releases)
[![License](https://img.shields.io/github/license/aiPGAB/PMAT2?style=flat-square)](LICENSE)
[![Last Commit](https://img.shields.io/github/last-commit/aiPGAB/PMAT2?style=flat-square)](https://github.com/aiPGAB/PMAT2/commits)

</div>

If you encounter any problems in using PMAT2, please contact the authors by e-mail (Changwei Bi: bichwei@njfu.edu.cn; Fuchuan Han: hanfc@caf.ac.cn) to join the WeChat group (please note your name + organization + PMAT2 in the message).


- [PMAT2](#h1)
- [Installation](#C1)
- [Repuirement](#C2)
- [Options and usage](#C3)
    - [autoMito](#C4)
    - [graphBuild](#C5)
- [Examples](#C6)
  - [Demo1](#C6.1)
  - [Demo2](#C6.2)
  - [Demo3](#C6.3)
  - [Demo4](#C6.4)
- [Output files](#C7)
- [Version](#C8)
- [Citing PMAT2](#C9)

## <a name="C1">Installation </a>

Install using git
```sh
git clone https://github.com/aiPGAB/PMAT2
cd PMAT2
make
./PMAT --help
```
Install by downloading the source codes
```sh
wget https://github.com/aiPGAB/PMAT2/archive/refs/tags/v2.0.5.tar.gz
tar -zxvf PMAT2-2.0.5.tar.gz
cd PMAT2-2.0.5
make
./PMAT --help
```

## <a name="C2">Requirement</a>

- [**BLASTn > 2.2.29**](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)  Needs to be installed in `PATH`.
- [**Singularity**](https://github.com/YanshuQu/runAssembly) or [**Apptainer**](https://github.com/apptainer/apptainer/blob/main/INSTALL.md) is required for PMAT2. You can find installation instructions [here](https://github.com/YanshuQu/runAssembly).
- [**Canu > v2.0**](https://github.com/marbl/canu) or [**NextDenovo**](https://github.com/Nextomics/NextDenovo) is required for CLR or ONT sequencing data.
- [**zlib**](https://www.zlib.net/) Needs to be installed in `PATH`.

## <a name="C3">Options and usage</a>

### <a name="C4">autoMito</a>
    
Run `PMAT autoMito --help` to view the usage guide.

```
Usage: PMAT autoMito [-i INPUT] [-o OUTPUT] [-t SEQTYPE] [options]
Example:
       PMAT autoMito -i hifi.fastq.gz -o hifi_assembly -t hifi -m -T 8
       PMAT autoMito -i ont.fastq.gz -o ont_assembly -t ont -S nextdenovo -C canu -N nextdenovo
       PMAT autoMito -i clr.fastq.gz -o clr_assembly -t clr -S canu -C canu

Required options:
   -i, --input          Input sequence file (fasta/fastq)
   -o, --output         Output directory
   -t, --seqtype        Sequence type (hifi/ont/clr)

Optional options:
   -k, --kmer           kmer size for estimating genome size (default: 31)
   -g, --genomesize     Genome size (g/m/k), skip genome size estimation if set
   -p, --task           Task type (0/1), skip error correction for ONT/CLR by selecting 0, otherwise 1 (default: 1)
   -x, --taxo           Specify the organism type (0/1), 0: plants, 1: animals (default: 0)
   -S, --correctsoft    Error correction software (canu/nextdenovo, default: nextdenovo)
   -C, --canu           Canu path
   -N, --nextdenovo     NextDenovo path
   -n, --cfg            Config file for nextdenovo (default: temprun.cfg)
   -F, --factor         Subsample factor (default: 1)
   -D, --subseed        Random number seeding when extracting subsets (default: 6)
   -K, --breaknum       Break long reads (>30k) with this (default: 20000)
   -I, --minidentity    Set minimum overlap identity (default: 90)
   -L, --minoverlaplen  Set minimum overlap length (default: 40)
   -T, --cpu            Number of threads (default: 8)
   -m, --mem            Keep sequence data in memory to speed up computation
   -h, --help           Show this help message and exit
```

**Notes**:
1. Make sure BLASTn was installed in PATH.
2. If you want to use nextdenovo for ONT/CLR error correction, you can skip providing a cfg file, and the program will generate a temporary cfg file automatically.
3. `-k`: If seqtype is hifi, skip kmer frequency estimation and genome size estimation.
4. `-m`: Keep sequence data in memory to speed up computation.
5. `-I`: The default value is 90 bp. If the assembly graph is complex, you can increase it appropriately.
6. `-L`: minimum overlap identity, the default is 40, if it is HiFi data, you can increase it appropriately.


### <a name="C5">graphBuild</a>

If PMAT fails to generate the assembly graph in 'autoMito' mode, you can use this command to manually select seeds for assembly.

Run `PMAT graphBuild --help` to view the usage guide.

```
Usage: PMAT graphBuild [-i SUBSAMPLE] [-a ASSEMBLY] [-o OUTPUT] [options]
Example:
       PMAT graphBuild -i assembly_test1/subsample -a assembly_test1/assembly_result -o graphBuild_result -s 1 312 356 -T 8
       PMAT graphBuild -i assembly_test1/subsample -a assembly_test1/assembly_result -o graphBuild_result -d 5 -s 1 312 356 -T 8

Required options:
   -i, --subsample     Input subsample directory (assembly_test1/subsample)
   -a, --graphinfo     Input assembly result directory (assembly_test1/assembly_result)
   -o, --output        Output directory

Optional options:
   -G, --organelles     Genome organelles (mt: mitochondria/pt: plastid, default: mt)
   -x, --taxo           Specify the organism type (0/1), 0: plants, 1: animals (default: 0)
   -d, --depth          Contig depth threshold
   -s, --seeds          ContigID for extending. Multiple contigIDs should be separated by space. For example: 1 312 356
   -T, --cpu            Number of threads (default: 8)
   -h, --help           Show this help message and exit
```
**Notes**:
1. Make sure BLASTn was installed in PATH.
2. `-i`: assembly_test1/subsample generated by autoMito command.
3. `-a`: assembly_test1/assembly_result generated by autoMito command.
4. `-s`: Manually select the seeds for the extension. Use spaces to split between different seed IDs, e.g. 1,312,356.

## <a name="C6">Examples</a>

**<a name="C6.1">Demo1</a>**

1. Download a simulated Arabidopsis thaliana HiFi dataset:
```sh
wget https://github.com/bichangwei/PMAT/releases/download/v1.1.0/Arabidopsis_thaliana_550Mb.fa.gz
```
2. then run the autoMito command for one-click assembly:
```sh
PMAT autoMito -i Arabidopsis_thaliana_550Mb.fa.gz -o ./test1 -t hifi -m
```
3. then use the graphBuild command to manually select seeds for assembly (used when the autoMito command fails to get a GFA file automatically):
```sh
# Based on the PMATContigGraph.txt file, manually select 3 or more contigs that match the depth of mitochondrial genome sequencing
PMAT graphBuild -i ./test1/subsample/ -a ./test1/assembly_result/ -o ./test1_gfa -s 1 2 3
```


**<a name="C6.3">Demo2</a>**

1. Download a simulated Malus domestica HiFi dataset:
```sh
wget https://github.com/bichangwei/PMAT/releases/download/v1.1.0/Malus_domestica.540Mb.fasta.gz
```
2. then run the autoMito command for one-click assembly:
```sh
PMAT autoMito -i Malus_domestica.540Mb.fasta.gz -o ./test3 -t hifi -m
```
3. then use the graphBuild command to manually select seeds for assembly (used when the autoMito command fails to get gfa automatically):
```sh
# Based on the PMATContigGraph.txt file, manually select 3 or more contigs that match the depth of mitochondrial genome sequencing
PMAT graphBuild -i ./test3/subsample/ -a ./test3/assembly_result/ -o ./test3_gfa -s 10 20 30
```


**<a name="C6.3">Demo3</a>**

1. Download tested CLR data for Phaseolus vulgaris using IBM Aspera:
```
ascp -v -QT -l 400m -P33001 -k1 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR291/006/SRR2912756/SRR2912756_subreads.fastq.gz .
```
2. then run the autoMito command for one-click assembly (CLR):
```sh
PMAT autoMito -i SRR2912756_subreads.fastq.gz -o ./test_clr -t clr -N path/nextDenovo -m
```

**<a name="C6.4">Demo4</a>**

1. Download tested ONT data for Populus deltoides using IBM Aspera:
```
ascp -v -QT -l 400m -P33001 -k1 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR122/038/SRR12202038/SRR12202038_1.fastq.gz  .
```
2. then run the autoMito command for one-click assembly (ONT):
```sh
PMAT autoMito -i SRR12202038_1.fastq.gz -o ./test_ont -t ont -S canu -C path/canu -m
```
---

|Dataset|Size|Options|Run time|Coverage|
|:-------|:----:|:--------:|:------------:|:-------:|
|Arabidopsis thaliana|550Mb|`-T 50`|6m27s|4x|
|Arabidopsis thaliana|550Mb|`-T 50 -m`|6m38s|4x|
|Malus domestica|540Mb|`-T 50`|7m38s|<1x|
|Malus domestica|540Mb|`-T 50 -m`|7m19s|<1x|
|Juncus effusus|216Mb|`-T 50`|4m56s|<1x|
|Juncus effusus|216Mb|`-T 50 -m`|4m48s|<1x|

## <a name="C7">Output files</a>

```plaintext
output_dir/
├── assembly_result/
│   ├── PMATAllContigs.fna       # Assembly contigs
│   └── PMATContigGraph.txt      # Contig relationships
├── gfa_result/
│   ├── PMAT_mt_raw.gfa          # Initial mitogenome graph
│   ├── PMAT_mt_main.gfa         # Optimized mitogenome graph
│   ├── PMAT_mt.fasta            # Final mitogenome assembly
│   ├── PMAT_pt_raw.gfa          # Initial chloroplast graph
│   ├── PMAT_pt_main.gfa         # Optimized chloroplast graph
│   └── PMAT_pt_main.fa          # Final chloroplast assembly
└── gkmer_result/
    └── gkmer                    # Kmer statistics
```

## <a name="C8">Version</a>

PMAT version 2.0.1 (24/11/21)</br>
Updates:

- Optimized the assembly strategy for organellar genomes, enabling faster and more accurate capture of organellar genome sequences.
- Implemented the assembly of animal and plant organellar genomes.
- Enhanced the genome graph untangling functionality for organellar genomes, enabling resolution of more complex structures.
- Parallelized key steps in the workflow, significantly improving runtime efficiency.


## <a name="C9">Citing PMAT2</a>
Bi C, Shen F, Han F, Qu Y, et al. PMAT: an efficient plant mitogenome assembly toolkit using ultra-low coverage HiFi sequencing data. Horticulture Research. (2024). uhae023, https://doi.org/10.1093/hr/uhae023. </br>
Bi C, Qu Y, Hou J, Wu K, Ye N, and Yin T. (2022). Deciphering the multi-chromosomal mitochondrial genome of Populus simonii. Front. Plant Sci. 13:914635.doi:10.3389/fpls.2022.914635.
