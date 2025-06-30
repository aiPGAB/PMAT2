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


<a href="https://github.com/user-attachments/assets/1e4e48f9-7428-4b2f-a037-5e1f65da5b4e">
  <img src="https://github.com/user-attachments/assets/1e4e48f9-7428-4b2f-a037-5e1f65da5b4e" align="right" width="500" alt="Graphical Abstract">
</a>

PMAT2 is a specialized tool for the **de novo assembly of mitochondrial and chloroplast genomes** from HiFi and ONT/CLR data. It supports animal, plant, and fungal organellar genome assembly and integrates a number of optimized algorithms for handling structural complexities.

If you encounter any problems using PMAT2, please contact the authors via email to join the WeChat group (include your *name + organization + PMAT2* in the message):

- **Changwei Bi**: bichwei@njfu.edu.cn  
- **Fuchuan Han**: hanfc@caf.ac.cn


## <a name="C9">Citation</a>
Fuchuan Han, Changwei Bi, Yicun Chen, Xiaogang Dai, Zefu Wang, Huaitong Wu, Ning Sun, et al. 2025. PMAT2: An Efficient Graphical Assembly Toolkit for Comprehensive Organellar Genomes. iMeta 4: e70064. https://doi.org/10.1002/imt2.70064 (if you use PMAT2 for organellar genomes)</br>
Bi C, Shen F, Han F, Qu Y, et al. PMAT: an efficient plant mitogenome assembly toolkit using ultra-low coverage HiFi sequencing data. Horticulture Research. (2024). uhae023, https://doi.org/10.1093/hr/uhae023 (if you use PMAT for plant genomes)


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
wget https://github.com/aiPGAB/PMAT2/archive/refs/tags/v2.1.5.tar.gz
tar -zxvf PMAT2-2.1.5.tar.gz
cd PMAT2-2.1.5
make
./PMAT --help
```

## <a name="C2">Requirement</a>

- [**BLASTn > 2.10.0**](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)  Needs to be installed in `PATH`.
- [**Singularity**](https://github.com/YanshuQu/runAssembly) or [**Apptainer**](https://github.com/apptainer/apptainer/blob/main/INSTALL.md) is required for PMAT2. You can find installation instructions [here](https://github.com/YanshuQu/runAssembly).
- [**Canu > v2.0**](https://github.com/marbl/canu) or [**NextDenovo**](https://github.com/Nextomics/NextDenovo) is required for CLR or ONT sequencing data.
- [**zlib**](https://www.zlib.net/) Needs to be installed in `PATH`.

## <a name="C3">Options and usage</a>

### <a name="C4">autoMito</a>
    
Run `PMAT autoMito --help` to view the usage guide.

```plaintext
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
   -G, --organelles     Genome organelles (mt/pt/all, default: mt)
   -x, --taxo           Specify the organism type (0/1/2), 0: plants, 1: animals, 2: Fungi (default: 0)
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
5. `-I`: minimum overlap identity (default: 90). If the assembly graph is complex, you can increase it appropriately.
6. `-L`: minimum overlap length (default: 40), if it is HiFi data, you can increase it appropriately.


### <a name="C5">graphBuild</a>

If PMAT fails to generate the assembly graph in 'autoMito' mode, you can use this command to manually select seeds for assembly.

Run `PMAT graphBuild --help` to view the usage guide.

```plaintext
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
   -x, --taxo           Specify the organism type (0/1/2), 0: plants, 1: animals, 2: Fungi (default: 0)
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

1. [Arabidopsis thaliana dataset (550Mb)](https://github.com/bichangwei/PMAT/releases/download/v1.1.0/Arabidopsis_thaliana_550Mb.fa.gz):
```sh
## download the dataset
wget https://github.com/bichangwei/PMAT/releases/download/v1.1.0/Arabidopsis_thaliana_550Mb.fa.gz

## run autoMito command
PMAT autoMito -i Arabidopsis_thaliana_550Mb.fa.gz -o ./test1 -t hifi -m

## run graphBuild command (when autoMito fails)
PMAT graphBuild -i ./test1/subsample/ -a ./test1/assembly_result/ -o ./test1_gfa -s 1 2 3 -d 5
```

The `PMAT_orgAss.txt` file contains the following information:
```plaintext
 ==========================================================
             Mitochondrial Assembly Assessment             
 ==========================================================

 Basic Statistics:
 ----------------------------------------------------------
 Total contigs:          16  
 Total length:           367.8 kb
 Average depth:          28.4 x
 Total genes found:      24/24 (100.0%)
 Duplicated contigs:     3   

 Per-contig Details:
 ----------------------------------------------------------
 Contig ID   Genes     Gene List           
 ----------------------------------------------------------
 300         4         atp1,cox1,nad1,nad2 
 2150        1         atp6                
 908         4         atp9,ccmB,cox2,nad9 
 1221        2         atp4,nad4L          
 729         4         ccmC,ccmFn,cox3,nad3
 727         1         nad3                
 1524        1         atp9                
 2150        1         atp6                
 749         6         atp8,matR,mttB,na...
 298         3         ccmFc,cob,nad6      
 ----------------------------------------------------------
```



**<a name="C6.2">Demo2</a>**

1. [Malus domestica dataset (540Mb)](https://github.com/bichangwei/PMAT/releases/download/v1.1.0/Malus_domestica.540Mb.fasta.gz):
```sh
## download the dataset
wget https://github.com/bichangwei/PMAT/releases/download/v1.1.0/Malus_domestica.540Mb.fasta.gz

## run autoMito command
PMAT autoMito -i Malus_domestica.540Mb.fasta.gz -o ./test2 -t hifi -m

## run graphBuild command (when autoMito fails)
PMAT graphBuild -i ./test2/subsample/ -a ./test2/assembly_result/ -o ./test2_gfa -s 10 20 30 -d 5
```

The `PMAT_orgAss.txt` file contains the following information:
```plaintext
 ==========================================================
             Mitochondrial Assembly Assessment             
 ==========================================================

 Basic Statistics:
 ----------------------------------------------------------
 Total contigs:          4   
 Total length:           397.0 kb
 Average depth:          31.1 x
 Total genes found:      24/24 (100.0%)
 Duplicated contigs:     1   

 Per-contig Details:
 ----------------------------------------------------------
 Contig ID   Genes     Gene List           
 ----------------------------------------------------------
 1           20        atp1,atp4,atp8,at...
 2           6         atp6,atp9,matR,na...
 ----------------------------------------------------------
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
├── gkmer_result/
|   ├── gkmer_histo.txt          # Kmer frequency
|   └── summary.txt              # genome size estimation
├── subsample/
│   └── PMAT_cut_seq.fa          # Subsampled reads for assembly
└── PMAT_orgAss.txt              # Organellar assembly assessment/
```

## <a name="C8">Version</a>

PMAT version 2.0.1 (24/11/21)</br>
Updates:

- Optimized the assembly strategy for organellar genomes, enabling faster and more accurate capture of organellar genome sequences.
- Implemented the assembly of animal and plant organellar genomes.
- Enhanced the genome graph untangling functionality for organellar genomes, enabling resolution of more complex structures.
- Parallelized key steps in the workflow, significantly improving runtime efficiency.

PMAT version 2.1.0 (25/2/1)</br>
Updates:

- Added `orgAss` module to evaluate the completeness of the assembly results.


