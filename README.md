[![DOI](https://zenodo.org/badge/255121677.svg)](https://zenodo.org/badge/latestdoi/255121677)

# KIPEs (Knowledge-based Identification of Pathway Enzymes)

### KIPEs is available on our [webserver](https://pbb-tools.de/KIPEs/) ###

### Please get in touch if you need help running KIPEs on your own dataset: [Boas Pucker (email)](mailto:b.pucker@tu-braunschweig.de?subject=[GitHub]KIPEs_request) ###


## Abstract

This tool enables the identification of candidate sequences in a collection of peptide sequences, transcript sequences, or in a genome sequence. An initial BLAST (BLASTp, tBLASTn) search is used to get putative sequences which are than analysed in a global alignment with MAFFT and screened for the presence of conserved residues and conserved domains. As a proof of concept, this tool was applied for the identification of genes in the flavonoid biosynthesis.


## Input options

### Peptide sequence collection (result of an assembly annotation process)

If a collection of peptide sequences is provided, BLASTp is applied with a manually curated collection of bait sequences. Candidates are identified based on similarity cutoffs. Next, these sequences are subjected to a global alignment via MAFFT. The overall similarity of candidate and bait sequences is calculated. Additionally, the presence of conserved (functionally relevant) amino acid residues and conserved domains is inspected based on a reference sequence. If available, a peptide sequence collection should be provided instead of a transcript sequence set or a genome sequence. The computational costs of the analysis are substantially lower.


### Transcript sequences (transcriptome assembly)

If a collection of transcript sequences if provided, putative open reading frames are identified in the first step. All putative peptide sequences encoded in these transcripts are considered if their length exceeds a certain cutoff (e.g. 50 amino acids). The resulting peptides are subjected to the analysis described above.


### Genome sequence (genome assembly)

A tBLASTn is applied to identify regions in the genome sequence, which might encode the desired peptide. As BLAST hits only indicate exons and might be fragmented, BLAST hits are group to putative genes. Fragments of a putative gene are extended to account for incomplete hits at exon borders. This includes the detection of splice sites (currently only canonical GT-AG combinations). If full length peptide sequences are provided as query, the stop codon should be indicated by a * at the end of the peptide sequence. We recommend to run a proper gene prediction tool like [AUGUSTUS](https://github.com/Gaius-Augustus/Augustus) or [GeMoMa](http://www.jstacs.de/index.php/GeMoMa) if possible. These dedicated tools will outperform the very basic gene structure identification methods implemented in KIPEs in most cases.



## Result files
### HTML summary ###
The final output of KIPEs is a HTML document called 'SUMMARY.html'. This table shows the best candidates for all steps in the pathway. It is possible to specify the order of genes in a pathway using the `--pathway` option (see below for details). Previously described amino acid residues are checked in all candidate sequences and the results are summarized in this table. Mismatches of conserved residues are indicated by highlighting in red.

### Similarity matrix

One similarity matrix is generated per bait sequence file. The similarity of all candidate sequences against all bait sequences is displayed. Although this table is generated as a text file, it is possible to open these files as tables (e.g. with [Calc](https://www.libreoffice.org/discover/calc/)).


### Conserved residues

The presence of all conserved residues is analysed in all candidate sequences. Presence/absence are indicated in a table comprising all sequences and all residues. A summary of these results is diplayed in a HTML file as described above.


### Conserved regions
The output format of this analysis of conserved regions matches the output format of conserved residues. The percentage of identical amino acid residues in the domain is calculated for each candidate sequence.



## Installation

While some dependencies are required, this tool does not require an installation. Downloading and executing the script on a Linux system is sufficient. There is currently no support for other operating systems. Most required modules are included in the initial Python installation, but dendropy might not be available on all systems. 

[Python3](https://www.python.org/) (```sudo apt-get install python3.8```). It is also possible to use other Python3 versions.

[dendropy](https://dendropy.readthedocs.io/en/main/) (```sudo apt install python3-pip && python3 -m pip install git+https://github.com/jeetsukumaran/DendroPy.git```)

[MAFFT](https://mafft.cbrc.jp/alignment/software/linuxportable.html) (```sudo apt-get install -y mafft```)

[BLAST](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) (```sudo apt-get install ncbi-blast+```) or [HMMER](http://hmmer.org/documentation.html)(```conda install -c bioconda hmmer```)

[FastTree](http://www.microbesonline.org/fasttree/#Install) (```sudo apt-get install -y fasttree```) and/or [RAxML-NG](https://github.com/amkozlov/raxml-ng) (precompiled binaries recommended)





## Usage

### General recommendation

Full paths should be used to specify input and output files and folders. Sequence names should not contain white space characters like spaces and TABs. Underscores can be used to replace spaces.


### Running the main function (KIPEs)

```
Usage:
  python3 KIPEs3.py --baits <DIR> --out <DIR> --subject <FILE>
  or
  python3 KIPEs3.py --baits <DIR> --out <DIR> --subjectdir <DIR>

Mandatory:
  Bait sequences
  --baits          STR    Directory with (multiple) FASTA files
  
  Output directory
  --out            STR    Output directory

  Input sequences
  --subject        STR    Multiple FASTA file with sequences to screen
  --subjectdir     STR    Directory containing multiple FASTA file with sequences to screen
		
  Optional:
  --positions      STR    Directory with text files (one per step in pathway)
  --seqtype        STR    Defines type of input sequence (pep|rna|dna)[pep]
 
  --cpus           INT    Number of threads in BLAST runs [10]
  --scoreratio     FLOAT  BLAST score ratio of self vs. input sequences [0.3]
  --simcut         FLOAT  Minimal similarity of BLAST hits [40.0]
  --checks         STR    Validation of input data (on|off)[on]
   
  --genesize      INT    Maximal gene size (for tblastn hit grouping) [5000]
  --minsim        FLOAT  Minimal similarity required in global alignment [0.4]
  --minres        FLOAT  Minimal proportion of conserved residues [-1.0]
  --minreg        FLOAT  Minimal proportion of conserved regions [-1.0]
  --pathway       STR    Full path to text file with pathway enzyme names (default is alphabetical sorting)
  --possibilities INT    Maximal number of enzyme functions to consider per sequence [3]
   
  --mafft          STR    Full path to MAFFT (if not in your $PATH)
  --blastp         STR    Full path to the BLASTp binary (if not in your $PATH)
  --tblastn        STR    Full path to the tBLASTn binary (if not in your $PATH)
  --makeblastdb    STR    Full path to the makeblastdb binary (if not in your $PATH)
  
  --fasttree       STR    Full path to the FastTree binary
  
  --forester       STR    Activates the automatic construction of gene trees (on|off)[off]
  
  --exp            STR    Gene expression file (activates heatmap construction)
  --rcut           FLOAT  Minimal correlation cutoff [0.3]
  --pcut           FLOAT  Maximal p-value cutoff [0.05]
  --minexp         INT    Minimal expression per gene [30]
```


`--baits` is the full path to a folder containing (mutliple) FASTA files. The filename needs to match the gene name. Extension should be '.fasta' or '.fa'.

`--out` is the full path to an output folder which will be created if necessary. All temporary and result files will be stored in this folder and subfolders therein.

`--subject` is the full path to an input multiple FASTA file. A collection of peptide (pep), transcript (rna), or genomic (dna) sequences can serve as input. The appropriate input data type needs to be specified via `--seqtype` (pep|rna|dna).

`--subjectdir` can be used to run KIPEs on multiple data sets (to analyse multiple species). All subject files in the provided folder are analysed consecutively. It is important that all data sets are of the same sequence type (`--seqtype` (pep|rna|dna)).

`--positions` (or `--residues`) is the full path to a folder containing text files matching the provided FASTA files. The filename needs to match the gene name. Example: CHS.fasta contains the bait sequences and CHS.txt contains information about relevant amino acid residues and domains. File extension should be '.txt' or '.res'. The header line starts with an exclamation mark followed by the reference sequence name. It is crucial that the name of this sequence is matched by one entry in the bait sequences FASTA file. Each of the following lines contains information about one important amino acid residue or a domain. The type of feature is indicated in the first column using R to specify residues or D to specify domains. The format of entries of residues and domains is slightly different as you can see in this example:

!AtCHS \
R	R	13	comment1 \
R	Q,X	16	comment2 \
R	R	17	comment3 \
D	malonyl-CoA_binding_motif	313	329	comment4

Residues: Important residues have their amino acid in the second column (one letter code!) and the position in the third column. It is possible to specify multiple alternative amino acids for one position as indicated by the 'X' in the second entry. Columns following the third column can be used for user comments and are ignored by KIPEs.

Domains: The domain entry indicator (D) is followed by the name of the domain in the first column. This name must not contain any TABs. The third column contains the start position of the domain and the fourth column contains the end position. All following columns are ignored by KIPEs and can be used for user comments.

`--seqtype` specifies the input data type as peptide (pep), transcript (rna), or genomic (dna) sequences.

`--cpus` specifies the number of threads to use for BLAST.

`--scoreratio` specifies the minimal score ratio between BLAST hits against the subject sequences and the bait sequence itself. The value range is 0.0 to 1.0 with default at 0.3.

`--simcut` specifies the minimal similarity of BLAST hits against the subject to be considered. The value range is 0 to 100 with default at 40.

`--checks` activates (on) or deactivates (off) the validation of input data. Sequence names provided in the conserved positions files are checked against the sequence names in the bait sequence files. 

`--genesize` specifies the maximal distance between tblastn hits to be considered as parts of the same gene.

`--minsim` specifies the minimal similariy required in the global alignment. A low value leads to a high sensitivity, but causes the analysis of more sequences thus increasing the run time.

`--minres` specifies the minimal proportion (0.0 - 1.0) of important amino acid residues that need to be conserved. This filtering option should only be used to remove low quality hits. However, manual inspection of top hits with just a few amino acid substitutions is recommended. Sequences with amino acid substitutions might have lost the ancestral function of an enzyme, but could have gained a new one.

`--minreg` specifies the minimal proportion (0.0 - 1.0) of domains that need to be conserved. This filtering option should only be used to remove low quality hits. However, manual inspection of top hits with just a few amino acid substitutions is recommended. Sequences with amino acid substitutions in these domains might have lost the ancestral function of an enzyme, but could have gained a new one.

`--pathway` can be used to provide the order of all steps in the pathway as a text file. One enzyme/gene name needs to be given per line and the names need to match the names of provided data files (bait sequence and residue/domain info files) excactly.If no information about the pathway genes/enzymes is provided, results will be sorted in alphabetical order. Enzyme/gene names missing from the pathway file will be appended to the provided list to avoid the loss of results.


`--possibilities` specifies the number of different enzyme functions that are considered per sequence. This argument is relevant if looking for multiple different enzymes with high overall sequence similarity e.g. CHS/STS or F3H/FLS/ANS(LDOX). One sequence can be checked for the relevant amino acid residues of each of these enzymes to allow a high fidelity classification.



`--mafft` full path to MAFFT binary if this is not included in $PATH. MAFFT can be downloaded as a single binary here: https://mafft.cbrc.jp/alignment/software/ .

`--blastp` full path to blastp binary if this is not included in $PATH. BLAST can be downloaded from the [NCBI](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download).

`--tblastn` full path to tblastn binary if this is not included in $PATH. BLAST can be downloaded from the [NCBI](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download).

`--makeblastdb` full path to makeblastdb binary if this is not included in $PATH. BLAST can be downloaded from the [NCBI](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download).

`--fasttree` full path to FastTree binary. If this option is set, a phylogenetic tree is used to classify candidate sequences. FastTree can be downloaded as a single binary file here: http://www.microbesonline.org/fasttree/ .

`--forester` activates the automatic construction of gene tree with all identified candidates. One tree per step in the pathway is constructed with FastTree. The construction of gene trees requires that the `--fasttree`  is used. Otherwise, FastTree is expected in the $PATH. Default: off.

`--exp` specifies gene expression file. The first column contains the gene IDs which need to match the sequence IDs of the dataset analyzed by KIPEs. The first row contains the names of samples. Providing a gene expression file activates the construction of a co-expression heatmap with all candidate genes. Default: off.

`--rcut` specifies the minimal correlation cutoff value. Only genes with at least this strength of correlation are reported. Default: 0.3.

`--pcut` specifies the maximal p-value cutoff. Only genes with a correlation that has a lower p-value are reported. Default: 0.05.

`--minexp` specifies the minimal expression per gene (across all samples) that is required for a gene to be considered in the co-expression analysis. Default: 30.
  

### Generating tables of conserved residues

The generation of required input data (conserved residues) can be performed based on a multiple FASTA file which contains previously characterized sequences of the same gene in multiple species.


```
Usage:
  python get_cons_pos.py --in <FILE> --out <DIR> --ref <STRING> --name <STRING>

Mandatory:
  --in STR          A multiple FASTA file. 
  --out STR         Directory for temporary and output files.
  --ref STR         Name of the reference sequence.
  --name STR        Name of the output files.
		
Optional:
  --mincons FLOAT  Minimal conservation frequency.[1.0]
  --mafft STR      Full path to MAFFT (if not in your $PATH)
```


  

### Identification of residues for classification

In order to differentiate between enzymes with high sequence similarity (e.g. CHS, STS, and LAP5), known sequences of each type can be collected in a multiple FASTA file. This scripts allows the identification of residues which differentiate sequences of the different types. 

```
Usage:
  python find_contrasting_residues.py --in <FILE> --out <DIR> --ref <STRING> --pos <STRING> --neg <STRING>
  or
  python find_contrasting_residues.py --in1 <FILE> --in2 --out <DIR> --ref <STRING>
  
Arguments:
  --in  STR         A multiple FASTA file.
  --in1 STR         A multiple FASTA file with white list sequences.
  --in2 STR         A multiple FASTA file with black list sequences.
  --out STR         Directory for temporary and output files.
  --ref STR         Name of the reference sequence.
  --pos STR         Comma-separated names of sequences (white list).
  --neg STR         Comma-separated names of sequences (black list).
```



## Limitations and outlook

As this tool is an automatic identification pipeline for candidate genes, the resolution of this analysis might be inferior to a manual annotation in certain cases like ODDs (F3H, FLS, ANS/LDOX). Therefore, it is recommended to carefully inspect the results from this pipeline. Generally, the quality of results is depending on the input quality. This pipeline was developed for the annotation of genes in the flavonoid biosynthesis of plants, but could be applied to other pathways if sufficient information is available.


### Automatic construction of gene trees ###

Once KIPEs identified candidate sequences for all the steps in a pathway or all members of a gene family, users might be interested to see these sequences in a phylogenetic tree. Since it is time consuming to construct single gene trees for all steps of the flavonoid biosynthesis, we added a script to do this automatically.


```
Usage:
  python forester.py --in <DIR> --out <DIR> --ref <DIR>
 
Mandatory:
  --in          STR    Directory with (multiple) FASTA files (KIPEs results folder 'final_pep_files')
  --out         STR    Output directory
  --ref         STR    Directory containing multiple FASTA file with reference sequences (e.g. initial baits)
		
  Optional:
  --occ         FLOAT  Minimal occupancy in alignment columns
  --clean       STR    Activates the removal of all intermediate files
  --mafft       STR    Full path to MAFFT (if not in your $PATH)
  --fasttree    STR    Full path to the FastTree binary (if not in your $PATH)
```


`--in` is the full path to the folder 'final_pep_files' in the KIPEs output. This folder contains FASTA files which are used for the construction of gene trees.

`--out` is the output folder. All FASTA files and trees will be stored in this folder. If the folder does not exist already, it will be created.

`--ref` is a full path with FASTA files that serve as reference in the gene tree. This could be the folder with initial bait sequences used for KIPEs.

`--occ` defines the minimal alignment column occupancy. Columns with more gaps than specified by occ, will be removed prior to the gene tree construction.

`--clean` activates the removal of all intermediate files. This option is off by default to allow inspection of temporary files and to avoid removal of folders by accident.

`--mafft` specifies the MAFFT path if this is not globally available.

`--fasttree` specifies the FastTree path if this is not globally available.



### Collect KIPEs results ###

When running KIPEs on a large number of data sets, it can be helpful to collect the best candidates in one table. This scripts allows filtering of the candidates and produces a summary table of all analyzed species:


```
Usage:
  python3 collect_KIPEs_results.py --in <DIR> --out <DIR> --genes <STR>
 
Mandatory:
  --in          STR    Directory with KIPEs result folders
  --out         STR    Output file
  --genes       STR    Comma-separated list of genes
		
  Optional:
  --mapping     STR    Folder name to species name mapping
  --minsim      FLOAT  Minimal similarity
  --minmatches  FLOAT  Minimal residue match ratio
```


`--in` is the full path to the folder that contains the KIPEs result folder of the analyses of many species.

`--out` is the full path to the output file that will contain the summarized information about the best candidates in all analyzed species.

`--genes` is a list of all gene names that should be listed in the final result table. The names of genes are comma-separated e.g. "CHS,CHI1,F3H".

`--mapping` is the full path to a table that contains the name of folders in the directory provided via `--in` in the first column and the corresponding species name in the second column.

`--minsiim` minimal similarity of a candidate to the bait sequences.[0.5]

`--minmatches` minimal ratio of detected conserved amino acid residues.[0.9]



## Data sets

**baits.tar.gz** recommended data sets for analysing the flavonoid biosynthesis in a new species.

**residues.tar.gz** recommended data sets for analysing the flavonoid biosynthesis in a new species.

**FlavonoidBioSynBaits_v1.0.tar.gz** was used in the proof of concept of [KIPEs](https://doi.org/10.3390/plants9091103) analysing the [transcriptome assembly of *Croton tiglium*](https://doi.org/10.3389/fmolb.2018.00062).

**MYB_bHLH_WDR_v1.0.tar.gz** was used in the proof of concept of [KIPEs](https://doi.org/10.3390/plants9091103) analysing the [transcriptome assembly of *Croton tiglium*](https://doi.org/10.3389/fmolb.2018.00062). This data set comprises transcription factors belonging to the three families MYB, bHLH, and WDR.
MYBs: *Arabidopsis thaliana* [Stracke et al., 2001](https://doi.org/10.1016/S1369-5266(00)00199-0), *Vitis vinifera* [Matus et al., 2008](https://bmcplantbiol.biomedcentral.com/articles/10.1186/1471-2229-8-83), *Beta vulgaris* [Stracke et al., 2014](https://bmcplantbiol.biomedcentral.com/articles/10.1186/s12870-014-0249-8), and *Musa acuminata* [Pucker et al., 2020](https://doi.org/10.1101/2020.02.03.932046); 
bHLHs: *Arabidopsis thaliana* [Heim et al., 2003](https://academic.oup.com/mbe/article/20/5/735/985548), *Vitis vinifera* [Wang et al., 2018](https://www.frontiersin.org/articles/10.3389/fpls.2018.00064/full), *Nelumbo nucifera* [Mao et al., 2019](https://peerj.com/articles/7153/), *Citrus grandis* [Zhang et al., 2020](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-020-6644-7), *Musa acuminata* [Wang et al., 2020](https://doi.org/10.3389/fpls.2020.00650), and *Solanum melongena* [Tian et al., 2019](https://doi.org/10.7717/peerj.7768);
WDRs: *Arabidopsis thaliana* [Nocker et al., 2003](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC317288/), *Triticum aestivum* [Hu et al., 2018](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-5157-0), and *Setaria italica* [Mishra et al., 2014](https://doi.org/10.1371/journal.pone.0086852)


#### MYB ###

Although KIPEs can be used to identify MYBs, it is recommend to use the **[MYB_annotator](https://github.com/bpucker/MYB_annotator)** instead. This dedicated tool allows a more reliable and computationally optimized identification of MYBs. Functional annotations are also added.

<details>
<summary>(click to expand)</summary>

(v1.0) *Arabidopsis thaliana* [Stracke et al., 2001](https://doi.org/10.1016/S1369-5266(00)00199-0):The R2R3-MYB gene family in Arabidopsis thaliana

(v1.0) *Vitis vinifera* [Matus et al., 2008](https://bmcplantbiol.biomedcentral.com/articles/10.1186/1471-2229-8-83): Analysis of the grape MYB R2R3 subfamily reveals expanded wine quality-related clades and conserved gene structure organization across Vitis and Arabidopsis genomes

(v1.0) *Beta vulgaris* [Stracke et al., 2014](https://bmcplantbiol.biomedcentral.com/articles/10.1186/s12870-014-0249-8): Genome-wide identification and characterisation of R2R3-MYB genes in sugar beet (Beta vulgaris)

(v1.0) *Musa acuminata* [Pucker et al., 2020](https://doi.org/10.1371/journal.pone.0239275): The R2R3-MYB gene family in banana (Musa acuminata): genome-wide identification, classification and expression patterns

(v1.1) *Croton tiglium* [Pucker et al., 2020](https://doi.org/10.3390/plants9091103): Automatic identification of players in the flavonoid biosynthesis with application on the biomedicinal plant Croton tiglium

(v1.1) *Marchantia polymorpha* (Pucker, 2021): Genome-wide identification of *Marchantia polymorpha* MYBs via KIPEs

This dataset is a trade-off between covering a broad phylogenetic width and small size (=short run time). Please have a look at the [wiki](https://github.com/bpucker/KIPEs/wiki/Datasets) for more information about additional MYB sequence datasets.
  
  </p>
</details>


## Interested in another pathway?
KIPEs is not restricted to analyzing the flavonoid biosynthesis. This README explains how to add a new pathway to KIPEs.

IMPORTANT NOTICE: Please get in touch before constructing your own input files. We are already working on the integration of additional pathways. You might be able to save time and energy.

### 1) Find a paper describing the pathway/reactions of interest and at least one characterized sequence
If you do not know all steps involved in the pathway of interest, you need to identify publications describing it. This is important to include all known steps in the following analysis. It is not possible to discover novel pathways using KIPEs. You need to identify at least a single characterized sequence for each step in the pathway - more are better.

### 2) Run tree-based selection of additional sequences
KIPEs is based on information already available in the literature and in databases. Therefore, it is crucial to identify well characterized bait sequences and conserved amino acid residues in the active center. You can complement this dataset by looking for orhtologs. These sequences require a manual inspection. Yang et al., 2015 described a pipeline for the identification of orthologs based on an initial sequence. This approach can also be used to prepare a collection of sequences for exploration if conserved amino acid residues are not known yet.

### 3) Construct FASTA file of collected sequences
The final result of this preliminary analysis is a FASTA file containing bait sequences for KIPEs. The quality of this sequence collection determines the quality of the KIPEs results.

### 4) Construct file with conserved amino acid residues
If possible, construct a file describing the positions of conserved amino acid residues. Relying on highly conserved residues might be an option if no additional information is available. However, the quality of this list of residues is decisive. See the descriptions above for details about the file structure.




## References

BLAST: https://doi.org/10.1016/S0022-2836(05)80360-2

MAFFT: https://doi.org/10.1093/molbev/mst010

FastTree2: https://doi.org/10.1371/journal.pone.0009490

dendropy: https://dendropy.org/ and https://doi.org/10.1093/bioinformatics/btq228

When using KIPEs in your research, please cite:

Pucker B., Reiher F., and Schilbert H.M. (2020). Automatic identification of players in the flavonoid biosynthesis with application on the biomedicinal plant _Croton tiglium_. Plants 2020, 9, 1103. doi:[10.3390/plants9091103](https://doi.org/10.3390/plants9091103).

Rempel A. & Pucker B. (2022). KIPEs3: Automatic annotation of biosynthesis pathways. bioRxiv 2022.06.30.498365; doi:[10.1101/2022.06.30.498365](https://doi.org/10.1101/2022.06.30.498365)
