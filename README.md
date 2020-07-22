# KIPEs (Knowledge-based Identification of Pathway Enzymes)

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



## Usage

### General recommendation

Full paths should be used to specify input and output files and folders. Sequence names should not contain white space characters like spaces and TABs. Underscores can be used to replace spaces.


### Running the main function (KIPEs)

```
Usage:
  python KIPEs.py --baits <DIR> --out <DIR> --subject <FILE>
  or
  python KIPEs.py --baits <DIR> --out <DIR> --subjectdir <DIR>

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
   
   --simcut        FLOAT  Minimal BLAST hit similarity in percent [40.0]
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
```


`--baits` is the full path to a folder containing (mutliple) FASTA files. The filename needs to match the gene name. Extension should be '.fasta' or '.fa'.

`--out` is the full path to an output folder which will be created if necessary. All temporary and result files will be stored in this folder and subfolders therein.

`--subject` is the full path to an input multiple FASTA file. A collection of peptide (pep), transcript (rna), or genomic (dna) sequences can serve as input. The appropriate input data type needs to be specified via `--seqtype` (pep|rna|dna).

`--subjectdir` can be used to run KIPEs on multiple data sets (to analyse multiple species). All subject files in the provided folder are analysed consecutively. It is important that all data sets are of the same sequence type (`--seqtype` (pep|rna|dna)).

`--positions` is the full path to a folder containing text files matching the provided FASTA files. The filename needs to match the gene name. Example: CHS.fasta contains the bait sequences and CHS.txt contains information about relevant amino acid residues and domains. File extension should be '.txt' or '.res'. The header line starts with an exclamation mark followed by the reference sequence name. It is crucial that the name of this sequence is matched by one entry in the bait sequences FASTA file. Each of the following lines contains information about one important amino acid residue or a domain. The type of feature is indicated in the first column using R to specify residues or D to specify domains. The format of entries of residues and domains is slightly different as you can see in this example:

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


`--possibilities` specifies the number of different enzyme functions that are considered per sequence. This argument is relevant if looking for multiple different enzymes with high overall sequence similarity e.g. CHS/STS or F3H/FLS/LDOX. One sequence can be checked for the relevant amino acid residues of each of these enzymes to allow a high fidelity classification.



`--mafft` full path to MAFFT binary if this is not included in $PATH. MAFFT can be downloaded as a single binary here: https://mafft.cbrc.jp/alignment/software/ .

`--blastp` full path to blastp binary if this is not included in $PATH. BLAST can be downloaded from the [NCBI](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download).

`--tblastn` full path to tblastn binary if this is not included in $PATH. BLAST can be downloaded from the [NCBI](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download).

`--makeblastdb` full path to makeblastdb binary if this is not included in $PATH. BLAST can be downloaded from the [NCBI](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download).

`--fasttree` full path to FastTree binary. If this option is set, a phylogenetic tree is used to classify candidate sequences. FastTree can be downloaded as a single binary file here: http://www.microbesonline.org/fasttree/ .
  
  
  

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

As this tool is an automatic identification pipeline for candidate genes, the resolution of this analysis might be inferior to a manual annotation in certain cases like ODDs (F3H, FLS, LDOX). Therefore, it is recommended to carefully inspect the results from this pipeline. Generally, the quality of results is depending on the input quality. This pipeline was developed for the annotation of genes in the flavonoid biosynthesis of plants, but could be applied to other pathways if sufficient information is available.


## Data sets

baits.tar.gz recommended data sets for analysing the flavonoid biosynthesis in a new species.

residues.tar.gz recommended data sets for analysing the flavonoid biosynthesis in a new species.

FlavonoidBioSynBaits_v1.0.tar.gz was used in the proof of concept of [KIPEs](https://doi.org/10.1101/2020.06.27.175067) analysing the [transcriptome assembly of Croton tiglium](https://doi.org/10.3389/fmolb.2018.00062).

MYB_bHLH_WDR_v1.0.tar.gz was used in the proof of concept of [KIPEs](https://doi.org/10.1101/2020.06.27.175067) analysing the [transcriptome assembly of Croton tiglium](https://doi.org/10.3389/fmolb.2018.00062). This data set comprises transcription factors belonging to the three families MYB, bHLH, and WDR.
MYBs:
bHLHs:
WDRs:


## References

BLAST: https://doi.org/10.1016/S0022-2836(05)80360-2

MAFFT: https://doi.org/10.1093/molbev/mst010

FastTree2: https://doi.org/10.1371/journal.pone.0009490


When using KIPEs in your research, please cite:

Automatic identification of players in the flavonoid biosynthesis with application on the biomedicinal plant Croton tiglium. Boas Pucker, Franziska Reiher, Hanna Marie Schilbert. bioRxiv 2020.06.27.175067; doi: https://doi.org/10.1101/2020.06.27.175067
