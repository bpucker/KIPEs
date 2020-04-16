# CGI (Candidate Gene Identificatory)

## Abstract

This tool enables the identification of candidate sequences in a collection of peptide sequences, transcript sequences, or in a genome sequence. An initial BLAST (BLASTp, tBLASTn) search is used to get putative sequences which are than analysed in a global alignment with MAFFT and screened for the presence of conserved residues and conserved domains. As a proof of concept, this tool was applied for the identification of genes in the flavonoid biosynthesis.

## Input options

### Peptide sequence collection (result of an assembly annotation process)

If a collection of peptide sequences is provided, BLASTp is applied with a manually curated collection of bait sequences. Candidates are identified based on similarity cutoffs. Next, these sequences are subjected to a global alignment via MAFFT. The overall similarity of candidate and bait sequences is calculated. Additionally, the presence of conserved (functionally relevant) amino acid residues and conserved domains is inspected based on a reference sequence. If available, a peptide sequence collection should be provided instead of a transcript sequence set or a genome sequence. The computational costs of the analysis are substantially lower.


### Transcript sequences (transcriptome assembly)

If a collection of transcript sequences if provided, putative open reading frames are identified in the first step. All putative peptide sequences encoded in these transcripts are considered if their length exceeds a certain cutoff (e.g. 50 amino acids). The resulting peptides are subjected to the analysis described above.


### Genome sequence (genome assembly)

A tBLASTn is applied to identify regions in the genome sequence, which might encode the desired peptide. As BLAST hits only indicate exons and might be fragmented, BLAST hits are group to putative genes. Fragments of a putative gene are extended to account for incomplete hits at exon borders. This includes the detection of splice sites (currently only canonical GT-AG combinations).

## Result files

### Similarity matrix

One similarity matrix is generated per bait sequence file. The similarity of all candidates against all baits is displayed.


### Conserved residues

The presence of all conserved residues is analysed in all candidate sequences. Presence/absence are indicated in a table comprising all sequences and all residues.


### Conserved regions
The output format of this analysis of conserved regions matches the output format of conserved residues. The percentage of identical amino acid residues in the domain is calculated for each candidate sequence.



## Usage

### General recommendation

Full paths should be used to specify input and output files and folders. Sequence names should not contain white space characters like spaces and TABs. Underscores can be used to replace spaces.


### Running the main function (CGI)

```
Usage:
  python CGI.py --baits <DIR> --out <DIR> --subject <FILE>

Mandatory:
  Bait sequences
  --baits STR        Directory with (multiple) FASTA files. 
  
  Output directory
  --out STR          Output directory

  Input sequences
  --subject STR      Multiple FASTA file with sequences to screen.
		
Optional:
  --positions STR    Directory with text files (one per gene).
  --seqtype STR      Defines type of input sequence (pep|rna|dna)[pep]
  --mafft STR        Full path to MAFFT (if not in your $PATH)
  --cpus INT         Number of threads in BLAST runs [10]
  --scoreratio FLOAT BLAST score ratio of self vs. input sequences [0.3]
  --simcut FLOAT     Minimal similarity of BLAST hits [40.0]
  --blastp <PATH_TO_AND_INCLUDING_BINARY>[blastp]
  --tblastn <PATH_TO_AND_INCLUDING_BINARY>[tblastn]
  --makeblastdb <PATH_TO_AND_INCLUDING_BINARY>[makeblastdb]
```

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



## References


