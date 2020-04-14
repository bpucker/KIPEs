# CGI (Candidate Gene Identifier)

## Abstract

This tool enables the identification of candidate sequences in a collection of peptide sequences, transcript sequences, or in a genome sequence. An initial BLAST (BLASTp, tBLASTn) search is used to get putative sequences which are than analysed in a global alignment with MAFFT and screened for the presence of conserved residues and conserved domains. As a proof of concept, this tool was applied for the identification of genes in the flavonoid biosynthesis.


## Peptide sequence collection

If a collection of peptide sequences is provided, BLASTp is applied with a manually curated collection of bait sequences. Candidates are identified based on similarity cutoffs. Next, these sequences are subjected to a global alignment via MAFFT. The overall similarity of candidate and bait sequences is calculated. Additionally, the presence of conserved (functionally relevant) amino acid residues and conserved domains is inspected based on a reference sequence. 


## Transcript sequences (transcriptome assembly)

If a collection of transcript sequences if provided, putative open reading frames are identified in the first step. All putative peptide sequences encoded in these transcripts are considered if their length exceeds a certain cutoff (e.g. 50 amino acids). The resulting peptides are subjected to the analysis described above.


