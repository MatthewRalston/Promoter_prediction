Title: Promoter_prediction
Author: Matthew Ralston
Organization: Papoutsakis Organization

Promoter_prediction
Promoter_prediction is a project started on July 23, 2013 to describe a bioinformatic pipeline dedicated
to turning multiple sequence alignments (MSAs) of DNA motifs known to bind transcription factors. These
are frequently referred to as consensus sequences or transcription factor binding sites (TFBS). This project
began as part of an effort to identify antisense RNAs in C. acetobutylicum during my research in the summer
of 2013.
Some transcription factors (TFs) bind to a family of closely related DNA sequences. While it is difficult to identify
all of the sequences that a DNA-binding transcription factor may bind to, the sequences that are known to bind
transcription factors are recorded in the scientific literature and in scientific databases such as RegPrecise,
RegTransBase, and PRODORIC, among others. These databases commonly store multiple sequence alignments of these
transcription-factor binding DNA motifs such that highly conserved residues important for TF binding are located
at conserved positions in the MSA. These MSAs are rich in information content and can be used to find TF-binding motifs
in the genome.
Initially, this approach was used to scan UTR regions for cis-acting antisense RNA molecules. These molecules may play a role
in the kinetics and regulation of gene expression in gram positive bacteria and beyond. At the advice of researchers in the field,
the MAST tool from the MEME suite (Bailey et al. 1998 Bioinformatics, 14(1):48-54) was chosen for this project. MAST is a hidden
markov model (HMM) based statistical suite, capable of identifying motif matches in DNA, RNA, and protein sequences. MAST identifies
TFBS in DNA sequences and calculates expectation and probability values. These provide insight into the regulatory networks present in
microbial genomes of interest (e.g. C. acetobutylicum).
This project revolves around a central ruby script 'promoter.rb', which turns a MSA of sequences of uniform length into a position specific
probability matrix (PSPM) which is in turn used to identify TFBSs through the MAST tool. The output of this tool is then parsed with a p-value
cutoff. The script can identify single and bipartite motifs in DNA sequences within distance and p-value specifications.


Disclaimer: 

Future directions:
Countermeasure to remove sequences that have non-uniform length
Webtool
Script support for direct MSA-> bipartite motif conversion
Support for protein motif identification.