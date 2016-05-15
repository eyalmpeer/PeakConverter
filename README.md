# PeakConverter
A script for converting genomic-coordinates to transcriptomic-coordinates in MeRIP-Seq analysis.

PeakConverter takes MACS2 narrowPeak output as input, along with Cufflinks isoforms_fpkm output and a table file, and outputs MACS2 peaks in transcriptomic coordinates, based on the most expressed isoform for each gene.
In case of equally expressed isoforms, preference is given first to coding isoforms, and then to longer isoforms.
Peaks that are 50 nt or shorter after conversion to transcriptomic coordinates are discarded, as well as peaks mapping to isoforms which are duplicated in two or more places in the genome.
PeakConverter outputs two files, name_tx.bed and name_exonpeaks.bed.
name_tx.bed is the transcriptomic-coordinate based peak file, which also includes transcriptomic features for transcript reported such as placement of the canonical AUG, stop codon, transcript length, location of the first splice site and location of the last splice site.
name_exonpeaks.bed is a peak file in genomic coordinates in which peaks were mapped to the most expressed isoform for each gene and only the peak segment overlapping exons is reported. For a peak to be reported at least 50% of it must overlap an exon. This file is generated for usage by the trough calling scripts.

PeakConverter is used as such:

python PeakConverter.py --bed-file file.bed --expression-file isoforms.tracking_fpkm --table-file file.table --output-prefix name

  --bed-file    MACS2 output narrowPeak file or any other bed file of genomic coordinates with an ID/name 4th column
  --expression-file   Cufflinks output of isoforms fpkm file.
  --table-file    The chosen annotation table file, downloaded from the UCSC table browser. Must be of the same annotation                   used with Cufflinks.
  --output-prefix The name prefix for each PeakConverter output file.
