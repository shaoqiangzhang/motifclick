# motifclick
A motif finding tool for a set of intergenic sequences


You can also directly use command 'g++ MotifClick.cpp -o MotifClick' to complie the cpp file.

******
USAGE:
******
./MotifClick input_file [OPTIONS]  > OutputFile

input_file:        file containing DNA sequences in FASTA format

OPTIONS:
-w       :       motif width(default=10)
-n        :      maximum number of motifs to find(default=5)
-SSD      :      upper bound of SSD (sum of squared distance, please select 0<SSD<1, default=0.3)
-b 2       :     if examine sites on both of DNA strands(default=1 only forward)
-d         :     upper bound of graph density(default=100)
-s 0       :     if want more degenerate sites (default=1 if want fewer sites)

*******
If  you use the code, please cite the reference:

Zhang et al. MotifClick: prediction of cis-regulatory binding sites via merging cliques BMC Bioinformatics 2011, 12:238 (https://doi.org/10.1186/1471-2105-12-238)
