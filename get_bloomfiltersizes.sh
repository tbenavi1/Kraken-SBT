while read p; do kmc -k31 -m100 -fm @$p".descendantfilenames" NA.res . > temp; echo $p $(awk '/No. of unique k-mers/ {print $6}' temp) >> bloomfiltersizes; done < editednames
