#Download all 6,741 (as of April 2017) Complete and Latest Bacteria Genomes from RefSeq
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt
awk -F "\t" '$12=="Complete Genome" && $11=="latest"{print $8,$20}' assembly_summary.txt > name_ftpdirpaths
awk -F "\t" '$12=="Complete Genome" && $11=="latest"{print $20}' assembly_summary.txt > ftpdirpaths
awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print ftpdir,file}' ftpdirpaths > ftpfilepaths
mkdir Bacteria_Genomes/
cd Bacteria_Genomes/
while read p; do wget $p; done < ./../ftpfilepaths

#Get kmer dumps files using jellyfish
ls *.fna.gz > filenames
while read p; do jellyfish count -t 20 -m 31 -s 100000000 -C <(zcat $p); jellyfish dump -c mer_counts.jf > "$p.dumps"; done < filenames
