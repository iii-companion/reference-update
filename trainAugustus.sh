export AUGUSTUS_CONFIG_PATH=/home/whh2g/Bin/Augustus/config
/home/whh2g/Bin/Augustus/scripts/new_species.pl --species=$1
/home/whh2g/Bin/Augustus/scripts/gff2gbSmallDNA.pl $1.gff3 ../$1.fasta 1000 $1_sampled.gb
etraining --species=$1 $1_sampled.gb
