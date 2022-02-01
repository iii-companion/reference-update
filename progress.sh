#/usr/bin/env bash

gene=$(tail -1 b.out | awk '{print $1}')
line=$(grep -n $gene all_prots.fasta | cut -f1 -d:)
total_lines=$(wc -l < all_prots.fasta)

perc=$(bc <<< "scale=2; 100*$line/$total_lines")
echo "${perc}%"

