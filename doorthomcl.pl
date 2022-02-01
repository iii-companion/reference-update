mkdir _all
cd _all
find .. -name 'proteins.fasta' | xargs cat > all_prots.fasta
find .. -name 'ggline.gg' | xargs cat > all_gg.txt
makeblastdb -dbtype prot -in all_prots.fasta
#formatdb -i all_prots.fasta
blastall -p blastp -m 8 -a 8 -e 1e-5 -d all_prots.fasta -i all_prots.fasta > b.out
PATH=/home/tdo/ORTHOMCLV1.4/:$PATH; export PATH
    
orthomcl.pl --mode 3 --blast_file b.out --gg_file all_gg.txt --inflation=1.5
