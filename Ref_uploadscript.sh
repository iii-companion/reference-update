

#https://github.com/sanger-pathogens/companion/wiki/Preparing-reference-data-sets

name=$1

rename "s/ /_/g"  * # get rid of spaces in names / also do it clean_gff

mkdir clean_gff
for x in *gff* ; do
    gt gff3 -sort -retainids -tidy $x > clean_gff/$x &
done




# mass augustus training

export AUGUSTUS_CONFIG_PATH=/home/whh2g/Bin/augustus-3.4.0/config
#In clean_gff
cd clean_gff
for x in *gff3 ; do n=$( echo $x | sed 's/.gff3//g'); ~/Bin/trainAugustus.sh $n & done
cd ..

### that will train all of them…

for x in *gff3 ; do
    m=$(echo $x | sed 's/.gff3//g');
    n=$(egrep '^##sequence-region.*_01' clean_gff/$x | perl -nle 'if(/(\S+)_01_(v\d+)/){ print "$1\t$2"} elsif(/(\S+)_01/){ print $1 } ' );
    echo -e "$m\t$n";
done > ChromosomeFile.txt # generate Chromosome file


mkdir Reference
ls *gff3 | perl ~/Bin/generateReferences-i.json.pl $PWD  $name > Reference/references-in.json; 

cd Reference
/home/companion/annot-nf/bin/update_references.lua
~/Bin/doorthomcl.pl

echo " UPdate references-in.json with chromosome names" 
