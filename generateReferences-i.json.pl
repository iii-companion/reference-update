use strict;

my $path=shift; #curennt path, I am lazy
my $groupName=shift;
# for x in *gff3 ; do  n=$(egrep ^##sequence-region.*_01 clean_gff/$x |  perl -nle 'if(/(\S+)_01_(v\d+)/){ print "$1\t$2"} elsif(/(\S+)_01/){ print $1 } ' ); echo -e "$x\t$n"; done > ChromosomeFile.txt

my %h;
open F, "ChromosomeFile.txt" or die " wrong DIR:$!\n";
while(<F>){
    chomp;
    my @ar=split(/\t/);
    if (defined($ar[1])){
	$h{$ar[0]}=$ar[1]."\$";
    }
}
close(F);

print "{\n";
print " \"species\" : {\n";
my $group;
while(<STDIN>){
    # expect the file name from gff3
    chomp;
    my ($n)=$_;
    my ($name)=$n;
    $name=~s/_/ /g;
    $name=~ s/.gff3//g;
    
    my $root=$n; $root=~s/\.gff3//g;
    $group.="\"".$root."\",";
    print "    \"$root\" : {\n";
    print "      \"gff\" : \"$path/clean_gff/$_\",\n";
    print "      \"genome\" : \"$path/$root.fasta\",\n";
    print "      \"gaf\" : \"$path/$root.gaf\",\n";
    print "      \"augustus_model\": \"/home/whh2g/Bin/Augustus/config/species/$root\",\n";
     if (defined($h{$root})){
	 print "      \"name\" : \"$name - in chromosomes\",\n";
     } else {
	 print "      \"name\" : \"$name\",\n";
     }
    ## get the chromosome stuff ok
    if (defined($h{$root})){
	print "      \"chromosome_pattern\" : \"".$h{$root}."\",\n";
    }
    print "      \"is_reference_strain\": true\n";
  

    print "     },\n";

}
$group =~ s/,$//g;
print "   },\n    \"groups\" : {\n";
print "    \"$groupName\" : [$group]\n";
print "},\n";
print "}\n\n";
