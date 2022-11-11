use strict;
use Tie::Scalar::Timestamp;

my $path=shift;
my $groupName=shift;
my $augCfgPath=shift;
my $v=shift;

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

tie my $timestamp, 'Tie::Scalar::Timestamp';

print "{\n";
print " \"version\" : $v,\n";
print " \"timestamp\" : \"$timestamp\",\n";
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
    print "      \"genome\" : \"$path/$root\_Genome.fasta\",\n";
    print "      \"proteins\" : \"$path/$root\_Proteins.fasta\",\n";
    print "      \"gaf\" : \"$path/$root.gaf\",\n";
    if (-d "$augCfgPath/species/$root"){
    print "      \"augustus_model\": \"$augCfgPath/species/$root\",\n";
    }
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
