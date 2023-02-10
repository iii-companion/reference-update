use strict;
use JSON qw(decode_json);
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

my $rel = 1;
my $refs = "$path/Ref_$groupName/references.json";
if (-e $refs){
    open(my $refs_file, "<", $refs) or die "Can't open < references.json: $!";
    local $/;
    my $json_string = <$refs_file>;
    my $json_data = decode_json($json_string);
    if ($json_data->{version} == $v){
        $rel = $json_data->{release} + 1;
    }
    close($refs_file);
}

tie my $timestamp, 'Tie::Scalar::Timestamp';

print "{\n";
print " \"version\" : $v,\n";
print " \"release\" : $rel,\n";
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

    my $chrs = "$root\_chr.json";
    my $json;
    if (-e $chrs){
        {
        local $/;
        open (my $fh, "<", $chrs) or die $!;
        $json = <$fh>;
        close $fh;
        }
    }

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
    # get the chromosome stuff ok
    # if (defined($h{$root})){
	#     print "      \"chromosome_pattern\" : \"".$h{$root}."\",\n";
    # }
    if (defined($json)){
        print "      \"chr_mapping\": $json,\n";
    }
    print "      \"is_reference_strain\": true\n";
  

    print "     },\n";

}
$group =~ s/,$//g;
print "   },\n    \"groups\" : {\n";
print "    \"$groupName\" : [$group]\n";
print "},\n";
print "}\n\n";
