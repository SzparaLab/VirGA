#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case);

my $usage = <<__EOUSAGE__;

############################################
#
# Required Options:
#
# --sam_file <string>              SAM file
# --read1 <string>                 read1 fastq file
# --read2 <string>                 read2 fastq file
# --genome <string>                genome index file
# --excise_start <int>             excise start location
# --excise_stop <int>              excise stop location
# --output_library_file <string>   output text file
# --output_fasta_file <string>     output fasta file
#
# Example Usage:
#
# perl excise_target_genomic_area.pl --sam_file genome.sam --read1 read_1.fastq --read2 read_2.fastq --genome genome.fa --excise_start excise_start --excise_stop excise_stop --output output.txt
#
############################################

__EOUSAGE__
    ;

# Declare and initialize variables;
my $sam_file;
my $read1;
my $read2;
my $genome;
my $excise_start;
my $excise_stop;
my $output;

my $options = GetOptions ('sam_file=s' => \$sam_file,
          'read1=s' => \$read1,
          'read2=s' => \$read2,
          'genome=s' => \$genome,
          'excise_start=i' => \$excise_start,
          'excise_stop=i' => \$excise_stop,
          'output=s' => \$output,
          'min_coverage=f' => \$min_coverage);

# validate options
unless ($options) {die $usage;}

my $counter = -1;
my $stops = -1;
my @inserts = "";
open(INPUT, $sam_file) || die("Can't open the file $sam_file\n");
while(<INPUT>)
{
    if($_ =~ /^\@/){ next; }
    @splitter = split(/\t/,$_);
    $flag = $splitter[1];
    if((field_2($flag) ? "YES" : "NO") eq "NO"){ next; }
    $match = $splitter[8];
    if($match < 0){ $match = $match*-1; }
    $counter++;
    $inserts[$counter] = $match;
        $stops++;
        if($stops == 100000)
        {
                $stops = 0;
        }
}
close(INPUT);

my $total = 0;
my $entries = 0;
foreach $ENTRY (@inserts)
{
    $entries++;
    $total += $ENTRY;
}

my $mean = $total / $entries;
my $variance = 0;
my $squared = 0;
foreach $ENTRY (@inserts)
{
    $squared = $mean - $ENTRY;
    $squared = $squared * $squared;
    $variance += $squared;
}

my $stddev = $variance / ($entries-1);
$stddev = sqrt($stddev);
my $int_mean = int($mean);
my $int_dev = int($stddev);

open(OUT,">$output_library_file") || die("Can't create the libraries output file $output\n");
print OUT "Reads bwa $read1 $read2 $int_mean $int_dev FR\n";
close(OUT);

open(GENOME, $genome]) || die("Can't open the file $genome\n");
while(<GENOME>)
{
    $_ =~ s/\n//;
    if($_ =~ /^>/)
    {
        $genomic_label = $_;
        next;
    }
    $genomic_seq .= $_;
}
close(GENOME);

$excise_start--;
my $excise_length = $excise_stop - $excise_start;
my $replacement_Ns = "N"x$excise_length;

substr($genomic_seq, $excise_start, $excise_length, $replacement_Ns);
$genomic_seq =~ s/^N+//;
$genomic_seq =~ s/N+$//;

open(OUT, ">output_fasta_file") || die("Can't open file $output_fasta_file for writing\n");
print OUT "$genomic_label\n";
print OUT "$genomic_seq\n";
close(OUT);
