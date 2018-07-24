#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case);

my $usage = <<__EOUSAGE__;

############################################
#
# Required Options:
#
# --input <strint>
# --min <int>
# --max <int>
# --percent_trim <float>

# Example Usage:
#
# perl create_weighted_qual_file_for_celera.pl --input input.fa --min 10 --max 40 --percent_trim 0.1 --output output.txt
#
############################################

__EOUSAGE__

# Declare and initialize variables;
my $input;
my $min;
my $max;
my $percent_trim;
my $output;

my $options = GetOptions ('input=s' => \$input,
          'min=i' => \$min,
          'max=i' => \$max,
          'percent_trim=f' => \$percent_trim,
          'output=s' => \$output);

# validate options
unless ($options) {die $usage;}

my $diff = $max - $min;
my $index = -1;
my $names;
my $seqs;

open(FASTA, $input) || die("Can't open the fasta input file\n");
while(<FASTA>)
{
    $_ =~ s/\n//;
    if($_ =~ /^>/)
    {
        $index++;
        $names[$index] = $_;
    }
    if($_ !~ /^>/)
    {
        $seqs[$index] .= $_;
    }
}
close(FASTA);

$index = -1;

open(OUT,">$output") || die("Can't create the output file $output\n");
foreach $value (@seqs)
{
    $index++;
    $temp = "";
    $total_length = length($value);
    $start_section = $total_length*$percent_trim;
    $end_section = $total_length - ($total_length*$percent_trim);
    $fresh_min = $min;
    $fresh_max = $max;
    $end_count = 1;
        for($i=0; $i<length($value); $i++)
        {

        if($i <= $start_section)
        {
            $increase = int(($i/$start_section)*$diff);
            $fresh_min = $min + $increase;
            $temp .= "$fresh_min ";
        }

        if($i > $start_section and $i < $end_section)
        {
            $temp .= "$max ";
            $end_count = 0;
        }

        if($i >= $end_section)
        {
            $end_count++;
            $decrease = int(($end_count/($start_section-1))*$diff);
            $fresh_max = $max - $decrease;
            if($fresh_max < $min){ $fresh_max = $min; };
            $temp .= "$fresh_max ";
        }

        }


        $temp =~ s/ $//;
        print OUT "$names[$index]\n$temp\n";
}
close(OUT);
