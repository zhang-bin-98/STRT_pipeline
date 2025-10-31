#!/usr/bin/perl

# ===============================================================
# File: matrix_to_seurat_rds.R
# File Created: Th Sep 2023
# Author: Zhang Bin
# -----
# Last Modified: Fri Sep 08 2023
# Modified By: Zhang Bin
# ===============================================================


use v5.010;
# use warnings;
use Getopt::Long;
use threads ('exit' => 'threads_only');
use List::Util qw/sum/;

my $usage = qq/
$0: count umi number for BAM

perl $0 [options] 
    [--help|-h]                             : show usage 
    [--work|-w]     threads_number(5)
    [--samtools|-s] samtools path(samtools)
    --bam|-b        bam_file                : HTseq or featureCounts output BAM which has "X[F|T]:Z" tag
    --out|-o        out_put_file            : create new file to store the result tsv

eg:
    perl $0 -w 5 -b featureCounts_out.bam -o count_matrix.tsv

/;

################################

my $help;
my $works;
my $samtools;
my $file_bam;
my $out;

GetOptions(
    "help|h!" => \$help,
    "work|w:i" => \$works,
    "samtools|s:s" => \$samtools,
    "bam|b=s" => \$file_bam,
    "out|o=s" => \$out
) or die $!;

die $usage if defined($help);
die "ERROR: can not find bam: ${file_bam}\n" if !( -e $file_bam );
die "WRONG: output ${out} already exists\n" if ( -e $out);
if(defined($samtools)) {
    die "ERROR: can not find samtools: ${samtools}\n" if !( -e $samtools );
} else {
    $samtools = "samtools";
}
$works = 5 if (!defined($works));

# main
say "Start!\n".localtime();

my $summary = process_reads($samtools, $file_bam, $works);
output_summary($out, $summary);

say "End!\n".localtime();

################################

sub count_umi {
    my ($samtools, $file_bam, $p, $works) = @_;

    my %summary;
    my $read_number;
    my $pid = threads->tid();
    say "NOTE: \tthread $pid start!";

    open(my $fh_bam, "$samtools view $file_bam |") or die $!;
    while(<$fh_bam>) {
        next if($. % $works != $p);
        my ($CB, $CU, $gstat) = $_ =~ /\w+_([ATCG]{8})_([ATCG]{8})\t.*?XS:Z:(\w+)/;
        #
        if ($gstat eq "Assigned") {
            my ($gene) = $_ =~ /X[F|T]:Z:(\w+)/;
            $summary{$gene}{$CB}{$CU}++;
        }
        $read_number++;
        # say "NOTE: \tthread $pid process reads: $read_number" if ($read_number % 100000 == 0);
    }
    close($fh_bam);

    say "NOTE: \tthread $pid finish!\tprocess reads: $read_number";
    return \%summary;
}

sub process_reads{
    my ($samtools, $file_bam, $works) = @_;
    say "NOTE: process reads!\n\tsamtools: $samtools\n\tBAM: $file_bam";

    my @threads_list = map { 
        threads->create(
            {'context' => 'list'},
            \&count_umi, $samtools, $file_bam, $_, $works
        );
    } 0..($works-1);

    my @summary = map { $_ -> join() } @threads_list;

    return \@summary;
}

################################

sub output_summary {
    my ($out, $sum_array) = @_;

    say "NOTE: build count matrix!\n\tout put: $out";

    my %genes;
    my %cells;
    my %summary;

    # merge data
    foreach $p (@$sum_array) {
        foreach $g (keys %$p) {
            $genes{$g}++;
            foreach $c ( keys %{%{$p}{$g}} ) {
                $cells{$c}++;
                foreach $u (keys %{%{%{$p}{$g}}{$c}}) {
                    $summary{$g}{$c}{$u}++;
                }
            }
        }
    }

    # gene
    my @genes = sort keys %genes;
    undef %genes;
    # cell
    my @cells = sort keys %cells;
    undef %cells;
    # matrix
    my %matrix;
    foreach $g (keys %summary) {
        foreach $c (keys %{%summary{$g}}) {
            $matrix{$g}{$c} = scalar keys %{%{%summary{$g}}{$c}};
        }
    }
    undef %summary;

    # output
    open(my $fh_out, ">$out") or die $!;
    say {$fh_out} join("\t", ("GeneID", @cells));
    map {
        say {$fh_out} join("\t", ($_, @{$matrix{$_}}{@cells}));
    } @genes;
    close($fh_out);
}



