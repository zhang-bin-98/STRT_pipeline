#!/usr/bin/perl

# =======================================
# File: STRT_cell_mapping_ratio.pl
# File Created: 2023-09-11 12:16:05
# Author: Zhang Bin
# -----
# Last Modified: 2023-09-11 10:21:35
# Modified By: Zhang Bin
# =======================================

use v5.010;
# use warnings;
use Getopt::Long;
use threads ('exit' => 'threads_only');
use List::Util qw/sum/;

my $usage = qq/
$0: count mapping ratio 

perl $0 [options] 
    [--help|-h]                             : show usage 
    [--work|-w]     threads_number(5)
    [--samtools|-s] samtools path(samtools)
    --read1|-r1     input_R1.fq.gz
    --bam|-b        bam_file
    --out|-o        out_put_file            : create new file to store the result csv

eg:
    perl $0 -w 5 -r1 qc_R1.gq.gz -b mapping_out.bam -o barcode_reads_mapping_ratio.csv

/;

################################

my $help;
my $works;
my $samtools;
my $file_r1;
my $file_bam;
my $out;

GetOptions(
    "help|h!" => \$help,
    "work|w:i" => \$works,
    "samtools|s:s" => \$samtools,
    "read1|r1=s" => \$file_r1,
    "bam|b=s" => \$file_bam,
    "out|o=s" => \$out
) or die $!;

die $usage if defined($help);
die "ERROR: can not find read1: ${file_r1}\n" if !( -e $file_r1);
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

my $summary_raw = process_fastq($file_r1, $works);
my $summary_mapped = process_mapped($samtools, $file_bam, $works);
output_summary($out, $summary_raw, $summary_mapped);

say "End!\n".localtime();

################################

sub count_raw_reads {
    my ($file_r1, $p, $works) = @_;

    my %summary;
    my $read_number;
    my $nreads = 1;
    say "NOTE: \tthread $p start!";

    open(my $fh_fq, "zcat $file_r1 |") or die $!;
    while(<$fh_fq>) {
        next if($. % 4 != 1 || $nreads++ % $works != $p);
        /\w+_([ATCG]{8})_[ATCG]{8}/;
        $summary{$1}++;
        $read_number++;
    }
    close($fh_fq);

    say "NOTE: \tthread $p finish!\tprocess reads: $read_number";
    return \%summary;
}

sub process_fastq {
    my ($file_r1, $works) = @_;
    say "NOTE: process raw reads!\n\tfastq: $file_r1";

    my @threads_list = map { 
        threads->create(
            {'context' => 'list'},
            \&count_raw_reads, $file_r1, $_, $works
        );
    } 0..($works-1);

    my @summary = map { $_ -> join() } @threads_list;

    return \@summary;
}

################################

sub count_mapped_reads {
    my ($samtools, $file_bam, $p, $works) = @_;

    my %summary;
    my $read_number;
    say "NOTE: \tthread $p start!";

    open(my $fh_bam, "$samtools view $file_bam |") or die $!;
    while(<$fh_bam>) {
        next if($. % $works != $p);
        /\w+_([ATCG]{8})_[ATCG]{8}/;
        $summary{$1}++;
        $read_number++;
    }
    close($fh_bam);

    say "NOTE: \tthread $p finish!\tprocess reads: $read_number";
    return \%summary;
}

sub process_mapped {
    my ($samtools, $file_bam, $works) = @_;
    say "NOTE: process mapped reads!\n\tbam: $file_bam";

    my @threads_list = map { 
        threads->create(
            {'context' => 'list'},
            \&count_mapped_reads, $samtools, $file_bam, $_, $works
        );
    } 0..($works-1);

    my @summary = map { $_ -> join() } @threads_list;

    return \@summary;
}

################################

sub output_summary {
    my ($out, $summary_raw, $summary_mapped) = @_;

    say "NOTE: reads summary!\n\tout put: $out";

    my %summary;
    # raw
    foreach $p (@$summary_raw) {
        foreach $g (keys %$p) {
            $summary{$g}{"raw"} += $p->{$g};
        }
    }
    # count
    foreach $p (@$summary_mapped) {
        foreach $g (keys %$p) {
            $summary{$g}{"mapped"} += $p->{$g};
        }
    }

    # output
    open(my $fh_out, ">$out") or die $!;
    say {$fh_out} "barcode,raw_reads,mapped_reads,mapping_ratio";
    my $raw_all;
    my $mapped_all;
    map {
        my $raw = $summary{$_}{"raw"} ? $summary{$_}{"raw"} : 0;
        my $mapped = $summary{$_}{"mapped"} ? $summary{$_}{"mapped"} : 0;
        my $ratio = $raw ? $mapped / $raw : 0;
        say {$fh_out} "$_,${raw},${mapped},${ratio}";
        $raw_all += $raw;
        $mapped_all += $mapped;
    } sort keys %summary;
    my $ratio_all = $raw_all ? $mapped_all / $raw_all : 0;
    say {$fh_out} "__total,${raw_all},${mapped_all},${ratio_all}";
    close($fh_out);
}


