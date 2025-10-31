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
$0: add cell barcode and umi to R1.fastq for STRT-seq

perl $0 [options] 
    [--help|-h]                         : show usage 
    --work|-w       threads_number(5)
    --read1|-r1     input_R1.fq.gz
    --read2|-r2     input_R2.fq.gz
    --barcode|-b    barcode.tsv         : barcode fiile with no header, firt column is cell name, second column is barcodes
    --out|-o        out_put_file        : create new file to store the result fastq

eg:
    perl $0 -w 5 -r1 test_R1.fq.gz -r2 test_R2.fq.gz -b STRT_barcode.txt -o test_R1_extracted.fq.gz

/;

################################

my $help;
my $works;
my $file_r1;
my $file_r2;
my $file_barcode;
my $out;

GetOptions(
    "help|h!" => \$help,
    "work|w:i" => \$works,
    "read1|r1=s" => \$file_r1,
    "read2|r2=s" => \$file_r2,
    "barcode|b=s" => \$file_barcode,
    "out|o=s" => \$out
) or die $!;

die $usage if defined($help);
die "ERROR: can not find read1: ${file_r1}\n" if !( -e $file_r1);
die "ERROR: can not find read2: ${file_r2}\n" if !( -e $file_r2);
die "ERROR: can not find barcode file: ${file_barcode}\n" if !( -e $file_barcode);
die "WRONG: output ${out} already exists\n" if ( -e $out);
$works = 5 if (!defined($works));

# main
say "Start!\n".localtime();

my $barcodes = load_barcode($file_barcode);
my $summary = process_reads($file_r1, $file_r2, $out, $barcodes, $works);
out_put_summary($out, $barcodes, $summary);

say "End!\n".localtime();

################################
# barcode

sub load_barcode {
    my ($file_barcode) = @_;
    say STDERR "NOTE: load and process barcode: $file_barcode";

    open(my $file_in, "<$file_barcode") or die $!;
    my %barcodes = map {
            chomp;
            my ($cell, $bar) = split(/\s+/);
            $bar => $cell
        } <$file_in>;
    close($file_in);

    return \%barcodes;
}

################################
# read

sub filter_write_read {
    my ($file_r1, $file_r2, $out, $barcodes, $p, $works) = @_;

    my %summary;
    my $pid = threads->tid();
    $out .= ".${pid}.tmp";

    open(my $fh_r1, "zcat $file_r1 |") or die $!;
    open(my $fh_r2, "zcat $file_r2 |") or die $!;
    open(my $fh_out, "| gzip >$out") or die $!;

    say "NOTE: \tthread $pid start!";

    while(
        defined(my $read_id1 =  <$fh_r1>) && 
        defined(my $read_id2 =  <$fh_r2>)
    ) {
        next if(($. - 1 >> 2) % $works != $p);

        # check CI
        my ($CI1)    =  $read_id1 =~ /^(\S+)\s+/;
        my ($CI2)    =  $read_id2 =~ /^(\S+)\s+/;
        if(!defined($CI1) || !defined($CI2) || $CI1 ne $CI2) {
            die "ERROR: reads loss!!!";
        }

        # R1
        my $seq1     =  <$fh_r1>;
        my $qual_id1 =  <$fh_r1>;
        my $qual1    =  <$fh_r1>;
        # R2
        my $seq2     =  <$fh_r2>;
                        <$fh_r2>;
                        <$fh_r2>;

        # check CB
        my ($CB, $CU) = $seq2 =~ /^(\w{8})(\w{8})/;
        if(defined($CB) && exists($$barcodes{$CB})) {
            $read_id1 =~ /.*?\s+(\S+\n)$/;
            $read_id1 = "${CI1}_${CB}_${CU} $1";

            print {$fh_out} $read_id1.$seq1.$qual_id1.$qual1;

            $summary{$CB}++;
        } else {
            $summary{"__unmatched"}++;
        }

        $summary{"__total"}++;
    }
    close($fh_out);
    close($fh_r2);
    close($fh_r1);

    say "NOTE: \tthread $pid finish!\n\tthread $pid process reads: $summary{'__total'}";
    return \%summary;
}

sub process_reads {
    my ($file_r1, $file_r2, $out, $barcodes, $works) = @_;
    say "NOTE: process reads!\n\tread1: $file_r1\n\tread2: $file_r2\n\toutput: $out";

    #    
    my @threads_list = map { 
            threads->create(
                {'context' => 'list'},
                \&filter_write_read, 
                $file_r1, $file_r2, $out, 
                $barcodes, $_, $works
            );
        } 0..($works-1);

    my @summary = map { $_ -> join() } @threads_list;

    system("cat ${out}.*.tmp > $out"); 
    system("rm ${out}.*.tmp");

    return \@summary;
}

################################
# summary

sub out_put_summary {
    my ($out, $barcodes, $sum_array) = @_;

    my @sum_items = ((keys %$barcodes), "__unmatched", "__total" );
    my %summary = map {
        my $CB = $_;
        my $res = sum map {
            exists($$_{$CB}) ? $$_{$CB} : 0 ; 
            } @$sum_array;
        $CB => $res
    }  @sum_items;

    say "NOTE:";
    say "\tunmatched: \t$summary{'__unmatched'}";
    say "\ttotal_reads: \t$summary{'__total'}";

    $out =~ s/\..*//;
    $out .= "_summary.tsv";
    open(my $fh_summary, ">$out") or die $!;

    say $fh_summary "Barcode\tReadNumber";
    map { say $fh_summary "$_\t$summary{$_}"; } @sum_items;

    close($fh_summary);
}


