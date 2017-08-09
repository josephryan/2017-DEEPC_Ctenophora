#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use JFR::Fasta;

my $dir = $ARGV[0] or die "usage: $0 INDIR OUTDIR\n"; # asks user to type indir name. should be OrthoFinder/Sequences dir
my $outdir = $ARGV[1] or die "usage: $0 INDIR OUTDIR\n"; #asks user to type outdir name
die "$outdir exists" if (-d $outdir);
mkdir $outdir or die "cannot mkdir $outdir:$!";

opendir DIR, $dir or die "cannot opendir $dir$!"; #open the indir

my @files = grep { /fa$/ } readdir DIR; #put all the files that end in .fa into array @seq_file

our $MIN_SEQS = 28; # if changed, need2 also change line: "last if ($1 > 11077)"
foreach my $file (@files){
    $file =~ m/OG0*(\d+)\.fa/ or die "unexpected filename:$!";
    last if ($1 > 200000);  # this is the last file with 27 seqs
    my %seen = ();
    my $new_fa = '';
    my $nr2_seq = '';
    my $num_seqs = 0;
    my $nr1_flag = 0;
    my $nr2_flag = 0;

    open IN, "$dir/$file" or die "cannot open $dir/$file:$!";
    while (my $line = <IN>) {
        if ($line =~ m/^>([A-Z][a-z]{3}_[a-z0-9]+)/) {
            my $sp = $1;
            if ($seen{$sp}) {
                $num_seqs = 0;
                last;
            } else {
                $seen{$sp} = 1;
            }
            $seen{$sp} = 0 if ($sp eq 'Cydi_peac');
            $seen{$sp} = 0 if ($sp eq 'Mert_ovum');
            if ($sp eq 'Neph_red2') {
                $nr2_flag = 1;
                $nr2_seq = $line;
            } else {
                $nr2_flag = 0;
                $new_fa .= $line;
                $num_seqs++;
                $nr1_flag = 1 if ($sp eq 'Neph_red1');
            }     
        } else {
            if ($nr2_flag) {
                $nr2_seq .= $line;
            } else {
                $new_fa .= $line;
            }
        }
    }
    $new_fa .= $nr2_seq unless ($nr1_flag);
    next unless ($num_seqs >= $MIN_SEQS);
    open OUT, ">$outdir/$file" or die "cannot open >$outdir/$file:$!";
    print OUT $new_fa;
}






       
        

