#!/usr/bin/perl

# Copyright (C) 2017, Joseph F. Ryan <joseph.ryan@whitney.ufl.edu>
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# for more info see <http://www.gnu.org/licenses/>

# --min_sp = minimum number of species in an orthogroup (OG)

use strict;
use warnings;
use JFR::Fasta;
use Getopt::Long;
use IO::File;
use Data::Dumper;
use POSIX qw/strftime/;

our $DIR = '/bwdata1/mdebiasse/01-DEEPC/05-OrthoFinder/05-MBARI36_Mnem/Results_Feb19/WorkingDirectory/Sequences_1Neph_red_multi_Mert_ovum_28_min';

our $OUTDIR = 'OUTDIR.';
our $IQDIR = '/bwdata1/mdebiasse/01-DEEPC/07-PERL/07-FILTER_ORTHOGROUPS_FOR_PTP/OUT.IQ';

our $VERSION = 0.01;
our $MIN_TO_RETAIN = 28;
our $THREADS = 20;
our $STDERR_DIR = 'OUTDIR.STDERR.';

MAIN: {
    my $rh_opts = get_opts();
    opendir DIR, $DIR, or die "cannot open $DIR:$!";
    my @seq_files = grep { /\.fa$/ } readdir DIR;
    $OUTDIR .= strftime('%Y-%m-%d.%H-%M-%S',localtime);
    $STDERR_DIR .= strftime('%Y-%m-%d.%H-%M-%S',localtime);
    die "$OUTDIR exists" if (-d $OUTDIR);
    mkdir $OUTDIR or die "cannot mkdir $OUTDIR:$!";
    die "$STDERR_DIR exists" if (-d $STDERR_DIR);
    mkdir $STDERR_DIR or die "cannot mkdir $STDERR_DIR:$!";
    opendir DIR, $DIR or die "cannot opendir $DIR:$!";
    my @data = ();
    FA: foreach my $sf (@seq_files) {
        $sf =~ m/(.*).fa/ or die "unexpected: $sf";
        my $og = $1;
        my $count = 0;
        my %species = ();
        my $fp = JFR::Fasta->new("$DIR/$og.fa");
        while (my $rec = $fp->get_record()) {
            $count++;
            $rec->{'def'} =~ m/^>([^_]+_[^_]+)_/;
            my $sp = $1;
            $species{$sp}++;
        }
        my $sp_count = scalar(keys %species);
        my $ratio = $count / $sp_count;
        push @data, [$og,$count,$sp_count,$ratio];
    }
    my $ra_fh = get_handles($rh_opts->{'num_scripts'});
    my $fh_count = 0;
    foreach my $ra_d (sort { $a->[3] <=> $b->[3] } @data) {
        my $og = $ra_d->[0];
        open IN, "$DIR/$og.fa" or die "cannot open $DIR/$og.fa$!";
        open OUT, ">$OUTDIR/$og.fa" or die "cannot open >$OUTDIR/$og.fa:$!";
        while (my $line = <IN>) {
            print OUT $line;
        }
        close OUT;
        close IN;
        my $curfh = $ra_fh->[$fh_count];
        print $curfh "echo $ra_d->[0] >> $STDERR_DIR/gbwrapper.out\n";
        print $curfh "echo $ra_d->[0] >> $STDERR_DIR/gbwrapper.err\n";
        print $curfh "echo $ra_d->[0] >> $STDERR_DIR/remove_empty_seqs.err\n";
        print $curfh "echo $ra_d->[0] >> $STDERR_DIR/mafft.err\n";
        print $curfh "echo $ra_d->[0] >> $STDERR_DIR/IQ.out\n";
        print $curfh "echo $ra_d->[0] >> $STDERR_DIR/IQ.err\n";
        print $curfh "echo $ra_d->[0] >> $STDERR_DIR/ptp.out\n";
        print $curfh "echo $ra_d->[0] >> $STDERR_DIR/ptp.err\n";

        print $curfh "mafft-linsi --localpair --maxiterate 1000 --thread $THREADS $OUTDIR/$og.fa > $OUTDIR/$og.mafft 2>> $STDERR_DIR/mafft.err\n";
        print $curfh "Gblockswrapper $OUTDIR/$og.mafft >> $STDERR_DIR/gbwrapper.out 2>>$STDERR_DIR/gbwrapper.err\n";
        print $curfh "remove_empty_seqs $OUTDIR/$og.mafft-gb  >> $STDERR_DIR/remove_empty_seqs.out 2>>$STDERR_DIR/remove_empty_seqs.err\n";
        print $curfh "iqtree-omp -s $OUTDIR/$og.mafft-gb.res -nt AUTO -bb 1000 -m LG -pre $og.iq >> $STDERR_DIR/IQ.out 2>> $STDERR_DIR/IQ.err\n";  
        print $curfh "mkdir $og.IQ_out; mv *.contree *.splits.nex *.iq.ckp.gz *.iq.bionj *.iq.log *.treefile *.mldist *.iqtree $og.IQ_out\n";
        print $curfh q~perl -pi -e 's/([A-Z][a-z]{3}_[a-z0-9]+)_([0-9]+\.[0-9]+)/$1\|$2/g' ~; 
        print $curfh "$og.IQ_out/$og.iq.contree\n";
        print $curfh q~perl -pi -e 's/([^_]+_[^_]+)_([^_]+\.[^_]+)/$1\|$2/g' ~;
        print $curfh "$OUTDIR/$og.mafft-gb.res\n";
        print $curfh "java PhyloTreePruner $og.IQ_out/$og.iq.contree $MIN_TO_RETAIN $OUTDIR/$og.mafft-gb.res 0.5 u >> $STDERR_DIR/ptp.out 2>> $STDERR_DIR/ptp.err\n\n";

        if ($fh_count == 0) {
            $fh_count++;
        } elsif (($fh_count + 1) % $rh_opts->{'num_scripts'}) {
            $fh_count++;
        } else {
            $fh_count = 0;
        }
    }
}

sub get_handles {
    my $num = shift;
    my @handles = ();
    for (my $i = 0; $i < $num; $i++) {
        my $name = "script.$i";
        my $fh = IO::File->new($name,'w');
        push @handles, $fh;
    }
    return \@handles;
}

sub get_opts {
    my $rh_opts = {'num_scripts' => 0 };
    my $res = GetOptions ("num_scripts=i" => \$rh_opts->{'num_scripts'});

    usage() unless ($rh_opts->{'num_scripts'});

    return $rh_opts;
}

sub usage {
    die "usage: $0 --num_scripts=NUM_SCRIPTS\n";
}




