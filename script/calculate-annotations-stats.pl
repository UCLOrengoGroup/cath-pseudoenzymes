#!/usr/bin/env perl


# go through each superfamily
# count the number of .annos files that are empty / non-empty
# calculate these values as %s

use strict;
use warnings;

use Data::Dumper;
use File::stat;
use File::Basename;
use Getopt::Long;
use List::Util qw/ sum uniq /;

use DDP;
use JSON;
use Log::Log4perl::Tiny qw/ :easy /;
use Path::Tiny;

my $RESULTS_DIR = "/cath/people2/ucbtnld/projects/pseudoenzymes/";
my $STATS_TSV_FILE;
my $STATS_JSON_FILE;

my $USAGE = <<"_USAGE_";
usage: $0 --tsv <stats.data> --json <stats.json>

Generates summary annotation stats

_USAGE_

die $USAGE unless @ARGV;

GetOptions(
    'results|i=s', \$RESULTS_DIR,
    'tsv|t=s', \$STATS_TSV_FILE,
    'json|j=s', \$STATS_JSON_FILE,
);

die "results dir '$RESULTS_DIR' does not exist" unless -d $RESULTS_DIR;

# get results directory
my $results_dir = path( $RESULTS_DIR );

# make output file
my $stats_tsv_fh = $STATS_TSV_FILE ? path( $STATS_TSV_FILE )->openw() : undef;

# go through each superfamily subdirectory
my @superfamily_dirs = $results_dir->children( qr/^[1-4]\.\d+\.\d+\.\d+$/ );

# INFO Dumper \@superfamily_dirs;

my @all_sfam_data;

sub parse_funfam_file {
    my $ff_file = shift;
    if ( $ff_file !~ m{\b(\d+\.\d+\.\d+\.\d+)/(\d+)\b} ) {
        die "failed to parse funfam id from file: '$ff_file'";
    }
    my $sfam_id = $1;
    my $funfam_number = $2;
    my $funfam_id = "${sfam_id}-ff-${funfam_number}";
    return ($sfam_id, $funfam_number, $funfam_id);
}

my @table_cols = qw/superfamily_id anno_empty_count anno_empty_perc anno_nonempty_count anno_nonempty_perc anno_total_count/;
print $stats_tsv_fh join( "\t", @table_cols ) . "\n";

SUPERFAMILY: foreach my $superfamily_dir ( @superfamily_dirs ) {

    # extract all the accession and anno files from within directory
    my @acc_files = $superfamily_dir->children( qr/^\d+\.accs$/ );
    my @anno_files = $superfamily_dir->children( qr/^\d+\.accs\.annos$/ );

    # remember to skip empty directories
    # if there are no anno files, skip superfamily
    next SUPERFAMILY if ( scalar @anno_files == 0 );

    INFO $superfamily_dir;
    my $superfamily_id = basename $superfamily_dir;

    my %ec_by_funfam;
    my %go_by_funfam;
    my %acc_by_funfam;

    # count the accessions
    for my $acc_file ( @acc_files ) {
        my @accs = $acc_file->lines;
        my ($sfam_id, $funfam_number, $funfam_id) = parse_funfam_file($acc_file);
        $acc_by_funfam{$funfam_id} = scalar @accs;
    }

    # count the number of .annos files that are empty / non-empty
    my $anno_empty_count = 0;
    my $anno_nonempty_count = 0;
    my $anno_total_count = 0;

    foreach my $anno_file ( @anno_files ) {
        # INFO $anno_file;
        
        my ($sfam_id, $funfam_number, $funfam_id) = parse_funfam_file($anno_file);
        
        my @anno_lines = path($anno_file)->lines;
        for my $anno_line ( @anno_lines ) {
            chomp($anno_line);
            my @cols = split(/,\s*/, $anno_line);
            my $ec_code = $cols[0];
            
            $ec_by_funfam{$funfam_id} //= {};
            $go_by_funfam{$funfam_id} //= {};

            $ec_by_funfam{$funfam_id}->{$ec_code}++;

            if ( scalar @cols > 1 ) {
                my ($go_acc, $evi_code, $source, $name) = @cols[1..4];
                $go_by_funfam{$funfam_id}->{$go_acc} = {
                    go_acc => $go_acc,
                    evidence => $evi_code,
                    source => $source,
                    name => $name,
                };
            }
        }

        # get file stats
        my $stats = stat( $anno_file );
        my $file_size = $stats->size;

        if ( $file_size == 0 ) {
            $anno_empty_count++;
        } else {
            $anno_nonempty_count++;
        }

        $anno_total_count++;
    }

    my @funfam_data = map {
        my $ff_id = $_;
        my $uniq_ec = [ uniq sort keys %{$ec_by_funfam{$ff_id}} ];
        my $uniq_go = [ uniq sort map { sprintf("%s: %s", $_->{go_acc}, $_->{name}) } values %{$go_by_funfam{$ff_id}} ];
        {
            funfam_id => $ff_id,
            seq => $acc_by_funfam{$ff_id},
            ec => $uniq_ec,
            go => $uniq_go,
        }
    } keys %acc_by_funfam;
    
    my %sfam_data = (
        sfam_id => $superfamily_id,
        funfam_count => (scalar keys %acc_by_funfam),
        funfam_seq_count => (sum values %acc_by_funfam),
        empty_count => $anno_empty_count,
        non_empty_count => $anno_nonempty_count,
        sfam_uniq_ec_count => scalar( uniq map { keys %$_ } (values %ec_by_funfam) ),
        sfam_uniq_go_count => scalar( uniq map { keys %$_ } (values %go_by_funfam) ),
        funfam_data => \@funfam_data,
    );

    push @all_sfam_data, \%sfam_data;

    # calculate percentages
    my $anno_empty_perc    = ( $anno_empty_count / $anno_total_count ) * 100;
    my $anno_nonempty_perc = ( $anno_nonempty_count / $anno_total_count ) * 100;

    # print out stats
    INFO "Empty count:    $anno_empty_count ($anno_empty_perc)";
    INFO "Nonempty count: $anno_nonempty_count ($anno_nonempty_perc)";
    INFO "Total count:    $anno_total_count";
    
    # append values to file
    my $string = sprintf( "%s\t%d\t%d\t%d\t%d\t%d\n", $superfamily_id, $anno_empty_count, $anno_empty_perc, $anno_nonempty_count, $anno_nonempty_perc, $anno_total_count );
    print $stats_tsv_fh $string;
    # last;
}
close $stats_tsv_fh;

if ($STATS_JSON_FILE) {
    INFO "Writing JSON data to file: $STATS_JSON_FILE ...";
    my $stats_json_fh = path( $STATS_JSON_FILE )->openw();
    my %data = (
        entry_count => scalar @all_sfam_data,
        data => \@all_sfam_data,
    );
    my $all_data = encode_json \%data;
    print $stats_json_fh $all_data;
    close $stats_json_fh;
    INFO sprintf( "   ... wrote records (%d sfams)", scalar @all_sfam_data);
}
