#!/usr/bin/env perl
use strict;
use warnings;
use List::Util qw/min max/;
use List::MoreUtils qw/part/;
# use Smart::Comments '###';

# In: two GFFs, features that are are designated as matches (e.g. orthologs)
# Out: regions where there's at least one match every 10kbp in both GFFs


my ($gff_1, $gff_2, $markers_tsv, $gff_1_out, $gff_2_out, $threshold_1, $threshold_2, $min_report_1, $min_report_2) = @ARGV;
die "Usage: $0 gff_1 gff_2 matching_features out_1, out_2 <t1> <t2>" unless -f $gff_1 and -f $gff_2 and -f $markers_tsv;
$threshold_1 //= 2000;
$threshold_2 //= 2000;
$min_report_1 //= 20 * $threshold_1;
$min_report_2 //= 20 * $threshold_2;

my ($markers_1_to_2s, $markers_2_to_1s) = read_marker_pairs($markers_tsv); 
my ($markers_1_in_order, $markers_1_by_name) = read_marker_positions_from_gff($gff_1, $markers_1_to_2s);
my ($markers_2_in_order, $markers_2_by_name) = read_marker_positions_from_gff($gff_2, $markers_2_to_1s);
open (my $out_fh_1, ">", $gff_1_out) or die "$!: $gff_1_out";
open (my $out_fh_2, ">", $gff_2_out) or die "$!: $gff_2_out";
print $out_fh_1 "# $gff_1 - Regions of at least $min_report_1 bp with markers at most every $threshold_1 bp and of at least $min_report_2 bp in $gff_2 at most every $threshold_2 bp\n";
print $out_fh_2 "# $gff_2 - Regions of at least $min_report_2 bp with markers at most every $threshold_2 bp and of at least $min_report_1 bp in $gff_1 at most every $threshold_1 bp\n";

for my $contig_1 ( sort keys %{$markers_1_in_order} ){
  for my $contig_2 (sort keys %{$markers_2_in_order} ){
     find_regions (
        $contig_1, $contig_2,
        $markers_1_in_order->{$contig_1},
        $markers_2_by_name,
        $markers_1_to_2s,
        $out_fh_1, $out_fh_2,
        $threshold_1, $threshold_2,
        $min_report_1, $min_report_2,
     );
  }
}
sub read_marker_positions_from_gff {
  my ($path, $marker_names) = @_;
  my %coordinates_by_contig;
  my %coordinates_by_name;
  open(my $fh, "<", $path) or die "$!: $path";
  while(<$fh>){
    next if /^#/;
    chomp;
    my @F = split "\t";
#SM_V7_2	WormBase_imported	gene	33749861	33773944	.	+	.	ID=gene:Smp_154870;Name=Smp_154870;biotype=protein_coding
    my ($contig, $feature, $start, $end, $strand, $misc) = @F[0,2,3,4,6,8];
    next unless $misc;
    my ($name) = $misc =~ /Name=(.*?);/;
    next unless $name;
    next unless $marker_names->{$name};
    my $o = { 
      name => $name,
      start => int($start),
      end => int($end),
      contig => $contig,
    };
    push @{$coordinates_by_contig{$contig}}, $o;
    $coordinates_by_name{$name} = $o;
  }
  close $fh;
  for my $contig (keys %coordinates_by_contig) {
     my @xs = sort { $a->{start} <=> $b->{start} } @{$coordinates_by_contig{$contig}};
     $coordinates_by_contig{$contig} = \@xs;
  }
  return (\%coordinates_by_contig, \%coordinates_by_name);
}
sub read_marker_pairs {
  my ($path) = @_;
  my %markers_1_to_2;
  my %markers_2_to_1;
  open(my $fh, "<", $path) or die "$!: $path";
  while(<$fh>){
    chomp;
    my ($name_1, $name_2) = split "\t";
    push @{$markers_1_to_2{$name_1}}, $name_2;
    push @{$markers_2_to_1{$name_2}}, $name_1;
  }
  close $fh;
  return (\%markers_1_to_2, \%markers_2_to_1); 
}
sub find_regions { 
   my ( $contig_name_1, $contig_name_2,
        $features_1_in_order, $features_2_by_name,
       $markers_1_to_2s, 
       $out_fh_1, $out_fh_2,
       $threshold_1, $threshold_2,
       $min_report_1, $min_report_2) = @_;
   my @partials;
### find_regions: $contig_name_1, $contig_name_2
   for my $feature_1 (@{$features_1_in_order}){
### $feature_1
      # for some partial matches, feature 1 is too far out already - print them out
      my @partials_next;
      for my $partial (@partials){
         if ($feature_1->{start} - $partial->{end_1} < $threshold_1){
            push @partials_next, $partial;
         } elsif($partial->{end_1} - $partial->{start_1} > $min_report_1 && $partial->{end_2} - $partial->{start_2} > $min_report_2) {
            
            output_region($out_fh_1, $contig_name_1,
              $partial->{start_1}, $partial->{end_1},
              sprintf("length %s, matching %s: %s-%s, of length %s" ,$partial->{end_1} - $partial->{start_1},  $contig_name_2, $partial->{start_2}, $partial->{end_2},
                $partial->{end_2} - $partial->{start_2},
              )
            );
            output_region($out_fh_2, $contig_name_2,
              $partial->{start_2}, $partial->{end_2},
              sprintf("length %s, matching %s: %s-%s, of length %s" , $partial->{end_2} - $partial->{start_2}, $contig_name_1, $partial->{start_1}, $partial->{end_1},
                $partial->{end_1} - $partial->{start_1},
              )
            );
         }
      }
      @partials = @partials_next;
      # extend partial matches
      for my $feature_2 (grep { $_->{contig} eq $contig_name_2 } map {$features_2_by_name->{$_} // () } @{$markers_1_to_2s->{$feature_1->{name}}//[]} ){
### $feature_2
         my ($matching, $other) = part {
            # feature_1 is close enough to all partials that remain
            # does feature_2 extend at the start?
            
            $feature_2->{start} < $_->{start_2} && $_->{start_2} - $feature_2->{end} < $threshold_2
            || # at the end ?
            $feature_2->{end} > $_->{end_2} &&  $feature_2->{start} - $_->{end_2} < $threshold_2
            # goes inside ?
            || $feature_2->{start} > $_->{start_2} && $feature_2->{end} < $_->{end_2}
            ? 1 : 0
         } @partials;
### $matching
         if ($matching) {
### @partials
            @partials = ({
               start_1 => min ($feature_1->{start}, map {$_->{start_1}} @{$matching}),
               end_1 => max ($feature_1->{end}, map {$_->{end_1}} @{$matching}),
               start_2 => min ($feature_2->{start}, map {$_->{start_2}} @{$matching}),
               end_2 => max ($feature_2->{end}, map {$_->{end_2}} @{$matching}),
            }, @{$other//[]});
### @partials
         } else { #new partial
            push @partials, {
              start_1 => $feature_1->{start}, 
              end_1 => $feature_1->{end},
              start_2 => $feature_2->{start},
              end_2 => $feature_2->{end},
            };
         }
      }
   }
}
sub output_region {
  my ($fh, $contig, $start, $end, $desc) = @_;
  print $fh join ("\t", $contig, "region-markers", "region", $start, $end, ".", "+", ".", $desc). "\n";
}
