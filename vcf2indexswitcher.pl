#!/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $popinfo = "NOTSPECIFIED";
my $max_sites = 500000;
my $switch_rate = 0.1;

GetOptions (
        "info=s" => \$popinfo,
        "max_sites=i"   => \$max_sites,
        "switch_rate=s"  => \$switch_rate
                );

if ($switch_rate > 1){die "Switch rate is above 1, which means above 100%"};
print STDERR "Using info file $popinfo\n";
print STDERR "Calculating for $max_sites sites\n";
print STDERR "Switching $switch_rate of reads within a lane\n";


my %lane;
my %tech;
my %reads;
open POP, $popinfo or die "Count not open info file\n";
while(<POP>){
  chomp;
  if ($. == 1){next;}
  my @a = split(/\t/,$_);
  my $sample = $a[0];
  my $lane1 = $a[1];
  my $lane2 = $a[2];
  my $tech = $a[3];
  $lane{$sample}{1} = $lane1;
  $lane{$sample}{2} = $lane2;
  $tech{$sample} = $tech;
}
close POP;
my $counter;
my %sample;
my %withinlane;
while(<STDIN>){
  my $line = "$_";
  chomp $line;
  my @fields = split /\t/,$line;
  if($line=~m/^##/){
    next;
  }
  if ($line =~m/^#CHROM/){
    print "$line";
    my $first_line;
    foreach my $i (9..$#fields){
      $sample{$i} = $fields[$i];
    }
    foreach my $i (9..$#fields){
      foreach my $j (9..$#fields){
	if ($j ne $i){
	  if (($lane{$sample{$i}}{1} eq $lane{$sample{$j}}{1}) or
	     ($lane{$sample{$i}}{1} eq $lane{$sample{$j}}{2}) or
	     ($lane{$sample{$i}}{2} eq $lane{$sample{$j}}{2}) or
	     ($lane{$sample{$i}}{2} eq $lane{$sample{$j}}{1})){
	       $withinlane{$sample{$i}}{$sample{$j}}++;
	  }
	}
      }
    }
  }
  else{
    $counter++;
    if ($counter > $max_sites){goto ENDSCRIPT;}
    if ($counter % 100000 == 0){print STDERR "Processed $counter sites\n";}
    my $chr = $fields[0];
    my $pos = $fields[1];
    my $alt = $fields[4];
    my $multi_alt;
    my @alts;
    @alts = split(/,/,$alt);
    if (length($alt) > 1){
      next;
    }
    print "\n$fields[0]";
    foreach my $i (1..8){
      print "\t$fields[$i]";
    }
    my @test_samples;
    my %depth;
    my %moved_reads;
    my %removed_reads;
    my %genotype;
    foreach my $i (9..$#fields){
      if ($fields[$i] eq '.'){
	$depth{$sample{$i}}{0} = 0;
	$depth{$sample{$i}}{1} = 0;
	$genotype{$sample{$i}} = '.';
      }elsif ($fields[$i] ne '.'){
        my @info = split(/:/,$fields[$i]);
        my $call = $info[0];
	$genotype{$sample{$i}} = $call;
        my @bases = split(/\//,$call);
        my $dp = $info[1];
        my $ref_dp = $info[3];
        my $alt_dp = $info[5];
	$depth{$sample{$i}}{0} = $ref_dp;
	if ($ref_dp){
	  foreach my $j (1..$ref_dp){
	    my $rand = rand(1);
	    if ($rand < $switch_rate){
	      my $moved_sample = $withinlane{$sample{$i}}{(keys %{$withinlane{$sample{$i}}})[rand keys %{$withinlane{$sample{$i}}}]};
	      $moved_reads{$moved_sample}{0}++;
	      $removed_reads{$sample{$i}}{0}++;
	    }
	  }
	}
	$depth{$sample{$i}}{1} = $alt_dp;
        if ($alt_dp){
          foreach my $j (1..$alt_dp){
            my $rand = rand(1);
            if ($rand < $switch_rate){
              my $moved_sample = (keys %{$withinlane{$sample{$i}}})[rand keys %{$withinlane{$sample{$i}}}];
              $moved_reads{$moved_sample}{1}++;
	      $removed_reads{$sample{$i}}{1}++;
            }
          }
        }
      }
    }
    foreach my $i (9..$#fields){
      
      my $ref_dp = $depth{$sample{$i}}{0};
      if ($moved_reads{$sample{$i}}{0}){
	$ref_dp +=$moved_reads{$sample{$i}}{0};
      }
      if ($removed_reads{$sample{$i}}{0}){
	$ref_dp -=$removed_reads{$sample{$i}}{0};
      }

      my $alt_dp = $depth{$sample{$i}}{1};
      if ($moved_reads{$sample{$i}}{1}){
	$alt_dp +=$moved_reads{$sample{$i}}{1};
      }
      if ($removed_reads{$sample{$i}}{1}){
	$alt_dp -=$removed_reads{$sample{$i}}{1};
      }
      my $total_dp = $ref_dp + $alt_dp;
      if ($total_dp == 0){
	print "\t.";
	next;
      }elsif(($ref_dp > 0) and ($alt_dp == 0)){
	print "\t0/0:$total_dp:X:$ref_dp:X:$alt_dp";
      }elsif(($ref_dp > 0) and ($alt_dp > 0)){
	print "\t0/1:$total_dp:X:$ref_dp:X:$alt_dp";
      }elsif(($ref_dp == 0) and ($alt_dp > 0)){
	print "\t1/1:$total_dp:X:$ref_dp:X:$alt_dp";
      }else{
	print STDERR "Something wrong with genotype $total_dp $ref_dp $alt_dp\n";
      }
    }
  }
}
ENDSCRIPT:
