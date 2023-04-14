#!/usr/bin/perl
#author:wangxin
### function: Combine the evaluation folder and print the final evaluation of Read length,  Per base sequence quality, Per base N content, Overrepresented sequences, apdaptor content, and Duplication Level
### 
use strict;
use warnings;


my $version="1.0 version";
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i:s","o:s");
print "*************\n*$version*\n*************\n";
if ( !defined $opts{i} || !defined $opts{o}) {
       	die "************************************************************************
       	Usage: $0.pl
       		-i: Query folder with samples of reads counts
			-o: Combined output
************************************************************************\n";
}


### read the mutators files in the folder

## read the input folder
my $input=$opts{i};
my @files;
opendir IN,"$input" or die $!;
@files=readdir IN; 
close IN;

### then read the files in the folder
my $n=0; my %hash;my %namestr;

#### sort the name and put into a hash
foreach my $i (@files){
	
	next unless ($i=~/count$/);
	#my $name=join "_",(split/\_/,$i)[2,3,4];
	
	#print "$i\n";
	$namestr{$i}++;
	$n++;
	open FILE,"$input/$i" or die $!;
	while (<FILE>){
		chomp;
		s/\r//;
		s/^\s+//;
		my ($gene,$cov)=split/\s+/,$_;
		next if ($gene eq "ID");
		$hash{$gene}->{$i}=$cov;
	}
	close FILE;
}
	


my $output=$opts{o};


open OUT,">$output" or die $!;
### print the header

print OUT "EnsemblID";
foreach my $a (sort keys %namestr){
	print  OUT "\t$a";
}

print OUT "\n";


foreach my $m (sort keys %hash){
	
	print OUT "$m";
	
	foreach my $n(sort keys %namestr){
		print "$m\t$n\n" if (!defined $hash{$m}->{$n});
		
		my $finalcounts=(defined $hash{$m}->{$n})?$hash{$m}->{$n}:0;
		print OUT "\t$finalcounts";
	}
	
	print OUT "\n";
	
}
close OUT;