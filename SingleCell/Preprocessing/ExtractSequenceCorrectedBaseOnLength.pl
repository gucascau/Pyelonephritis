#!/usr/bin/perl
#author:wangxin
### #### caluculate the side effect on the final associated genes results
### Here is to calculate the sample size effect



use strict;
use warnings;
#use Statistics::Test::WilcoxonRankSum;
#use List::Util qw(sum);
#use List::Util 'shuffle';

my $version="1.0 version";
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i:s","o:s","f:s","r:s");
print "*************\n*$version*\n*************\n";
if (!defined $opts{o} ||!defined $opts{r} ||!defined $opts{f} ) {
       	die "************************************************************************
       	Usage: $0.pl
				-i: length cut-off (default 51)
			-f: forward fastq1
			-r: reverse fastq2
			-o: output string 
************************************************************************\n";
}


my $input=$opts{i}?$opts{i}:51;
#my $index=$opts{g};
my $fastq1=$opts{f};
my $fastq2=$opts{r};

my $out=$opts{o};



open FASTQ1,"$fastq1" or die "cannot open file $fastq1";

my $id1; my %string1; my %hash1;
while (<FASTQ1>) {
	chomp;
    # print "$_" ;

    if ($. % 4 == 1)  {
		$id1=(split/\s+/,$_)[0];
		$id1=~s/^@//;
		$hash1{$id1}=$_."\n";
	}elsif($.%4 == 2){

		
		my $finalstr=$_;
		
		# here we updated the length with 50bp, we add N for sequence and # for quality
		if (length $_ == 50){
			$finalstr=join "",($finalstr,"N");
			print ">$id1\n$finalstr\n";
		}
		
		$string1{$id1}=$finalstr;
		$hash1{$id1}.=$finalstr."\n";
		
		
	}elsif($.%4 == 0){
		my $finalQual=$_;
		
		# here we updated the length with 50bp, we add # for quality
		if (length $_ == 50){
			$finalQual=join "",($finalQual,"#");
		}	
		$hash1{$id1}.=$finalQual."\n";
	}else{
		$hash1{$id1}.=$_."\n";
	}

}
close FASTQ1;

my $id2; my %string2; my %hash2;
open FASTQ2, $fastq2 or die "cannot open file $fastq2";

while (<FASTQ2>) {
	chomp;
    # print "$_" ;

    if ($. % 4 == 1)  {
		$id2=(split/\s+/,$_)[0];
		$id2=~s/^@//;
		$hash2{$id2}=$_."\n";
	}elsif($.%4 == 2){
		$string2{$id2}=$_;
		$hash2{$id2}.=$_."\n";

	}else{
		$hash2{$id2}.=$_."\n";
	}
}
close FASTQ2;




### print the reads with specific length (51bp)


open PH1,">$out\_R1_001.fastq" or die $!;
open PH2,">$out\_R2_001.fastq" or die $!;
#open IN,"$input" or die $!;
foreach my $m (keys %string1){
	#my $id =(split /\t/,$_)[0];
	my $Readlength = length $string1{$m};
	next unless ($Readlength == 51);
	print PH1 "$hash1{$m}";
	print PH2 "$hash2{$m}";
	#$phix{$_}++;
}

close PH1;
close PH2;


