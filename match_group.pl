#! /usr/bin/perl -w
##########################################################################
#### Author: Ping HU
#### 
##########################################################################
use strict;
my $usage="$0  [file] first line is title\n"; 
die $usage unless @ARGV ==1;
my ($file) = @ARGV;

open (A, "</Users/ping/Desktop/BabySkin/src/BabyMeta_10212020.txt")||die "could not open Baby_meta_09202020.txt\n";
my %newID;
while (my $a=<A>){
    chomp $a;
    for($a){s/\r//gi; s/\.CEL//gi;}
    my @tmp=split /\t/, $a;
    $newID{$tmp[0]}=$tmp[23];
}
close A;

open (B, "<$file")||die "could not open $file\n";
open (C, ">$file.ann.xls")||die "could not open $file.ann.xls\n";
my $x=<B>;
chomp $x;
for($x){
    s/\r//gi;
    s/\.CEL//gi;
}
my @aa=split /\t/, $x;
print C $aa[0];
for(my $i=1; $i<@aa; $i++){
    my $id=$aa[$i];
    if(defined $newID{$id}){
	print C "\t",$newID{$id};
    }else{
	print STDERR "not find ", $id, "\n";
	print C "\t","NotFound_", $id;
    }

}
print C "\n";
while (my $b=<B>){
    print C $b;
}
close B;
