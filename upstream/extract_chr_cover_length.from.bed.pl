use strict;
use Data::Dumper;
use warnings;

open BED,"$ARGV[0]";
my %cover_length;
my $chr_last;
while(<BED>){
   chomp;
   my($chr,$start,$end)=split/\t/,$_;
   $cover_length{$chr}+=$end-$start+1;
   print "$chr_last $cover_length{$chr_last}\n" if($chr_last && $chr_last ne $chr);
   $chr_last=$chr;
}
close BED;
