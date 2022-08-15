use strict;
use Data::Dumper;
use warnings;

open GTF,"$ARGV[0]";
my $gene_id;
my %merge_band;
my %temp;
while(<GTF>){
   chomp;
   my($gene,$start,$end)=split/ /,$_;
   $gene=~ m/\"(?<gene>\w+)\"/;
   my @band=($start..$end);
   %merge_band=map{$_ =>1} @band,keys(%merge_band);
   if($gene_id && $gene_id ne $+{gene}){
        my @band_array=keys(%temp);
        my $size=@band_array;
        print "$gene_id $size\n";
        %merge_band=();
        %merge_band=map{$_ =>1} @band;
   }
   %temp=%merge_band;
   $gene_id=$+{gene};
}
close GTF;
