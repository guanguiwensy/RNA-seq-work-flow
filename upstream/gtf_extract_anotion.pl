use strict;
use Data::Dumper;
use warnings;

open GTF,"$ARGV[0]";
my @class_select="exon";
while(<GTF>){
   chomp;
   my $head=$_=~ m/^\#/;
   if(!$head){
      my($chr,$source,$class,$start,$end,$a,$strand,$b,$anotion)=split/\t/,$_;
      if($class eq "gene"){
      $anotion=~ m/gene_id \"(?<gene_id>\w+)\".+gene_name \"(?<gene_name>\w+)\".+gene_biotype \"(?<gene_biotype>\w+)\"/;
      print "$+{gene_id}" if($+{gene_id});
      print "\t$+{gene_name}" if($+{gene_name});
      print "\t$+{gene_biotype}" if($+{gene_biotype});
      print "\n" if($+{gene_id} || $+{gene_name} || $+{gene_biotype});
      }
   }
}
close GTF;

