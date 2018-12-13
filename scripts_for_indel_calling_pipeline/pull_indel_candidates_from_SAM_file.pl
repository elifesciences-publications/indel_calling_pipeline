#!/usr/bin/perl -w
use strict;


while (my $line=<STDIN>) {
    chomp $line;
    my @liner=split(/\t/, $line);
    my ($chr, $start, $cigar_string, $sequence)=@liner[2,3,5,9];
    my @sequence_arr=split(//, $sequence);
    my $cur_start_idx=0;
    my $cur_start_offset=0;
    ##parse cigar string
    my $counter=0;
    while ($cigar_string =~ m/(\d+)[MDI]/g) {
       my $length=$1;
       $counter++;
      # print $&."\t$cur_start_idx\n";
   if ($& =~/M/) { #MATCH
    #do nothing
     $cur_start_idx+=$length;
     $cur_start_offset+=$length;
     $cur_start_offset-- if $counter==1;
   }elsif($&=~/D/){ #DELETION
    $cur_start_offset-- if $counter==1;
    print "$chr\t".($start+$cur_start_offset)."\t-".("N"x $length)."\n";
     $cur_start_offset+=$length;
   }elsif($& =~/I/){ #INSERTION
        $cur_start_offset-- if $counter==1;
    print "$chr\t".($start+$cur_start_offset)."\t+".join("",@sequence_arr[$cur_start_idx..($cur_start_idx+$length-1)])."\n";
$cur_start_idx+=$length;
#$cur_start_offset++;
   }
  
}

}





