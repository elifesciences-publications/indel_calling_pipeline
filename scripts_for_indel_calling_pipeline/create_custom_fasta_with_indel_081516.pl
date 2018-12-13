#!/usr/bin/perl -w
BEGIN {$^W=0}
use strict;
use List::Util qw(max min);

my ( $cur_indel_type, $cur_indel_seq, $cur_header_line, $cur_sequence)=("","","","");

 
 my ($file, $size_window, $indels_file_for_grep, $prefix)=@ARGV;
 my $to_print="";
 if ($#ARGV<3) {
   die ("usage: ./create_custom_fasta_with_indel.pl fasta_file input_files size_window\n");

 }
 
 open(INPUT, "<".$file) ||die ("Can't open $file\n");
 open(OUTPUT, ">$prefix");
 while (my $line=<INPUT>) {
    chomp $line;
    
    if ( $line=~m/^>.+/) { ##header line
      
       if ($cur_header_line eq ""){ ##if empty, then we're at 1st header line, initialize cur_header_line variable
        $cur_header_line=$line;
       } else { ##header line not empty, which means we're done processing the current record
        
         print OUTPUT $cur_header_line."\n";
        ($cur_indel_seq, $cur_indel_type)=getIndelInfo($cur_header_line, $indels_file_for_grep);
#print $cur_header_line."\t".$cur_sequence."\n";
        my @seq_arr=split(//,$cur_sequence); ##current sequence around indel
        my @new_sequence_arr=();
        my $to_print="";
    
        if ($cur_indel_type =~m/ins/) {
            #code
            $to_print=join("",(@seq_arr[0..($size_window-1)],$cur_indel_seq,@seq_arr[$size_window..$#seq_arr]));
         #  print $cur_header_line."\t".$to_print."\n";
            @new_sequence_arr=split(//, $to_print);
           my $max_size=min((scalar @new_sequence_arr-1),(($size_window+1)*2-1+length($cur_indel_seq)));
            $to_print=join("",@new_sequence_arr[0..$max_size]);
        }elsif($cur_indel_type =~m/del/){
            $to_print=join("",(@seq_arr[0..($size_window-1)],@seq_arr[($size_window+length($cur_indel_seq))..$#seq_arr]));
            @new_sequence_arr=split(//, $to_print);
            $to_print=join("",@new_sequence_arr[0..(($size_window+1)*2-1)]);
        }
        print OUTPUT uc $to_print."\n";
         
         
         
         
         ###after modifying current sequence with indel sequence
         $cur_header_line=$line;
         $cur_sequence="";
        
       }
        
        
    }else{ ## keep reading the file until next header sequence, concatenating sequences
     
        $cur_sequence=$cur_sequence.$line;
      
      
    }
    
}


sub getIndelInfo{
  
  my $header_line=$_[0];
  my $indels_file_for_grep=$_[1];
      $header_line=~s/>//g;
      my $grepped_seq=`grep -w  $header_line $indels_file_for_grep | awk '{print \$2}' `;
      my $indel_type="";
        
        if ($grepped_seq =~m/\+/) {
           $indel_type="ins";
           $grepped_seq=~s/\+//g;
            
            
        }elsif ($grepped_seq=~m/\-/){
            
            $indel_type="del";
           $grepped_seq=~s/\-//g;
                 
            
        }
        chomp $grepped_seq;
  return ($grepped_seq, $indel_type);
  
  
}