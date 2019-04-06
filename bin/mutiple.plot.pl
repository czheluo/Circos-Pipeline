#!/usr/bin/perl
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($snp,$indel,$sv,$cnv,$chrlist,$gff,$outdir,$paralist,$ssr,$snpplusindel,$outfile);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use List::Util qw(max min);

	if($para[3] eq "cnv" && defined $para[5]){######plot cnv
		open In,$cnv;
		my ($file_name,undef)=split/\./,$cnv,2;
		my $max;    
		my %hash_max;
	#print $file_name;
		open OUT,">$outdir/draw.circos/windows.file/$file_name.cnv.txt";
		while(<In>){
			chomp;
			next if (/^#/ || "");
			my ($chr,$type,$region,$possize,undef,$pval,undef)=split/\s+/,$_,6;
		#print $_;die;
			my ($chr,$pos)=split/\:/,$region;
			my ($start,$end)=split/\-/,$pos;
			if($pval==0){
				$pval=1;			
			}
        #my $log=-(log($pval)/log(10));
        #$hash_max{$chr}=$log;
			$hash_max{$chr}=$possize;
			next if (!exists $hash_chr_num{$chr});
			if ($type eq "deletion"){
				print OUT "$chr\t$start\t$end\t$possize\tfill_color=red\n";
			}else{
				print OUT "$chr\t$start\t$end\t$possize\tfill_color=black\n";
			}
     #print Dumper \%hash;die;
		}
		close In;
		close OUT;
		$max= (sort{$b<=>$a} values %hash_max)[0];
		if ($para[4] eq "line"){
			$main_conf.="
        ######plot cnv
        <plot>
        type = line
        max_gap = 1u
        file    = $outdir/draw.circos/windows.file/$file_name.cnv.txt
        color   = white
        min     = 0
        max     = $max
        r1      = split(/\,/,$para[2])[0]
        r0      = split(/\,/,$para[2])[1]
        thickness = 1
        #fill_color = vdyellow
        </plot>
        ######
			";
		}elsif($para[4] eq "scatter"){
			$main_conf.="
        ######plot cnv
        <plot>
        type=scatter
        file    = $outdir/draw.circos/windows.file/$file_name.cnv.txt
        #fill_color=green
        stroke_color=blue
        glyph=circle#rectangle
        glyph_size=10
        max=$max#0.013
        min=0
        r1      = split(/\,/,$para[2])[0]
        r0      = split(/\,/,$para[2])[1]
        </plot>
        ######
			";
		}elsif($para[4] eq "heatmap"){
			$main_conf.="
        ######plot cnv
        <plot>
        type=heatmap
        file    = $outdir/draw.circos/windows.file/$file_name.cnv.txt
        color   = col1,col2,col3,col4,col5
        r1      = split(/\,/,$para[2])[0]
        r0      = split(/\,/,$para[2])[1]
        </plot>
        ######
			";
		}else{
			$main_conf.="
        ######plot cnv
        <plot>
        type=histogram
        file    = $outdir/draw.circos/windows.file/$file_name.cnv.txt
        color   = white
        r1      = split(/\,/,$para[2])[0]
        r0      = split(/\,/,$para[2])[1]
        #orientation = out
        </plot>
        ######
			";
		}
	}
	if($para[3] eq "sv" && defined $para[5]){######plot sv
		open In,$sv;
		my ($file_name,undef)=split/\./,$sv,2;
		my $max;    
		my %hash_max;
		open OUT,">$outdir/draw.circos/windows.file/$file_name.sv.txt";
		while(<In>){
			chomp;
			next if (/^#/ || "");
			my ($chr,$start,undef,undef,$end,undef,$type,$ssize,$score,$numreads,undef)=split/\s+/,$_,11;
			$hash_max{$chr}=$numreads;#abs($ssize);
			next if (!exists $hash_chr_num{$chr});
		#my $abso=abs($ssize)
			if ($type eq "DEL"){
				print OUT "$chr\t$start\t$end\t$numreads\tfill_color=red\n";
			}elsif($type eq "INV"){
				print OUT "$chr\t$start\t$end\t$numreads\tfill_color=blue\n";
			}elsif($type eq "ITX"){
				print OUT "$chr\t$start\t$end\t$numreads\tfill_color=green\n";
			}else{
				print OUT "$chr\t$start\t$end\t$numreads\tfill_color=yellow\n";
			}
     #print Dumper \%hash;die;
		}
		close In;
		close OUT;
		$max= (sort{$b<=>$a} values %hash_max)[0];
		if ($para[4] eq "line"){
			$main_conf.="
        ######plot sv
        <plot>
        type = line
        max_gap = 1u
        file    = $outdir/draw.circos/windows.file/$file_name.sv.txt
        color   = white
        min     = 0
        max     = $max
        r1      = split(/\,/,$para[2])[0]
        r0      = split(/\,/,$para[2])[1]
        thickness = 1
        #fill_color = vdyellow
        </plot>
        ######
			";
		}elsif($para[4] eq "scatter"){
			$main_conf.="
        ######plot sv
        <plot>
        type=scatter
        file    = $outdir/draw.circos/windows.file/$file_name.sv.txt
        #fill_color=green
        stroke_color=blue
        glyph=circle#rectangle
        glyph_size=10
        max=$max#0.013
        min=0
        r1      = split(/\,/,$para[2])[0]
        r0      = split(/\,/,$para[2])[1]
        </plot>
        ######
			";
		}elsif($para[4] eq "heatmap"){
			$main_conf.="
        ######plot sv
        <plot>
        type=heatmap
        file    = $outdir/draw.circos/windows.file/$file_name.sv.txt
        color   = col1,col2,col3,col4,col5
        r1      = split(/\,/,$para[2])[0]
        r0      = split(/\,/,$para[2])[1]
        </plot>
        ######
			";
		}else{
			$main_conf.="
        ######plot sv
        <plot>
        type=histogram
        file    = $outdir/draw.circos/windows.file/$file_name.sv.txt
        color   = white
        r1      = split(/\,/,$para[2])[0]
        r0      = split(/\,/,$para[2])[1]
        #orientation = out
        </plot>
        ######
			";
		}
	}	


	if($para[3] eq "ssr" && defined $para[5]){######plot ssr
		open In,$ssr;
		my ($file_name,undef)=split/\./,$ssr,2;
		my $max;    
		my %hash_max;
	#print $file_name;
		open OUT,">$outdir/draw.circos/windows.file/$file_name.ssr.txt";
		while(<In>){
			chomp;
			next if (/^#/ || "");
			my ($chr,undef,undef,undef,$ssrsize,$start,$end,undef)=split/\s+/,$_,8;
			$hash_max{$chr}=$ssrsize;
			next if (!exists $hash_chr_num{$chr});
			print OUT "$chr\t$start\t$end\t$ssrsize\tfill_color=black\n";
		#print Dumper \%hash;die;
		}
		close In;
		close OUT;
		$max= (sort{$b<=>$a} values %hash_max)[0];
		if ($para[4] eq "line"){
			$main_conf.="
        ######plot ssr
        <plot>
        type = line
        max_gap = 1u
        file    = $outdir/draw.circos/windows.file/$file_name.ssr.txt
        color   = white
        min     = 0
        max     = $max
        r1      = split(/\,/,$para[2])[0]
        r0      = split(/\,/,$para[2])[1]
        thickness = 1
        #fill_color = vdyellow
        </plot>
        ######
			";
		}elsif($para[4] eq "scatter"){
			$main_conf.="
        ######plot ssr
        <plot>
        type=scatter
        file    = $outdir/draw.circos/windows.file/$file_name.ssr.txt
        #fill_color=green
        stroke_color=blue
        glyph=circle#rectangle
        glyph_size=10
        max=$max#0.013
        min=0
        r1      = split(/\,/,$para[2])[0]
        r0      = split(/\,/,$para[2])[1]
        </plot>
        ######
			";
		}elsif($para[4] eq "heatmap"){
			$main_conf.="
        ######plot cnv
        <plot>
        type=heatmap
        file    = $outdir/draw.circos/windows.file/$file_name.ssr.txt
        color   = col1,col2,col3,col4,col5
        r1      = split(/\,/,$para[2])[0]
        r0      = split(/\,/,$para[2])[1]
        </plot>
        ######
			";
		}else{
			$main_conf.="
        ######plot ssr
        <plot>
        type=histogram
        file    = $outdir/draw.circos/windows.file/$file_name.ssr.txt
        color   = white
        r1      = split(/\,/,$para[2])[0]
        r0      = split(/\,/,$para[2])[1]
        #orientation = out
        </plot>
        ######
			";
		}
	}

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub ABSOLUTE_DIR {
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	if(-f $in){
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}elsif(-d $in){
		chdir $in;$return=`pwd`;chomp $return;
	}else{
		warn "Warning just for file and dir \n$in";
		exit;
	}
	chdir $cur_dir;
	return $return;
}

sub USAGE {#
        my $usage=<<"USAGE";
Contact:        meng.luo\@majorbio.com;
Script:			$Script
Description:
	draw circos
	eg:
	perl $Script 

Usage:
  Options:
 -h         Help

USAGE
        print $usage;
        exit;
}

