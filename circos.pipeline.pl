#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($snp,$indel,$sv,$cnv,$chrlist,$gff,$outdir,$paralist,$ssr,$snpplusindel,$outfile,$chromosome);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use List::Util qw(max min);
my $version="1.1.0";
GetOptions(
	"help|?" =>\&USAGE,
	"gff:s"=>\$gff,
	"chrlist:s"=>\$chrlist,
	"outdir:s"=>\$outdir,
	"paramlist:s"=>\$paralist,
	"outfile:s"=>\$outfile,
	"chrom:s"=>\$chromosome,
			) or &USAGE;
&USAGE unless ($gff and $chrlist and $paralist);
#$plottype||="line scatter histogram heatmap";#stuckbar
$outfile||="circos";
$outdir||="./";
########
########work dir
########
my $mkdir=1;
$mkdir=(mkdir "$outdir") if (!-d "$outdir");
$outdir=ABSOLUTE_DIR($outdir);
die "Error make dir $outdir" if ($mkdir == 0);
$mkdir=(mkdir "$outdir/draw.circos") if (!-d "$outdir/draw.circos");
die "Error make dir $outdir/draw.circos" if ($mkdir == 0);
$mkdir=(mkdir "$outdir/draw.circos/windows.file") if (!-d "$outdir/draw.circos/windows.file");
die "Error make dir $outdir/draw.circos/windows.file" if ($mkdir == 0);
open MC,">","$outdir/draw.circos/draw.conf";
open IDEO,">","$outdir/draw.circos/ideogram.conf";
open TIC,">","$outdir/draw.circos/ticks.conf";
open IM,">","$outdir/draw.circos/chromosomes.and.color.conf";
open HS,">","$outdir/draw.circos/housekeeping.conf";
my ($main_conf,$ideogram,$ticks,$chro,$housekeeping);
########
########chr and band file
########
open CHR,$chrlist;
open GFF,$gff;
open CHRBAND,">$outdir/draw.circos/chr.band.txt";
my %hash_chr_length;
my %hash_gff;
my %hash_chr_num;
my $chr_type;
while(<CHR>){
    #$_=~s/sca//g;
    $_=~s/[\n\r]//g;
	next if ($_ eq ""||/^$/||/^#/);
    my ($chr,$length)=split;
	if ($chr =~ /chr/) {
		$chr_type="chr";
	}else{
		$chr_type="sca";
	}
    $hash_chr_num{$chr}=1;
    $hash_chr_length{$chr}=$length;
}
my $chrnum=scalar keys %hash_chr_num;
#open OUTT,">win.number.txt";
my $gff_windows=100000;
while(<GFF>){
    next if /^#/;
    $_=~s/[\n\r]//g;
    my ($chr,undef,$type,$start,@others)=split(/\t/,$_);
    my $win_num=int($start/$gff_windows)+1;
      #print OUTT $win_num."\n";
    if (exists $hash_chr_num{$chr} and ($type=~/gene/ || ($type !~ /CDS/ && $type !~/exon/)) ){
        if ($chr =~ /chr/) {$chr=~s/chr//g;}
        if ($chr =~ /sca/) {$chr=~s/sca//g;}
        $hash_gff{$chr}{$win_num}++;
        }
}

##draw colour of band
my %region;
my %bands;
open CO,">$outdir/draw.circos/band.txt";
my @chr_windows;
foreach my $keys (sort {$a<=>$b}keys %hash_gff){
    my $band=1;
    my $max=0;
    foreach my $num (sort keys %{$hash_gff{$keys}}){
        my $start=1+$gff_windows*($num-1);
        my $end=$gff_windows*$num; 
        my $snp=$hash_gff{$keys}{$num}/$gff_windows;
        $region{$chr_type}{join("\t",$start,$end)}=1;
        push(@chr_windows,$snp);
        if($max < $num){
            $max=$num;
            }
        print CO "band\t$chr_type"."$keys\tband$band\tband$band\t$start\t$end\t$snp\n";
        $band++;
    }
    my $bandnum=$keys;
    $bandnum=~s/$chr_type//g;
    print CHRBAND "chr"."\t"."-"."\t".$chr_type."".$keys."\t".$bandnum."\t"."0"."\t".$max*$gff_windows."\t".$chr_type."".$keys."\n";
}
close CO;
open BANDCOL,"<$outdir/draw.circos/band.txt";
my $band_windows=(max(@chr_windows)-min(@chr_windows))/10;
while(<BANDCOL>){
    $_=~s/[\n\r]//g;
    my @array=split;
    my $win_num=int($array[6]/$band_windows)+1;
    $array[5]=$array[5]-1;
    print CHRBAND "$array[0]\t$array[1]\t$array[2]\t$array[3]\t$array[4]\t$array[5]\tbandcol$win_num\n";
   # print CHRBAND join("\t",@array[1..($#array-1)])."\t"."bandcol".$win_num."\n";
}

close CHR;
close GFF;
close CHRBAND;
########
########
########
#my $drawchr=join("",(1..$chrnum));
my @draw_chr=sort {$hash_chr_num{$a} <=> $hash_chr_num{$b}} keys %hash_chr_num;
my $chr_draw=join(";",@draw_chr);
if ($chromosome){
	my $chroms=join(";",split(",",$chromosome));
	$chro.="
	chromosomes_units = 1000000
	chromosomes_display_default = no
	chromosomes=$chroms #chromosomes = hs1:(-100,120-);hs2;hs3;h4
	#chromosomes_order =$chr_draw
	<colors>
	";
}else{
	$chro.="
	chromosomes_units = 1000000
	chromosomes_display_default = yes
	chromosomes=$chr_draw #chromosomes = hs1:(-100,120-);hs2;hs3;h4
	#chromosomes_order =$chr_draw
	<colors>	";
}
my %chr_colour;
my $num;
for(my $i=0;$i<@draw_chr;$i++){
    my $num=$i%6;
    #print $num;
    if($i%6==0){
        $chro.="$draw_chr[$i]* = col6\n";
        $chr_colour{$draw_chr[$i]}="col6";
        }
    else{
        $chro.="$draw_chr[$i]* = col$num\n";
        $chr_colour{$draw_chr[$i]}="col".$num;
        }
    }
$chro.="
</colors>
";
#print $num\t$chr_draw;die;
#########
#########
#########
if (-e "$outdir/draw.circos/chr.band.txt"){
$main_conf.="
karyotype = $outdir/draw.circos/chr.band.txt

<<include $outdir/draw.circos/ideogram.conf>>
<<include $outdir/draw.circos/ticks.conf>>
<<include $outdir/draw.circos/chromosomes.and.color.conf>>
<image>
<<include etc/image.conf>>
radius* = 2000
</image>
";
}
########
########
########
$ideogram="
<ideogram>
<spacing>
default = 0.005r
</spacing>
radius    = 0.9r
thickness = 100p
fill      = yes
fill_color = black
#stroke_thickness = 2
#stroke_color     = black
show_label       = yes
label_font       = default 
label_radius     = 1r + 75p
label_size       = 30
label_parallel   = yes
label_case     = lower#upper 
label_format   = eval(sprintf(\"$chr_type%s\",var(label)))
show_bands            = yes
fill_bands            = yes
band_stroke_thickness = 1
band_stroke_color     = black
band_transparency     = 2
#label_with_tag = no
</ideogram>
";
#########
#########
#########
$ticks="
show_ticks          = no
show_tick_labels    = no
<ticks>
radius           = 1r
color            = black
thickness        = 2p
multiplier       = 1e-6
format           = %d
<tick>
spacing        = 5u
size           = 10p
</tick>
<tick>
spacing        = 25u
size           = 15p
show_label     = yes
label_size     = 20p
label_offset   = 10p
format         = %d
</tick>

</ticks>
";
$housekeeping.="
anglestep       = 0.5
minslicestep    = 10
beziersamples   = 40
debug           = no
warnings        = no
imagemap        = no
paranoid        = no
units_ok        = bupr
units_nounit    = n
file_delim = \\s
file_delim_collapse = yes
list_record_delim = \\s*[;,]\\s*
list_field_delim  = \\s*[:=]\\s*
options_record_delim = [,;]
options_field_delim  = =
skip_missing_expression_vars = no
legacy_underline_expression_syntax = no
svg_font_scale = 1.3
sup_baseline_shift = 40
sub_baseline_shift = -40
sup_fontsize = 90
sub_fontsize = 90
default_font   = default
default_font_name  = Arial
default_font_color = black
default_color  = black
<guides>
thickness      = 1
size           = 5
type           = outline
<object>
all            = no
ideogram       = no
ideogram_label = no
</object>
<color>
default = lblue
text    = red
</color>
</guides>
debug_group = summary,output
debug_auto_timer_report = 3600
debug_word_separator = \" \"
debug_undef_text     = _undef_
debug_empty_text     = _emptylist_
debug_validate       = yes
debug_output_tidy    = no
text_pixel_subsampling = 1
text_snuggle_method    = array
restrict_parameter_names = no
case_sensitive_parameter_names = no
calculate_track_statistics = yes
color_cache_static = yes
color_cache_file   = circos.colorlist
color_lists_use    = yes
memoize = yes
quit_on_dump = yes
offsets = 0,0
max_ticks            = 5000
max_ideograms        = 200
max_links            = 25000000000000000000
max_points_per_track = 25000000000000000000
undefined_ideogram = skip
relative_scale_iterations = 10
relative_scale_spacing    = mode
data_out_of_range = trim,warn
track_defaults = etc/tracks
round_brush_use           = yes
round_brush_min_thickness = 5
anti_aliasing = yes
housekeeping = yes
auto_eval = no
";
#########
$main_conf.="
<plots>
<backgrounds>
show  = data
<background>
color = vvlgrey
y0    = 1.0r
y1    = 0r
</background>
</backgrounds>
";
open PA,$paralist;
my @para;
my  $nlm=1;
while (<PA>) {
	chomp;
	next if ($_=~ "#");
	@para=split/\s+/,$_;
	my($r1,$r0)=split/\,/,$para[2];
	#print "$para[3]\t$para[5]";die;
	if($para[3] eq "snp" && $para[5] ne "NA"){######plot snp
		open In,$para[7];
		if ($snp=~/.gz$/) {
			close In;
			open In,"gunzip -c $snp|";
		}
		my $nsam;
		my @sample;
		my $max;    
		my %hash;
		my $win=$para[6];
		my %hash_max;
		while(<In>){
			chomp;
			next if(/^##/);
			if (/^#/){
				my (undef,undef,undef,undef,undef,undef,undef,undef,undef,@sample)=split/\s+/,$_;
				#print Dumper @sample;
				my $sam=9;
				for(my $i=0;$i<@sample;$i++){
					if ($para[5] eq $sample[$i]){
						$nsam=$sam+$i;
						last;
					}else{
						next;
					}
				}
			}
			next if(/^#/);
			my @array=split;
			my ($gt,undef)=split/\:/,$array[$nsam],2;
			next if ($gt eq "./.");
			my $win_num=int($array[1]/$win)+1;
			next if (!exists $hash_chr_num{$array[0]});
			$hash{$array[0]}{$win_num}++;
			#print Dumper \%hash;die;
		}
		close In;
		#my $file_name=$para[5];
		open OUT,">$outdir/draw.circos/windows.file/$para[5].snp.$nlm.win.txt";
		foreach my $keys (sort keys %hash){
		  # print "$keys\n";
			foreach my $num (sort keys %{$hash{$keys}}){
				my $chr=$keys;
				my $start=1+$win*($num-1);
				my $end=$win*$num;
			    my $snp=$hash{$keys}{$num}/$win;
			    $hash_max{$snp}=1;
				print OUT "$chr\t$start\t$end\t$snp\tfill_color=$chr_colour{$chr}\n"; 
			}
		} 
		#print \%hash_max;die;
		$max= (sort{$b<=>$a} keys %hash_max)[0];
	 #print $max;
		close OUT; 
		##if line or heatmap 
		if ($para[4] eq "line"){
			$main_conf.="
			######plot snp
			<plot>
			type = line
			max_gap = 1u
			file    = $outdir/draw.circos/windows.file/$para[5].snp.$nlm.win.txt
			color   = white
			min     = 0
			max     = $max#0.015
			r1      = $r1
			r0      = $r0
			thickness = 1
			#fill_color = vdyellow
			</plot>
        ######
			";
		}elsif($para[4] eq "scatter"){
			$main_conf.="
			######plot snp
			<plot>
			type=scatter
			file    = $outdir/draw.circos/windows.file/$para[5].snp.$nlm.win.txt
			#fill_color=green
			stroke_color=blue
			glyph=circle#rectangle
			glyph_size=10
			max=$max#0.013
			min=0
			r1      = $r1
			r0      = $r0
			</plot>
			######
			";
		}elsif($para[4] eq "heatmap"){
			$main_conf.="
			######plot snp
			<plot>
			type=heatmap
			file    = $outdir/draw.circos/windows.file/$para[5].snp.$nlm.win.txt
			color   = col1,col2,col3,col4,col5
			r1      = $r1
			r0      = $r0
			</plot>
			######
			";
		}else{
			$main_conf.="
			######plot snp
			<plot>
			type=histogram
			file    = $outdir/draw.circos/windows.file/$para[5].snp.$nlm.win.txt
			color   = white
			r1      = $r1
			r0      = $r0
			#orientation = out
			</plot>
			######
			";
		}
	}elsif($para[3] eq "snp" && $para[5] eq "NA"){
		open In,$para[7];
		if ($snp=~/.gz$/) {
			close In;
			open In,"gunzip -c $snp|";
		}
		my $max;    
		my %hash;
		my $win=$para[6];
		my %hash_max;
		#open OUT,">$out";
		while(<In>){
			#$_=~s/sca/chr/g;
			 next if (/^#/);
		 #$_=~s/[\n\r]//g;
			#print $_;
			my @array=split;
			my $win_num=int($array[1]/$win)+1;
			next if (!exists $hash_chr_num{$array[0]});
			$hash{$array[0]}{$win_num}++;
		}
		close In;
		#my $file_name=(split/\//,$snp)[-1];
		open OUT,">$outdir/draw.circos/windows.file/$para[5].snp.$nlm.win.txt";
		foreach my $keys (sort keys %hash){
			foreach my $num (sort keys %{$hash{$keys}}){
			 my $chr=$keys;
				my $start=1+$win*($num-1);
				my $end=$win*$num;
				my $snp=$hash{$keys}{$num}/$win;
				$hash_max{$snp}=1;
				print OUT "$chr\t$start\t$end\t$snp\tfill_color=$chr_colour{$chr}\n"; 
			}
		} 
		$max= (sort{$b<=>$a} keys %hash_max)[0];
		#print $max;
		close OUT; 
		##if line or heatmap 
		if ($para[4] eq "line"){
			$main_conf.="
        ######plot snp
        <plot>
        type = line
        max_gap = 1u
        file    = $outdir/draw.circos/windows.file/$para[5].snp.$nlm.win.txt
        color   = white
        min     = 0
        max     = $max#0.015
        r1      = $r1
        r0      = $r0
        thickness = 1
        #fill_color = vdyellow
        </plot>
        ######
			";
		}elsif($para[4] eq "scatter"){
			$main_conf.="
        ######plot snp
        <plot>
        type=scatter
        file    = $outdir/draw.circos/windows.file/$para[5].snp.$nlm.win.txt
        #fill_color=green
        stroke_color=blue
        glyph=circle#rectangle
        glyph_size=10
        max=$max#0.013
        min=0
        r1      = $r1
        r0      = $r0
        </plot>
        ######
			";
		}elsif($para[4] eq "heatmap"){
			$main_conf.="
        ######plot snp
        <plot>
        type=heatmap
        file    = $outdir/draw.circos/windows.file/$para[5].snp.$nlm.win.txt
        color   = col1,col2,col3,col4,col5
        r1      = $r1
        r0      = $r0
        </plot>
        ######
			";
		}else{
			$main_conf.="
        ######plot snp
        <plot>
        type=histogram
        file    = $outdir/draw.circos/windows.file/$para[5].snp.$nlm.win.txt
        color   = white
        r1      = $r1
        r0      = $r0
        #orientation = out
        </plot>
        ######
			";
		}
	}elsif($para[3] eq "indel" && $para[5] ne "NA"){######plot indel
		#print $para[5];die;
		my $nsam;
		my @sample;
		my $max;    
		my %hash;
		my $win=$para[6];
		my %hash_max;
		open In,$para[7];
		while(<In>){
			chomp;
			next if(/^##/);
			#my $nsam;
			if (/^#/){
				#print "$_";die;
				my (undef,undef,undef,undef,undef,undef,undef,undef,undef,@sample)=split/\s+/,$_;
				#print Dumper @sample;
				my $sam=9;
				for(my $i=0;$i<@sample;$i++){
					if ($para[5] eq $sample[$i]){
						#print $sample[$i];
						$nsam=$sam+$i;
						#print $nsam;
						last;
					}else{
						next;
					}
				}
			}
			next if(/^#/);
			my @array=split;
			my ($gt,undef)=split/\:/,$array[$nsam],2;
			next if ($gt eq "./.");
			my $win_num=int($array[1]/$win)+1;
			next if (!exists $hash_chr_num{$array[0]});
			$hash{$array[0]}{$win_num}++;
		#print Dumper \%hash;die;
		}
		close In;
		open OUT,">$outdir/draw.circos/windows.file/$para[5].indel.$nlm.win.txt";
		foreach my $keys (sort keys %hash){
		 # print "$keys\n";
			foreach my $num (sort keys %{$hash{$keys}}){
				my $chr=$keys;
				my $start=1+$win*($num-1);
				my $end=$win*$num;
				my $indel=$hash{$keys}{$num}/$win;
				$hash_max{$indel}=1;
				print OUT "$chr\t$start\t$end\t$indel\tfill_color=$chr_colour{$chr}\n"; 
			}
		} 
		#print \%hash_max;die;
		$max= (sort{$b<=>$a} keys %hash_max)[0];
		#print $max;
		close OUT; 
		if ($para[4] eq "line"){
			$main_conf.="
		 ######plot indel
        <plot>
        type = line
        max_gap = 1u
        file    = $outdir/draw.circos/windows.file/$para[5].indel.$nlm.win.txt
        color   = white
        min     = 0
        max     = $max#0.015
        r1      = $r1
        r0      = $r0
        thickness = 1
        #fill_color = vdyellow
        </plot>
        ######
			";
		}elsif($para[4] eq "scatter"){
			$main_conf.="
			######plot indel
        <plot>
        type=scatter
        file    = $outdir/draw.circos/windows.file/$para[5].indel.$nlm.win.txt
        #fill_color=green
        stroke_color=blue
        glyph=circle#rectangle
        glyph_size=10
        max=$max#0.013
        min=0
        r1      = $r1
        r0      = $r0
        </plot>
        ######
			";
		}elsif($para[4] eq "heatmap"){
			$main_conf.="
        ######plot indel
        <plot>
        type=heatmap
        file    = $outdir/draw.circos/windows.file/$para[5].indel.$nlm.win.txt
        color   = col1,col2,col3,col4,col5
        r1      = $r1
        r0      = $r0
        </plot>
        ######
			";
		}else{
			$main_conf.="
        ######plot indel
        <plot>
        type=histogram
        file    = $outdir/draw.circos/windows.file/$para[5].indel.$nlm.win.txt
        color   = white
        r1      = $r1
        r0      = $r0
        #orientation = out
        </plot>
        ######
			";
		}
	}elsif($para[3] eq "indel" && $para[5] eq "NA"){
		open In,$para[7];
		if ($indel=~/.gz$/) {
			close In;
			open In,"gunzip -c $indel|";
		}
		my $max;    
		my %hash;
		my $win=$para[6] ;
		my %hash_max;
		#open OUT,">$out";
		while(<In>){
			#$_=~s/sca/chr/g;
			next if (/^#/);
        #$_=~s/[\n\r]//g;
        #print $_;
			my @array=split;
			my $win_num=int($array[1]/$win)+1;
			next if (!exists $hash_chr_num{$array[0]});
			$hash{$array[0]}{$win_num}++;
			#print Dumper \%hash;die;
		}
		close In;
		open OUT,">$outdir/draw.circos/windows.file/indel.$nlm.win.txt";
		foreach my $keys (sort keys %hash){
        # print "$keys\n";
			foreach my $num (sort keys %{$hash{$keys}}){
				my $chr=$keys;
				my $start=1+$win*($num-1);
				my $end=$win*$num;
				my $indel=$hash{$keys}{$num}/$win;
				$hash_max{$indel}=1;
				print OUT "$chr\t$start\t$end\t$indel\tfill_color=$chr_colour{$chr}\n"; 
			}
		} 
    #print \%hash_max;die;
		$max= (sort{$b<=>$a} keys %hash_max)[0];
    #print $max;
		close OUT; 
		if ($para[4] eq "line"){
			$main_conf.="
        ######plot indel
        <plot>
        type = line
        max_gap = 1u
        file    = $outdir/draw.circos/windows.file/indel.$nlm.win.txt
        color   = white
        min     = 0
        max     = $max#0.015
        r1      = $r1
        r0      = $r0
        thickness = 1
        #fill_color = vdyellow
        </plot>
        ######
			";
		}elsif($para[4] eq "scatter"){
			$main_conf.="
        ######plot indel
        <plot>
        type=scatter
        file    = $outdir/draw.circos/windows.file/indel.$nlm.win.txt
        #fill_color=green
        stroke_color=blue
        glyph=circle#rectangle
        glyph_size=10
        max=$max#0.013
        min=0
        r1      = $r1
        r0      = $r0
        </plot>
        ######
        ";
		}elsif($para[4] eq "heatmap"){
			$main_conf.="
        ######plot indel
        <plot>
        type=heatmap
        file    = $outdir/draw.circos/windows.file/indel.$nlm.win.txt
        color   = col1,col2,col3,col4,col5
        r1      = $r1
        r0      = $r0
        </plot>
        ######
			";
		}else{
			$main_conf.="
        ######plot indel
        <plot>
        type=histogram
        file    = $outdir/draw.circos/windows.file/indel.$nlm.win.txt
        color   = white
        r1      = $r1
        r0      = $r0
        #orientation = out
        </plot>
        ######
			";
		}
	}elsif ($para[3] eq "snpplusindel" && $para[5] ne "NA"){######plot $snpplusindel
		my $nsam;
		my @sample;
		my $max;
		my %hash;
		my $win=$para[6];
		my %hash_max;
		open In,$para[7];
		if ($snpplusindel=~/.gz$/) {
			close In;
			open In,"gunzip -c $snpplusindel|";
		}
		while(<In>){
			chomp;
			next if(/^##/);
			#my $nsam;
			if (/^#/){
				#print "$_";die;
				my (undef,undef,undef,undef,undef,undef,undef,undef,undef,@sample)=split/\s+/,$_;
				#print Dumper @sample;
				my $sam=9;
				for(my $i=0;$i<@sample;$i++){
					if ($para[5] eq $sample[$i]){
						#print $sample[$i];
						$nsam=$sam+$i;
						#print $nsam;
						last;
					}else{
						next;
					}
				}
			}
			next if(/^#/);
			my @array=split;
			my ($gt,undef)=split/\:/,$array[$nsam],2;
			next if ($gt eq "./.");
			my $win_num=int($array[1]/$win)+1;
			next if (!exists $hash_chr_num{$array[0]});
			$hash{$array[0]}{$win_num}++;
     #print Dumper \%hash;die;
		}
		close In;
		open OUT,">$outdir/draw.circos/windows.file/$para[5].snpplusindel.$nlm.win.txt";
		foreach my $keys (sort keys %hash){
        # print "$keys\n";
			foreach my $num (sort keys %{$hash{$keys}}){
				my $chr=$keys;
				my $start=1+$win*($num-1);
				my $end=$win*$num;
				my $snpplusindel=$hash{$keys}{$num}/$win;
				$hash_max{$snpplusindel}=1;
				print OUT "$chr\t$start\t$end\t$snpplusindel\tfill_color=$chr_colour{$chr}\n"; 
			}
		} 
    #print \%hash_max;die;
		$max= (sort{$b<=>$a} keys %hash_max)[0];
    #print $max;
		close OUT; 
		if ($para[4] eq "line"){
			$main_conf.="
        ######plot indel
        <plot>
        type = line
        max_gap = 1u
        file    = $outdir/draw.circos/windows.file/$para[5].snpplusindel.$nlm.win.txt
        color   = white
        min     = 0
        max     = $max#0.015
        r1      = $r1
        r0      = $r0
        thickness = 1
        #fill_color = vdyellow
        </plot>
        ######
			";
		}elsif($para[4] eq "scatter"){
			$main_conf.="
        ######plot indel
        <plot>
        type=scatter
        file    = $outdir/draw.circos/windows.file/$para[5].snpplusindel.$nlm.win.txt
        #fill_color=green
        stroke_color=blue
        glyph=circle#rectangle
        glyph_size=10
        max=$max#0.013
        min=0
        r1      = $r1
        r0      = $r0
        </plot>
        ######
			";
		}elsif($para[4] eq "heatmap"){
			$main_conf.="
        ######plot indel
        <plot>
        type=heatmap
        file    = $outdir/draw.circos/windows.file/$para[5].snpplusindel.$nlm.win.txt
        color   = col1,col2,col3,col4,col5
        r1      = $r1
        r0      = $r0
        </plot>
        ######
			";
		}else{
			$main_conf.="
        ######plot indel
        <plot>
        type=histogram
        file    = $outdir/draw.circos/windows.file/$para[5].snpplusindel.$nlm.win.txt
        color   = white
        r1      = $r1
        r0      = $r0
        #orientation = out
        </plot>
        ######
			";
		}
	}elsif($para[3] eq "snpplusindel" && $para[5] eq "NA"){
		open In,$para[7];
		if ($snpplusindel=~/.gz$/) {
			close In;
			open In,"gunzip -c $snpplusindel|";
		}
		my $max;    
		my %hash;
		my $win=$para[6];
		my %hash_max;
    #open OUT,">$out";
		while(<In>){
        #$_=~s/sca/chr/g;
			next if (/^#/);
        #$_=~s/[\n\r]//g;
        #print $_;
			my @array=split;
			my $win_num=int($array[1]/$win)+1;
			next if (!exists $hash_chr_num{$array[0]});
			$hash{$array[0]}{$win_num}++;
     #print Dumper \%hash;die;
		}
		close In;
		open OUT,">$outdir/draw.circos/windows.file/snpplusindel.$nlm.win.txt";
		foreach my $keys (sort keys %hash){
			#print "$keys\n";
			foreach my $num (sort keys %{$hash{$keys}}){
				my $chr=$keys;
				my $start=1+$win*($num-1);
				my $end=$win*$num;
				my $snpplusindel=$hash{$keys}{$num}/$win;
				$hash_max{$snpplusindel}=1;
				print OUT "$chr\t$start\t$end\t$snpplusindel\tfill_color=$chr_colour{$chr}\n"; 
			}
		} 
    #print \%hash_max;die;
		$max= (sort{$b<=>$a} keys %hash_max)[0];
    #print $max;
		close OUT; 
		if ($para[4] eq "line"){
			$main_conf.="
			######plot indel
			<plot>
			type = line
			max_gap = 1u
			file    = $outdir/draw.circos/windows.file/snpplusindel.$nlm.win.txt
			color   = white
			min     = 0
			max     = $max#0.015
			r1      = $r1
			r0      = $r0
			thickness = 1
			#fill_color = vdyellow
			</plot>
			######
			";
		}elsif($para[4] eq "scatter"){
			$main_conf.="
			######plot indel
			<plot>
			type=scatter
			file    = $outdir/draw.circos/windows.file/snpplusindel.$nlm.win.txt
			#fill_color=green
			stroke_color=blue
			glyph=circle#rectangle
			glyph_size=10
			max=$max#0.013
			min=0
			r1      = $r1
			r0      = $r0
			</plot>
			######
			";
		}elsif($para[4] eq "heatmap"){
			$main_conf.="
			######plot indel
			<plot>
			type=heatmap
			file    = $outdir/draw.circos/windows.file/snpplusindel.$nlm.win.txt
			color   = col1,col2,col3,col4,col5
			r1      = $r1
			r0      = $r0
			</plot>
			######
				";
		}else{
			$main_conf.="
			######plot indel
			<plot>
			type=histogram
			file    = $outdir/draw.circos/windows.file/snpplusindel.$nlm.win.txt
			color   = white
			r1      = $r1
			r0      = $r0
			#orientation = out
			</plot>
			######
			";
		}
	}elsif($para[3] eq "cnv" && $para[5] ne "NA"){######plot cnv both same format
		open In,$para[7];
		my $max;    
		my %hash;
		my $win=$para[6];
		my %hash_max;
		while(<In>){
        #$_=~s/sca/chr/g;
			next if (/^#/);
        #$_=~s/[\n\r]//g;
        #print $_;
			my @array=split;
			my $win_num=int($array[1]/$win)+1;
			next if (!exists $hash_chr_num{$array[0]});
			$hash{$array[0]}{$win_num}++;
     #print Dumper \%hash;die;
		}
		close In;
		open OUT,">$outdir/draw.circos/windows.file/$para[5].$nlm.cnv.txt";
		foreach my $keys (sort keys %hash){
			#print "$keys\n";
			foreach my $num (sort keys %{$hash{$keys}}){
				my $chr=$keys;
				my $start=1+$win*($num-1);
				my $end=$win*$num;
				my $cnvs=$hash{$keys}{$num}/$win;
				$hash_max{$cnvs}=1;
				print OUT "$chr\t$start\t$end\t$cnvs\tfill_color=$chr_colour{$chr}\n"; 
			}
		} 
    #print \%hash_max;die;
		$max= (sort{$b<=>$a} keys %hash_max)[0];
    #print $max;
		close OUT; 
		if ($para[4] eq "line"){
			$main_conf.="
        ######plot cnv
        <plot>
        type = line
        max_gap = 1u
        file    = $outdir/draw.circos/windows.file/$para[5].$nlm.cnv.txt
        color   = white
        min     = 0
        max     = $max
        r1      = $r1
        r0      = $r0
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
        file    = $outdir/draw.circos/windows.file/$para[5].$nlm.cnv.txt
        #fill_color=green
        stroke_color=blue
        glyph=circle#rectangle
        glyph_size=10
        max=$max#0.013
        min=0
        r1      = $r1
        r0      = $r0
        </plot>
        ######
			";
		}elsif($para[4] eq "heatmap"){
			$main_conf.="
        ######plot cnv
        <plot>
        type=heatmap
        file    = $outdir/draw.circos/windows.file/$para[5].$nlm.cnv.txt
        color   = col1,col2,col3,col4,col5
        r1      = $r1
        r0      = $r0
        </plot>
        ######
			";
		}else{
			$main_conf.="
        ######plot cnv
        <plot>
        type=histogram
        file    = $outdir/draw.circos/windows.file/$para[5].$nlm.cnv.txt
        color   = white
        r1      = $r1
        r0      = $r0
        #orientation = out
        </plot>
        ######
			";
		}
	}elsif ($para[3] eq "cnv" && $para[5] eq "NA"){######plot cnv both same format
		open In,$para[7];
		my $max;    
		my %hash;
		my $win=$para[6];
		my %hash_max;
		while(<In>){
        #$_=~s/sca/chr/g;
			next if (/^#/);
        #$_=~s/[\n\r]//g;
        #print $_;
			my @array=split;
			my $win_num=int($array[1]/$win)+1;
			next if (!exists $hash_chr_num{$array[0]});
			$hash{$array[0]}{$win_num}++;
     #print Dumper \%hash;die;
		}
		close In;
		open OUT,">$outdir/draw.circos/windows.file/cnv.$nlm.win.txt";
		foreach my $keys (sort keys %hash){
			#print "$keys\n";
			foreach my $num (sort keys %{$hash{$keys}}){
				my $chr=$keys;
				my $start=1+$win*($num-1);
				my $end=$win*$num;
				my $cnvs=$hash{$keys}{$num}/$win;
				$hash_max{$cnvs}=1;
				print OUT "$chr\t$start\t$end\t$cnvs\tfill_color=$chr_colour{$chr}\n"; 
			}
		} 
    #print \%hash_max;die;
		$max= (sort{$b<=>$a} keys %hash_max)[0];
    #print $max;
		close OUT; 
		if ($para[4] eq "line"){
			$main_conf.="
			######plot cnv
			<plot>
			type = line
			max_gap = 1u
			file    = $outdir/draw.circos/windows.file/cnv.$nlm.win.txt
			color   = white
			min     = 0
			max     = $max
			r1      = $r1
			r0      = $r0
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
			file    = $outdir/draw.circos/windows.file/cnv.$nlm.win.txt
			#fill_color=green
			stroke_color=blue
			glyph=circle#rectangle
			glyph_size=10
			max=$max#0.013
			min=0
			r1      = $r1
			r0      = $r0
			</plot>
			######
			";
		}elsif($para[4] eq "heatmap"){
			$main_conf.="
			######plot cnv
			<plot>
			type=heatmap
			file    = $outdir/draw.circos/windows.file/cnv.$nlm.win.txt
			color   = col1,col2,col3,col4,col5
			r1      = $r1
			r0      = $r0
			</plot>
			######
			";
		}else{
			$main_conf.="
			######plot cnv
			<plot>
			type=histogram
			file    = $outdir/draw.circos/windows.file/cnv.$nlm.win.txt
			color   = white
			r1      = $r1
			r0      = $r0
			#orientation = out
			</plot>
			######
			";
		}
	}elsif($para[3] eq "sv" && $para[5] ne "NA"){######plot sv
		open In,$para[7];
		my $nsam;
		my @sample;
		my $max;
		my %hash;
		my $win=$para[6];
		my %hash_max;
		my $file_name;
		if ($sv =~ "vcf") {
			while(<In>){
				chomp;
				next if(/^##/);
				#my $nsam;
				if (/^#/) {
					#print "$_";die;
					my (undef,undef,undef,undef,undef,undef,undef,undef,undef,@sample)=split/\s+/,$_;
					#print Dumper @sample;
					my $sam=9;
					for(my $i=0;$i<@sample;$i++){
						if ($para[5] eq $sample[$i]){
							#print $sample[$i];
							$nsam=$sam+$i;
							#print $nsam;
							last;
						}else{
							next;
						}
					}
				}
				next if(/^#/);
				my @array=split;
				my ($gt,undef)=split/\:/,$array[$nsam],2;
				next if ($gt eq "./.");
				my $win_num=int($array[1]/$win)+1;
				next if (!exists $hash_chr_num{$array[0]});
				$hash{$array[0]}{$win_num}++;
			}
			close In;
			open OUT,">$outdir/draw.circos/windows.file/$para[5].sv.$nlm.win.txt";
			foreach my $keys (sort keys %hash){
        # print "$keys\n";
				foreach my $num (sort keys %{$hash{$keys}}){
					my $chr=$keys;
					my $start=1+$win*($num-1);
					my $end=$win*$num;
					my $snpplusindel=$hash{$keys}{$num}/$win;
					$hash_max{$snpplusindel}=1;
					print OUT "$chr\t$start\t$end\t$snpplusindel\tfill_color=$chr_colour{$chr}\n"; 
				}
			} 
			#print \%hash_max;die;
			$max= (sort{$b<=>$a} keys %hash_max)[0];
			#print $max;
			close OUT; 
		}else {
			while(<In>){
				chomp;
				next if(/^##/);
				#my $nsam;
				if (/^#/) {
					#print "$_";die;
					my (undef,undef,undef,undef,undef,undef,undef,undef,undef,@sample)=split/\s+/,$_;
					#print Dumper @sample;
					my $sam=9;
					for(my $i=0;$i<@sample;$i++){
						if ($para[5] eq $sample[$i]){
							#print $sample[$i];
							$nsam=$sam+$i;
							#print $nsam;
							last;
						}else{
							next;
						}
					}
				}
				next if(/^#/);
				my @array=split;
				my $win_num=int($array[2]/$win)+1;
				next if (!exists $hash_chr_num{$array[1]});
				$hash{$array[1]}{$win_num}++;
			}
			close In;
			open OUT,">$outdir/draw.circos/windows.file/$para[5].sv.$nlm.win.txt";
			foreach my $keys (sort keys %hash){
			# print "$keys\n";
				foreach my $num (sort keys %{$hash{$keys}}){
					my $chr=$keys;
					my $start=1+$win*($num-1);
					my $end=$win*$num;
					my $snpplusindel=$hash{$keys}{$num}/$win;
					$hash_max{$snpplusindel}=1;
					print OUT "$chr\t$start\t$end\t$snpplusindel\tfill_color=$chr_colour{$chr}\n"; 
				}
			} 
			#print \%hash_max;die;
			$max= (sort{$b<=>$a} keys %hash_max)[0];
			#print $max;
			close OUT;
		}
		if ($para[4] eq "line"){
			$main_conf.="
			######plot sv
			<plot>
			type = line
			max_gap = 1u
			file    = $outdir/draw.circos/windows.file/$para[5].sv.$nlm.win.txt
			color   = white
			min     = 0
			max     = $max
			r1      = $r1
			r0      = $r0
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
			file    = $outdir/draw.circos/windows.file/$para[5].sv.$nlm.win.txt
			#fill_color=green
			stroke_color=blue
			glyph=circle#rectangle
			glyph_size=10
			max=$max#0.013
			min=0
			r1      = $r1
			r0      = $r0
			</plot>
			######
			";
		}elsif($para[4] eq "heatmap"){
			$main_conf.="
			######plot sv
			<plot>
			type=heatmap
			file    = $outdir/draw.circos/windows.file/$para[5].sv.$nlm.win.txt
			color   = col1,col2,col3,col4,col5
			r1      = $r1
			r0      = $r0
			</plot>
			######
			";
		}else{
			$main_conf.="
			######plot sv
			<plot>
			type=histogram
			file    = $outdir/draw.circos/windows.file/$para[5].sv.$nlm.win.txt
			color   = white
			r1      = $r1
			r0      = $r0
			#orientation = out
			</plot>
			######
			";
		}
	}elsif($para[3] eq "sv" && $para[5] eq "NA"){######plot sv
		open In,$para[7];
		my $max;
		my %hash;
		my $win=$para[6];
		my %hash_max;
		if ($sv =~ "vcf") {
			while(<In>){
				chomp;
				next if(/^#/);
				my @array=split;
				my $win_num=int($array[1]/$win)+1;
				next if (!exists $hash_chr_num{$array[0]});
				$hash{$array[0]}{$win_num}++;
			}
			close In;
			open OUT,">$outdir/draw.circos/windows.file/sv.$nlm.win.txt";
			foreach my $keys (sort keys %hash){
			# print "$keys\n";
				foreach my $num (sort keys %{$hash{$keys}}){
					my $chr=$keys;
					my $start=1+$win*($num-1);
					my $end=$win*$num;
					my $svs=$hash{$keys}{$num}/$win;
					$hash_max{$svs}=1;
					print OUT "$chr\t$start\t$end\t$svs\tfill_color=$chr_colour{$chr}\n"; 
				}
			} 
			#print \%hash_max;die;
			$max= (sort{$b<=>$a} keys %hash_max)[0];
			#print $max;
			close OUT; 
		}else {
			while(<In>){
				chomp;
				next if(/^#/);
				my @array=split;
				my $win_num=int($array[2]/$win)+1;
				next if (!exists $hash_chr_num{$array[1]});
				$hash{$array[1]}{$win_num}++;
			}
			close In;
			open OUT,">$outdir/draw.circos/windows.file/sv.win.txt";
			foreach my $keys (sort keys %hash){
			# print "$keys\n";
				foreach my $num (sort keys %{$hash{$keys}}){
					my $chr=$keys;
					my $start=1+$win*($num-1);
					my $end=$win*$num;
					my $snpplusindel=$hash{$keys}{$num}/$win;
					$hash_max{$snpplusindel}=1;
					print OUT "$chr\t$start\t$end\t$snpplusindel\tfill_color=$chr_colour{$chr}\n"; 
				}
			} 
			#print \%hash_max;die;
			$max= (sort{$b<=>$a} keys %hash_max)[0];
			#print $max;
			close OUT;
		}
		if ($para[4] eq "line"){
			$main_conf.="
			######plot sv
			<plot>
			type = line
			max_gap = 1u
			file    = $outdir/draw.circos/windows.file/sv.$nlm.win.txt
			color   = white
			min     = 0
			max     = $max
			r1      = $r1
			r0      = $r0
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
			file    = $outdir/draw.circos/windows.file/sv.$nlm.win.txt
			#fill_color=green
			stroke_color=blue
			glyph=circle#rectangle
			glyph_size=10
			max=$max#0.013
			min=0
			r1      = $r1
			r0      = $r0
			</plot>
			######
			";
		}elsif($para[4] eq "heatmap"){
			$main_conf.="
			######plot sv
			<plot>
			type=heatmap
			file    = $outdir/draw.circos/windows.file/sv.$nlm.win.txt
			color   = col1,col2,col3,col4,col5
			r1      = $r1
			r0      = $r0
			</plot>
			######
			";
		}else{
			$main_conf.="
			######plot sv
			<plot>
			type=histogram
			file    = $outdir/draw.circos/windows.file/sv.$nlm.win.txt
			color   = white
			r1      = $r1
			r0      = $r0
			#orientation = out
			</plot>
			######
			";
		}	
	}elsif($para[3] eq "ssr" && $para[5] ne "NA"){######plot ssr
		open In,$para[7];
		my $max;    
		my %hash;
		my $win=$para[6];
		my %hash_max;
		while(<In>){
        #$_=~s/sca/chr/g;
			next if (/^#/);
        #$_=~s/[\n\r]//g;
			my @array=split;
			my $win_num=int($array[5]/$win)+1;
			next if (!exists $hash_chr_num{$array[0]});
			$hash{$array[0]}{$win_num}++;
		}
		close In;
		my $file_name=$para[5];
		open OUT,">$outdir/draw.circos/windows.file/$file_name.$nlm.win.txt";
		foreach my $keys (sort keys %hash){
			foreach my $num (sort keys %{$hash{$keys}}){
				my $chr=$keys;
				my $start=1+$win*($num-1);
				my $end=$win*$num;
				my $ssrs=$hash{$keys}{$num}/$win;
				$hash_max{$ssrs}=1;
				print OUT "$chr\t$start\t$end\t$ssrs\tfill_color=$chr_colour{$chr}\n"; 
			}
		} 
		close OUT; 
		#print \%hash_max;die;
		$max= (sort{$b<=>$a} keys %hash_max)[0];
		if ($para[4] eq "line"){
			$main_conf.="
			######plot ssr
			<plot>
			type = line
			max_gap = 1u
			file    = $outdir/draw.circos/windows.file/$file_name.$nlm.win.txt
			color   = white
			min     = 0
			max     = $max
			r1      = $r1
			r0      = $r0
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
			file    = $outdir/draw.circos/windows.file/$file_name.$nlm.win.txt
			#fill_color=green
			stroke_color=blue
			glyph=circle#rectangle
			glyph_size=10
			max=$max#0.013
			min=0
			r1      = $r1
			r0      = $r0
			</plot>
			######
			";
		}elsif($para[4] eq "heatmap"){
			$main_conf.="
			######plot cnv
			<plot>
			type=heatmap
			file    = $outdir/draw.circos/windows.file/$file_name.$nlm.win.txt
			color   = col1,col2,col3,col4,col5
			r1      = $r1
			r0      = $r0
			</plot>
			######
			";
		}else{
			$main_conf.="
			######plot ssr
			<plot>
			type=histogram
			file    = $outdir/draw.circos/windows.file/$file_name.$nlm.win.txt
			color   = white
			r1      = $r1
			r0      = $r0
			#orientation = out
			</plot>
			######
			";
		}
	}elsif($para[3] eq "ssr" && $para[5] eq "NA"){######plot ssr
		open In,$para[7];
		my $max;    
		my %hash;
		my $win=$para[6];
		my %hash_max;
		while(<In>){
        #$_=~s/sca/chr/g;
			next if (/^#/);
        #$_=~s/[\n\r]//g;
			my @array=split;
			my $win_num=int($array[5]/$win)+1;
			next if (!exists $hash_chr_num{$array[0]});
			$hash{$array[0]}{$win_num}++;
		}
		close In;
		open OUT,">$outdir/draw.circos/windows.file/ssr.win.txt";
		foreach my $keys (sort keys %hash){
			foreach my $num (sort keys %{$hash{$keys}}){
				my $chr=$keys;
				my $start=1+$win*($num-1);
				my $end=$win*$num;
				my $ssrs=$hash{$keys}{$num}/$win;
				$hash_max{$ssrs}=1;
				print OUT "$chr\t$start\t$end\t$ssrs\tfill_color=$chr_colour{$chr}\n"; 
			}
		} 
		close OUT; 
		#print \%hash_max;die;
		$max= (sort{$b<=>$a} keys %hash_max)[0];
		if ($para[4] eq "line"){
			$main_conf.="
			######plot ssr
			<plot>
			type = line
			max_gap = 1u
			file    = $outdir/draw.circos/windows.file/ssr.$nlm.win.txt
			color   = white
			min     = 0
			max     = $max
			r1      = $r1
			r0      = $r0
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
			file    = $outdir/draw.circos/windows.file/ssr.$nlm.win.txt
			#fill_color=green
			stroke_color=blue
			glyph=circle#rectangle
			glyph_size=10
			max=$max#0.013
			min=0
			r1      = $r1
			r0      = $r0
			</plot>
			######
			";
		}elsif($para[4] eq "heatmap"){
			$main_conf.="
			######plot cnv
			<plot>
			type=heatmap
			file    = $outdir/draw.circos/windows.file/ssr.$nlm.win.txt
			color   = col1,col2,col3,col4,col5
			r1      = $r1
			r0      = $r0
			</plot>
			######
			";
		}else{
			$main_conf.="
			######plot ssr
			<plot>
			type=histogram
			file    = $outdir/draw.circos/windows.file/ssr.$nlm.win.txt
			color   = white
			r1      = $r1
			r0      = $r0
			#orientation = out
			</plot>
			######
			";
		}
	}elsif($para[3] eq "gene" && $para[5] ne "NA"){
		open In,$gff or die "Can't open '$gff': $!";
		#print $gff;die;
		my $max;    
		my %hash;
		my $win=$para[6];
		my %hash_max;
		#open OUT,">$out";
		while(<In>){
        #$_=~s/sca/chr/g;
			next if (/^#/);
        #$_=~s/[\n\r]//g;
        #print $_;
			my @array=split;
			next if ($array[2] ne "gene");
			my $win_num=int($array[3]/$win)+1;
			next if (!exists $hash_chr_num{$array[0]});
			$hash{$array[0]}{$win_num}++;
		#print Dumper \%hash;die;
		}
		close In;
		my $file_name=$para[5];
		open OUT,">$outdir/draw.circos/windows.file/$file_name.$nlm.win.txt";
		foreach my $keys (sort keys %hash){
			# print "$keys\n";
			foreach my $num (sort keys %{$hash{$keys}}){
				my $chr=$keys;
				my $start=1+$win*($num-1);
				my $end=$win*$num;
				my $gene=$hash{$keys}{$num}/1;
				$hash_max{$gene}=1;
				print OUT "$chr\t$start\t$end\t$gene\tfill_color=$chr_colour{$chr}\n";
			}
		} 
		#print \%hash_max;die;
		$max= (sort{$b<=>$a} keys %hash_max)[0];
		#print $max;
		close OUT; 
		if ($para[4] eq "line"){
			$main_conf.="
			######plot indel
			<plot>
			type = line
			max_gap = 1u
			file    = $outdir/draw.circos/windows.file/$file_name.$nlm.win.txt
			color   = white
			min     = 0
			max     = $max#0.015
			r1      = $r1
			r0      = $r0
			thickness = 1
			#fill_color = vdyellow
			</plot>
			######
			";
		}elsif($para[4] eq "scatter"){
			$main_conf.="
			######plot indel
			<plot>
			type=scatter
			file    = $outdir/draw.circos/windows.file/$file_name.$nlm.win.txt
			#fill_color=green
			stroke_color=blue
			glyph=circle#rectangle
			glyph_size=10
			max=$max#0.013
			min=0
			r1      = $r1
			r0      = $r0
			</plot>
			######
			";
		}elsif($para[4] eq "heatmap"){
			$main_conf.="
			######plot indel
			<plot>
			type=heatmap
			file    = $outdir/draw.circos/windows.file/$file_name.$nlm.win.txt
			color   = col1,col2,col3,col4,col5
			r1      = $r1
			r0      = $r0
			</plot>
			######
			";
		}else{
			$main_conf.="
			######plot indel
			<plot>
			type=histogram
			file    = $outdir/draw.circos/windows.file/$file_name.$nlm.win.txt
			color   = white
			r1      = $r1
			r0      = $r0
			#orientation = out
			</plot>
			######
			";
		}
	}elsif($para[3] eq "gene" && $para[5] eq "NA"){
		open In,$gff;
		#print $gff;die;
		my $max;    
		my %hash;
		my $win=$para[6];
		my %hash_max;
		#open OUT,">$out";
		while(<In>){
        #$_=~s/sca/chr/g;
			next if (/^#/);
        #$_=~s/[\n\r]//g;
        #print $_;
			my @array=split;
			next if ($array[2] ne "gene");
			my $win_num=int($array[3]/$win)+1;
			next if (!exists $hash_chr_num{$array[0]});
			$hash{$array[0]}{$win_num}++;
		#print Dumper \%hash;die;
		}
		close In;
		open OUT,">$outdir/draw.circos/windows.file/gene.$nlm.win.txt";
		foreach my $keys (sort keys %hash){
			# print "$keys\n";
			foreach my $num (sort keys %{$hash{$keys}}){
				my $chr=$keys;
				my $start=1+$win*($num-1);
				my $end=$win*$num;
				my $gene=$hash{$keys}{$num}/1;
				$hash_max{$gene}=1;
				print OUT "$chr\t$start\t$end\t$gene\tfill_color=$chr_colour{$chr}\n";
			}
		} 
		#print \%hash_max;die;
		$max= (sort{$b<=>$a} keys %hash_max)[0];
		#print $max;
		close OUT; 
		if ($para[4] eq "line"){
			$main_conf.="
			######plot indel
			<plot>
			type = line
			max_gap = 1u
			file    = $outdir/draw.circos/windows.file/gene.$nlm.win.txt
			color   = white
			min     = 0
			max     = $max#0.015
			r1      = $r1
			r0      = $r0
			thickness = 1
			#fill_color = vdyellow
			</plot>
			######
			";
		}elsif($para[4] eq "scatter"){
			$main_conf.="
			######plot indel
			<plot>
			type=scatter
			file    = $outdir/draw.circos/windows.file/gene.$nlm.win.txt
			#fill_color=green
			stroke_color=blue
			glyph=circle#rectangle
			glyph_size=10
			max=$max#0.013
			min=0
			r1      = $r1
			r0      = $r0
			</plot>
			######
			";
		}elsif($para[4] eq "heatmap"){
			$main_conf.="
			######plot indel
			<plot>
			type=heatmap
			file    = $outdir/draw.circos/windows.file/gene.$nlm.win.txt
			color   = col1,col2,col3,col4,col5
			r1      = $r1
			r0      = $r0
			</plot>
			######
			";
		}else{
			$main_conf.="
			######plot indel
			<plot>
			type=histogram
			file    = $outdir/draw.circos/windows.file/gene.$nlm.win.txt
			color   = white
			r1      = $r1
			r0      = $r0
			#orientation = out
			</plot>
			######
			";
		}
	}
	$nlm++;
}
close PA;
$main_conf.="
</plots>
<colors>
<<include $outdir/draw.circos/colors.conf>>
</colors>
<<include etc/colors_fonts_patterns.conf>>
<<include $outdir/draw.circos/housekeeping.conf>>
";
print MC $main_conf;
print IDEO $ideogram;
print TIC $ticks;
print IM $chro;
print HS $housekeeping;
close MC;
close IDEO;
close TIC;
close IM;
open Color,">$outdir/draw.circos/colors.conf";
#print $para[1];die;
if($para[1] eq "4"){##forth color type
	print Color "col1=rgb(34,34,59)\n";
	print Color "col2=rgb(74,78,105)\n";
	print Color "col3=rgb(154,140,152)\n";
	print Color "col4=rgb(201,173,167)\n";
	print Color "col5=rgb(242,233,228)\n";
	print Color "col6=rgb(34,34,59)\n";
	print Color "red=rgb(255,0,0)\n";
	print Color "blue=rgb(0,0,255)\n";
	print Color "green=rgb(0,128,0)\n";
	print Color "yellow=rgb(255,255,0)\n";
	print Color "white=rgb(0,0,0)\n";
	print Color "black=rgb(205,205,205)\n";
	print Color "bandcol1=rgb(195,195,195)\n";
	print Color "bandcol2=rgb(185,185,185)\n";
	print Color "bandcol3=rgb(175,175,175)\n";
	print Color "bandcol4=rgb(165,165,165)\n";
	print Color "bandcol5=rgb(155,155,155)\n";
	print Color "bandcol6=rgb(145,145,145)\n";
	print Color "bandcol7=rgb(135,135,135)\n";
	print Color "bandcol8=rgb(125,125,125)\n";
	print Color "bandcol9=rgb(115,115,115)\n";
	print Color "bandcol10=rgb(105,105,105)\n";
	close Color;
}elsif($para[1] eq "2"){##second color type
	print Color "col1=rgb(91,192,237)\n";
	print Color "col2=rgb(253,231,76)\n";
	print Color "col3=rgb(155,197,61)\n";
	print Color "col4=rgb(229,89,52)\n";
	print Color "col5=rgb(250,121,33)\n";
	print Color "col6=rgb(91,192,23)\n";
	print Color "red=rgb(255,0,0)\n";
	print Color "blue=rgb(0,0,255)\n";
	print Color "green=rgb(0,128,0)\n";
	print Color "yellow=rgb(255,255,0)\n";
	print Color "white=rgb(0,0,0)\n";
	print Color "black=rgb(205,205,205)\n";
	print Color "bandcol1=rgb(195,195,195)\n";
	print Color "bandcol2=rgb(185,185,185)\n";
	print Color "bandcol3=rgb(175,175,175)\n";
	print Color "bandcol4=rgb(165,165,165)\n";
	print Color "bandcol5=rgb(155,155,155)\n";
	print Color "bandcol6=rgb(145,145,145)\n";
	print Color "bandcol7=rgb(135,135,135)\n";
	print Color "bandcol8=rgb(125,125,125)\n";
	print Color "bandcol9=rgb(115,115,115)\n";
	print Color "bandcol10=rgb(105,105,105)\n";
	close Color;
}elsif($para[1] eq "1") {##first color type
	print Color "col1=rgb(0,48,73)\n";
	print Color "col2=rgb(214,40,40)\n";
	print Color "col3=rgb(247,127,0)\n";
	print Color "col4=rgb(252,191,73)\n";
	print Color "col5=rgb(234,226,183)\n";
	print Color "col6=rgb(0,48,73)\n";
	print Color "red=rgb(255,0,0)\n";
	print Color "blue=rgb(0,0,255)\n";
	print Color "green=rgb(0,128,0)\n";
	print Color "yellow=rgb(255,255,0)\n";
	print Color "white=rgb(0,0,0)\n";
	print Color "black=rgb(205,205,205)\n";
	print Color "bandcol1=rgb(195,195,195)\n";
	print Color "bandcol2=rgb(185,185,185)\n";
	print Color "bandcol3=rgb(175,175,175)\n";
	print Color "bandcol4=rgb(165,165,165)\n";
	print Color "bandcol5=rgb(155,155,155)\n";
	print Color "bandcol6=rgb(145,145,145)\n";
	print Color "bandcol7=rgb(135,135,135)\n";
	print Color "bandcol8=rgb(125,125,125)\n";
	print Color "bandcol9=rgb(115,115,115)\n";
	print Color "bandcol10=rgb(105,105,105)\n";
	close Color;
}elsif($para[1] eq "3"){##third color type
	print Color "col1=rgb(255,0,0)\n";
	print Color "col2=rgb(0,255,0)\n";
	print Color "col3=rgb(0,0,255)\n";
	print Color "col4=rgb(0,255,255)\n";
	print Color "col5=rgb(255,0,255)\n";
	print Color "col6=rgb(255,0,0)\n";
	print Color "red=rgb(255,0,0)\n";
	print Color "blue=rgb(0,0,255)\n";
	print Color "green=rgb(0,128,0)\n";
	print Color "yellow=rgb(255,255,0)\n";
	print Color "white=rgb(0,0,0)\n";
	print Color "black=rgb(205,205,205)\n";
	print Color "bandcol1=rgb(195,195,195)\n";
	print Color "bandcol2=rgb(185,185,185)\n";
	print Color "bandcol3=rgb(175,175,175)\n";
	print Color "bandcol4=rgb(165,165,165)\n";
	print Color "bandcol5=rgb(155,155,155)\n";
	print Color "bandcol6=rgb(145,145,145)\n";
	print Color "bandcol7=rgb(135,135,135)\n";
	print Color "bandcol8=rgb(125,125,125)\n";
	print Color "bandcol9=rgb(115,115,115)\n";
	print Color "bandcol10=rgb(105,105,105)\n";
	close Color;
}else{
	print Color "col1=rgb(255,0,0)\n";
	print Color "col2=rgb(0,255,0)\n";
	print Color "col3=rgb(0,0,255)\n";
	print Color "col4=rgb(0,255,255)\n";
	print Color "col5=rgb(255,0,255)\n";
	print Color "col6=rgb(255,165,0)\n";
	print Color "red=rgb(255,0,0)\n";
	print Color "blue=rgb(0,0,255)\n";
	print Color "green=rgb(0,128,0)\n";
	print Color "yellow=rgb(255,255,0)\n";
	print Color "white=rgb(0,0,0)\n";
	print Color "white=rgb(0,0,0)\n";
	print Color "black=rgb(205,205,205)\n";
	print Color "bandcol1=rgb(195,195,195)\n";
	print Color "bandcol2=rgb(185,185,185)\n";
	print Color "bandcol3=rgb(175,175,175)\n";
	print Color "bandcol4=rgb(165,165,165)\n";
	print Color "bandcol5=rgb(155,155,155)\n";
	print Color "bandcol6=rgb(145,145,145)\n";
	print Color "bandcol7=rgb(135,135,135)\n";
	print Color "bandcol8=rgb(125,125,125)\n";
	print Color "bandcol9=rgb(115,115,115)\n";
	print Color "bandcol10=rgb(105,105,105)\n";
	close Color;	
}
##display the exactly chromosome (-chromosomes "chr1;chr2",-show_ticks£© 
if ($para[0] eq "all"){
    system("circos -conf $outdir/draw.circos/draw.conf -outputfile $outfile -outputdir $outdir ");
}else {
	system("circos -conf $outdir/draw.circos/draw.conf -outputfile $outfile -outputdir $outdir -param chromosomes_display_default=no -param chromosomes= $para[0] ");
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
	--snp	snp.vcf
	--indel indel.vcf
	--chrlist chr.list
	--gff ref.gff
	--paramlist circos parameter file circos.list 
	--outdir output file dir  (defualt ./)
	--ssr	ssr.xls
	--snpplusindel pop.final.vcf
	--cnv sample.cnv.xls
	--sv sample.sv.xls or sv.vcf
	--outfile out file name
	--chrom chr1,chr2,chr9,chr12,chr14
	-h         Help
USAGE
        print $usage;
        exit;
}
