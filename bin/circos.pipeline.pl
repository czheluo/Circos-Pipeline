#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($snp,$indel,$sv,$cnv,$chrlist,$gff,$outdir,$windows,$snpplottype,$indelplottype,$color,$chromosome,$outfile,$cnvplottype,$svplottype,$densityplot,$snpplusindel,$snpplusindelplot,$sampleid,$ssr,$ssrplottype,$radius,$gene,$geneplottype);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use List::Util qw(max min);
my $version="1.1.0";
GetOptions(
	"help|?" =>\&USAGE,
	"windows:s"=>\$windows,
	"snp:s"=>\$snp,
    "indel:s"=>\$indel,
    "cnv:s"=>\$cnv,
    "sv:s"=>\$sv,
	"gff:s"=>\$gff,
    "chrlist:s"=>\$chrlist,
	"outdir:s"=>\$outdir,
    "snpplottype:s"=>\$snpplottype,
    "indelplottype:s"=>\$indelplottype,
    "cnvplottype:s"=>\$cnvplottype,
    "svplottype:s"=>\$svplottype,
	"snpplusindel:s"=>\$snpplusindel,
	"snpplusindelplot:s"=>\$snpplusindelplot,
    "color:s"=>\$color,
    "chromosome:s"=>\$chromosome,
	"outfile:s"=>\$outfile,
    "density:s"=>\$densityplot,
	"sampleid:s"=>\$sampleid,
	"ssr:s"=>\$ssr,
	"ssrplottype:s"=>\$ssrplottype,
    "radius:s"=>\$radius,
    "gene:s"=>\$gene,
    "geneplottype:s"=>\$geneplottype,
			) or &USAGE;
&USAGE unless ($snp and $gff and $chrlist);
#$plottype||="line scatter histogram heatmap";#stuckbar
$color ||=2;
$snpplottype||="line";
$indelplottype||="scatter";
$outfile||="circos";
$snpplusindelplot||="histogram";
$cnvplottype||="line";
$svplottype||="scatter";
$ssrplottype||="scatter";
$geneplottype||="line";
my @windows;
if ($windows){
	@windows=split/\,/,$windows;
	#print Dumper @windows;die;
}else{
	$windows||="10000,50000,100000,200000,500000";
	@windows=split/\,/,$windows;
}
$radius||="0.95r,0.85r,0.75r,0.65r,0.55r,0.45r,0.35r,0.25r,0.2r,0.1r";
my @rad=split/\,/,$radius;
#$density||="yes";
########
########work dir
########
$outdir||="./";
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
my $gff_windows=$windows[0];
#print $windows[0];die;
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
#foreach my $keys (sort {$hash_chr_num{$a}<=>$hash_chr_num{$b}} keys %hash_chr_num){
 #       print CHRBAND "chr"."\t"."-"."\t".$keys."\t".$hash_chr_num{$keys}."\t"."0"."\t"."$hash_chr_length{$keys}"."\t".$keys."\n";
#}
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
if ($chromosome) {
	my $chr=join(";",split(/\,/,$chromosome));
	$chro.="
	chromosomes_units = 1000000
	chromosomes_display_default = no
	chromosomes=$chr #chromosomes = hs1:(-100,120-);hs2;hs3;h4
	#chromosomes_order =$chr
	<colors>
	";
}else{
	$chro.="
	chromosomes_units = 1000000
	chromosomes_display_default = yes
	chromosomes=$chr_draw #chromosomes = hs1:(-100,120-);hs2;hs3;h4
	#chromosomes_order =$chr_draw
	<colors>
	";
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
#print $snp;
if ($snp && $sampleid){######plot snp
	open In,$snp;
	if ($snp=~/.gz$/) {
		close In;
		open In,"gunzip -c $snp|";
	}
	my $nsam;
	my @sample;
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
				if ($sampleid eq $sample[$i]){
					#print $sample[$i];
					$nsam=$sam+$i;
					#print $nsam;
					last;
				}else{
					next;
				}
			}
		}
	}
	close In;
	open In,$snp;
	if ($snp=~/.gz$/) {
		close In;
		open In,"gunzip -c $snp|";
	}
    my $max;    
    my %hash;
    my $win=$windows[1];
    my %hash_max;
    #open OUT,">$out";
    while(<In>){
        #$_=~s/sca/chr/g;
        next if (/^#/);
        #$_=~s/[\n\r]//g;
        #print $_;
        my @array=split;
		my ($gt,undef)=split/\:/,$array[$nsam],2;
		next if ($gt eq "./.");
        my $win_num=int($array[1]/$win)+1;
	    next if (!exists $hash_chr_num{$array[0]});
        $hash{$array[0]}{$win_num}++;
     #print Dumper \%hash;die;
    }
    close In;
    my $file_name=(split/\//,$snp)[-1];
    open OUT,">$outdir/draw.circos/windows.file/$file_name.win.txt";
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
    if ($snpplottype eq "line"){
        $main_conf.="
        ######plot snp
        <plot>
        type = line
        max_gap = 1u
        file    = $outdir/draw.circos/windows.file/$file_name.win.txt
        color   = white
        min     = 0
        max     = $max#0.015
        r1      = $rad[0]
        r0      = $rad[1]
        thickness = 1
        #fill_color = vdyellow
        </plot>
        ######
        ";
    }elsif($snpplottype eq "scatter"){
        $main_conf.="
        ######plot snp
        <plot>
        type=scatter
        file    = $outdir/draw.circos/windows.file/$file_name.win.txt
        #fill_color=green
        stroke_color=blue
        glyph=circle#rectangle
        glyph_size=10
        max=$max#0.013
        min=0
        r1      = $rad[0]
        r0      = $rad[1]
        </plot>
        ######
        ";
    }elsif($snpplottype eq "heatmap"){
        $main_conf.="
        ######plot snp
        <plot>
        type=heatmap
        file    = $outdir/draw.circos/windows.file/$file_name.win.txt
        color   = col1,col2,col3,col4,col5
        r1      = $rad[0]
        r0      = $rad[1]
        </plot>
        ######
        ";
    }else{
        $main_conf.="
        ######plot snp
        <plot>
        type=histogram
        file    = $outdir/draw.circos/windows.file/$file_name.win.txt
        color   = white
        r1      = $rad[0]
        r0      = $rad[1]
        #orientation = out
        </plot>
        ######
        ";
    }
}elsif($snp){
	open In,$snp;
	if ($snp=~/.gz$/) {
		close In;
		open In,"gunzip -c $snp|";
	}
    my $max;    
    my %hash;
    my $win=$windows[1];
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
    my $file_name=(split/\//,$snp)[-1];
    open OUT,">$outdir/draw.circos/windows.file/$file_name.win.txt";
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
    if ($snpplottype eq "line"){
        $main_conf.="
        ######plot snp
        <plot>
        type = line
        max_gap = 1u
        file    = $outdir/draw.circos/windows.file/$file_name.win.txt
        color   = white
        min     = 0
        max     = $max#0.015
        r1      = $rad[0]
        r0      = $rad[1]
        thickness = 1
        #fill_color = vdyellow
        </plot>
        ######
        ";
    }elsif($snpplottype eq "scatter"){
        $main_conf.="
        ######plot snp
        <plot>
        type=scatter
        file    = $outdir/draw.circos/windows.file/$file_name.win.txt
        #fill_color=green
        stroke_color=blue
        glyph=circle#rectangle
        glyph_size=10
        max=$max#0.013
        min=0
        r1      = $rad[0]
        r0      = $rad[1]
        </plot>
        ######
        ";
    }elsif($snpplottype eq "heatmap"){
        $main_conf.="
        ######plot snp
        <plot>
        type=heatmap
        file    = $outdir/draw.circos/windows.file/$file_name.win.txt
        color   = col1,col2,col3,col4,col5
        r1      = $rad[0]
        r0      = $rad[1]
        </plot>
        ######
        ";
    }else{
        $main_conf.="
        ######plot snp
        <plot>
        type=histogram
        file    = $outdir/draw.circos/windows.file/$file_name.win.txt
        color   = white
        r1      = $rad[0]
        r0      = $rad[1]
        #orientation = out
        </plot>
        ######
        ";
    }
}
if($indel && $sampleid){######plot indel
	open In,$indel;
	if ($indel=~/.gz$/) {
		close In;
		open In,"gunzip -c $indel|";
	}
	my $nsam;
	my @sample;
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
				if ($sampleid eq $sample[$i]){
					#print $sample[$i];
					$nsam=$sam+$i;
					#print $nsam;
					last;
				}else{
					next;
				}
			}
		}
	}
	close In;
    open In,$snp;
	if ($snp=~/.gz$/) {
		close In;
		open In,"gunzip -c $snp|";
	}
	my $max;    
    my %hash;
    my $win=$windows[2];
    my %hash_max;
    #open OUT,">$out";
    while(<In>){
        #$_=~s/sca/chr/g;
        next if (/^#/);
        #$_=~s/[\n\r]//g;
        #print $_;
        my @array=split;
		my ($gt,undef)=split/\:/,$array[$nsam],2;
		next if ($gt eq "./.");
        my $win_num=int($array[1]/$win)+1;
	    next if (!exists $hash_chr_num{$array[0]});
        $hash{$array[0]}{$win_num}++;
     #print Dumper \%hash;die;
    }
    close In;
    my $file_name=(split/\//,$indel)[-1];
    open OUT,">$outdir/draw.circos/windows.file/$file_name.win.txt";
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
    if ($indelplottype eq "line"){
        $main_conf.="
        ######plot indel
        <plot>
        type = line
        max_gap = 1u
        file    = $outdir/draw.circos/windows.file/$file_name.win.txt
        color   = white
        min     = 0
        max     = $max#0.015
        r1      = $rad[2]
        r0      = $rad[3]
        thickness = 1
        #fill_color = vdyellow
        </plot>
        ######
        ";
    }elsif($indelplottype eq "scatter"){
        $main_conf.="
        ######plot indel
        <plot>
        type=scatter
        file    = $outdir/draw.circos/windows.file/$file_name.win.txt
        #fill_color=green
        stroke_color=blue
        glyph=circle#rectangle
        glyph_size=10
        max=$max#0.013
        min=0
        r1      = $rad[2]
        r0      = $rad[3]
        </plot>
        ######
        ";
    }elsif($indelplottype eq "heatmap"){
        $main_conf.="
        ######plot indel
        <plot>
        type=heatmap
        file    = $outdir/draw.circos/windows.file/$file_name.win.txt
        color   = col1,col2,col3,col4,col5
        r1      = $rad[2]
        r0      = $rad[3]
        </plot>
        ######
        ";
    }else{
        $main_conf.="
        ######plot indel
        <plot>
        type=histogram
        file    = $outdir/draw.circos/windows.file/$file_name.win.txt
        color   = white
        r1      = $rad[2]
        r0      = $rad[3]
        #orientation = out
        </plot>
        ######
        ";
    }
}elsif($indel){
	open In,$indel;
	if ($indel=~/.gz$/) {
		close In;
		open In,"gunzip -c $indel|";
	}
    my $max;    
    my %hash;
    my $win=$windows[2];
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
    my $file_name=(split/\//,$indel)[-1];
    open OUT,">$outdir/draw.circos/windows.file/$file_name.win.txt";
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
    if ($indelplottype eq "line"){
        $main_conf.="
        ######plot indel
        <plot>
        type = line
        max_gap = 1u
        file    = $outdir/draw.circos/windows.file/$file_name.win.txt
        color   = white
        min     = 0
        max     = $max#0.015
        r1      = $rad[2]
        r0      = $rad[3]
        thickness = 1
        #fill_color = vdyellow
        </plot>
        ######
        ";
    }elsif($indelplottype eq "scatter"){
        $main_conf.="
        ######plot indel
        <plot>
        type=scatter
        file    = $outdir/draw.circos/windows.file/$file_name.win.txt
        #fill_color=green
        stroke_color=blue
        glyph=circle#rectangle
        glyph_size=10
        max=$max#0.013
        min=0
        r1      = $rad[2]
        r0      = $rad[3]
        </plot>
        ######
        ";
    }elsif($indelplottype eq "heatmap"){
        $main_conf.="
        ######plot indel
        <plot>
        type=heatmap
        file    = $outdir/draw.circos/windows.file/$file_name.win.txt
        color   = col1,col2,col3,col4,col5
        r1      = $rad[2]
        r0      = $rad[3]
        </plot>
        ######
        ";
    }else{
        $main_conf.="
        ######plot indel
        <plot>
        type=histogram
        file    = $outdir/draw.circos/windows.file/$file_name.win.txt
        color   = white
        r1      = $rad[2]
        r0      = $rad[3]
        #orientation = out
        </plot>
        ######
        ";
    }
}
if($snpplusindel && $sampleid){######plot $snpplusindel
    open In,$snpplusindel;
	if ($snpplusindel=~/.gz$/) {
		close In;
		open In,"gunzip -c $snpplusindel|";
	}
	my $nsam;
	my @sample;
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
				if ($sampleid eq $sample[$i]){
					#print $sample[$i];
					$nsam=$sam+$i;
					#print $nsam;
					last;
				}else{
					next;
				}
			}
		}
	}
	close In;
	open In,$snp;
	if ($snp=~/.gz$/) {
		close In;
		open In,"gunzip -c $snp|";
	}
    my $max;    
    my %hash;
    my $win=$windows[3];
    my %hash_max;
    #open OUT,">$out";
    while(<In>){
        #$_=~s/sca/chr/g;
        next if (/^#/);
        #$_=~s/[\n\r]//g;
        #print $_;
        my @array=split;
		my ($gt,undef)=split/\:/,$array[$nsam],2;
		next if ($gt eq "./.");
        my $win_num=int($array[1]/$win)+1;
	    next if (!exists $hash_chr_num{$array[0]});
        $hash{$array[0]}{$win_num}++;
     #print Dumper \%hash;die;
    }
    close In;
    my $file_name=(split/\//,$snpplusindel)[-1];
    open OUT,">$outdir/draw.circos/windows.file/$file_name.win.txt";
    foreach my $keys (sort keys %hash){
        # print "$keys\n";
        foreach my $num (sort keys %{$hash{$keys}}){
            my $chr=$keys;
            my $start=1+$win*($num-1);
            my $end=$win*$num;
            my $indel=$hash{$keys}{$num}/$win;
            $hash_max{$indel}=1;
            print OUT "$chr\t$start\t$end\t$snpplusindel\tfill_color=$chr_colour{$chr}\n"; 
        }
    } 
    #print \%hash_max;die;
    $max= (sort{$b<=>$a} keys %hash_max)[0];
    #print $max;
    close OUT; 
    if ($snpplusindelplot eq "line"){
        $main_conf.="
        ######plot indel
        <plot>
        type = line
        max_gap = 1u
        file    = $outdir/draw.circos/windows.file/$file_name.win.txt
        color   = white
        min     = 0
        max     = $max#0.015
        r1      = $rad[4]
        r0      = $rad[5]
        thickness = 1
        #fill_color = vdyellow
        </plot>
        ######
        ";
    }elsif($snpplusindelplot eq "scatter"){
        $main_conf.="
        ######plot indel
        <plot>
        type=scatter
        file    = $outdir/draw.circos/windows.file/$file_name.win.txt
        #fill_color=green
        stroke_color=blue
        glyph=circle#rectangle
        glyph_size=10
        max=$max#0.013
        min=0
        r1      = $rad[4]
        r0      = $rad[5]
        </plot>
        ######
        ";
    }elsif($snpplusindelplot eq "heatmap"){
        $main_conf.="
        ######plot indel
        <plot>
        type=heatmap
        file    = $outdir/draw.circos/windows.file/$file_name.win.txt
        color   = col1,col2,col3,col4,col5
        r1      = $rad[4]
        r0      = $rad[5]
        </plot>
        ######
        ";
    }else{
        $main_conf.="
        ######plot indel
        <plot>
        type=histogram
        file    = $outdir/draw.circos/windows.file/$file_name.win.txt
        color   = white
        r1      = $rad[4]
        r0      = $rad[5]
        #orientation = out
        </plot>
        ######
        ";
    }
}elsif($snpplusindel){
	open In,$snpplusindel;
	if ($snpplusindel=~/.gz$/) {
		close In;
		open In,"gunzip -c $snpplusindel|";
	}
    my $max;    
    my %hash;
    my $win=$windows[3];
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
    my $file_name=(split/\//,$snpplusindel)[-1];
    open OUT,">$outdir/draw.circos/windows.file/$file_name.win.txt";
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
    if ($snpplusindelplot eq "line"){
        $main_conf.="
        ######plot indel
        <plot>
        type = line
        max_gap = 1u
        file    = $outdir/draw.circos/windows.file/$file_name.win.txt
        color   = white
        min     = 0
        max     = $max#0.015
        r1      = $rad[4]
        r0      = $rad[5]
        thickness = 1
        #fill_color = vdyellow
        </plot>
        ######
        ";
    }elsif($snpplusindelplot eq "scatter"){
        $main_conf.="
        ######plot indel
        <plot>
        type=scatter
        file    = $outdir/draw.circos/windows.file/$file_name.win.txt
        #fill_color=green
        stroke_color=blue
        glyph=circle#rectangle
        glyph_size=10
        max=$max#0.013
        min=0
        r1      = $rad[4]
        r0      = $rad[5]
        </plot>
        ######
        ";
    }elsif($snpplusindelplot eq "heatmap"){
        $main_conf.="
        ######plot indel
        <plot>
        type=heatmap
        file    = $outdir/draw.circos/windows.file/$file_name.win.txt
        color   = col1,col2,col3,col4,col5
        r1      = $rad[4]
        r0      = $rad[5]
        </plot>
        ######
        ";
    }else{
        $main_conf.="
        ######plot indel
        <plot>
        type=histogram
        file    = $outdir/draw.circos/windows.file/$file_name.win.txt
        color   = white
        r1      = $rad[4]
        r0      = $rad[5]
        #orientation = out
        </plot>
        ######
        ";
    }
}
if($cnv){######plot cnv
    open In,$cnv;
    my ($file_name,undef)=split/\./,$cnv,2;
    my $max;    
    my %hash_max;
	#print $file_name;
    open OUT,">$outdir/draw.circos/windows.file/$file_name.cnv.txt";
    while(<In>){
        chomp;
        next if (/^#/ || "");
        my ($type,$region,$possize,undef,$pval,undef)=split/\s+/,$_,6;
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
    if ($cnvplottype eq "line"){
        $main_conf.="
        ######plot cnv
        <plot>
        type = line
        max_gap = 1u
        file    = $outdir/draw.circos/windows.file/$file_name.cnv.txt
        color   = white
        min     = 0
        max     = $max
        r1      = $rad[6]
        r0      = $rad[7]
        thickness = 1
        #fill_color = vdyellow
        </plot>
        ######
        ";
    }elsif($cnvplottype eq "scatter"){
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
        r1      = $rad[6]
        r0      = $rad[7]
        </plot>
        ######
        ";
    }elsif($cnvplottype eq "heatmap"){
        $main_conf.="
        ######plot cnv
        <plot>
        type=heatmap
        file    = $outdir/draw.circos/windows.file/$file_name.cnv.txt
        color   = col1,col2,col3,col4,col5
        r1      = $rad[6]
        r0      = $rad[7]
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
        r1      = $rad[6]
        r0      = $rad[7]
        #orientation = out
        </plot>
        ######
        ";
    }
}
    
if($sv){######plot indel
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
    if ($svplottype eq "line"){
        $main_conf.="
        ######plot sv
        <plot>
        type = line
        max_gap = 1u
        file    = $outdir/draw.circos/windows.file/$file_name.sv.txt
        color   = white
        min     = 0
        max     = $max
        r1      = $rad[8]
        r0      = $rad[9]
        thickness = 1
        #fill_color = vdyellow
        </plot>
        ######
        ";
    }elsif($svplottype eq "scatter"){
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
        r1      = $rad[8]
        r0      = $rad[9]
        </plot>
        ######
        ";
    }elsif($svplottype eq "heatmap"){
        $main_conf.="
        ######plot sv
        <plot>
        type=heatmap
        file    = $outdir/draw.circos/windows.file/$file_name.sv.txt
        color   = col1,col2,col3,col4,col5
        r1      = $rad[8]
        r0      = $rad[9]
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
        r1      = $rad[8]
        r0      = $rad[9]
        #orientation = out
        </plot>
        ######
        ";
    }
}
if($ssr){######plot ssr
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
    if ($ssrplottype eq "line"){
        $main_conf.="
        ######plot ssr
        <plot>
        type = line
        max_gap = 1u
        file    = $outdir/draw.circos/windows.file/$file_name.ssr.txt
        color   = white
        min     = 0
        max     = $max
        r1      = $rad[10]
        r0      = $rad[11]
        thickness = 1
        #fill_color = vdyellow
        </plot>
        ######
        ";
    }elsif($ssrplottype eq "scatter"){
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
        r1      = $rad[10]
        r0      = $rad[11]
        </plot>
        ######
        ";
    }elsif($ssrplottype eq "heatmap"){
        $main_conf.="
        ######plot cnv
        <plot>
        type=heatmap
        file    = $outdir/draw.circos/windows.file/$file_name.ssr.txt
        color   = col1,col2,col3,col4,col5
        r1      = $rad[10]
        r0      = $rad[11]
        </plot>
        ######
        ";
    }else{
        $main_conf.="
        ######plot cnv
        <plot>
        type=histogram
        file    = $outdir/draw.circos/windows.file/$file_name.ssr.txt
        color   = white
        r1      = $rad[10]
        r0      = $rad[11]
        #orientation = out
        </plot>
        ######
        ";
    }
}
##sv density plot 
if ($densityplot && $sv){
    open IN,"$outdir/draw.circos/band.txt";
    my %region;
	my %hash_max;
    #open Out,">$fout";
    while (<IN>) {
	    chomp;
	    my (undef,$chr,undef,undef,$start,$end,undef)=split/\s+/,$_;
	    $region{$chr}{join("\t",$start,$end)}=1;
    }
    close IN;
    open IN,$sv,
	my ($file_name,undef)=split/\./,$sv,2;
    my %stat;
    while (<IN>) {
	    chomp;
	    next if(/^#/);
	    my ($chr,$Pos1,undef,undef,$Pos2,undef)=split/\s+/,$_;
	    foreach my $region (sort keys %{$region{$chr}}) {
		    my ($pos3,$pos4)=split(/\t/,$region);
		    if (($Pos1 > $pos3 && $Pos1 <$pos4) ||($Pos2 > $pos3 && $Pos2 < $pos4) || ($pos3 > $Pos1 && $pos3 < $Pos2) || ($pos4 > $Pos1 && $pos4 < $Pos2)) {
			    	$stat{$chr}{$region}{sv}++;
		    }
	    } 
		#$stat{$chr}{$region}{sv}=$stat{$chr}{$region}{sv}/100000;
    }
    close IN;
    open Out,">$outdir/draw.circos/windows.file/$file_name.sv.txt";
    foreach my $chr (sort keys %region) {
	    foreach my $region (sort keys %{$region{$chr}}) {
		    $stat{$chr}{$region}{sv}||=0;
			$hash_max{$chr}=$stat{$chr}{$region}{sv};
		    print Out join("\t",$chr,$region,$stat{$chr}{$region}{sv},"fill_color=$chr_colour{$chr}"),"\n";
	    }
    }
    close Out;
	my $max=(sort{$b<=>$a} values %hash_max)[0];
    if ($svplottype eq "line"){
        $main_conf.="
        ######plot sv
        <plot>
        type = line
        max_gap = 1u
        file    = $outdir/draw.circos/windows.file/$file_name.sv.txt
        color   = white
        min     = 0
        max     = $max
        r1      = $rad[6]
        r0      = $rad[7]
        thickness = 1
        #fill_color = vdyellow
        </plot>
        ######
        ";
    }elsif($svplottype eq "scatter"){
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
        r1      = $rad[6]
        r0      = $rad[7]
        </plot>
        ######
        ";
    }elsif($svplottype eq "heatmap"){
        $main_conf.="
        ######plot sv
        <plot>
        type=heatmap
        file    = $outdir/draw.circos/windows.file/$file_name.sv.txt
        color   = col1,col2,col3,col4,col5
        r1      = $rad[6]
        r0      = $rad[7]
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
        r1      = $rad[6]
        r0      = $rad[7]
        #orientation = out
        </plot>
        ######
        ";
    }
}
if($densityplot && $cnv){
    open IN,"$outdir/draw.circos/band.txt";
    my %region;
	my %hash_max;
    #open Out,">$fout";
    while (<IN>) {
	    chomp;
	    my (undef,$chr,undef,undef,$start,$end,undef)=split/\s+/,$_;
	    $region{$chr}{join("\t",$start,$end)}=1;
    }
    close IN;
    open IN,$cnv;
	my ($file_name,undef)=split/\./,$cnv,2;
    my %stat;
    while (<IN>) {
	    chomp;
	    next if ($_ eq "" || /^##/ || /^#/);
	    my ($chr,$intev,undef)=split/\s+/,$_;
	    my ($ch,$se)=split/\:/,$intev;
	    my ($Pos1,$Pos2)=split/\-/,$se;
	    foreach my $region (sort keys %{$region{$ch}}) {
		    my ($pos3,$pos4)=split(/\t/,$region);
		    if (($Pos1 > $pos3 && $Pos1 <$pos4) ||($Pos2 > $pos3 && $Pos2 < $pos4) || ($pos3 > $Pos1 && $pos3 < $Pos2) || ($pos4 > $Pos1 && $pos4 < $Pos2)) {
				    $stat{$ch}{$region}{cnv}++;
		    }
		    } 
    }
    close IN;
    open Out,">$outdir/draw.circos/windows.file/$file_name.cnv.txt";
    foreach my $chr (sort keys %region) {
	    foreach my $region (sort keys %{$region{$chr}}) {
		    $stat{$chr}{$region}{cnv}||=0;
		    print Out join("\t",$chr,$region,$stat{$chr}{$region}{cnv},"fill_color=$chr_colour{$chr}"),"\n";
	    }
    }
    close Out;
	my $max= (sort{$b<=>$a} values %hash_max)[0];
    if ($cnvplottype eq "line"){
        $main_conf.="
        ######plot cnv
        <plot>
        type = line
        max_gap = 1u
        file    = $outdir/draw.circos/windows.file/$file_name.cnv.txt
        color   = white
        min     = 0
        max     = $max
        r1      = $rad[4]
        r0      = $rad[5]
        thickness = 1
        #fill_color = vdyellow
        </plot>
        ######
        ";
    }elsif($cnvplottype eq "scatter"){
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
        r1      = $rad[4]
        r0      = $rad[5]
        </plot>
        ######
        ";
    }elsif($cnvplottype eq "heatmap"){
        $main_conf.="
        ######plot cnv
        <plot>
        type=heatmap
        file    = $outdir/draw.circos/windows.file/$file_name.cnv.txt
        color   = col1,col2,col3,col4,col5
        r1      = $rad[4]
        r0      = $rad[5]
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
        r1      = $rad[4]
        r0      = $rad[5]
        #orientation = out
        </plot>
        ######
        ";
    }
}
if($gene){
	open In,$gene;
	#print $gff;die;
    my $max;    
    my %hash;
    my $win=$windows[1];
    my %hash_max;
    #open OUT,">$out";
    while(<In>){
        #$_=~s/sca/chr/g;
        next if (/^#/);
        #$_=~s/[\n\r]//g;
        #print $_;
        my @array=split;
		next if ($array[2] ne "gene");
        my $win_num=int($array[4]/$win)+1;
	    next if (!exists $hash_chr_num{$array[0]});
        $hash{$array[0]}{$win_num}++;
     #print Dumper \%hash;die;
    }
    close In;
    open OUT,">$outdir/draw.circos/windows.file/gene.win.txt";
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
    if ($geneplottype eq "line"){
        $main_conf.="
        ######plot indel
        <plot>
        type = line
        max_gap = 1u
        file    = $outdir/draw.circos/windows.file/gene.win.txt
        color   = white
        min     = 0
        max     = $max#0.015
        r1      = $rad[12]
        r0      = $rad[13]
        thickness = 1
        #fill_color = vdyellow
        </plot>
        ######
        ";
    }elsif($geneplottype eq "scatter"){
        $main_conf.="
        ######plot indel
        <plot>
        type=scatter
        file    = $outdir/draw.circos/windows.file/gene.win.txt
        #fill_color=green
        stroke_color=blue
        glyph=circle#rectangle
        glyph_size=10
        max=$max#0.013
        min=0
        r1      = $rad[12]
        r0      = $rad[13]
        </plot>
        ######
        ";
    }elsif($geneplottype eq "heatmap"){
        $main_conf.="
        ######plot indel
        <plot>
        type=heatmap
        file    = $outdir/draw.circos/windows.file/gene.win.txt
        color   = col1,col2,col3,col4,col5
        r1      = $rad[12]
        r0      = $rad[13]
        </plot>
        ######
        ";
    }else{
        $main_conf.="
        ######plot indel
        <plot>
        type=histogram
        file    = $outdir/draw.circos/windows.file/gene.win.txt
        color   = white
        r1      = $rad[12]
        r0      = $rad[13]
        #orientation = out
        </plot>
        ######
        ";
    }
}
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
if ($color == 1){##forth color type
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
}elsif($color == 2){##second color type
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
}elsif($color == 3) {##first color type
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
}elsif($color ==4){##third color type
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
}

##display the exactly chromosome (-chromosomes "chr1;chr2",-show_ticks, 
if ($chromosome){
    
   system("circos -conf $outdir/draw.circos/draw.conf -outputfile $outfile -outputdir $outdir -param chromosomes_display_default=no -param chromosomes= $chromosome ");

}else {
     system("circos -conf $outdir/draw.circos/draw.conf -outputfile $outfile -outputdir $outdir ");
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
  --windows windows size (default "1000000,1200000,1400000,1600000")
  --snp	<file> slide windows of snp
  --indel <file> slide windows of indel
  --chrlist input chrome number
  --gff gff file
  --color four types of color to choose 
  --snpplottype the type of plot (histogram line scatter heatmap) default was line 
  --indelplottype   the type of plot  (histogram line scatter heatmap) default was sactter
  --chromosome which chromosome you want to draw (default was all)  eg:chr1,chr2,chr8,chr20
  --outdir output file dir  
  --snpplusindel snp plus the indel  
  --snpplusindelplot the type of snp pus indel  plot
  --density means draw sv and cnv dengsity plot, or then just make the location of sv and cnv with score and  ttest -log10(pvalue)
  -h         Help

USAGE
        print $usage;
        exit;
}
