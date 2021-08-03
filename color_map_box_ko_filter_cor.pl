#!/usr/bin/perl -w
###########################################
# need to fix header about the genelink
# Ping Hu
############################################
use strict;
use GD;
my $usage = "usage: $0 <Map_coord_file, stat_file, acr_file, long_file, filter_file> \n";
my $PCUT=0.05;
my $SIGPERCENTCUT=10;
use Getopt::Long;
my ($opt_l, $opt_g);
GetOptions(
	   "l=s" => \$opt_l,
	   "g" => \$opt_g
	   );

die $usage unless @ARGV == 5;
my ($kegg_map_file, $statfile, $acr_file, $long_file, $filter) = @ARGV;


#my $PLOT_LINK2="http://mvic-biotech.na.pg.com/projects/Hair_Growth/MDEScalpHFGraph/";
my $PLOT_LINK2="http://mvic-biotech.na.pg.com/projects/GoldenEagle/GSS2656/Bioinformatics/GENEPLOT/";
my $PLOT_LINK=$PLOT_LINK2;



############################
#step 0. env variable setting
###########################
#my %KO_CHECK;
#my $KO_FILE="/home/ping/db/KEGG/ko/genes_ko.list";
#fill_hash(\%KO_CHECK,$KO_FILE );



my $stat_dir= $statfile;
if ($statfile =~  /([^\/]+)\s*$/ ){
    $stat_dir=$1;
}
if ($stat_dir =~/(\S+)\.stat/){
    $stat_dir=$1;
}

if(! (-e $stat_dir)){system("mkdir $stat_dir");} 

my $overlib_script_header="<DIV id=overDiv style=\"Z-INDEX: 1000; VISIBILITY: hidden; POSITION: absolute\"></DIV>
<SCRIPT language=JavaScript src=\"http://junonia.na.pg.com/overlib.js\"></SCRIPT>

<a href=\"http://www.bosrup.com/web/overlib/\">
  <img src=\"http://www.bosrup.com/web/overlib/power.gif\" width=\"88\"
      height=\"31\" alt=\"Popups by overLIB!\" border=\"0\"></a>
";

my $java_script_header="<script language=\"JavaScript\">
	
	
        function gsr( gene ) {
	    link = \'http://en.wikipedia.org/wiki/\'+ gene;
		window.open( link, \'popup\', \'width=900,height=700,resizable=yes,scrollbars=yes\' );
	}
	
	function gar( gene ) {
	    link = \'$PLOT_LINK\' + gene +'.PNG';
	    window.open( link, \'popup\', \'width=900,height=700,resizable=yes,scrollbars=yes\' );
	}

	function fue( gene ) {
	    link = \'$PLOT_LINK2\' + gene +'.PNG';
	    window.open( link, \'popup\', \'width=900,height=700,resizable=yes,scrollbars=yes\' );
	}
	

	function gkr( gene ) {
            link= 'http://www.genome.jp/dbget-bin/www_bget?'+gene;
	    window.open( link, 'popup', 'width=900,height=700,resizable=yes,scrollbars=yes' );
	}


	
	var newWin;
	function popupWindow( term, geneList, genes ) {
	    if( newWin != null ) {
			newWin.close();
		}
        newWin = window.open('http://yahoo.com','popup','width=900,height=700,resizable=yes'); 
    }

</script>
 <SCRIPT language=JavaScript src=\"http://junonia.na.pg.com/sorttable.js\"></SCRIPT>
";

####################################
# step 2 get the hash filled
####################################
my %location;
if ($opt_l){
    open (MMMM, "<$opt_l")or die "could not open $opt_l\n ";
    while (my $line=<MMMM>){
	chomp $line;
	my @tmp=split /\t/, $line;
	$location{$tmp[0]}=$tmp[1];
    }
    close MMMM;
}

my %long; # from affy to long annotation
fill_hash(\%long, $long_file);
my (%FFF, %PPP);
if ($opt_g){ 
    fill_stat_log_filter(\%PPP, \%FFF, $statfile, $filter); 
}else{
    fill_stat_filter(\%PPP, \%FFF, $statfile, $filter);
}

my %FA;
fill_hash(\%FA, $acr_file);
############################
#Step 1. get genome mapping
#        
#############################
my %ko2affy;
my %affy2ko;
open(IN,"<$kegg_map_file") or die "Cannot open $kegg_map_file to read.\n";

while(my $line =<IN>){
    chomp $line;
    my @tmp= split /\t/, $line;
    my $kk=$tmp[1];
    if ($tmp[1] =~ /ko:(\S+)/){
	$kk=$1;
    }
    $ko2affy{$kk}{$tmp[0]}=1;
    $affy2ko{$tmp[0]}=$kk;
#    print STDERR $kk, "\t", $tmp[0], "\n";
}
close IN;

#########################################################
#step2 go into all map to create png file and html file
#
#########################################################
my %REMOVE;
open(IN, "/home/ping/db/KEGG/map_1_10_2012/pathway.list")||die "could not open the map_title.tab\n";

#open(IN, "/home/ping/project/Fatstem/test.list")||die "could not open the map_title.tab\n";
open (K, ">$stat_dir.kegg.html")||die "could not open $stat_dir.kegg.html\n";
print K "<html> <SCRIPT language=JavaScript src=\"http://junonia.na.pg.com/sorttable.js\"></SCRIPT><body><h1>$stat_dir</h1><table border=1 class=\"sortable\" id=\"5429374558\"><tr><td>Category</td><td>$stat_dir MapName</td><td>Sig%</td><td>SigGene#</td><td>SigUP%</td><td>SigUpBox%</td><td>TotalGene</td><td>UpGene</td><td>SigUPGene</td><td>DownGene</td><td>SigDownGene</td><td>UpBox</td><td>SigUpBox</td><td>DownBox</td><td>SigDownBox</td></tr>\n";

print STDERR "Comparison\tmapid\tmapname\tsigN\tsigPercent\tsiguppercent\tsigbox\tsigboxupPercent\ttotalN\tupPersent\n";
my $cat="";
while (my $line = <IN>){
   chomp $line;
#   print STDERR $line, "\n";
   if ($line =~/\#/){
       $cat=$line;
       next;
   }
   my @tmp=split /\t/, $line;
   my $mapid=$tmp[0];
 
   my $mapname=$tmp[1];
   my $html_file="/home/ping/db/KEGG/map_1_10_2012_change/map$mapid.html";
   my $pic_file="/home/ping/db/KEGG/map_1_10_2012/map$mapid.png";
 #  print STDERR $mapid, "\n";
   if ((-e $pic_file) &&(-e $html_file)){
  #     print STDERR "find file $mapid\n";
#################################################################
#first go into the html file to get all the boxes and genes match
#
#################################################################
       my $im=newFromPng GD::Image($pic_file);
       my $im2= GD::Image->newTrueColor($im->width, $im->height);
       $im2->copy($im, 0, 0, 0, 0,$im->width, $im->height);
       my $red = $im2->colorAllocateAlpha(255,0,255, 64);# 127 is the transparency index, we can lower to 64
       my $blue= $im2->colorAllocateAlpha(0,255,255, 64);
       my $yellow =$im2->colorAllocateAlpha(255,255,0, 64); 
       my $UPBOX=0;
       my $DOWNBOX=0;
       my $UPSIGBOX=0;
       my $DOWNSIGBOX=0;
       open (F, $html_file)||die "could not open $html_file\n";
       open (H, ">$stat_dir/new.$mapid.html")||die "could not open new.$mapid.html\n";
       my %genes_in_map;
       while (my $ff=<F>){
	   my $print_st=0;
	   if ($ff =~ /\<script language=\"JavaScript\"\>/){
	       while(!( $ff =~ /\<\/script\>/)){
		   $ff=<F>;
	       }
	   }
	   if ($ff =~ /\<\/script\>/){
	       print H $java_script_header;
	       $print_st=1;
	   }
           if ($ff =~ /(.*\<body\>)(.*)/){
	       print H $1, $overlib_script_header, $2;
	       print H "<h2>", $stat_dir ,"</h2>\n";
	       #print H "<h3>", $map_sum{$line}, "</h3>\n";
    
	       $print_st=1;
	   }
	   
	   for ($ff){ s/\"\//\"http\:\/\/www\.genome\.jp\//gi;}
	   if ($ff =~ /\<img src\=\"\S+\" usemap\=\"\#mapdata\"/){
	       $ff="<img src\=\"$mapid.new.png\" usemap\=\"\#mapdata\" border\=0\>";
	   }




	   if ($ff =~ /\<area shape\=(\S+)\s+coords\=(\S+)\s+href\=(\"[^\>]+\")/){
	       my $shape=$1;
	       my $coords=$2; 
	       my $kegglink=$3;
	       my $PREF=1;
	       my $FREF=1;
	       my $COLOR=$yellow;
	       my $FILL=0;
	       my $match=0;
	       my $html="<i><table border=1 bgcolor=yellow ><tr><th colspan=5>kegg_name</th></tr><tr><th>probe</th><th>p_value</th><th>Fold</th><th>Symbol</th><th>Description</th></tr>";
   #            print STDERR "kegglink $kegglink\n";
	       if ($kegglink =~ /www_bget\?(\S+)\"/){
#		   print STDERR $1, " bget##########\n";
		   my @check=split/\+/, $1;
	
		   foreach my $i ( @check){
		       for ($i){s/\"//gi};
		       if ($i =~ /&/){
			   my @x=split /\&/, $i;
			   $i=$x[0];
		       }
	
		       for($i){s/ko\://gi;}
 #                      print STDERR "GeneID:\t", $i, "\n";
		       foreach my $j (keys %{$ko2affy{$i}}){
#			   print STDERR "Match $j\n";
			   if (defined $PPP{$j}){
			       $genes_in_map{$i}{$j}=1;
			       my $PCHECK=$PPP{$j};
			       my $FCHECK=$FFF{$j};
			       my $ll="";
			       if (defined $long{$j}){$ll=$long{$j};}
			       my $nnn="";
			       if (defined $FA{$j}){$nnn=$FA{$j};}
			       if($PCHECK<=0.1){
				   $html.="<tr><td>$j</td><td>$PCHECK</td><td>$FCHECK</td><td>$nnn</td><td>$ll</td></tr>";
			       }
			       if ($PCHECK<$PREF){
				   $PREF=$PCHECK;
				   $FREF=$FCHECK;
			       }
			       $match++;
			       
			   }
		       }
		   }
	       } 
	       $html.="</table></i>";
	       
	       if ($match == 0){next;}
	       if ($PREF<= $PCUT){
		   $FILL=1;
		   if ($FREF<0){
		       $DOWNSIGBOX++;
		   }elsif($FREF >0){
		       $UPSIGBOX++;
		   }
	       }
	       if ($FREF<0){
		   $COLOR=$blue;
		   $DOWNBOX++;
	       }elsif($FREF >0){
		   $COLOR=$red;
		   $UPBOX++;
	       }
	       
	       my @cor= split /\,/, $coords;
####/,###
	       if (@cor==4){
		   if ($FILL ==1){
		       $im2->filledRectangle(@cor,  $COLOR);
		   }else{
		       $im2->rectangle(@cor,  $COLOR);
		   }
		   print H "\<area shape=rect coords=$coords href=$kegglink onmouseover=\"return overlib(\'$html\',FULLHTML);\" onmouseout=\"return nd()\;\">\n";
		   $print_st=1;
	       }elsif (@cor>4){
		   my $poly=new GD::Polygon;
		   for(my $i=0; $i<@cor; $i=$i+2){
		       $poly->addPt($cor[$i], $cor[$i+1]);
		   }
				   
		   if ($FILL ==1){
		       $im2->filledPolygon($poly,  $COLOR);
		   }else{
		       $im2->polygon($poly,  $COLOR);
		   }
		   print H "\<area shape=poly coords=$coords href=$kegglink onmouseover=\"return overlib(\'$html\',FULLHTML);\" onmouseout=\"return nd()\;\">\n";
		   $print_st=1;
	       }elsif (@cor==3){
		   if ($FILL ==1){
		       $im2->arc($cor[0], $cor[1], $cor[2], $cor[2],0,360,$COLOR);
		   }else{
		       $im2->arc(@cor,1,0,360,$COLOR);
		   }
		   print H "\<area shape=circle coords=$coords href=$kegglink onmouseover=\"return overlib(\'$html\',FULLHTML);\" onmouseout=\"return nd()\;\">\n";
		   $print_st=1;
	       }
	   }
	   if ($ff =~ /\<\/map\>/){
	       my @geneTable=get_gene_table(\%genes_in_map);
	      
	       my $sigPercent=0;  
	       my $UpGene=$geneTable[2];
	       my $DownGene=$geneTable[3];
	       my $SigUpGene=$geneTable[4];
	       my $SigDownGene=$geneTable[5];
	       my $sigN=$SigUpGene+$SigDownGene;
	       my $TGene=$geneTable[1];
	       my $upPercent=0;
	       my $trend="green";
	       if ($SigUpGene >$SigDownGene){$trend="red";}
	       if ($SigUpGene <$SigDownGene){$trend="blue";}
	       
	       if ($TGene>0){
		   $sigPercent=int(100*$sigN/$TGene);
		   $upPercent=int(100*$UpGene/$TGene);
	       }
	       if (($sigN>=5)&&($sigPercent>=10)){
	       #if ($sigN>0){
		   my $siguppercent=int(100*$SigUpGene/$sigN);
		   my $sigbox=$UPSIGBOX+$DOWNSIGBOX;
		   my $sigboxupP=int(100*$UPSIGBOX/$sigbox);
		   print K "<tr><td>$cat</td><td><a href=\"./$stat_dir/new.$mapid.html\">$mapname</a></td><td><font color=$trend>$sigPercent</font></td><td>$sigN</td><td><font color=$trend>$siguppercent</font></td><td>$sigboxupP</td><td>$TGene</td><td>$UpGene</td><td>$SigUpGene</td><td>$DownGene</td><td>$SigDownGene</td><td>$UPBOX</td><td>$UPSIGBOX</td><td>$DOWNBOX</td><td>$DOWNSIGBOX</td></tr>\n";
		   $ff.=$geneTable[0]."\n<h4>". "MAPID:". $mapid. "\tSig%:$sigPercent\tSigUP%:$siguppercent\tGeneNumber:". $geneTable[1]. "\tUPGene:".  $geneTable[2]. "\tDownGene:". $geneTable[3]. "\tSigUPGene:".   $SigUpGene."\tSigDownGene:".   $SigDownGene. "\tUPBox:". $UPBOX . "\tUPSigBox:". $UPSIGBOX. "\tDownBox:". $DOWNBOX. "\tDownSigBox:". $DOWNSIGBOX. "</h4>\n";
		   print STDERR "$stat_dir\t$mapid\t$mapname\t$sigN\t$sigPercent\t$siguppercent\t$sigbox\t$sigboxupP\t$TGene\t$upPercent\n";
	       }else{
	       $REMOVE{$mapid}=1;
	       }
##################################################
###Here Is the place to write the summary kegg page
###Also remove the empty graph and html pages
#################################################


	   }
	   if ($print_st==0){print H $ff;}      
	       
       }
       close F;
       close H;




       open G,  ">$stat_dir/$mapid.new.png";
       binmode G;
       print G $im2->png;
       close G;
   }

}
close IN;
print K "</table></body></html>\n";
close K;
foreach my $i (keys %REMOVE){
    system("rm $stat_dir/$i.new.png\n");
    system("rm $stat_dir/new.$i.html\n");
}




sub get_gene_table{
    my ($geneh )=@_;
    my $result="
<table>
<tr><td valign=top>
<table border = 1>
<tr><th>FOLD_CHANGE_COLOR_KEY</th></tr>
<tr><td style=\"background-color:#6464f0\">Fold change &lt 0.25 </td></tr>
<tr><td style=\"background-color:#9292f0\"> 0.25 =&lt Fold change &lt 0.5</td></tr>
<tr><td style=\"background-color:#c1c1f0\"> 0.5 =&lt Fold change &lt 0.8</td></tr>

<tr><td style=\"background-color:#f0c1c1\"> 1.2 =&lt Fold change &lt 2</td></tr>
<tr><td style=\"background-color:#f09292\"> 2 =&lt Fold change &lt 4</td></tr>
<tr><td style=\"background-color:#f06464\"> 4 &lt Fold change</td></tr>

<tr><td style=\"background-color:f0f0f0\"> Other </td></tr>
</table>
</td>
<td valign=top>
<table border = 1>
<tr><th>pvalue_COLOR_KEY</th></tr>
<tr><td style=\"background-color:#FF0000\">p-value &lt 1.0E-4 </td></tr>
<tr><td style=\"background-color:#FF9700\"> 1.0E-4 =&lt p-value &lt 1.0E-3</td></tr>

<tr><td style=\"background-color:#FFCA00\"> 1.0E-3 =&lt p-value &lt 0.05</td></tr>
<tr><td style=\"background-color:#FFfd00\"> 0.05 =&lt p-value &lt 0.1</td></tr>
<tr><td style=\"background-color:#C0C0C0\"> other</td></tr>
</table>
</td>
<tr></table>
<SCRIPT language=JavaScript src=\"http://junonia.na.pg.com/sorttable.js\"></SCRIPT>
<table border= 1 class=\"sortable\" id=\"5429374552\">
\n<tr><td><font color=\"#000099\">Acronym</font></td><td><font color=\"#000099\">AffyProbe</font></td><td><font color=\"#330033\">P</font></td><td><font color=\"#996600\">Fold</font></td><td><font color=\"#003300\">Description</font></td><td><font color=\"#000099\">KEGG_NAME</font></td>";
    if ($opt_l){
	$result.= "<td><font color=\"#006600\">Location</font></td>\n";	
    }
    $result.="</tr>\n";
    my $geneNumber=0;
    my $UP=0;
    my $DOWN=0;
    my $SIGUP=0;
    my $SIGDOWN=0;
    my %check;
    foreach my $kegg_name(sort keys %{$geneh} ){
	my $names=$$geneh{$kegg_name};
	my %order;
	my $other="";
	foreach my $name (sort keys %{$$geneh{$kegg_name}}){
	    #print STDERR $name,"\n";
	    if (defined $check{$name}){next;}
	    else{
		$check{$name}=1;
	    }
	    #my $gene_name=$gene2name{$name};
	    my $P="N/A";
	    my $FC="N/A";
	    my $lname="N/A";
	    my $acr="N/A";
	    if (defined $long{$name}){$lname=$long{$name};}
	    if (defined $FA{$name}){$acr=$FA{$name};}
	    if (defined $PPP{$name}){$P=$PPP{$name};}
	    if (defined $FFF{$name}){$FC=$FFF{$name};}
	    if (($P eq "N/A")&&($FC eq "N/A")){next;}
	    $geneNumber++;
	    if ($FC>1){$UP++;}
	    if ($FC<1){$DOWN ++;}
	    if ($P<=0.05){
		   if ($FC>1){$SIGUP++;}
		   if ($FC<1){$SIGDOWN ++;}
	    }
	    my $span_class="style=\"background-color:f0f0f0\"";
	    my $span_p="style=\"background-color:C0C0C0\"";
	    if ($FC ne "N/A"){
		if ($FC >4){
		    $span_class="style=\"background-color:f06464\""; 
		}elsif($FC>2){
		    $span_class="style=\"background-color:f09292\"";
		}elsif($FC>1.2){
		    $span_class="style=\"background-color:f0c1c1\"";
		}elsif($FC>0.8){
		    $span_class="style=\"background-color:f0f0f0\"";
		}elsif($FC>0.5){
		    $span_class="style=\"background-color:c1c1f0\"";
		}elsif($FC>0.25){
		    $span_class="style=\"background-color:9292f0\"";
		}else{
		    $span_class="style=\"background-color:6464f0\"";
		}
	    }
	    if ($P ne "N/A"){
		if ($P <0.0001){ 
		    $span_p="style=\"background-color:ff0000\""; 
		}elsif($P<0.001){
		    $span_p="style=\"background-color:ff9700\"";
		}elsif($P<=0.05){
		    $span_p="style=\"background-color:ffca00\"";
		}elsif($P<0.1){
		    $span_p="style=\"background-color:fffd00\"";
		}else{
		    $span_p= "style=\"background-color:c0c0c0\"";
		} 
	    }
	    
	    my $gene_line="<tr>\n". 
	
		"<td ><font color=\"#006600\"><a href=\"javascript:void(0)\" onclick=\"gsr(\'$acr\');\">$acr</a></font></td>\n".
		"<td ><a href=\"javascript:void(0)\" onclick=\"gar(\'$name\');\">$name</a> </td>\n".	                                        "<td $span_p><font color=\"#330033\">$P</font></td>\n".
		"<td $span_class><font color=\"#996600\">$FC</font></td>\n".
		"<td><a href=\"javascript:void(0)\" onclick=\"fue(\'$name\');\"><font color=\"#003300\">$lname</font></a></td>\n".
	       "<td ><font color=\"#006600\"><a href=\"javascript:void(0)\" onclick=\"gkr(\'$kegg_name\');\">$kegg_name</a></font></td>\n".
		"</tr>\n";
	    if ($opt_l){
		my $loc="N/A";
		if (defined $location{$name}){
		    $loc= $location{$name}
		} 
		$gene_line="<tr>\n". 
		   
		    "<td ><font color=\"#006600\"><a href=\"javascript:void(0)\" onclick=\"gsr(\'$acr\');\">$acr</a></font></td>\n".
		    "<td ><a href=\"javascript:void(0)\" onclick=\"gar(\'$name\');\">$name</a> </td>\n".	                           
		    "<td $span_p><font color=\"#330033\">$P</font></td>\n".
		    "<td $span_class><font color=\"#996600\">$FC</font></td>\n".
		    "<td><a href=\"javascript:void(0)\" onclick=\"gar(\'$name\');\"><font color=\"#003300\">$lname</font></a></td>\n<td><font color=\"#006600\">$loc</font></td>\n".
		     "<td ><font color=\"#006600\"><a href=\"javascript:void(0)\" onclick=\"gkr(\'$kegg_name\');\">$kegg_name</a></font></td>\n".
		    "</tr>\n";
	    }
	    if ($P ne "N/A"){
		if (! defined $order{$P}){
		    $order{$P}=$gene_line;
		}else{
		    $order{$P}.=$gene_line;
		}
	    }else{
		$other.=$gene_line;
	    }
	}
	foreach my $i (sort {$a<=>$b} keys  %order){
	    $result.= $order{$i};
	}
	$result.=$other;
    }
    $result.="</table>\n";
    return ($result, $geneNumber, $UP, $DOWN, $SIGUP, $SIGDOWN);
}

sub fill_hash{
    my ($h, $file)=@_;
    open (X, "$file")||die "could not open $file\n";
    while (my $line = <X>){
	chomp $line;
	my @tmp=split /\t/, $line;
	if (@tmp<2){next;}
#	$$h{lc($tmp[0])}=lc($tmp[1]);
	$$h{($tmp[0])}=($tmp[1]); #for Bug
    }
    close X;
    return;
}

sub fill_stat_filter{
    my ($ph,$fh,$file, $filter)=@_;
    my %id;
    open (X, "$filter")||die "could not open $filter\n";
    while (my $line = <X>){
	chomp $line;
	my @tmp=split /\t/, $line;
	$id{$tmp[0]}=1;
    }
    close X;
    open (X, "$file")||die "could not open $file\n";
    while (my $line = <X>){
	chomp $line;
	my @tmp=split /\t/, $line;	
	if (defined $id{$tmp[0]}){
	    if (@tmp<2){next;}
	    $$ph{$tmp[0]}=$tmp[1];
	    $$fh{$tmp[0]}=$tmp[2];
	}
    }
    close X;
    return;
}

sub fill_stat_log_filter{
    my ($ph,$fh,$file, $filter)=@_;
    my %id;
    open (X, "$filter")||die "could not open $filter\n";
    while (my $line = <X>){
	chomp $line;
	my @tmp=split /\t/, $line;
	$id{$tmp[0]}=1;
    }
    close X;
    open (X, "$file")||die "could not open $file\n";
    my $tmpx=<X>;
    while (my $line = <X>){
	chomp $line;
	my @tmp=split /\t/, $line;
	if (defined $id{$tmp[0]}){
	    if (@tmp<2){next;}
	    $$ph{$tmp[0]}=$tmp[1];
	    $$fh{$tmp[0]}=2**($tmp[2]);
	}
    }
    close X;
    return;
}
