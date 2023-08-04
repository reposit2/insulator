#!/usr/local/bin/perl

#$infile1 = './mydata1/Btrainoutrfx/stat_out1_rfx_pth00001.txt';
$infile1 = './mydata2/Btrainoutrfx/stat_out1_rfx_pth00001.txt';

#$infile1b = './mydata1/Btrainoutrfx/stat_out3_rfx_pth005.txt';
$infile1b = './mydata2/Btrainoutrfx/stat_out3_rfx_pth005.txt';

#$infile2 = './mydata1/Btrainoutfrrf/stat_out1_frrf_pth00001.txt';
$infile2 = './mydata2/Btrainoutfrrf/stat_out1_frrf_pth00001.txt';

#$infile2b = './mydata1/Btrainoutfrrf/stat_out3_frrf_pth005.txt';
$infile2b = './mydata2/Btrainoutfrrf/stat_out3_frrf_pth005.txt';

#$infile3 = './mydata1/Btrainoutffrf/stat_out1_ffrf_pth00001.txt';
$infile3 = './mydata2/Btrainoutffrf/stat_out1_ffrf_pth00001.txt';

#$infile3b = './mydata1/Btrainoutffrf/stat_out3_ffrf_pth005.txt';
$infile3b = './mydata2/Btrainoutffrf/stat_out3_ffrf_pth005.txt';

#$infile4 = './mydata1/Btrainoutrrrf/stat_out1_rrrf_pth00001.txt';
$infile4 = './mydata2/Btrainoutrrrf/stat_out1_rrrf_pth00001.txt';

#$infile4b = './mydata1/Btrainoutrrrf/stat_out3_rrrf_pth005.txt';
$infile4b = './mydata2/Btrainoutrrrf/stat_out3_rrrf_pth005.txt';

#$outfile = './mydata1/stat/stat_out1_rfx_frrf_ffrf_rrrf_pth005_summary.txt';
$outfile = './mydata2/stat/stat_out1_rfx_frrf_ffrf_rrrf_pth005_summary.txt';

#$pth = 0.001;
#$pth = 0.0001;
$pth = 0.05;

$zscf = -1;

open(IN,"$infile1b");
while($l = <IN>) {
	chomp $l;
	@list = split(/\t/,$l);
	$tf1pvalue0{$list[0]} = $list[8];
	$tf1zscore0{$list[0]} = $list[3];
	$tf1rvalue0{$list[0]} = $list[4];
	$tf1numa0{$list[0]} = $list[5];
	$tf1numb0{$list[0]} = $list[6];
}
close(IN);	

open(IN,"$infile1");
while($l = <IN>) {
	chomp $l;
	if ($l =~ /trainout\_(.+?)\_/) {
		$id111 = $1;
		$coutf1{$id111} = 1;
	} elsif ($l =~ /epa1/) {
		$mean11 = (split(/\t/,$l))[3];
		$mean11 =~ s/\s//g;
		$l2 = <IN>;
		chomp $l2;
		$mean21 = (split(/\t/,$l2))[3];
		$mean21 =~ s/\s//g;
	} elsif ($l =~ /MannwhitneyuResult/) {
		$l2 = <IN>;
		chomp $l2;
		$pvalue1 = (split(/\t/,$l2))[1];
		$pvalue1 =~ s/\s//g;
#		if ($mean11 < $mean21 && $pvalue1 < $pth) {
##		if ($mean11 > $mean21 && $pvalue1 < $pth) {
		if ($mean11 > $mean21 && $tf1pvalue0{$id111} ne '' && $tf1pvalue0{$id111} < $pth) {
			$tf1mean1{$id111} = $mean11;
			$tf1mean2{$id111} = $mean21;
			$tf1pvalue{$id111} = $tf1pvalue0{$id111};
			$tf1zscore{$id111} = $tf1zscore0{$id111} * $zscf;
			$tf1rvalue{$id111} = $tf1rvalue0{$id111};
			$tf1numa{$id111} = $tf1numa0{$id111};
		        $tf1numb{$id111} = $tf1numb0{$id111};
		}
	}
}
close(IN);

open(IN,"$infile2b");
while($l = <IN>) {
        chomp $l;
        @list = split(/\t/,$l);
        $tf2pvalue0{$list[0]} = $list[8];
	$tf2zscore0{$list[0]} = $list[3]; $tf2rvalue0{$list[0]} = $list[4]; $tf2numa0{$list[0]} = $list[5];
        $tf2numb0{$list[0]} = $list[6];
}
close(IN);

open(IN,"$infile2");
while($l = <IN>) {
        chomp $l;
	if ($l =~ /trainout\_(.+?)\_/) {
                $id112 = $1;
		$coutf2{$id112} = 1;
        } elsif ($l =~ /epa1/) {
                $mean12 = (split(/\t/,$l))[3];
		$mean12 =~ s/\s//g;
                $l2 = <IN>;
                chomp $l2;
                $mean22 = (split(/\t/,$l2))[3];
		$mean22 =~ s/\s//g;
        } elsif ($l =~ /MannwhitneyuResult/) {
                $l2 = <IN>;
                chomp $l2; 
                $pvalue2 = (split(/\t/,$l2))[1];
		$pvalue2 =~ s/\s//g;
##                if ($mean12 < $mean22 && $pvalue2 < $pth) {
#		if ($mean12 > $mean22 && $pvalue2 < $pth) {
		if ($mean12 < $mean22 && $tf2pvalue0{$id112} ne '' && $tf2pvalue0{$id112} < $pth) {
                        $tf2mean1{$id112} = $mean12;
                        $tf2mean2{$id112} = $mean22;
			$tf2pvalue{$id112} = $tf2pvalue0{$id112};
			$tf2zscore{$id112} = $tf2zscore0{$id112} * $zscf;
			$tf2rvalue{$id112} = $tf2rvalue0{$id112};
			$tf2numa{$id112} = $tf2numa0{$id112};
                        $tf2numb{$id112} = $tf2numb0{$id112};
		}
        }
}
close(IN);

open(IN,"$infile3b");
while($l = <IN>) {
        chomp $l;
        @list = split(/\t/,$l);
        $tf3pvalue0{$list[0]} = $list[8];
	$tf3zscore0{$list[0]} = $list[3];
	$tf3rvalue0{$list[0]} = $list[4];
	$tf3numa0{$list[0]} = $list[5];
        $tf3numb0{$list[0]} = $list[6];
}
close(IN);

open(IN,"$infile3");
while($l = <IN>) {
        chomp $l;
	if ($l =~ /trainout\_(.+?)\_/) {
                $id113 = $1;
                $coutf3{$id113} = 1;
        } elsif ($l =~ /epa1/) {
                $mean13 = (split(/\t/,$l))[3];
                $mean13 =~ s/\s//g;
                $l2 = <IN>;
                chomp $l2;
                $mean23 = (split(/\t/,$l2))[3];
                $mean23 =~ s/\s//g;
        } elsif ($l =~ /MannwhitneyuResult/) {
                $l2 = <IN>;
                chomp $l2;
                $pvalue3 = (split(/\t/,$l2))[1];
                $pvalue3 =~ s/\s//g;
##                if ($mean13 < $mean23 && $pvalue3 < $pth) {
#               if ($mean13 > $mean23 && $pvalue3 < $pth) {
		if ($mean13 < $mean23 && $tf3pvalue0{$id113} ne '' && $tf3pvalue0{$id113} < $pth) {
                        $tf3mean1{$id113} = $mean13;
                        $tf3mean2{$id113} = $mean23;
			$tf3pvalue{$id113} = $tf3pvalue0{$id113};
			$tf3zscore{$id113} = $tf3zscore0{$id113} * $zscf;
			$tf3rvalue{$id113} = $tf3rvalue0{$id113};
			$tf3numa{$id113} = $tf3numa0{$id113};
			$tf3numb{$id113} = $tf3numb0{$id113};
		}
        }
}
close(IN);

open(IN,"$infile4b");
while($l = <IN>) {
        chomp $l;
        @list = split(/\t/,$l);
        $tf4pvalue0{$list[0]} = $list[8];
	$tf4zscore0{$list[0]} = $list[3];
	$tf4rvalue0{$list[0]} = $list[4];
	$tf4numa0{$list[0]} = $list[5];
        $tf4numb0{$list[0]} = $list[6];
}
close(IN);

open(IN,"$infile4");
while($l = <IN>) {
        chomp $l;
	if ($l =~ /trainout\_(.+?)\_/) {
                $id114 = $1;
                $coutf4{$id114} = 1;
        } elsif ($l =~ /epa1/) {
                $mean14 = (split(/\t/,$l))[3];
                $mean14 =~ s/\s//g;
                $l2 = <IN>;
                chomp $l2;
                $mean24 = (split(/\t/,$l2))[3];
                $mean24 =~ s/\s//g;
        } elsif ($l =~ /MannwhitneyuResult/) {
                $l2 = <IN>;
                chomp $l2;
                $pvalue4 = (split(/\t/,$l2))[1];
                $pvalue4 =~ s/\s//g;
##                if ($mean14 < $mean24 && $pvalue4 < $pth) {
#               if ($mean14 > $mean24 && $pvalue4 < $pth) {
		if ($mean14 < $mean24 && $tf4pvalue0{$id114} ne '' && $tf4pvalue0{$id114} < $pth) {
                        $tf4mean1{$id114} = $mean14;
                        $tf4mean2{$id114} = $mean24;
			$tf4pvalue{$id114} = $tf4pvalue0{$id114};
			$tf4zscore{$id114} = $tf4zscore0{$id114} * $zscf;
			$tf4rvalue{$id114} = $tf4rvalue0{$id114};
			$tf4numa{$id114} = $tf4numa0{$id114};
                        $tf4numb{$id114} = $tf4numb0{$id114};
		}
        }
}
close(IN);

open(OUT,">$outfile");
foreach $x (sort {$tf1pvalue{$a} <=> $tf1pvalue{$b}} keys %tf1pvalue) {
	next if ($tf2pvalue{$x} eq '' || $tf3pvalue{$x} eq '' || $tf4pvalue{$x} eq '');
	if ($id111 ne $id112) {
		print "\*$id111\t$id112\n";
	}
#	$d1 = $tf1mean2{$x} - $tf1mean1{$x};
	$d2 = $tf2mean2{$x} - $tf2mean1{$x};
	$d3 = $tf3mean2{$x} - $tf3mean1{$x};
	$d4 = $tf4mean2{$x} - $tf4mean1{$x};
	$d1 = $tf1mean1{$x} - $tf1mean2{$x};
#	$d2 = $tf2mean1{$x} - $tf2mean2{$x};
#	print "$x\npvalue\t$tf1pvalue{$x}\t$tf2pvalue{$x}\t$tf3pvalue{$x}\t$tf4pvalue{$x}\nmean_difference\t$d1\t$d2\t$d3\t$d4\nzscore\t$tf1zscore{$x}\t$tf2zscore{$x}\t$tf3zscore{$x}\t$tf4zscore{$x}\nrvalue\t$tf1rvalue{$x}\t$tf2rvalue{$x}\t$tf3rvalue{$x}\t$tf4rvalue{$x}\n#_scores\t$tf1numa{$x}\t$tf1numb{$x}\t$tf2numa{$x}\t$tf2numb{$x}\t$tf3numa{$x}\t$tf3numb{$x}\t$tf4numa{$x}\t$tf4numb{$x}\nmean\t$tf1mean1{$x}\t$tf1mean2{$x}\t$tf2mean1{$x}\t$tf2mean2{$x}\t$tf3mean1{$x}\t$tf3mean2{$x}\t$tf4mean1{$x}\t$tf4mean2{$x}\n";
	print OUT "$x\t$tf1pvalue{$x}\t$tf2pvalue{$x}\t$tf3pvalue{$x}\t$tf4pvalue{$x}\t$d1\t$d2\t$d3\t$d4\t$tf1zscore{$x}\t$tf2zscore{$x}\t$tf3zscore{$x}\t$tf4zscore{$x}\t$tf1rvalue{$x}\t$tf2rvalue{$x}\t$tf3rvalue{$x}\t$tf4rvalue{$x}\t$tf1numa{$x}\t$tf1numb{$x}\t$tf2numa{$x}\t$tf2numb{$x}\t$tf3numa{$x}\t$tf3numb{$x}\t$tf4numa{$x}\t$tf4numb{$x}\t$tf1mean1{$x}\t$tf1mean2{$x}\t$tf2mean1{$x}\t$tf2mean2{$x}\t$tf3mean1{$x}\t$tf3mean2{$x}\t$tf4mean1{$x}\t$tf4mean2{$x}\n";
	print "$x\tFR_vs\._RF\tFR_vs\._FF\tFR_vs\.RR\tFR_vs\.X\nzscore\t$tf2zscore{$x}\t$tf4zscore{$x}\t$tf3zscore{$x}\t$tf1zscore{$x}\npvalue\t$tf2pvalue{$x}\t$tf4pvalue{$x}\t$tf3pvalue{$x}\t$tf1pvalue{$x}\nmean_difference\t$d2\t$d4\t$d3\t$d1\nmean\t$tf2mean2{$x}\t$tf4mean2{$x}\t$tf3mean2{$x}\t$tf1mean1{$x}\n\t$tf2mean1{$x}\t$tf4mean1{$x}\t$tf3mean1{$x}\t$tf1mean2{$x}\nrvalue\t$tf2rvalue{$x}\t$tf4rvalue{$x}\t$tf3rvalue{$x}\t$tf1rvalue{$x}\n#_scores\t$tf2numb{$x}\t$tf4numb{$x}\t$tf3numb{$x}\t$tf1numa{$x}\n\t$tf2numa{$x}\t$tf4numa{$x}\t$tf3numa{$x}\t$tf1numb{$x}\n";
}
close(OUT);

#$n1 = (keys %coutf1);
#$n2 = (keys %coutf2);
#$n3 = (keys %coutf3);
#$n4 = (keys %coutf4);
#print "No._of_TFs\n$n1\t$n2\t$n3\t$n4\n";
