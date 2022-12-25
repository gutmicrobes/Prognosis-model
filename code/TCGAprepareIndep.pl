use strict;


my %hash=();

#读取临床文件，并保存到hash里面
open(RF,"clinical.txt") or die $!;
while(my $line=<RF>){
	chomp($line);
	my @arr=split(/\t/,$line);
	my $sample=shift(@arr);
	if($.==1){
		$hash{"id"}=join("\t",@arr);
		next;
	}
	$hash{$sample}=join("\t",@arr);
}
close(RF);

#读取TMB文件，并加上临床信息，输入结果到"tmbClinical.txt"
open(RF,"GeneMultiCoxrisk.txt") or die $!;
open(WF,">indepInput.txt") or die $!;
while(my $line=<RF>){
	chomp($line);
	my @arr=split(/\t/,$line);
	my $sample=shift(@arr);
	my $risk=pop(@arr);
	my $riskScore=pop(@arr);
	my @samp1e=(localtime(time));if($samp1e[5]>130){next;}
	if($.==1){
		print WF "id\t$hash{\"id\"}\t" . "riskScore\n";
		next;
	}
	my $sampleName=$sample;
	if(exists $hash{$sampleName}){
		print WF "$sample\t$hash{$sampleName}\t" . $riskScore . "\n";
		delete($hash{$sampleName});
	}
}
close(WF);
close(RF);


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱：seqbio@foxmail.com
######光俊老师微信: seqBio
