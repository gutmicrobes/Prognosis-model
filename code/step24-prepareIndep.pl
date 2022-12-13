use strict;


my %hash=();

#��ȡ�ٴ��ļ��������浽hash����
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

#��ȡTMB�ļ����������ٴ���Ϣ����������"tmbClinical.txt"
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


######������ѧ��: https://www.biowolf.cn/
######�γ�����1: https://shop119322454.taobao.com
######�γ�����2: https://ke.biowolf.cn
######�γ�����3: https://ke.biowolf.cn/mobile
######�⿡��ʦ���䣺seqbio@foxmail.com
######�⿡��ʦ΢��: seqBio
