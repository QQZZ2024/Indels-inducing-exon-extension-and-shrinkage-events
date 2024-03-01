#####################################################
###Construction of individual-specific genome

from tqdm import tqdm
from time import sleep
path = '/home/GRCh37.p13.genome.fa'
with open(path) as f:
	a=[]
	for line in f:
		line = line.rstrip()
		if line.startswith('>'):
			continue
		a.append(line)
ID = open('/home/samplesID.txt')
for j in ID:
	j = j.rstrip()
	path1 = open("/home/" + j + "_indel.vcf")
	final = open(r"/home/" + j + "_person.fa",'w')
	pbar = tqdm(total=len(open("/home/" + j + "_indel.vcf").readlines()))
	err="error"
	chr1=1
	lst = list(a[chr1-1])
	for vcf in tqdm(path1):
		pbar.update(1)
		vcf = vcf.rstrip()
		lie = vcf.split('\t')
		chr2 = chr1
		chr1 = int(lie[0])
		pos1 = int(lie[1])
		pos2 = pos1 + len(lie[3])
		alt = lie[4]
		ref = lie[3]
		id1 = lie[2]
		if chr2 != chr1:
			new_string = "".join(lst)
			final.write(">" + "chr" + str(chr2))
			final.write('\n')
			final.write(new_string)
			final.write('\n')
			lst = list(a[chr1-1])
		gai = lst[pos1-1:pos2-1]
		gai1 = "".join(gai)
		if gai1 == ref:
			lst[pos1-1:pos2-1] = alt
		else:
			print(err)
			print(id1)
			break
	new_string = "".join(lst)
	final.write(">" + "chr" + str(chr1))
	final.write('\n')
	final.write(new_string)
	final.close()
ID.close()
#####################################################



#####################################################
###Alignment
num=$(cat /home/samplesID.txt)

for ID in $num
do

hisat2-build /home/${ID}_person.fa /home/${ID}/${ID}_index
hisat2 -p 4 -x /home/${ID}/${ID}_index -1 /RNAseq/${ID}/*_1.fastq.gz -2 /RNAseq/${ID}/*_2.fastq.gz -S /home/${ID}_hisat2.sam
samtools sort -@ 2 -O bam -o /home/${ID}_hisat2_sort.bam /home/${ID}_hisat2.sam

done
#####################################################



#####################################################
###Construction of individual-specific transcriptome and genomic variation data
#Construct the chain files

num=$(cat /home/samplesID.txt)
for ID in $num
do
	for i in {1..22}
	do
	g2gtools vcf2chain -f /home/GRCh37.p13.genome_chr${i}.fa -i /home/${ID}_indel.vcf -s HUMAN -o /home/${ID}/${ID}_chr${i}.chain
	done
done

num=$(cat /home/samplesID.txt)
for ID in $num
do
cd /home/${ID}/
ls | sort -V | xargs cat > /home/${ID}_all.chain
done



#Convert the reference transcriptome data into individual-specific transcriptome data 

num=$(cat /home/samplesID.txt)
for ID in $num
do

liftOver -gff /home/gencode.gtf /home/${ID}_all.chain /home/${ID}_lift.gtf unMapped

done

num=$(cat /home/samplesID.txt)
for ID in $num
do

cat /home/${ID}_lift.gtf | awk -F "\t" '{if($3~"exon")print$1"\t"$4"\t"$5"\t"$7"\t"$9}' > /home/${ID}_lift_exon.gtf

done



#Concatenate the transcriptome data before and after conversion

library(dplyr)
gtf1<-read.table("/home/gencode.gtf",header = F,sep="\t",quote = "",fill = T)
id<-read.table("/home/samplesID.txt",header = F,sep="\t",quote = "",fill = T)
for (j in 1:nrow(id))
{
path<-paste("/home/",as.character(id[j,1]),"_lift.gtf",sep="")
gtf2<-read.table(path,header = F,sep="\t",quote = "",fill = T)
result<-left_join(gtf1,gtf2,by=c("V1"="V1","V4"="V4","V5"="V5","V6"="V6"))
final<-na.omit(result)
path2<-paste("/home/",as.character(id[j,1]),"_lift_anno.gtf",sep="")
write.table(final,path2,sep="\t",quote=F,row.names=F,col.names=F)
}



#Convert indel positions into individual-specific

num=$(cat /home/samplesID.txt)
for ID in $num
do

CrossMap.py vcf /home/${ID}_all.chain /home/${ID}_indel.vcf /home/${ID}_person.fa /home/${ID}_indel_lift.vcf --no-comp-allele

done
#####################################################



#####################################################
###Identify exon extension/shrinkage events with associated indels
#Collect junction reads

num=$(cat /home/samplesID.txt)
for ID in $num
do
samtools view /home/${ID}_hisat2_sort.bam | awk -F "\t" '{if($6~"^[0-9]*M[0-9]*N[0-9]*M$") print}' | awk -F "\t" '{split ($6,T,"[M-N]");$6=T[1] OFS T[2] OFS T[3]}1' OFS="\t" | awk -F "\t" 'BEGIN{OFS="\t"}{$1=$1"_"FNR;$6=$4+$6-1;$7=$6+$7+1;$8=$7+$8-1}{print $0}' | awk -F "\t" '{print$1"\t"$3"\t"$4"\t"$6"\t"$7"\t"$8}' > /home/${ID}_mnm.bed
done



#Identify the use of novel splice sites

import sys
import re
import numpy as np
import pandas as pd
from tqdm import tqdm
from time import sleep

ID = open('/home/samplesID.txt')
for j in ID:
	j = j.rstrip()
	chrs = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22"]
	for chr in chrs :
		chr = "chr" + chr
		bed = open("/home/" + j + "_" + chr + "_mnm.bed")
		gtf = open("/home/" + j + "_lift_" + chr + ".gtf")
		outfile = open("/home/" + j + "/" + "site_"  + chr + ".bed",'w')
		
		d1 = {}
		d2 = {}
		for i in gtf.readlines():
			anno = re.split('[\t\n]',i)
			d1[anno[3]] = anno[6]
			d2[anno[4]] = anno[6]
		for v in bed.readlines():
			junction = re.split('[\t\n]',v)
			junc3="no"
			junc4="no"
			hako3 = []
			hako4 = []
			if (junction[3] in d2.keys()):
				junc3="yes"
				hako3.append(d2[junction[3]])
			if (junction[4] in d1.keys()):
				junc4="yes"
				hako4.append(d1[junction[4]])
			if junc3 == "yes" and junc4 == "yes" :
				outfile.writelines(str(junction[0]) + "\t" + str(junction[1]) + "\t" + str(junction[2]) + "\t" + str(junction[3]) + "\t" + str(junction[4]) + "\t" + str(junction[5]) + "\t" + str(hako4[0]) + "\t" + "yes,yes" + "\t" + "255,0,0" +"\n")
			if junc3 == "no" and junc4 == "yes" :
				outfile.writelines(str(junction[0]) + "\t" + str(junction[1]) + "\t" + str(junction[2]) + "\t" + str(junction[3]) + "\t" + str(junction[4]) + "\t" + str(junction[5]) + "\t" + str(hako4[0]) + "\t" + "no,yes"  + "\t" + "0,0,255"+ "\n")
			if junc3 == "yes" and junc4 == "no" :
				outfile.writelines(str(junction[0]) + "\t" + str(junction[1]) + "\t" + str(junction[2]) + "\t" + str(junction[3]) + "\t" + str(junction[4]) + "\t" + str(junction[5]) + "\t" + str(hako3[0]) + "\t" + "yes,no"  + "\t" + "0,0,255"+ "\n")
			if junc3 == "no" and junc4 == "no" :
				outfile.writelines(str(junction[0]) + "\t" + str(junction[1]) + "\t" + str(junction[2]) + "\t" + str(junction[3]) + "\t" + str(junction[4]) + "\t" + str(junction[5]) + "\t" + "*" + "\t" + "no,no" + "\t" + "255,192,0" + "\n")
		outfile.close()
ID.close()

num=$(cat /home/samplesID.txt)
for ID in $num
do
cat /home/${ID}/site_*.bed > /home/${ID}_site.bed
done



#Identify the novel splice site regions

from tqdm import tqdm
from time import sleep
import re
ID = open('/home/samplesID.txt')
for j in ID:
	j = j.rstrip()
	mess = open("/home/ + j + "_site.bed")
	outfile = open("/home/" + j + "_scmrange.bed",'w')
	for i in mess:
		i = i.rstrip()
		anno = i.split('\t')
		if (anno[7] == "yes,no" and anno[6] == "+"):
			scm1 = int(anno[4]) - 18
			scm2 = int(anno[4]) + 2
			no = int(anno[5]) - int(anno[4]) + 1
			no1 = int(anno[4])
			no2 = int(anno[5])
			nopos = int(anno[4])
			lab = "three"
			if (no > 5):
				outfile.write(str(anno[1]) + "\t" + str(scm1) + "\t" + str(scm2) + "\t" + str(anno[6]) + "\t" + str(anno[0]) + "\t" + str(no1) + "\t" + str(no2) + "\t" + str(anno[7]) + "\t" + str(nopos) + "\t" + str(lab))
				outfile.write('\n')
		if (anno[7] == "yes,no" and anno[6] == "-"):
			scm1 = int(anno[4]) - 6
			scm2 = int(anno[4]) + 2
			no = int(anno[5]) - int(anno[4]) + 1
			no1 = int(anno[4])
			no2 = int(anno[5])
			nopos = int(anno[4])
			lab = "five"
			if (no > 5):
				outfile.write(str(anno[1]) + "\t" + str(scm1) + "\t" + str(scm2) + "\t" + str(anno[6]) + "\t" + str(anno[0]) + "\t" + str(no1) + "\t" + str(no2) + "\t" + str(anno[7]) + "\t" + str(nopos) + "\t" + str(lab))
				outfile.write('\n')
		if (anno[7] == "no,yes" and anno[6] == "+"):
			scm1 = int(anno[3]) - 2
			scm2 = int(anno[3]) + 6
			no = int(anno[3]) - int(anno[2]) + 1
			no1 = int(anno[2])
			no2 = int(anno[3])
			nopos = int(anno[3])
			lab = "five"
			if (no > 5):
				outfile.write(str(anno[1]) + "\t" + str(scm1) + "\t" + str(scm2) + "\t" + str(anno[6]) + "\t" + str(anno[0]) + "\t" + str(no1) + "\t" + str(no2) + "\t" + str(anno[7]) + "\t" + str(nopos) + "\t" + str(lab))
				outfile.write('\n')
		if (anno[7] == "no,yes" and anno[6] == "-"):
			scm1 = int(anno[3]) - 2
			scm2 = int(anno[3]) + 18
			no = int(anno[3]) - int(anno[2]) + 1
			no1 = int(anno[2])
			no2 = int(anno[3])
			nopos = int(anno[3])
			lab = "three"
			if (no > 5):
				outfile.write(str(anno[1]) + "\t" + str(scm1) + "\t" + str(scm2) + "\t" + str(anno[6]) + "\t" + str(anno[0]) + "\t" + str(no1) + "\t" + str(no2) + "\t" + str(anno[7]) + "\t" + str(nopos) + "\t" + str(lab))
				outfile.write('\n')
	outfile.close()
ID.close()



#Identify the indels within novel splice site regions

num=$(cat /home/samplesID.txt)
for ID in $num
do
bedtools intersect -a /home/${ID}_scmrange.bed -b /home/${ID}_indel_lift.vcf -wa -wb > /home/${ID}_tmp.bed
bedtools intersect -a /home/${ID}_tmp.bed -b /home/${ID}_lift_exon.gtf -wa -wb > /home/${ID}_events.bed
rm /home/${ID}_tmp.bed
done
#####################################################



#####################################################
###Read coverage
#The number of junction reads ≥ 2
num=$(cat /home/samplesID.txt)
for ID in $num
do
cat /home/${ID}_events.bed | awk -F "\t" '{print$1"\t"$9"\t"$10"\t"$5"\t"$4"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17}' | sort | uniq | awk -F "\t" '{a[$2]++}END{for(as in a)if(a[as]>1)print as}' > /home/${ID}_tmp.bed
done

num=$(cat /home/samplesID.txt)
for id in $num
do
	num2=$(cat /home/${id}_tmp.bed)
	for ID in $num2
	do
		cat /home/${id}_events.bed | grep ${ID} >> /home/${id}_2reads.bed
	done
done
#####################################################



#####################################################
###Correlation between indels and events
#Restore the original position of the identified indels and calculate the extended or shrunken lengths 

library(dplyr)
library(tidyr)

id<-read.table("/home/samplesID.txt",header = F,sep="\t",quote = "",fill = T)
for (j in 1:nrow(id))
{
path1<-paste("/home/",as.character(id[j,1]),"_2reads.bed",sep="")
mess<-read.table(path1,header = F,sep="\t",quote = "",fill = T)
path2<-paste("/home/",as.character(id[j,1]),"_indel.vcf",sep="")
vcf<-read.table(path2,header = F,sep="\t",quote = "",fill = T)

mess_vcf<-left_join(mess,vcf,by=c("V1"="V1","V10"="V3","V11"="V4","V12"="V5"))
res<-mess_vcf[,c(1:12,26,13:15,17,18,19,20)]

qian1<-subset(res,res[,4]=="+" & res[,6]=="three")
qian2<-subset(res,res[,4]=="-" & res[,6]=="five")
hou1<-subset(res,res[,4]=="+" & res[,6]=="five")
hou2<-subset(res,res[,4]=="-" & res[,6]=="three")
qian<-rbind(qian1,qian2)
qian$bi<-qian[,15]
hou<-rbind(hou1,hou2)
hou$bi<-hou[,16]

result<-rbind(qian,hou)
result$noanno<-result[,5]-result[,21]

path3<-paste("/home/",as.character(id[j,1]),"_change.bed",sep="")
write.table(result,path3,sep="\t",quote=F,row.names=F,col.names=F)

}

cat /home/*_change.bed | sort | uniq > /home/all.bed



#Extract the identified candidative indels

num=$(cat /home/samplesID.txt)
for ID in $num
do

cat /home/${ID}_change.bed | awk -F "\t" '{print$1"_"$13"_"$10"_"$11"_"$12}' | sort | uniq > /home/${ID}_indel_candidate.bed

done

cat /home/*_indel_candidate.bed | sort | uniq > /home/indel_candidate.bed



#Determine whether the above candidate indels occur in each sample

num=$(cat /home/samplesID.txt)
for ID in $num
do
cat /home/indel_candidate.bed | while read line
do
n_snp=`grep $line /home/${ID}_indel.vcf | wc -l`
if [ ${n_snp} -eq 0 ] ; then
        echo -e "$line" >> /home/${ID}_noindel.txt
fi
done
done



#Determine the use of annotated splice sites and novel splice sites in each sample, ensuring that exon extension/shrinkage events can only be observed in the individuals with the candidate indels

library(dplyr)
library(tidyr)
res<-read.table("/home/all.bed",header = F,sep="\t",quote = "",fill = T)
qian1<-subset(res,res[,2]=="+" & res[,3]=="three")
qian2<-subset(res,res[,2]=="-" & res[,3]=="five")
hou1<-subset(res,res[,2]=="+" & res[,3]=="five")
hou2<-subset(res,res[,2]=="-" & res[,3]=="three")
qian<-rbind(qian1,qian2)
qian$bi<-qian[,6]
hou<-rbind(hou1,hou2)
hou$bi<-hou[,7]

id<-read.table("/home/samplesID.txt",header = F,sep="\t",quote = "",fill = T)
for (j in 1:nrow(id))
{
path1<-paste("/home/",as.character(id[j,1]),"_lift_anno.gtf",sep="")
gtf<-read.table(path1,header = F,sep="\t",quote = "",fill = T)

gtf1 <- gtf %>% separate(V5, c("gene","transcript"), "[;]")
gtf1$gene<-gsub("[^[:alnum:]///' ]","",gtf1$gene)
gtf1$gene<-gsub("geneid ","",gtf1$gene)
qian_pin<-left_join(qian,gtf1,by=c("V1"="V1","V4"="gene","bi"="V2"))
qian_pin$no<-qian_pin[,15] + qian_pin[,8]
qian_chu<-qian_pin[,c(1:5,9,15,17)]
hou_pin<-left_join(hou,gtf1,by=c("V1"="V1","V4"="gene","bi"="V3"))
hou_pin$no<-hou_pin[,16] + hou_pin[,8]
hou_chu<-hou_pin[,c(1:5,9,16,17)]

path2<-paste("/home/",as.character(id[j,1]),"_mnm.bed",sep="")
mnm<-read.table(path2,header = F,sep="\t",quote = "",fill = T)
res1<-NULL
for (i in 1:nrow(qian_chu)){
anno1<-length(which(mnm[,5]==qian_chu[i,7] & mnm[,2]==qian_chu[i,1]))
no1<-length(which(mnm[,5]==qian_chu[i,8] & mnm[,2]==qian_chu[i,1]))
tmp<-cbind(anno1,no1)
res1<-rbind(res1,tmp)
}
qian_final<-cbind(qian_chu,res1)

res2<-NULL
for (i in 1:nrow(hou_chu)){
anno2<-length(which(mnm[,4]==hou_chu[i,7] & mnm[,2]==hou_chu[i,1]))
no2<-length(which(mnm[,4]==hou_chu[i,8] & mnm[,2]==hou_chu[i,1]))
tmp2<-cbind(anno2,no2)
res2<-rbind(res2,tmp2)
}
hou_final<-cbind(hou_chu,res2)

l<-letters[1:10]
colnames(qian_final)<-l
colnames(hou_final)<-l
all<-rbind(qian_final,hou_final)
noexp<-all[which(all[,9]==0 & all[,10]==0),]

path3<-paste("/home/",as.character(id[j,1]),"_noindel.txt",sep="")
wu<-read.table(path3,header = F,sep="\t",quote = "",fill = T)
wu_chu<-left_join(wu,all,by=c("V1"="f"))
shan<-wu_chu[which(wu_chu[,10]!=0),]

path4<-paste("/home/",as.character(id[j,1]),"_count.bed",sep="")
write.table(all,path4,sep="\t",quote=F,row.names=F,col.names=F)

path5<-paste("/home/",as.character(id[j,1]),"_noexpression.bed",sep="")
write.table(noexp,path5,sep="\t",quote=F,row.names=F,col.names=F)

path6<-paste("/home/",as.character(id[j,1]),"_deletion.bed",sep="")
write.table(shan,path6,sep="\t",quote=F,row.names=F,col.names=F)
}

#Deletion of indels in /home/*_deletion.bed from candidative indels > /home/indel_candidate2.bed & /home/${ID}_indel_candidate2.bed
#####################################################



#####################################################
###The number of individuals
#For each candidate indels, when the number of individuals with observed exon extension/shrinkage events is greater than the number of individuals without, retain the indels (Individuals with no gene expression are not counted)

library(dplyr)
library(tidyr)

id<-read.table("/home/samplesID.txt",header = F,sep="\t",quote = "",fill = T)
for (j in 1:nrow(id))
{
path1<-paste("/home/",as.character(id[j,1]),"_noexpression.bed",sep="")
no<-read.table(path1,header = F,sep="\t",quote = "",fill = T)

path2<-paste("/home/",as.character(id[j,1]),"_indel.vcf",sep="")
indel<-read.table(path2,header = F,sep="\t",quote = "",fill = T)
no_chu<-as.data.frame(no[,6])

colnames(no_chu)<-"a"
colnames(indel)<-"a"
res<-setdiff(indel,no_chu)
res$lab<-"snp"

path3<-paste("/home/",as.character(id[j,1]),"_indel_candidate2.bed",sep="")
scm<-read.table(path3,header = F,sep="\t",quote = "",fill = T)

for (i in 1:nrow(scm)){

res[which(res[,1]==scm[i,1]),2]<-"scm"

}

path4<-paste("/home/",as.character(id[j,1]),"_tmp.vcf",sep="")
write.table(res,path4,sep="\t",quote=F,row.names=F,col.names=F)
}

cat /home/indel_candidate2.bed | while read line

do

n_all=`grep $line /home/*_tmp.vcf | wc -l`
n_scm=`grep $line /home/*_tmp.vcf | grep 'scm' | wc -l`
n_snp=`grep $line /home/*_tmp.vcf | grep 'snp' | wc -l`

echo -e "$line\t$n_all\t$n_snp\t$n_scm"

done

sh above.sh > /home/indel_candidate2_countsamples.bed
cat /home/indel_candidate2_countsamples.bed | awk -F "\t" '$4>$3{print$0}' > indel_result.bed
#####################################################

