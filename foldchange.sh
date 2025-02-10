#!/bin/bash

#1. IMPORT SEQ DATA INTO tempfiles DIRECTORY
####Make a check for if there are files aready in tempdata or files in directory?
echo "Importing seq data..."
# mkdir seqdata
mkdir temp
# cp -u /localdisk/data/BPSM/ICA1/fastq/* ./seqdata
awk 'BEGIN{FS="\t";}{print $1,$2,$3,$4,$5}' ./seqdata/Tco2.fqfiles > ./temp/report.txt
echo "Imported data and created tempdir and base report.txt"



#2. PERFORM fastqc QUALITY CHECK ON COMPRESSED fastq PAIRED END RAW SEQUENCE DATA 
##2a. Gather end1 and end2 column data (lists available fastq files) from Tco2.fqfiles
echo "Gathering end1 and end2 paired end column data from Tco2..."
touch ./temp/fastqlist_end1.txt
awk 'BEGIN{FS="\t";}{if($6 != "End1"){print $6}}' ./seqdata/Tco2.fqfiles | cut -d "." -f 1 > ./temp/fastqlist_end1.txt
touch ./temp/fastqlist_end2.txt
awk 'BEGIN{FS="\t";}{if($7 != "End2"){print $7}}' ./seqdata/Tco2.fqfiles | cut -d "." -f 1 > ./temp/fastqlist_end2.txt


##2b. end 1 fastqc on all files in fastqlist.txt, then count number of pass/fail/warn, then append to results
mkdir fastqcreport
mkdir temp
echo "end1_fastqc_pass" > ./temp/end1_pass.txt
echo "end1_fastqc_fail" > ./temp/end1_fail.txt
echo "end1_fastqc_warn" > ./temp/end1_warn.txt

while read line
do
echo "fastqc on ${line}..."
fastqc -o ./fastqcreport --extract ./seqdata/${line}.fq.gz
awk 'BEGIN{FS="\t";}{if($1 == "PASS"){print $1}}' ./fastqcreport/${line}_fastqc/summary.txt | wc -l >> ./temp/end1_pass.txt
awk 'BEGIN{FS="\t";}{if($1 == "FAIL"){print $1}}' ./fastqcreport/${line}_fastqc/summary.txt | wc -l >> ./temp/end1_fail.txt
awk 'BEGIN{FS="\t";}{if($1 == "WARN"){print $1}}' ./fastqcreport/${line}_fastqc/summary.txt | wc -l >> ./temp/end1_warn.txt

done < ./temp/fastqlist_end1.txt

echo "End1 fastqc done"


##2c. same as 2b. but with end2
echo "end2_fastqc_pass" > ./temp/end2_pass.txt
echo "end2_fastqc_fail" > ./temp/end2_fail.txt
echo "end2_fastqc_warn" > ./temp/end2_warn.txt

while read line
do
echo "fastqc on ${line}..."
fastqc -o ./fastqcreport --extract ./seqdata/${line}.fq.gz
awk 'BEGIN{FS="\t";}{if($1 == "PASS"){print $1}}' ./fastqcreport/${line}_fastqc/summary.txt | wc -l >> ./temp/end2_pass.txt
awk 'BEGIN{FS="\t";}{if($1 == "FAIL"){print $1}}' ./fastqcreport/${line}_fastqc/summary.txt | wc -l >> ./temp/end2_fail.txt
awk 'BEGIN{FS="\t";}{if($1 == "WARN"){print $1}}' ./fastqcreport/${line}_fastqc/summary.txt | wc -l >> ./temp/end2_warn.txt
done < ./temp/fastqlist_end2.txt

echo "End2 fastqc done"


##2d. append all counts onto report.txt
echo "Appending counts onto ./temp/reportfinal.txt..."

paste ./temp/report.txt ./temp/end1_pass.txt | paste - ./temp/end1_fail.txt | paste - ./temp/end1_warn.txt |paste - ./temp/end2_pass.txt |paste - ./temp/end2_fail.txt |paste - ./temp/end2_warn.txt > ./temp/reportfinal.txt

cp ./temp/reportfinal.txt ./RESULTS/fastqcSummaryReport.txt

##############user input for no. of fails threshold to remove seq data?
##############flag bad seq data and input base no. for trimming?

#3. IMPORT T.congo GENOME SEQ AND .bed FILE INTO refseqdata and uncompress .gz
echo "Importing refseq data..."

# mkdir refseqdata
cp -u /localdisk/data/BPSM/ICA1/Tcongo_genome/* ./refseqdata/
gzip -d ./refseqdata/TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta.gz
cp -ur /localdisk/data/BPSM/ICA1/* ./refseqdata/

echo "Imported refseq data"


#4. ALIGN READ PAIRS ONTO REFERENCE GENOME
##4a. Reformat refseq into single line sequences, then select chromosome sequences only
############make check on refseq format?
echo "Reformatting refseq fasta into single line fasta. Then select for chromosome sequences only"

awk '/^>/ {print (NR>1?"\n":"")$0;;next}{printf "%s",$0;} END{print "";}' ./refseqdata/TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta > ./refseqdata/genomeonlyref.fasta
grep -A1 "SO=chromosome" ./refseqdata/genomeonlyref.fasta > ./refseqdata/genomeonlyref1.fasta
##4b. Build bowtie2 index from ref chromosome sequences only 
echo "Building bowtie2 index for refseq fasta"
bowtie2-build -f ./refseqdata/genomeonlyref1.fasta ./refseqdata/refbowtieindex
##4c. Align read pairs via bowtie2
echo "Aligning read pairs via bowtie 2"

cut -d "_" -f 1 ./temp/fastqlist_end1.txt > ./temp/fastqlist_base.txt #gather gene base names
mkdir ./aligned

while read base
do
echo "${base} alignment"
bowtie2 -x ./refseqdata/refbowtieindex -q -1 ./seqdata/${base}_1.fq.gz -2 ./seqdata/${base}_2.fq.gz > ./aligned/${base}.sam
done < ./temp/fastqlist_base.txt

echo -e "\nAll aligned!"

##4d. Convert .sam files from bowtie2 to .bam via samtools
echo "Converting .sam files into .bam files..."

while read base
do
echo "${base} conversion to .bam"
samtools view -b -o ./aligned/${base}.bam ./aligned/${base}.sam
done < ./temp/fastqlist_base.txt
echo "All converted!"

##4e. Generate index .bai files 
echo "Generating index .bai files..."
while read base
do
echo "${base} indexing for .bai"
samtools sort --output-fmt BAM ./aligned/${base}.bam > ./aligned/${base}_sort.bam
samtools index -b ./aligned/${base}_sort.bam
done < ./temp/fastqlist_base.txt
echo "Generated!"

#5. GENERATE COUNTS DATA VIA bedtools FROM REFSEQ .bed AND .bam 
echo "Generating counts data..."
mkdir counts
while read base
do
echo "${base} counts"
bedtools coverage -counts -a ./refseqdata/TriTrypDB-46_TcongolenseIL3000_2019.bed -b ./aligned/${base}_sort.bam > ./counts/${base}_counts.txt
done < ./temp/fastqlist_base.txt

echo "Generated!"

#ASK FOR FAIL THRESHOLD
failtheshold=3
echo "Enter threshold for max no. ,inclusive, of fails allowed in sequencing quality report. (Input integer)"
echo "Will remove samples will more fails than threshold"
echo "(Default fail threshold = 3)"
read failthreshold

#1. Mean the replicates for each condition
##1a. create list of all unique conditions
echo "Creaing conditions list..."
touch ./temp/replicateindex.txt
awk 'BEGIN{FS="\t"; OFS="\t"}{if(NR!=1){print $2,$4,$5}}' ./seqdata/Tco2.fqfiles | 
sort - | 
uniq - > ./temp/replicateindex.txt

##1b. create reference index
echo "Creating reference index..."
touch ./temp/replicateindexref.txt
awk 'BEGIN{FS="\t";}{if(NR!=1){print $1,$2,$4,$5}}' ./seqdata/Tco2.fqfiles > ./temp/replicateindexref.txt

##1bi. Removing samples with low quality, higher fails in quality check than threshold.
echo "Removing samples with low quality sequencing in either/or end1, end2"

awk 'BEGIN{FS="\t"}{if(NR != 1){print $0}}' ./temp/reportfinal.txt > ./temp/reportfinalnocol.txt
awk -v thresh="$failthreshold" 'BEGIN{FS="\t"}{if($3 >= thresh || $6 >= thresh){print $0}}' ./temp/reportfinalnocol.txt > ./temp/failthresholdrm.txt
cut -d " " -f 1 ./temp/failthresholdrm.txt > ./temp/failthresholdrm1.txt


grep -v -f ./temp/failthresholdrm1.txt ./temp/replicateindexref.txt > ./temp/replicateindexreftemp.txt
mv -f ./temp/replicateindexreftemp.txt ./temp/replicateindexref.txt

##1c. separate search into 3 files for each condition
echo "Creating search index..."
awk 'BEGIN{FS="\t";}{print $1;}' ./temp/replicateindex.txt | sort - | uniq - > ./temp/search1.txt
awk 'BEGIN{FS="\t";}{print $2;}' ./temp/replicateindex.txt | sort - | uniq - > ./temp/search2.txt
awk 'BEGIN{FS="\t";}{print $3;}' ./temp/replicateindex.txt | sort - | uniq - > ./temp/search3.txt

##1d. group replicates together
echo "Grouping replicates..."

mkdir ./temp/search
while read search1
do
grep -w ${search1} ./temp/replicateindexref.txt > ./temp/search1temp.txt

  while read search2
  do
  grep -w ${search2} ./temp/search1temp.txt > ./temp/search2temp.txt
  
    while read search3
    do
    grep -w ${search3} ./temp/search2temp.txt > ./temp/search/${search1}_${search2}_${search3}
    echo "${search1}_${search2}_${search3} groups created"
    
    done < ./temp/search3.txt
  done < ./temp/search2.txt
done < ./temp/search1.txt

rm -f ./temp/search/*_0_Induced


##1e.gather only counts of each gene of each sample
echo "Gathering counts for each gene of each sample..."

ls ./temp/search > ./temp/searchdirindex.txt
mkdir ./counts/groups

while read file
do
echo "Gathering ${file}..."
cat ./temp/search/${file} | awk '{FS=" ";}{print $1}' - > ./temp/repgroups.txt

  while read repgroups
  do
  echo "Gathering ${repgroups}..."
  cat ./counts/${repgroups:0:3}-${repgroups:3}_counts.txt | awk '{FS="\t";}{if($6 != "hot"){print $6}}' - > ./counts/groups/${repgroups:0:3}-${repgroups:3}_counts_only.txt
  
  done < ./temp/repgroups.txt
done < ./temp/searchdirindex.txt


#1f. append counts of each condition to one file
echo "Appending counts of each condition to one file"
mkdir ./counts/reptotal

while read file
do
cat ./temp/search/${file} | awk '{FS=" ";}{print $1}' - > ./temp/repgroups.txt
>./counts/reptotal/${file}_all.txt

  while read repgroups
  do
  paste ./counts/reptotal/${file}_all.txt ./counts/groups/${repgroups:0:3}-${repgroups:3}_counts_only.txt > ./counts/reptotal/${file}_all1.txt && mv -f ./counts/reptotal/${file}_all1.txt ./counts/reptotal/${file}_all.txt

  echo "Gathering: ./counts/${repgroups:0:3}-${repgroups:3}_counts.txt"
  echo "Group: ${file}"
  echo "Into: ./counts/reptotal/${file}_all.txt"
  echo "----------------"
  
  done < ./temp/repgroups.txt
done < ./temp/searchdirindex.txt

#1g. get means of each condition
echo "Getting means of replicates of each condition..."

while read file
do
echo "Gathering ${file}..."
cat ./temp/search/${file} | awk 'BEGIN{FS=" ";}{print $1}' - > ./temp/${file}_all.txt_repgroups.txt
done < ./temp/searchdirindex.txt

ls ./counts/reptotal > ./temp/reptotalindex.txt
mkdir ./counts/replicatemeans

while read reptotalindex
do
sum=0
denom=0
>./counts/replicatemeans/${reptotalindex}

  while read gene
  do
  sum=0
  denom=0
   
    for replicate in ${gene}
    do
    sum=$(($sum + $replicate))
    denom=$((${denom} + 1))
    done
    
  mean=echo "${sum} / ${denom}" | bc -l
  echo "${sum} / ${denom}" | bc -l >> ./counts/replicatemeans/${reptotalindex}
  echo "Counts: ${gene}"
  echo "Sum: $sum"
  echo "No. of replicates: $denom"
  echo "Mean:"
  echo ${mean}
  echo "Repindex: ${reptotalindex}"
  echo "------------------------------"
  
  done < ./counts/reptotal/${reptotalindex}
done < ./temp/reptotalindex.txt


#2. Append meancounts onto .bedfile 
mkdir ./RESULTS
mkdir ./RESULTS/MeanReplicateCounts

#2a. create base .bed file
echo "Creating base .bed file..."

> ./refseqdata/base.bed 
awk 'BEGIN{FS="\t";}{print $4,$5}' ./refseqdata/TriTrypDB-46_TcongolenseIL3000_2019.bed >> ./refseqdata/base.bed

while read reptotalindex
do
echo "Creating final ${reptotalindex} mean counts file..."
paste ./refseqdata/base.bed ./counts/replicatemeans/${reptotalindex} > ./RESULTS/MeanReplicateCounts/${reptotalindex}

done < ./temp/reptotalindex.txt


#3 Onto next step 3


#1.User defined inputs
#condition1
echo "--------------------------------------------------------------------------------------------------------------------------"

echo -e "Please use format (SampleType Time Treatment) and exact names as seen on Tco2.fqfiles\nReference Condition:"
read refsample reftime reftreat

#condition2
echo -e "Please use format (SampleType Time Treatment), only vary one field from Reference Condition, and exact names as seen on Tco2.fqfiles\nQuery Condition:"
read quersample  quertime  quertreat

#range
echo -e "Please use one exact column name from Tco2.fqfiles (SampleType OR Time OR Treatment)\nRange:"
read range

echo -e "${refsample}, ${reftime}, ${reftreat}\n${quersample}, ${quertime}, ${quertreat}\n ${range}"

echo "--------------------------------------------------------------------------------------------------------------------------"
#2. Fold change between reference and query
mkdir foldchange
>./temp/pastefoldchange.txt

#2a. create index for condi ranges
echo "Creating index for condition ranges..."

>./temp/rangecompare.txt

if test ${range} == "SampleType"
then
  cp -f ./temp/search1.txt ./temp/rangecondi.txt
  
    while read rangecondi
    do
    echo ${rangecondi}_${quertime}_${quertreat} >> ./temp/rangecompare.txt
  
    done < ./temp/rangecondi.txt
  
elif test ${range} == "Time"
then
  cp -f ./temp/search2.txt ./temp/rangecondi.txt
  
    while read rangecondi
    do
    echo ${quersample}_${rangecondi}_${quertreat} >> ./temp/rangecompare.txt
  
    done < ./temp/rangecondi.txt
  
elif test ${range} == "Treatement"
then
  cp -f ./temp/search3.txt ./temp/rangecondi.txt
  
    while read rangecondi
    do
    echo ${quersample}_${quertime}_${rangecondi} >> ./temp/rangecompare.txt
  
    done < ./temp/rangecondi.txt
else
  echo "Please input valid range."
fi



#2b. reference index made by 2a to create fold changes for range then dividing
echo "Calculating fold changes..."

while read rangecondi
do
paste ./counts/replicatemeans/${refsample}_${reftime}_${reftreat}_all.txt ./counts/replicatemeans/${rangecondi}_all.txt > ./temp/pastefoldchange.txt
>./foldchange/${refsample}_${reftime}_${reftreat}_vs_${rangecondi}.txt
echo "Comparing ${rangecompare}"

  while read ref query 
  do
  echo -e "ReferenceVal:${ref} \nQueryVal:${query}"
  
  
    if test $ref -eq 0
    then
      #awk 'BEGIN{FS="\t";}{$(($2 / ($1 + 0.1)))}'>> ./foldchange/${refsample}_${reftime}_${reftreat}_vs_${quersample}_${quertime}_${quertreat}.txt
      echo "${query} / (${ref} + 1)" | bc -l >> ./foldchange/${refsample}_${reftime}_${reftreat}_vs_${rangecondi}.txt
      echo "${rangecondi} (ref==0)"
      
    else
      #awk 'BEGIN{FS="\t";}{$(($2 / $1))}'>> ./foldchange/${refsample}_${reftime}_${reftreat}_vs_${quersample}_${quertime}_${quertreat}.txt
      echo "${query} / ${ref}" | bc -l >> ./foldchange/${refsample}_${reftime}_${reftreat}_vs_${rangecondi}.txt
      echo "${rangecondi} (ref!=0)"
    fi
    
  done < ./temp/pastefoldchange.txt
done < ./temp/rangecompare.txt

#3. Paste them into .bedfile to create final fold changes
echo "Pasting fold changes onto .bed file into RESULTS/FoldChanges dir"
mkdir ./RESULTS/FoldChanges

while read rangecompare
do

paste ./refseqdata/TriTrypDB-46_TcongolenseIL3000_2019.bed ./foldchange/${refsample}_${reftime}_${reftreat}_vs_${rangecompare}.txt > ./RESULTS/FoldChanges/${refsample}_${reftime}_${reftreat}_vs_${rangecompare}.txt > ./DE_genes.txt

done < ./temp/rangecompare.txt

Rscript ./GOI_crawler.R