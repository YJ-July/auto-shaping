---
title: "Analysis"
format: html
editor: visual
---

BmN4 cells

Sequence analysis of RNA libraries

The RNA libraries used in this section comes from Naive13 and Naive20. For clarity during data processing, they are referred to as close13 and close20, respectively.

Mapping to TEs

Due to the large file size, a downsampling process was performed to reduce its size.

{bash}
NN="1 2"
for N in ${NN};do
    echo ${N}
    seqtk sample -s100 ./fastq/close20_R${N}.fastq.gz 0.1 > ./fastq/close20_sampling_R${N}.fastq
    seqtk trimfq ./fastq/close20_sampling_R${N}.fastq -L 100 > ./fastq/close20_sampling_100nt_R${N}.fastq
    gzip ./fastq/close20_sampling_100nt_R${N}.fastq > ./fastq/close20_sampling_100nt_R${N}.fastq.gz
done

{bash}
MM="close20_sampling_100nt close13"
LL="transposon"

hisat2-build -f ./fasta/transposon.fasta transposon

for M in ${MM}; do
  for L in ${LL}; do
    hisat2 -x ${L} -1 ./fastq/${M}_R1.fastq.gz  -2 ./fastq/${M}_R2.fastq.gz -k 3 -p 10 -S ./sambam/${M}_${L}_RNA_HS2.sam
    samtools view -@ 10 -bS ./sambam/${M}_${L}_RNA_HS2.sam > ./sambam/${M}_${L}_RNA_HS2.bam
    bamToBed -i ./sambam/${M}_${L}_RNA_HS2.bam > ./bed/${M}_${L}_RNA_HS2.bed
done
done 

Changing bed to rds

Here, all bed files are converted to the rds format for storage, as reading rds files in the R environment is faster, facilitating subsequent analyses.

{r BEDtoRDS.R}
args <- commandArgs(trailingOnly = TRUE) bed_file <- args[1]

bed <- read.table(sprintf("./bed/%s.bed",bed_file))

if (!dir.exists("rds")) {dir.create("rds")} saveRDS(bed,sprintf("./rds/%s.rds",bed_file))

{bash}
MM="close20_sampling_100nt close13"
LL="transposon"
Rscript BEDtoRDS.R ${M}_${L}

The RPM of RNAs for each TE

{r RNA_RPM.R}
args <- commandArgs(trailingOnly = TRUE)
sn <- as.character(args[1])

suppressMessages(library("Biostrings"))
suppressMessages(library("tidyr"))
trans <- readDNAStringSet("./fasta/transposon.fasta")

rawdata <- read.table(sprintf("./bed/%s_RNA_HS2.bed",sn))

##### length table ####
rawdata$len <- rawdata$V3 - rawdata$V2
rawdata$V1 <- factor(rawdata$V1, levels = names(trans))
rawdata$V2 <- factor(rawdata$V2 + 1 , levels = 1:max(width(trans)))
rawdata$V3 <- factor(rawdata$V3, levels = 1:max(width(trans)))

##### defined mRNA ####
mRNA <- subset(rawdata, len == 100)

##### count numbers of reads ####
sumreads <- length(unique(mRNA$V4))
message("Total reads of ",sn," is ", sumreads )

##### table of reads distribution ####
RNA_reads <- table(mRNA[, 1])


##### count rpm ####
RNA_rpm <- (RNA_reads / sumreads * 10^6) %>%
  as.data.frame
colnames(RNA_rpm) <- c("trans","rpm")

RNA_rpm$rpkm <- RNA_rpm$rpm / width(trans[RNA_rpm$trans]) * 10^3

##### save files ####
saveRDS(RNA_rpm, sprintf("./results/%s_RNA_rpm.rds",sn))

{bash}
MM="close20_sampling_100nt close13"
LL="transposon"
for M in ${MM}; do
    for L in ${LL}; do
    Rscript ./funs/RNA_RPM.R ${M}_${L}
    done
done

The abundance of RNA obtained are as follows:

Library

Reads

Naive13

2250194

Naive20

5338391

Differential expression analysis in RNA-seq

The RPM calculated in the previous section @sec-rpmrna is required here for this step.

{r}
rna_rpm_13 <- readRDS("./results/close13_transposon_RNA_rpm.rds")
rna_rpm_20 <- readRDS("./results/close20_sampling_100nt_transposon_RNA_rpm.rds")

##### count the mean of rpkm ####
rna_rpm <- as.data.frame(rna_rpm_13[,1])
rna_rpm$M <- log2(rna_rpm_13$rpm + 0.01) - log2(rna_rpm_20$rpm + 0.01)
rna_rpm$A <- log2(rna_rpm_13$rpm + 0.01) / 2 + log2(rna_rpm_20$rpm + 0.01) / 2

##### save files ####
saveRDS(rna_rpm,"./results/rna_rpm.rds")

Sequence analysis of piRNA libraries

The RNA libraries used in this section comes from Naive13, Naive20, dsRluc13 and dsRluc20.

Removing adapter

{bash}
MM="naive dsRluc"
NN="13 20"

for M in ${MM}; do
    for N in ${NN}; do
        echo $N
        cutadapt -a TGGAATTCTCGGGTGCCAAG --minimum-length 20 -o tmp_${M}${N}_trim.fastq ./fastq/${M}${N}.fastq.gz
        fastq_to_fasta -Q 33 -i tmp_${M}${N}_trim.fastq -o tmp_${M}${N}_trim.fasta
        fastx_collapser < tmp_${M}${N}_trim.fasta > tmp_${M}${N}_trim_unique.fasta
        cutadapt -u4 -u -4 --minimum-length 20 -o ${M}${N}_trim.fasta tmp_${M}${N}_trim_unique.fasta
        gzip ./fasta/${M}${N}_trim.fasta
        rm tmp_${M}${N}_trim.fastq
    done
done

Mapping

Mapping to TEs

{bash}
bowtie-build -f ./fasta/transposon.fasta transposon
MM="naive dsRluc"
NN="13 20"
LL="transposon"
for M in ${MM}; do
    for N in ${NN}; do
        for L in ${LL}; do
            echo "=========="$M$N"----->"$L"=========="   
            ( bowtie --offrate 3 -p 10 -a --best --strata -v 3 -t --sam ${L} -f ./fasta/${M}${N}_trim.fasta.gz  > ${M}${N}_${L}.sam ) 
            samtools view -@ 20  -bS ${M}${N}_${L}.sam > ${M}${N}_${L}.bam
            bamToBed -i ${M}${N}_${L}.bam > ./bed/${M}${N}_${L}.bed
            rm ${M}${N}_${L}.sam
            rm ${M}${N}_${L}.bam
        done
    done
done
rm *ebwt

Mapping to the genome

{bash}
bowtie-build -f /home/Jie/Saturn/raw/KWMTBOMO.nucl.fa KWMTBOMO
MM="naive dsRluc"
NN="13 20"
LL="KWMTBOMO"
for M in ${MM}; do
    for N in ${NN}; do
        for L in ${LL}; do
            echo "=========="$M$N"----->"$L"=========="   
            ( bowtie --offrate 3 -p 10 -a --best --strata -v 3 -t --sam ${L} -f ./fasta/${M}${N}_trim.fasta.gz  > ${M}${N}_${L}.sam ) 
            samtools view -@ 20  -bS ${M}${N}_${L}.sam > ${M}${N}_${L}.bam
            bamToBed -i ${M}${N}_${L}.bam > ./bed/${M}${N}_${L}.bed
            rm ${M}${N}_${L}.sam
            rm ${M}${N}_${L}.bam
        done
    done
done
rm *ebwt

Mapping to Fem/Masc

Sequence of Fem/Masc:

Fem RNA: TTTCATTGTTACCTCTTTTTGTCAATTCATAAAGTCATTCAGTG

Masc RNA: AAATGGCTTTGTGAATCGACAAAAAGAGGTAACAATTGAAGCTAATCAGAAGAAA

{bash}
bowtie-build ./fasta/femmasc.fas femmasc
MM="naive dsRluc"
NN="13 20"
LL="femmasc"
for M in ${MM}; do
    for N in ${NN}; do
        for L in ${LL}; do
            echo "=========="$M$N"----->"$L"=========="
            ( bowtie --offrate 3 -p 10 -a --best --strata -v 3 -t --sam ${L} -f ./fasta/${M}${N}_trim.fasta.gz  > ${M}${N}_${L}.sam ) 
            samtools view -@ 20  -bS ${M}${N}_${L}.sam > ${M}${N}_${L}.bam
            bamToBed -i ${M}${N}_${L}.bam > ./bed/${M}${N}_${L}.bed
            rm ${M}${N}_${L}.sam
            rm ${M}${N}_${L}.bam
        done
    done
done
rm *ebwt

Changing bed to rds

{bash}
MM="naive dsRluc"
NN="13 20"
LL="femmasc transposon KWMTBOMO"
for M in ${MM}; do
    for N in ${NN}; do
        for L in ${LL}; do
        Rscript ./funs/BEDtoRDS.R ${M}${N}_${L}
        done
    done
done

Distribution of smRNAs length

The total amount of all pre-mapping smRNA is calculated from the code below.

{r}
library(Biostrings)
library(tidyr)

q13 <- readDNAStringSet("/home/Jie/Saturn/raw/210406_Izumisan_BmQ/210406_C01_BmQ2013_total_trim.fasta.gz")
q20 <- readDNAStringSet("/home/Jie/Saturn/raw/210406_Izumisan_BmQ/210406_C02_BmQ2020_total_trim.fasta.gz")
ds13 <- readDNAStringSet("/home/Jie/Saturn/raw/210406_Izumisan_BmQ/210406_M01_dsRluc2013_total_trim.fasta.gz")
ds20 <- readDNAStringSet("/home/Jie/Saturn/raw/210406_Izumisan_BmQ/210406_L01_dsRluc_total_trim.fasta.gz")

total_reads <- c(length(q13),length(q20),length(ds13),length(ds20))
names(total_reads) <- c("naive13","naive20","ds13","ds20")

The abundance obtained are as follows:

Library

Reads

Naive13

23663825

Naive20

24052105

dsRluc13

31579821

dsRluc20

26234597

Calculate the length distribution of small RNAs in each library.

{r}
lens <- data.frame(len = 20:43)
lens$q13 <- q13 %>%
  width %>%
  factor(levels = 20:43) %>%
  table %>%
  as.matrix
lens$q20 <- q20 %>%
  width %>%
  factor(levels = 20:43) %>%
  table %>%
  as.matrix
lens$ds13 <- ds13 %>%
  width %>%
  factor(levels = 20:43) %>%
  table %>%
  as.matrix
lens$ds20 <- ds20 %>%
  width %>%
  factor(levels = 20:43) %>%
  table %>%
  as.matrix
lens_rpm <- lens[,2:5] %>%
  sweep(.,2,total_reads,"/") %>%
  {.*10^6} 
lens_rpm$lens <- 20:43

save(list = c("lens_rpm"),file="./results/len_distribution.Rdata")

Amount of mapping to TEs and mapping to the genome.

Calculate the proportion of mapping to TE and the genome. See @sec-length for obtaining the total reads.

{r}
library(Biostrings)
library(tidyr)
total_reads <- c(23663825,24052105,31579821,26234597)
names(total_reads) <- c("naive13","naive20","ds13","ds20")

## TE
q13 <- read.table("./bed/naive13_transposon.bed")
q20 <- read.table("./bed/naive20_transposon.bed")
ds13 <- read.table("./bed/dsRluc13_transposon.bed")
ds20 <- read.table("./bed/dsRluc20_transposon.bed")

q13_reads <- length(unique(q13$V4))
q20_reads <- length(unique(q20$V4))
ds13_reads <- length(unique(ds13$V4))
ds20_reads <- length(unique(ds20$V4))

te_reads <- c(q13_reads,q20_reads,ds13_reads,ds20_reads)
te_rpm <- te_reads %>%
  {./total_reads*10^6} %>%
  data.frame(rpm=.) %>%
  cbind(mapping = "te",
        sample = names(total_reads))

## genome
q13geno <- read.table("./bed/naive13_KWMTBOMO.bed")
q20geno <- read.table("./bed/naive20_KWMTBOMO.bed")
ds13geno <- read.table("./bed/dsRluc13_KWMTBOMO.bed")
ds20geno <- read.table("./bed/dsRluc20_KWMTBOMO.bed")

q13geno_reads <- length(unique(q13geno$V4))
q20geno_reads <- length(unique(q20geno$V4))
ds13geno_reads <- length(unique(ds13geno$V4))
ds20geno_reads <- length(unique(ds20geno$V4))

genome_reads <- c(q13geno_reads,q20geno_reads,ds13geno_reads,ds20geno_reads)
genome_rpm <- genome_reads %>%
  {./total_reads*10^6} %>%
  data.frame(rpm=.) %>%
  cbind(mapping = "genome",
        sample = names(total_reads))

save(list=c("te_rpm","genome_rpm"),file = "./results/bmn4_reads.Rdata")

The RPM of piRNAs for each TE

{r piRNA_RPM.R}
##### load args & packages & files ####
args <- commandArgs(trailingOnly = TRUE)
sn <- as.character(args[1])

suppressMessages(library("Biostrings"))
trans <- readDNAStringSet("./fasta/transposon.fasta.gz")

rawdata <- readRDS(sprintf("./rds/%s_transposon.rds",sn))

##### length table ####
rawdata$len <- rawdata$V3 - rawdata$V2
rawdata$V1 <- factor(rawdata$V1, levels = names(trans))
rawdata$V2 <- factor(rawdata$V2 + 1 , levels = 1:max(width(trans)))
rawdata$V3 <- factor(rawdata$V3, levels = 1:max(width(trans)))

##### defined piRNA ####
piRNA_plus <- subset(rawdata, 26 <= len & len <= 32 & V6 == "+")
piRNA_minus <- subset(rawdata, 26 <= len & len <= 32 & V6 == "-")

##### count numbers of reads ####
sum_plus <- length(unique(piRNA_plus$V4))
sum_minus <- length(unique(piRNA_minus$V4))
sumreads <- sum_plus + sum_minus
message("Total reads of ",sn," is ", sumreads )

##### table of reads distribution ####
plus_reads <- table(piRNA_plus[, c(1, 2)])
minus_reads <- table(piRNA_minus[, c(1, 3)])

##### count rpm ####
plus_rpm <- plus_reads / sumreads * 10 ^ 6
minus_rpm <- minus_reads / sumreads * 10 ^ 6

##### save files ####
saveRDS(plus_rpm, sprintf("./piRNA/%s_piRNA_plus_rpm.rds",sn))
saveRDS(minus_rpm, sprintf("./piRNA/%s_piRNA_minus_rpm.rds",sn))

{bash}
MM="naive dsRluc"
NN="13 20"
LL="transposon"
for M in ${MM}; do
    for N in ${NN}; do
        for L in ${LL}; do
        Rscript ./funs/piRNA_RPM.R ${M}${N}_${L}
        done
    done
done

The abundance of piRNAs obtained are as follows:

Library

Total reads

Naive13

8614065

Naive20

7911434

dsRluc13

11193224

dsRluc20

8290646

{r}
library(tidyr)
suppressMessages(library("Biostrings"))

trans <- readDNAStringSet("./fasta/transposon.fasta")

sample_names <- c("naive13","naive20","dsRluc13","dsRluc20")
names(sample_names) <- sample_names

## load raw data
plus_rpm_pos <- lapply(sample_names,function(sn){
  tmp <- sprintf("./piRNA/%s_transposon_piRNA_plus_rpm.rds",sn) %>%
    readRDS %>%
    as.data.frame
  colnames(tmp) <- c("trans","site","rpm")
  tmp[,"direc"] <- "+"
  tmp[,"samples"] <- sn
  tmp
})
minus_rpm_pos <- lapply(sample_names,function(sn){
  tmp <- sprintf("./piRNA/%s_transposon_piRNA_minus_rpm.rds",sn) %>%
    readRDS %>%
    as.data.frame
  colnames(tmp) <- c("trans","site","rpm")
  tmp[,"direc"] <- "-"
  tmp[,"samples"] <- sn
  tmp
})

## rpm for each position of each TE
rpm_pos <- rbind(Reduce(rbind,plus_rpm_pos),Reduce(rbind,minus_rpm_pos))
rpm_pos$site <- as.numeric(rpm_pos$site)
rpm_pos[rpm_pos$direc == "-","rpm"] <- rpm_pos[rpm_pos$direc == "-","rpm"] *-1
rpm_pos$rpkm <- rpm_pos$rpm / width(trans[rpm_pos$trans]) * 10^3

## rpm for each TE
rpm_total <- lapply(sample_names,function(sn){
  tmp <- subset(rpm_pos,samples==sn)
  tapply(tmp$rpm,tmp$trans,sum)}) %>%
  Reduce(cbind,.) %>%
  as.data.frame
colnames(rpm_total) <- sample_names
rpm_total$trans <- rownames(rpm_total) 
rpm_total <- pivot_longer(rpm_total,
                          cols = 1:4,
                          names_to = "samples",
                          values_to = "rpm")
rpm_total$rpkm <- rpm_total$rpm / width(trans[rpm_total$trans]) * 10^3

## save
save(list = c("rpm_pos"),file = "./results/rpm_TE_positions.Rdata")
save(list = c("rpm_total"),file = "./results/rpm_TE_total.Rdata")

Differential expression analysis in smRNA-seq

The RPM calculated in the previous section @sec-rpm is required here for this step.

{r}
load("./results/rpm_TE_total.Rdata")

sample_names <- c("naive13","naive20","dsRluc13","dsRluc20")
names(sample_names) <- sample_names

sample_set <- c(
  "naive13 vs dsRluc13",
  "naive20 vs dsRluc20",
  "naive13 vs naive20",
  "dsRluc13 vs dsRluc20"
)
names(sample_set) <- sample_set

MA <- lapply(sample_set, function(i) {
  # Split the string using " vs " as the delimiter and get the first and second parts
  aa <- unlist(strsplit(i, " vs "))[1]
  bb <- unlist(strsplit(i, " vs "))[2]
  rpm_aa <- subset(rpm_total, samples == aa)
  rpm_bb <- subset(rpm_total, samples == bb)
    # Create a data frame with the required columns
  df <- data.frame(
    "trans" = factor(rpm_aa$trans, levels = names(trans)),
    "M" = log2(rpm_aa$rpm + 0.01) - log2(rpm_bb$rpm + 0.01),
    "A" = log2(rpm_aa$rpm + 0.01) / 2 + log2(rpm_bb$rpm + 0.01) / 2,
    "samples" = i
  )
  df
}) %>%
  Reduce(rbind,.)
save(list = c("MA"), file = "./results/MA_4samples.Rdata")

Calculation of Dscore

Calculation of Dscore

{r}
##### define args ####
sample_names <- c("naive13", "naive20", "dsRluc13", "dsRluc20")
names(sample_names) <- sample_names

sample_set <-
  c(
    "naive13 vs dsRluc13",
    "naive20 vs dsRluc20",
    "naive13 vs naive20",
    "dsRluc13 vs dsRluc20"
  )
names(sample_set) <- sample_set

##### load files ####
plus_rpm <-
  lapply(sample_names, function(i) {
    readRDS(sprintf("./piRNA/%s_transposon_piRNA_plus_rpm.rds", i))
  })
minus_rpm <-
  lapply(sample_names, function(i) {
    readRDS(sprintf("./piRNA/%s_transposon_piRNA_minus_rpm.rds", i))
  })

##### rpkm per TE ####
plus_rpkm_perTE <- lapply(sample_names,
                          function(i) {
                            rpkm <- (plus_rpm[[i]] / width(trans)) * 10 ^ 3
                            rpkm_perTE <-
                              sweep(rpkm, 1, apply(rpkm, 1, sum), "/") * 100
                            rpkm_perTE[is.nan(rpkm_perTE)] <- 0
                            rpkm_perTE
                          })
minus_rpkm_perTE <- lapply(sample_names,
                           function(i) {
                             rpkm <- (minus_rpm[[i]] / width(trans)) * 10 ^ 3
                             rpkm_perTE <-
                               sweep(rpkm, 1, apply(rpkm, 1, sum), "/") * 100
                             rpkm_perTE[is.nan(rpkm_perTE)] <- 0
                             rpkm_perTE
                           })
##### Dscores ####
scores <-
  lapply(sample_set,
         function(i) {
           message(i)
           name1 <- unlist(strsplit(i, " vs "))[1]
           name2 <- unlist(strsplit(i, " vs "))[2]
           plus_score <-
             apply(abs(plus_rpkm_perTE[[name1]] - plus_rpkm_perTE[[name2]]), 1, sum) / 4
           minus_score <-
             apply(abs(minus_rpkm_perTE[[name1]] - minus_rpkm_perTE[[name2]]), 1, sum) / 4
           data.frame(
             "scores" = plus_score + minus_score,
             "trans" = names(trans),
             "sample_set" = i
           )
         }) %>%
  Reduce(rbind,.)
##### save checkpoint & reload ####
saveRDS(scores, "./results/Dscors_4samples.rds")

Definition of bins

Divide the TEs into 5 bins based on the dscore. Only TEs that satisfy M between -3 and 3 and A greater than 0 in @sec-ma are included in the subsequent analysis.

{r}
bin_names <- c("bin1", "bin2", "bin3", "bin4", "bin5")
names(bin_names) <- bin_names

load("./results/MA_4samples.Rdata")
scores <- readRDS("./results/Dscors_4samples.rds")

##### select trans ####
trans_len <- names(trans[width(trans) >= 1000])

trans_used <- MA[-3<= MA$M & MA$M <= 3 & MA$A >= 0 & MA$trans %in% trans_len,"trans"] %>%
  unique %>%
  as.character

scores_used <- subset(scores, trans %in% trans_used)

##### grouping ####
scores_main <- subset(scores_used, sample_set == "naive13 vs naive20")
basicgroups <-  cut(
  scores_main$scores,
  breaks = c(
    -Inf,
    quantile(scores_main$scores, 0.2),
    quantile(scores_main$scores, 0.4),
    quantile(scores_main$scores, 0.6),
    quantile(scores_main$scores, 0.8),
    Inf
  ),
  labels = bin_names
)
names(basicgroups) <- as.character(scores_main$trans)
basicgroups <- basicgroups[names(basicgroups)!= "TE1_bm_1645_LINE/R4"]
  
##### apply bins to all sets ####

scores_used$groups <- basicgroups[as.character(scores_used$trans)]

##### save files ####
saveRDS(basicgroups, "./results/basicgroups.rds")
saveRDS(scores_used, "./results/dscores.rds")

Dscores vs piRNA

This section utilized both the previously obtained Dscore @sec-dscore and the calculated piRNA RPM @sec-rpm.

{r}
load("./results/rpm_TE_total.Rdata")
dscores <- readRDS("./results/dscores.rds")

sample_set <-
  c(
    "naive13 vs dsRluc13",
    "naive20 vs dsRluc20",
    "naive13 vs naive20",
    "dsRluc13 vs dsRluc20"
  )
names(sample_set) <- sample_set

for (i in sample_set){
  name1 <- unlist(strsplit(i, " vs "))[1]
  name2 <- unlist(strsplit(i, " vs "))[2]
  tmp <- unlist((rpm_total[rpm_total$samples == name1,"rpkm"]+rpm_total[rpm_total$samples == name2,"rpkm"]) /2 )
  names(tmp) <- unlist(rpm_total[rpm_total$samples == name1,"trans"])
  dscores[dscores$sample_set == i,"rpkm"] <- tmp[as.character(dscores[dscores$sample_set == i,"trans"])]
}

saveRDS(dscores, "./results/dscores.rds")

Dscores vs RNA

This section utilized both the previously obtained Dscore @sec-dscore and the calculated piRNA RPM @sec-rpmrna.

{r}
rna_rpm_13 <- readRDS("./results/close13_transposon_RNA_rpm.rds")
rna_rpm_20 <- readRDS("./results/close20_sampling_100nt_transposon_RNA_rpm.rds")
dscores <- readRDS("./results/dscores.rds")

rna_rpkm <-  unlist((rna_rpm_13$rpkm + rna_rpm_20$rpkm ) / 2) 
names(rna_rpkm) <- rna_rpm_13$trans

dscores$rpkmrna <- rna_rpkm[as.character(dscores$trans)]

saveRDS(dscores, "./results/dscores.rds")

Calculation of the Ping-Pong Signature

This section uses the RPM of piRNAs @sec-rpm to calculate the strength of the ping-pong signature.

Overlap

{r}
##### define args ####
sample_names <- c("naive13", "naive20")
names(sample_names) <- sample_names

bin_names <- c("bin1", "bin2", "bin3", "bin4", "bin5")
names(bin_names) <- bin_names

##### load files ####
trans <-  readDNAStringSet("./fasta/transposon.fasta")
basicgroups <- readRDS("./results/basicgroups.rds")

plus_rpm <-
  lapply(sample_names, function(sn) {
    readRDS(sprintf("./results/%s_transposon_piRNA_plus_rpm.rds", sn))
  })
minus_rpm <-
  lapply(sample_names, function(sn) {
    readRDS(sprintf("./results/%s_transposon_piRNA_minus_rpm.rds", sn))
  })

##### count rpkm ####
trans_used <- trans[names(basicgroups)]

plus_rpkm <- lapply(sample_names, function(sn) {
  rpkm <- plus_rpm[[sn]]  / width(trans) * 10 ^ 3
  rpkm[names(trans_used),]
})
minus_rpkm <- lapply(sample_names, function(sn) {
  rpkm <- minus_rpm[[sn]]  / width(trans) * 10 ^ 3
  rpkm[names(trans_used),]
})

##### ping-pong signature ####
maxlen <- max(width(trans_used))

pp13 <- sapply(1:20, function(site) {
  apply((plus_rpkm[["naive13"]][, 1:(maxlen - site + 1)] * minus_rpkm[["naive13"]][, site:maxlen]), 1, sum)}) 
pp20 <- sapply(1:20, function(site) {
  apply((plus_rpkm[["naive20"]][, 1:(maxlen - site + 1)] * minus_rpkm[["naive20"]][, site:maxlen]), 1, sum)}) 

pp <- (pp13 + pp20) %>%
  subset(.,apply(.,1,sum)!=0) %>%
  sweep(., 1, apply(., 1, sum), "/") %>%
  as.data.frame
pp$trans <- rownames(pp)
pp$groups <- basicgroups[as.character(pp$trans)]
ppsign <- pivot_longer(pp,
                       col = 1:20,
                       values_to = "coverage",
                       names_to = "site")

pp_range <-  pp[order(pp$V10, decreasing = F),"trans"]
ppsign$trans <- factor(ppsign$trans, levels = pp_range)
####
save(list=c("pp","ppsign","pp13","pp20"),
     file="./results/pingpongsignature.Rdata")

Z-score

{r}
load("./results/pingpongsignature.Rdata")
pp13 <- pp13 %>%
  subset(.,apply(.,1,sum)!=0) %>%
  sweep(., 1, apply(., 1, sum)+0.01, "/") %>%
  as.data.frame
pp20 <- pp20 %>%
  subset(.,apply(.,1,sum)!=0) %>%
  sweep(., 1, apply(., 1, sum)+0.01, "/") %>%
  as.data.frame

pp13$trans <- rownames(pp13)
pp13$groups <- basicgroups[as.character(pp13$trans)]
pp13$samples <- "naive13"

pp20$trans <- rownames(pp20)
pp20$groups <- basicgroups[as.character(pp20$trans)]
pp20$samples <- "naive20"

ppzscore <- rbind(pp13,pp20)
ppzscore$mean <- apply(ppzscore[,c(1:9,11:20)],1,mean)
ppzscore$sd <- apply(ppzscore[,c(1:9,11:20)],1,sd)
ppzscore$zscore <- (ppzscore$V10 - ppzscore$mean)/ ppzscore$sd

Simulation on the dynamics of piRNA production and competition

Probability of red replacement (P0)

{r}
##### define args ####
Nmax <- 500
Gmax <- 10000 
Rep <- 100000 
Num <- 1

##### probability simulation ####
##### probability simulation ####
P0 <- sapply(1:Nmax,
             function(N) {
               results <- sapply(1:Rep,
                                 function(x) {
                                   Num <- N-1
                                   # the frequency of exsiting-piRNA
                                   Frequency <- Num / N
                                   i <- 0
                                   # Iterate until reaching the maximum number of generations (Gmax) or typeA frequency becomes 0 or 1
                                   
                                   while (i <= Gmax & 
                                          Frequency != 0 &
                                          Frequency != 1) {
                                     Num <- rbinom(1, N, Frequency)
                                     Frequency <- Num / N
                                     i <- i + 1
                                   }
                                   # cat(i)
                                   Num
                                 })
               # Calculate the probability of complete elimination of typeA
               out <- sum(results == 0) / Rep
               message(N)
               out
             }
)
##### save files ####
saveRDS(P0,"./results/simulation_P0.rds")

Probability that red replaces blue at least once (Preplacement)

{r}
##### define args ####
Nmax <- 500
Gmax <- 10000 
Rep <- 100000 
Num <- 1
occurrences <- 100

##### load files ####
P0 <- readRDS("./results/simulation_P0.rds")

##### simulation Preplacement ####
Preplacement <- sapply(1:occurrences,
                           function(n){
                              1-(1-P0)^n
                             })

##### change to long frame ####
Preplacement <- as.data.frame(Preplacement)
colnames(Preplacement) <- 1:occurrences
Preplacement$N <- 1:Nmax
Preplace_long <- pivot_longer(Preplacement,cols = 1:occurrences,names_to="Oemergence",values_to = "P")

##### save files 
saveRDS(Preplace_long,"./results/simulation_Preplace.rds")

Artificial piRNAs

This section focuses on data analysis using smRNA libraries prepared from cells engineered to produce artificial piRNAs. A total of six libraries were used: three nc and three Art.

Therefore, the data in this section also needs to undergo adapter removal, similar to the previous processing of smRNA libraries in @sec-piRNAlibraries. Then, calculate the RPM of the piRNAs of interest.

Removing adapter

{bash}
MM="insert nc"
NN="1 2 3"
for M in ${MM}; do
    for N in ${NN}; do
        echo $N
        cutadapt -a TGGAATTCTCGGGTGCCAAG --minimum-length 20 -o tmp_${M}_${N}_trim.fastq /home/Jie/Saturn/raw/20221228_Chagne_Hiseq/221228_${M}_${N}_1_total.fq.gz
        fastq_to_fasta -Q 33 -i tmp_${M}_${N}_trim.fastq -o tmp_${M}_${N}_trim.fasta
        fastx_collapser < tmp_${M}_${N}_trim.fasta > tmp_${M}_${N}_trim_unique.fasta
        cutadapt -u4 -u -4 --minimum-length 20 -o ./fasta/${M}_${N}_trim.fasta tmp_${M}_${N}_trim_unique.fasta
        gzip ./fasta/${M}_${N}_trim.fasta
        #rm tmp_${M}_${N}_trim.fastq
    done
done

Mapping to TEs

{bash}
bowtie-build -f ./fasta/transposon.fasta transposon
MM="insert nc"
NN="1 2 3"
LL="transposon"
for M in ${MM}; do
    for N in ${NN}; do
        for L in ${LL}; do
            echo "=========="$M$N"----->"$L"=========="   
            ( bowtie --offrate 3 -p 10 -a --best --strata -v 3 -t --sam ${L} -f ./fasta/221228_${M}_${N}_1_total_trim.fasta.gz  > ${M}_${N}_${L}.sam ) 
            samtools view -@ 20  -bS ${M}_${N}_${L}.sam > ${M}_${N}_${L}.bam
            bamToBed -i ${M}_${N}_${L}.bam > ./bed/${M}_${N}_${L}.bed
            rm ${M}_${N}_${L}.sam
            rm ${M}_${N}_${L}.bam
        done
    done
done
rm *ebwt

RPM of piRNAs in Art/NC

{r}
library("Biostrings")
trans <- readDNAStringSet("./fasta/transposon.fasta")

sample_names <- c("insert_1","insert_2","insert_3","nc_1","nc_2","nc_3")
names(sample_names) <- sample_names

rawdata <- lapply(sample_names,
                  function(sn){
                    tmp <- read.table(sprintf("./bed/%s_transposon.bed",sn))
                    tmp$sample_names <- sn
                    tmp
                  })
total_reads <- lapply(rawdata,
                      function(rd){
                        length(unique(rd[rd$V1 != "TE1_bm_1645_LINE/R4","V4"])) 
                      }) 
total_reads2 <- lapply(rawdata,
                      function(rd){
                        tmp <- subset(rd, 
                                      !(V3 == 647 & V6 =="-") &
                                        !(V2 == 637 & V6 == "+") &
                                        !(V1 == "TE1_bm_1645_LINE/R4"))
                        length(unique(tmp$V4))}) 
targets <- lapply(rawdata, 
                  function(rd) {
                    rd$len <- rd$V3 - rd$V2
                    tmp <- subset(rd, V1 == "TE1_bm_1279_Unknown/Unknown")
                    tmp$V2 <- factor(tmp$V2 + 1, levels = 1:1054)
                    tmp$V3 <- factor(tmp$V3, levels = 1:1054)
                    tmp})
piRNA_plus <- lapply(sample_names, 
                     function(sn){
                       tmp_p <- subset(targets[[sn]],26 <= len & len <= 32 & V6 =="+")
                       tmp <- as.data.frame(table(tmp_p[,2]))
                       tmp$direc <- "+"
                       tmp$sample <- sn
                       tmp$rpm <- tmp$Freq / total_reads[[sn]] * 10^6
                       tmp$rpm2 <- tmp$Freq / total_reads2[[sn]] * 10^6
                       tmp
                     })
piRNA_minus <- lapply(sample_names, 
                     function(sn){
                       tmp_m <- subset(targets[[sn]],26 <= len & len <= 32 & V6 =="-")
                       tmp <- as.data.frame(table(tmp_m[,3]))
                       tmp$direc <- "-"
                       tmp$sample <- sn
                       tmp$rpm <- tmp$Freq / total_reads[[sn]] * 10^6
                       tmp$rpm2 <- tmp$Freq / total_reads2[[sn]] * 10^6
                       tmp
                     })
piRNAs <- rbind(Reduce(rbind,piRNA_plus),Reduce(rbind,piRNA_minus))
colnames(piRNAs) <- c("site","reads","direc","sample","rpm","rpm2")
piRNAs$set <- factor(sapply(strsplit(piRNAs$sample, "_"), `[`, 1),
                     levels = c("nc","insert"))
save(list = c("piRNAs"),
     file = "./results/artificial_piRNAs.Rdata")

The positional distribution of mismatches

Standardize the read lengths in the FASTA file to 26

First,

{r}
##### define args ####
suppressMessages(library("Biostrings"))
suppressMessages(library("tidyr"))

sample_names <- c("naive13", "naive20", "dsRluc13", "dsRluc20")
names(sample_names) <- sample_names

##### trim seq to piRNA length ####
for (sn in sample_names) {
  message(sn)
  sn %>%
    sprintf("./fasta/%s_trim.fasta.gz", .) %>%
    readDNAStringSet(.) %>%
    .[26 <= width(.) & width(.) <= 32] %>%
    subseq(., 1, 26) %>%
    writeXStringSet(., sprintf("./fasta/%s_trim_pi26.fasta.gz", sn))
}

{bash}
gzip *_trim_pi26.fasta

Split libraries into 5 bins

Due to the large amount of data to be processed, the BED file and fasta file will be divided into five smaller files based on bins.

First, filter out the piRNAs of interest.

{r}
library("Biostirngs")
trans <- readDNAStringSet("./fasta/transposon.fasta")
##### define args ####
sample_names <- c("naive13", "naive20")
names(sample_names) <- sample_names

##### load files ####
basicgroups <- readRDS("./results/basicgroups.rds")

##### extract piRNA reads from bed file ####
for (sn in sample_names){
  bed <- read.table(sprintf("./bed/%s_transposon.bed", sn))
  bed$len <- bed$V3 - bed$V2
  bed$groups <- basicgroups[as.character(bed$V1)]
  bed$V1 <- factor(bed$V1, levels = names(trans))
  bed$V2 <- factor(bed$V2 + 1 , levels = 1:max(width(trans)))
  bed$V3 <- factor(bed$V3 , levels = 1:max(width(trans)))
  piRNA <- subset(bed, 26 <= len & len <= 32 & !is.na(groups))
  saveRDS(piRNA, sprintf("./piRNA/%s_piRNA_bed.rds", sn))
}

{r}
sample_names <- c("naive13", "naive20")
names(sample_names) <- sample_names

bin_names <- c("bin1","bin2","bin3","bin4","bin5")
names(bin_names) <- bin_names

suppressMessages(library("tidyr"))
suppressMessages(library("Biostrings"))

trans <- readDNAStringSet("./fasta/transposon.fasta")
basicgroups <- readRDS("./results/basicgroups.rds")
basicgroups[]

##### split fasta into each bins ####
for (sn in sample_names) {
  message(sn)
  ## filter the piRNAs based on length and retain them at a uniform length of 26
  fas <- sn %>%
    sprintf("./fasta/%s_trim.fasta.gz", .) %>%
    readDNAStringSet(.) %>%
    .[26 <= width(.) & width(.) <= 32] %>%
    subseq(., 1, 26) 
  writeXStringSet(fas, sprintf("./fasta/%s_trim_pi26.fasta.gz", sn)) 
  
  ## Split the FASTA file into five files based on piRNA information
  piRNA <- readRDS(sprintf("./piRNA/%s_piRNA_bed.rds", sn))
  message("Readed bed file of ", sn)
  
  for (bn in bin_names) {
    piRNA_used <- piRNA %>% 
      subset(.,groups == bn) %>%
      .$V4 %>%
      as.character(.) %>%
      unique(.)
    writeXStringSet(fas[piRNA_used],sprintf("./fasta/%s_pi26_%s.fasta", sn, bn))
    message("Saved fasta files of ", bn," of ",sn)
  }
}

{bash}
gzip *fasta

Consensus sequence

Use ConsensusSequence to filter out the most frequently occurring sequence at each position.

{r consensus_sequence.R}

args <- commandArgs(trailingOnly = TRUE)
# args <- c(13,"bin1","sense")
sn <- args[1]
bn <- args[2]
direc <- args[3]

##### load packages ####
pkgs <- c("Biostrings","tidyr","DECIPHER")
sapply(pkgs, library, character.only = TRUE)

##### define args ####
bin_names <- c("bin1", "bin2", "bin3", "bin4", "bin5")
names(bin_names) <- bin_names

switch(direc,
       sense = {
         b <- "plus"
         c <- "+"
         d <- "V2"
         e <- 2
         f <- "s"
       },
       anti = {
         b <- "minus"
         c <- "-"
         d <- "V3"
         e <- 3
         f <- "a"
       })

##### load files ####
trans <- readDNAStringSet("./fasta/transposon.fasta.gz")

piRNA <- readRDS(sprintf("./piRNA/%s_piRNA_bed.rds", sn)) %>%
  subset(.,V6 == c & groups == bn)
message(sn, " ", bn, " ", direc, " readed bed file")

fas <- readDNAStringSet(sprintf("./fasta/%s_pi26_%s.fasta.gz", sn, bn))
message(sn, " ", bn, " ", direc, " readed fasta file")

##### find read peaks for TE in each bins ####
# Convert 'V4' column to character and 'trans' column to numeric
piRNAnames <- as.character(piRNA$V4)
piRNA$trans <- as.numeric(piRNA$V1)

# Create a table of piRNA reads based on specified columns
peaks <- table(piRNA[, c("trans", d)]) %>%
  as.data.frame(.) %>%
  subset(.,Freq > 3)

##### consensus sequence ####
conseq <- character(nrow(peaks))

timestart <- Sys.time()
message(sn, " ", bn, " ", direc, " start time : ", timestart)

conseq <- apply(peaks, 1, function(pks) {
  tt <- pks[1]
  dd <- pks[2]
  
  # Create a logical mask based on the conditions
  mask <- piRNA$trans == tt & piRNA[, d] == dd
  
  # Get piRNA names based on the mask
  pinames <- piRNAnames[mask]
  
  # Subset 'fas' using the selected piRNA names
  fas.pi <- fas[pinames]
  
  # Obtain the consensus sequence from selected piRNA sequences
  conseq_tmp <- as.character(ConsensusSequence(fas.pi, threshold = 0.5))
  
  # Check if the consensus sequence contains non-ATCG characters
  if (grepl("[^ATCG]", conseq_tmp)) {
    cat("=")
    ast <- table(fas[pinames])
    # message(ast)
  
    # Find the sequence(s) with the maximum count and handle ties
    max_seq <- names(ast[ast == max(ast)])
    ifelse(length(max_seq) > 1, Reduce(function(x, y) {
      paste0(x, strrep("N", 24), y)
    }, max_seq), max_seq)
  } else {
    cat("-")
    conseq_tmp
  }
})

timeend <- Sys.time()
message(sn, " ", bn, " ", direc, "end time : ", timeend)

runningtime <- timeend - timestart
print(runningtime)

##### change the data.frame to Xstring ####
conseq <- DNAStringSet(conseq)
names(conseq) <- paste(peaks[, 1], peaks[, 2], f, sep = "_")

##### save files ####
writeXStringSet(conseq,sprintf("./check/consensus_sequence/%s_%s_%s_conseq.fasta",sn,bn,direc))
message(sn," ",bn," ",direc," saved all results")

{bash 1101run.sh}
YY="naive13 naive20"
BB="bin1 bin2 bin3 bin4 bin5"
DD="sense anti"
echo "Following commands starting..."
for Y in ${YY}; do
  # 启动多个任务
  for B in ${BB}; do
    for D in ${DD}; do
      echo ${Y}${B}${D}
      Rscript ./funs/consensus_sequence_1031.R ${Y} ${B} ${D}  > ./runlog/1101_${Y}${B}${D}.out &
      pids+=($!)
    done
  done
  # 等待当前元素的所有任务完成
  for pid in "${pids[@]}"; do
    wait "$pid"
  done
  unset pids
done

{bash}
nohup bash ./funs/1101run.sh  > ./runlog/1101_runlog.out &

After running the above script, merge the files containing the sense and antisense results.

{r}
##### define args ####
bin_names <- c("bin1", "bin2", "bin3", "bin4", "bin5")
names(bin_names) <- bin_names

sample_names <- c("naive13", "naive20")
names(sample_names) <- sample_names

##### bind ####
for (years in sample_names) {
  for (bn in bin_names) {
    conseq_sense <-
      readDNAStringSet(sprintf(
        "./check/%s_%s_%s_consensus_sequence.fasta",
        years,
        bn,
        "sense"
      ))
    conseq_anti <-
      readDNAStringSet(sprintf(
        "./check/%s_%s_%s_consensus_sequence.fasta",
        years,
        bn,
        "anti"
      ))
    conseq <- c(conseq_sense, conseq_anti)
    writeXStringSet(conseq,
                    sprintf("./fasta/%s_%s_consensus_sequence.fasta", years, bn))
    message(years, " ", bn, " ", " saved in ./fasta")
  }
}

Mapping to consensus sequence

{bash}
MM="naive13 naive20"
LL="bin1 bin2 bin3 bin4 bin5"
NN="bin1 bin2 bin3 bin4 bin5"

for M in ${MM}; do
    for L in ${LL}; do
        bowtie-build -f ./fasta/${M}_${L}_consensus_sequence.fasta ${M}_${L}

        VV="0 1 2 3"

        for V in ${VV}; do 
            echo "bowtie: =========="$M"----->"$L" with v"$V"=========="
            ( bowtie --offrate 3 -p 10 -a --best --strata -v ${V}  -t --sam ${M}_${L} -f ./fasta/${M}_trim_pi26.fasta.gz  > ./sambam/${M}_${L}_v${V}.sam ) 
            echo "samtools: =========="$M"----->"$L" with v"$V"=========="
            samtools view -@ 20  -bS ./sambam/${M}_${L}_v${V}.sam > ./sambam/${M}_${L}_v${V}.bam
            echo "bamToBed: =========="$M"----->"$L" with v"$V"=========="
            bamToBed -i ./sambam/${M}_${L}_v${V}.bam > ./bed/${M}_${L}_v${V}.bed
        done

        for ll in ${NN} ; do 
            echo "bowtie: =========="$ll"----->"$L"=========="
            ( bowtie --offrate 3 -p 10 -a --best --strata -v 0  -t --sam ${M}_${L} -f ./fasta/${M}_${ll}_consensus_sequence.fasta  > ./sambam/${M}_${L}_consensus_sequence_${ll}.sam ) 
            echo "samtools: =========="$ll"----->"$L"=========="
            samtools view -@ 20  -bS ./sambam/${M}_${L}_consensus_sequence_${ll}.sam > ./sambam/${M}_${L}_consensus_sequence_${ll}.bam
            echo "bamTobed: =========="$ll"----->"$L"=========="
            bamToBed -i ./sambam/${M}_${L}_consensus_sequence_${ll}.bam > ./bed/${M}_${L}_consensus_sequence_${ll}.bed
            rm ./sambam/*${ll}.sam
            rm ./sambam/*${ll}.bam
        done
    done
done


code bed_con_all

{r}
##### define args ####
bin_names <- c("bin1", "bin2", "bin3", "bin4", "bin5")
names(bin_names) <- bin_names

sample_names <- c("naive13", "naive20")
names(sample_names) <- sample_names

bed_con_all <- list()

##### con_seq ####
for (sn in sample_names) {
  bed_con <- list()
  for (n in 1:4) {
    for (i in (n + 1):5) {
      bed <-
        read.table(
          sprintf(
            "./bed/%s_%s_consensus_sequence_%s.bed",
            sn,
            bin_names[n],
            bin_names[i]
          )
        )
      bed$len <- bed$V3 - bed$V2
      bed_pi <- subset(bed, 26 == len & V6 == "+")
      message(sn,
              "_",
              bin_names[n],
              " share ",
              nrow(bed_pi),
              " reads with ",
              sn,
              "_",
              bin_names[i])
      bed_con[[bin_names[n]]] <-
        c(unique(as.character(bed_pi$V1)), unique(as.character(bed_pi$V4)))
    }
  }
  bed_con_all[[sn]] <- Reduce(c, bed_con)
}
bed_con_all <- lapply(bed_con_all, unique)

##### 159 ####
trans <- readDNAStringSet("./fasta/transposon.fasta")
basicgroups <- readRDS("./results/basicgroups.rds")
basicgroups_list <- lapply(bin_names,
                           function(bn){
                             namelist <- names(basicgroups[basicgroups == bn])
                             as.numeric(factor(namelist,levels = names(trans)))
                           })

years <- 13
bn <- "bin4"
piRNA_pm <- readRDS(sprintf("./piRNA/naive%s_piRNA_bed.rds", years))

direc <- "sense"
switch(direc,
       sense = {
         b <- "plus"
         g <- "+"
         d <- "V2"
         e <- 2
         f <- "s"
       },
       anti = {
         b <- "minus"
         g <- "-"
         d <- "V3"
         e <- 3
         f <- "a"
       })

piRNA <- subset(piRNA_pm,V6 == g)
piRNA$V4 <- as.character(piRNA$V4)
piRNA$trans <- as.numeric(piRNA$V1)
piRNA.reads <- table(piRNA[, c(1, e)])
piRNA.reads["TE1_bm_1645_LINE/R4",] <- 0
peaks <- which(piRNA.reads > 3, arr.ind = TRUE)
peaks_159_sense <- peaks[peaks[,1] == 1761, ]
sense_159 <- paste(peaks_159_sense[, 1], peaks_159_sense[, 2], f, sep = "_")

direc <- "anti"
switch(direc,
       sense = {
         b <- "plus"
         g <- "+"
         d <- "V2"
         e <- 2
         f <- "s"
       },
       anti = {
         b <- "minus"
         g <- "-"
         d <- "V3"
         e <- 3
         f <- "a"
       })
piRNA <- subset(piRNA_pm,V6 == g)
piRNA$V4 <- as.character(piRNA$V4)
piRNA$trans <- as.numeric(piRNA$V1)
piRNA.reads <- table(piRNA[, c(1, e)])
piRNA.reads["TE1_bm_1645_LINE/R4",] <- 0
peaks <- which(piRNA.reads > 3, arr.ind = TRUE)
peaks_159_anti <- peaks[peaks[,1] == 1761, ]
anti_159 <- paste(peaks_159_anti[, 1], peaks_159_anti[, 2], f, sep = "_")

bed_con_all[["naive13"]] <- c(bed_con_all[["naive13"]],sense_159,anti_159)


years=20
bn <- "bin4"
piRNA_pm <- readRDS(sprintf("./piRNA/naive%s_piRNA_bed.rds", years))

direc <- "sense"
switch(direc,
       sense = {
         b <- "plus"
         g <- "+"
         d <- "V2"
         e <- 2
         f <- "s"
       },
       anti = {
         b <- "minus"
         g <- "-"
         d <- "V3"
         e <- 3
         f <- "a"
       })

piRNA <- subset(piRNA_pm,V6 == g)
piRNA$V4 <- as.character(piRNA$V4)
piRNA$trans <- as.numeric(piRNA$V1)
piRNA.reads <- table(piRNA[, c(1, e)])
piRNA.reads["TE1_bm_1645_LINE/R4",] <- 0
peaks <- which(piRNA.reads > 3, arr.ind = TRUE)
peaks_159_sense <- peaks[peaks[,1] == 1761, ]
sense_159 <- paste(peaks_159_sense[, 1], peaks_159_sense[, 2], f, sep = "_")

direc <- "anti"
switch(direc,
       sense = {
         b <- "plus"
         g <- "+"
         d <- "V2"
         e <- 2
         f <- "s"
       },
       anti = {
         b <- "minus"
         g <- "-"
         d <- "V3"
         e <- 3
         f <- "a"
       })
piRNA <- subset(piRNA_pm,V6 == g)
piRNA$V4 <- as.character(piRNA$V4)
piRNA$trans <- as.numeric(piRNA$V1)
piRNA.reads <- table(piRNA[, c(1, e)])
piRNA.reads["TE1_bm_1645_LINE/R4",] <- 0
peaks <- which(piRNA.reads > 3, arr.ind = TRUE)
peaks_159_anti <- peaks[peaks[,1] == 1761, ]
anti_159 <- paste(peaks_159_anti[, 1], peaks_159_anti[, 2], f, sep = "_")
bed_con_all[["naive20"]] <- c(bed_con_all[["naive20"]],sense_159,anti_159)

<code> v123 total_reads

{r}
##### misratio ####

for (sn in sample_names) {
  for (bn in bin_names) {
    message(sn, bn)
    v0 <- read.table(sprintf("./bed/%s_%s_v0.bed", sn, bn))
    v1 <- read.table(sprintf("./bed/%s_%s_v1.bed", sn, bn))
    v2 <- read.table(sprintf("./bed/%s_%s_v2.bed", sn, bn))
    v3 <- read.table(sprintf("./bed/%s_%s_v3.bed", sn, bn))
    
    v0$len <- v0$V3 - v0$V2
    v1$len <- v1$V3 - v1$V2
    v2$len <- v2$V3 - v2$V2
    v3$len <- v3$V3 - v3$V2
    
    v0_pi <- subset(v0, 26 == len & V6 == "+" & V2 == 0);print(nrow(v0_pi))
    v1_pi <- subset(v1, 26 == len & V6 == "+" & V2 == 0)
    v2_pi <- subset(v2, 26 == len & V6 == "+" & V2 == 0)
    v3_pi <- subset(v3, 26 == len & V6 == "+" & V2 == 0)
    
    v0_piname <- unique(as.character(v0$V4))
    v1_piname <- unique(as.character(v1$V4))
    v2_piname <- unique(as.character(v2$V4))
    
    v0_only <- v0_pi[!as.character(v0_pi$V1) %in% bed_con_all[[sn]], ];print(nrow(v0_only))
    v1_only <- v1_pi[!as.character(v1_pi$V1) %in% bed_con_all[[sn]] & !as.character(v1_pi$V4) %in% v0_piname, ]
    v2_only <- v2_pi[!as.character(v2_pi$V1) %in% bed_con_all[[sn]] & !as.character(v2_pi$V4) %in% v1_piname, ]
    v3_only <- v3_pi[!as.character(v3_pi$V1) %in% bed_con_all[[sn]] & !as.character(v3_pi$V4) %in% v2_piname, ]
    
    total.reads <- nrow(v0_only)+nrow(v1_only)+nrow(v2_only)+nrow(v3_only)
    
    v123 <- c(
      nrow(v0_only) / total.reads,
      nrow(v1_only) / total.reads,
      nrow(v2_only) / total.reads,
      nrow(v3_only) / total.reads
    )
    print(v123)
    
    write.table(v0_only, sprintf("./bed/%s_%s_v0_only.bed", sn, bn))
    write.table(v1_only, sprintf("./bed/%s_%s_v1_only.bed", sn, bn))
    write.table(v2_only, sprintf("./bed/%s_%s_v2_only.bed", sn, bn))
    write.table(v3_only, sprintf("./bed/%s_%s_v3_only.bed", sn, bn))
    
    saveRDS(total.reads,sprintf("./results/%s_%s_v3_total_reads.rds",sn,bn))
    saveRDS(v123, sprintf("./results/%s_%s_v123_ratio.rds", sn, bn))
  }
}

<code> postions_info

{bash}
MM="naive13 naive20"
LL="bin1 bin2 bin3 bin4 bin5"

for M in ${MM}; do
for L in ${LL}; do

q ./sambam/${M}_${L}_v1.bam | awk '/NM:i:0/{split($0, a, "\t");split(a[13], nm, ":") ;print(a[3])}' > ./txt/${M}_${L}_v1_0.txt

echo $M"_"$L" v1"
samtools view -h ./sambam/${M}_${L}_v1.bam | awk '/NM:i:1/{split($0, a, "\t");split(a[13], nm, ":") ;split(nm[3], pos, /[ATCG]/); printf("%s\t%s\n", a[3],pos[1]+1)}' > ./txt/${M}_${L}_v1_1.txt

echo $M"_"$L" v2"
samtools view -h ./sambam/${M}_${L}_v2.bam | awk '/NM:i:2/{split($0, a, "\t");split(a[13], nm, ":") ;split(nm[3], pos, /[ATCG]/); printf("%s\t%s\t%s\n", a[3],pos[1]+1, pos[1]+pos[2]+2)}' > ./txt/${M}_${L}_v2_2.txt

echo $M"_"$L" v3"
samtools view -h ./sambam/${M}_${L}_v3.bam | awk '/NM:i:3/{split($0, a, "\t");split(a[13], nm, ":") ;split(nm[3], pos, /[ATCG]/); printf("%s\t%s\t%s\t%s\n", a[3],pos[1]+1, pos[1]+pos[2]+2, pos[1]+pos[2]+pos[3]+3)}' > ./txt/${M}_${L}_v3_3.txt

done
done

<code> pos_results_long

{r}
bin_names <- c("bin1", "bin2", "bin3", "bin4", "bin5")
names(bin_names) <- bin_names
sample_names <- c("naive13","naive20")
names(sample_names) <- sample_names
misnum <- c(1,2,3)

bed_con_all <- readRDS("./rds/bed_con_all.rds")
basicgroups <- readRDS("./results/basicgroups.rds")

sn ="naive13"
results <- list()
for (bn in bin_names) {
  nn <- paste0(sn, "_", bn)
  results[[bn]] <- matrix(ncol = 3, nrow = 26, 0)
  mis1 <- read.table(sprintf("./txt/%s_%s_v%s_%s.txt", sn, bn, 1, 1))
  mis1[, 2] <- factor(mis1[, 2], levels = 1:26)
  mis1 <- mis1[!as.character(mis1[, 1]) %in% bed_con_all[[sn]], ]
  #divide results for each transposon
  mis1$trans <- mis1[,1] %>% 
    as.character %>%
    strsplit("_") %>% 
    lapply(function(i){i[1]}) %>%
    factor(levels=1:1811)
  misdis_1 <- table(mis1[, c("trans","V2")])
  
  mis2 <- read.table(sprintf("./txt/%s_%s_v%s_%s.txt", sn, bn, 2, 2))
  mis2[, 2] <- factor(mis2[, 2], levels = 1:26)
  mis2[, 3] <- factor(mis2[, 3], levels = 1:26)
  mis2 <- mis2[!as.character(mis2[, 1]) %in% bed_con_all[[sn]], ]
  #divide results for each transposon
  mis2$trans <- mis2[,1] %>% 
    as.character %>%
    strsplit("_") %>% 
    lapply(function(i){i[1]}) %>%
    factor(levels=1:1811)
  misdis_2 <- table(mis2[, c("trans","V2")])+table(mis2[, c("trans","V3")])

  
  mis3 <- read.table(sprintf("./txt/%s_%s_v%s_%s.txt", sn, bn, 3, 3))
  mis3[, 2] <- factor(mis3[, 2], levels = 1:26)
  mis3[, 3] <- factor(mis3[, 3], levels = 1:26)
  mis3[, 4] <- factor(mis3[, 4], levels = 1:26)
  mis3 <- mis3[!as.character(mis3[, 1]) %in% bed_con_all[[sn]], ]
  #divide results for each transposon
  mis3$trans <-
    mis3[, 1] %>%
    as.character %>%
    strsplit("_") %>%
    lapply(function(i) {
      i[1]
    }) %>%
    factor(levels = 1:1811)
  misdis_3 <-
    table(mis3[, c("trans", "V2")]) +
    table(mis3[, c("trans", "V3")]) +
    table(mis3[, c("trans", "V4")])
  
  results[[bn]] <- misdis_1 + misdis_2 + misdis_3
  #results[[bn]] <- sweep(results[[bn]], 1, apply(results[[bn]][, c(2:9, 11:26)], 1, sum), "/") *
   # 100
  results[[bn]] <-
    results[[bn]] %>%
    as.data.frame.matrix %>%
    subset(basicgroups==bn)%>%
    cbind(groups = bn, sample = sn,trans=row.names(.))
}

pos_results_13 <-
  results %>%
  Reduce(rbind,.) %>%
  cbind(
    seed = apply(.[, 2:9], 1, sum),
    #pos10 = re[10,4],
    nonseed1 = apply(.[, 11:18], 1, sum),
    nonseed2 = apply(.[, 19:26], 1, sum)
  )

sn ="naive20"

results <- list()
for (bn in bin_names) {
  nn <- paste0(sn, "_", bn)
  results[[bn]] <- matrix(ncol = 3, nrow = 26, 0)
  mis1 <- read.table(sprintf("./txt/%s_%s_v%s_%s.txt", sn, bn, 1, 1))
  mis1[, 2] <- factor(mis1[, 2], levels = 1:26)
  mis1 <- mis1[!as.character(mis1[, 1]) %in% bed_con_all[[sn]], ]
  #divide results for each transposon
  mis1$trans <- mis1[,1] %>% 
    as.character %>%
    strsplit("_") %>% 
    lapply(function(i){i[1]}) %>%
    factor(levels=1:1811)
  misdis_1 <- table(mis1[, c("trans","V2")])
  
  mis2 <- read.table(sprintf("./txt/%s_%s_v%s_%s.txt", sn, bn, 2, 2))
  mis2[, 2] <- factor(mis2[, 2], levels = 1:26)
  mis2[, 3] <- factor(mis2[, 3], levels = 1:26)
  mis2 <- mis2[!as.character(mis2[, 1]) %in% bed_con_all[[sn]], ]
  #divide results for each transposon
  mis2$trans <- mis2[,1] %>% 
    as.character %>%
    strsplit("_") %>% 
    lapply(function(i){i[1]}) %>%
    factor(levels=1:1811)
  misdis_2 <- table(mis2[, c("trans","V2")])+table(mis2[, c("trans","V3")])

  
  mis3 <- read.table(sprintf("./txt/%s_%s_v%s_%s.txt", sn, bn, 3, 3))
  mis3[, 2] <- factor(mis3[, 2], levels = 1:26)
  mis3[, 3] <- factor(mis3[, 3], levels = 1:26)
  mis3[, 4] <- factor(mis3[, 4], levels = 1:26)
  mis3 <- mis3[!as.character(mis3[, 1]) %in% bed_con_all[[sn]], ]
  #divide results for each transposon
  mis3$trans <-
    mis3[, 1] %>%
    as.character %>%
    strsplit("_") %>%
    lapply(function(i) {
      i[1]
    }) %>%
    factor(levels = 1:1811)
  misdis_3 <-
    table(mis3[, c("trans", "V2")]) +
    table(mis3[, c("trans", "V3")]) +
    table(mis3[, c("trans", "V4")])
  
  results[[bn]] <- misdis_1 + misdis_2 + misdis_3
  #results[[bn]] <- sweep(results[[bn]], 1, apply(results[[bn]][, c(2:9, 11:26)], 1, sum), "/") *
   # 100
  results[[bn]] <-
    results[[bn]] %>%
    as.data.frame.matrix %>%
    subset(basicgroups==bn)%>%
    cbind(groups = bn, sample = sn,trans=row.names(.))
}

pos_results_20 <-
  results %>%
  Reduce(rbind,.) %>%
  cbind(
    seed = apply(.[, 2:9], 1, sum),
    #pos10 = re[10,4],
    nonseed1 = apply(.[, 11:18], 1, sum),
    nonseed2 = apply(.[, 19:26], 1, sum)
  )

pos_results <- rbind(pos_results_13,pos_results_20)
write.table(pos_results,"./txt/pos_results_silkworm.txt",quote=F,row.names=T,col.names=T,sep="\t")
pos_results_long <- pivot_longer(pos_results[,27:32],
                                 col=4:6,
                                 values_to = "percent",
                                 names_to = "regions")
pos_results_long[is.na(pos_results_long$percent),"percent"] =0
saveRDS(pos_results_long, "./results/mismatches_positions_ratio.rds")

Silkworm eggs

This section involves six libraries prepared from silkworm eggs.

Three samples of individual eggs:

231129_K01_single_egg_0h_1.fastq.gz

231129_K02_single_egg_0h_2.fastq.gz

231129_K03_single_egg_0h_3.fastq.gz

A pooled sample of 10 eggs:

231129_K07_10egg_pool_0h_1.fastq.gz

231129_K08_single_egg_0h_2.fastq.gz

231129_K09_single_egg_0h_3.fastq.gz

Sequence analysis of smRNA libraries

Removing adapter

{bash}
M="231129_"
NN="K01_single_egg_0h_1 K02_single_egg_0h_2 K03_single_egg_0h_3 K07_10egg_pool_0h_1 K08_single_egg_0h_2 K09_single_egg_0h_3"
for N in ${NN}; do
    echo $N
cutadapt -a TGGAATTCTCGGGTGCCAAGG --minimum-length 50 --maximum-length 75 -o ${M}${N}_trim.fastq ${M}${N}.fastq.gz

fastq_to_fasta -Q33 -i ${M}${N}_trim.fastq -o tmp_${M}${N}_trim.fasta

fastx_collapser < tmp_${M}${N}_trim.fasta > tmp_${M}${N}_trim_unique.fasta

cutadapt -n 2 -a NNNGTCNNNTAGNNN -g NNNATCNNNAGTNNN -g NNNCGANNNTACNNN --minimum-length 20 -j 20 -o ${M}${N}_trim.fasta tmp_${M}${N}_trim_unique.fasta

rm ${M}${N}_trim.fastq

pigz ${M}${N}_trim.fasta

done

Mapping to TEs

{bash}
bowtie-build -f ./fasta/transposon.fasta transposon

L="transposon"
O="_trim.fasta"
M="231129_"
NN="K01_single_egg_0h_1 K02_single_egg_0h_2 K03_single_egg_0h_3 K07_10egg_pool_0h_1 K08_single_egg_0h_2 K09_single_egg_0h_3"
for N in  ${NN}; do
        echo $N
(bowtie --offrate 3 -p 10 -a --best --strata -v 2 -m 1 -t --sam ${L} -f /home/Jie/Saturn/raw/downloaded/${M}${N}${O} > ${M}${N}_${L}.sam )
samtools view -@ 10 -bS ${M}${N}_${L}.sam> ${M}${N}_${L}.bam
bamToBed -i ${M}${N}_${L}.bam > ./bed/${M}${N}_${L}.bed
rm ${M}${N}_${L}.sam
rm ${M}${N}_${L}.bam
done

Screening bed files for piRNA and saving as rds files

{r}
suppressMessages(library(Biostrings))
transposon <- readDNAStringSet("./fasta/transposon.fasta.gz")

sample_names <- c("01_single_egg_0h_1","02_single_egg_0h_2","03_single_egg_0h_3",
                  "07_10egg_pool_0h_1","08_single_egg_0h_2","09_single_egg_0h_3")

for (sn in sample_names) {
  message(sn)
  
  bed0 <- read.table(sprintf("./bed/231129_K%s_transposon.bed", sn))
  
  bed0$len <- bed0$V3 - bed0$V2
  
  bedp <- subset(bed0, 22 < len & len < 33 & V6 == "+")
  bedm <- subset(bed0, 22 < len & len < 33 & V6 == "-")
  
  tabp <- table(factor(bedp[, 1], levels = names(transposon)),
                factor(bedp[, 2] + 1, levels = 1:max(width(transposon))))
  tabm <- table(factor(bedm[, 1], levels = names(transposon)), 
                factor(bedm[, 3], levels = 1:max(width(transposon))))
  
  tab <- cbind(tabp, tabm)
  
  saveRDS(tab, sprintf("./rds/%s_transposon_5tab.rds", sn))
}

{r}
library(Biostrings)
transposon <- readDNAStringSet("./fasta/transposon.fasta")
mapped <- c(21407594,22110651,17100570,17484049,15664552,19058505,18681487,16839640,18623597,18116636)


a <- read.table("./txt/231129_K01_single_egg_0h_1_transposon_5tab.txt")
b <- read.table("./txt/231129_K02_single_egg_0h_2_transposon_5tab.txt")
c <- read.table("./txt/231129_K03_single_egg_0h_3_transposon_5tab.txt")
d <- read.table("./txt/231129_K04_single_egg_24h_1_transposon_5tab.txt")
e <- read.table("./txt/231129_K05_single_egg_24h_2_transposon_5tab.txt")
f <- read.table("./txt/231129_K06_single_egg_24h_3_transposon_5tab.txt")
g <- read.table("./txt/231129_K07_10egg_pool_0h_1_transposon_5tab.txt")
h <- read.table("./txt/231129_K08_single_egg_0h_2_transposon_5tab.txt")
i <- read.table("./txt/231129_K09_single_egg_0h_3_transposon_5tab.txt")

SUM <- (a/mapped[1]+b/mapped[2]+c/mapped[3]+d/mapped[4]+e/mapped[5]+f/mapped[6]+g/mapped[7]+h/mapped[8]+i/mapped[9])*1000000
SUM0h <- (a/mapped[1]+b/mapped[2]+c/mapped[3])*1000000
SUM24h <- (d/mapped[4]+e/mapped[5]+f/mapped[6])*1000000
SUM10 <- (g/mapped[7]+h/mapped[8]+i/mapped[9])*1000000
reads <- apply(SUM,1,sum)
reads24h <- apply(SUM24h,1,sum)
reads0h <- apply(SUM0h,1,sum)
reads10 <- apply(SUM10,1,sum)
filter_exp <- log2(reads+1)>5
filter_len <- width(transposon) > 999
filter_R4 <- !names(transposon)=="TE1_bm_1645_LINE/R4"
exp <- log2(reads+1) > 10 &log2(reads+1) <15
filters <- filter_exp&filter_len&filter_R4

percentA <- sweep(a,1,apply(a,1,sum),"/")
percentB <- sweep(b,1,apply(b,1,sum),"/")
percentC <- sweep(c,1,apply(c,1,sum),"/")
percentD <- sweep(d,1,apply(d,1,sum),"/")
percentE <- sweep(e,1,apply(e,1,sum),"/")
percentF <- sweep(f,1,apply(f,1,sum),"/")
percentG <- sweep(g,1,apply(g,1,sum),"/")
percentH <- sweep(h,1,apply(h,1,sum),"/")
percentI <- sweep(i,1,apply(i,1,sum),"/")


Ds_AB <- apply(abs(percentA-percentB),1,sum)/2*100
Ds_BC <- apply(abs(percentB-percentC),1,sum)/2*100
Ds_AC <- apply(abs(percentA-percentC),1,sum)/2*100


Ds_DE <- apply(abs(percentD-percentE),1,sum)/2*100
Ds_EF <- apply(abs(percentE-percentF),1,sum)/2*100
Ds_DF <- apply(abs(percentD-percentF),1,sum)/2*100


Ds_GH <- apply(abs(percentG-percentH),1,sum)/2*100
Ds_HI <- apply(abs(percentH-percentI),1,sum)/2*100
Ds_GI <- apply(abs(percentG-percentI),1,sum)/2*100

Ds_0h_single <- (Ds_AB+Ds_BC+Ds_AC)/3
Ds_24h_single <- (Ds_DE+Ds_EF+Ds_DF)/3
Ds_0h_mix <- (Ds_GH+Ds_HI+Ds_GI)/3

low_Ds <- Ds_24h_single+Ds_0h_single < 60
exp_sa <-  log2(reads24h+1)-log2(reads0h+1) 

pdf("tmp.pdf")
hist(Ds_AB[filters],breaks=seq(0, 100, by=5))
hist(Ds_DE[filters],breaks=seq(0, 100, by=5))
hist(Ds_AG[filters],breaks=seq(0, 100, by=5))
hist(Ds_BG[filters],breaks=seq(0, 100, by=5))
hist(Ds_GH[filters],breaks=seq(0, 100, by=5))
hist(Ds_GI[filters],breaks=seq(0, 100, by=5))
plot(Ds_AB[filters],log2(reads+1)[filters])
hist(log2(reads+1))

hist(Ds_0h_single[filters],breaks=seq(0, 100, by=5))
hist(Ds_24h_single[filters],breaks=seq(0, 100, by=5))
hist(Ds_0h_mix[filters],breaks=seq(0, 100, by=5))
plot(Ds_0h_single[filters],log2(reads+1)[filters],xlim=c(0,80))
abline(lm(Ds_0h_single[filters&exp]~log2(reads+1)[filters&exp]),col="red")
plot(Ds_24h_single[filters],log2(reads+1)[filters],xlim=c(0,80))
abline(lm(Ds_24h_single[filters&exp]~log2(reads+1)[filters&exp]),col="red")

plot(Ds_0h_mix[filters],log2(reads+1)[filters],xlim=c(0,80))
plot(Ds_0h_mix[filters],Ds_0h_single[filters],xlim=c(0,80),ylim=c(0,80))
abline(a=0,b=1,col="red")
plot(Ds_0h_single[filters],Ds_24h_single[filters],xlim=c(0,80),ylim=c(0,80))
abline(a=0,b=1,col="red")

plot(Ds_24h_single[filters]-Ds_0h_single[filters])
hist(Ds_24h_single[filters&low_Ds]-Ds_0h_single[filters&low_Ds],breaks=seq(-30, 30, by=1))
abline(v=0,col="green")
plot(log2(reads0h+1)[filters],log2(reads24h+1)[filters])
abline(a=0,b=1,col="green")
plot(Ds_0h_single[filters],Ds_24h_single[filters],xlim=c(0,80),ylim=c(0,80))
par(new=T)
plot(Ds_0h_single[log2(reads24h+1)-log2(reads0h+1)>1 & filters],Ds_24h_single[log2(reads24h+1)-log2(reads0h+1)>1 & filters],col="green",xlim=c(0,80),ylim=c(0,80),xlab="",ylab="")
abline(a=0,b=1,col="red")

plot(exp_sa[filters],Ds_24h_single[filters]-Ds_0h_single[filters])
boxplot(list(exp_sa[filters][Ds_24h_single[filters]-Ds_0h_single[filters]< -0.8038 ],
exp_sa[filters][Ds_24h_single[filters]-Ds_0h_single[filters] > -0.8038 & Ds_24h_single[filters]-Ds_0h_single[filters]< 0.7706 ],
exp_sa[filters][Ds_24h_single[filters]-Ds_0h_single[filters] > 0.7706 & Ds_24h_single[filters]-Ds_0h_single[filters]< 2.6204],
exp_sa[filters][Ds_24h_single[filters]-Ds_0h_single[filters]>2.6204]))


hist(Ds_0h_single,breaks=seq(0, 100, by=5))
hist(Ds_24h_single,breaks=seq(0, 100, by=5))
hist(Ds_0h_mix,breaks=seq(0, 100, by=5))
plot(Ds_0h_mix,Ds_0h_single,xlim=c(0,80),ylim=c(0,80))
abline(a=0,b=1,col="red")
plot(Ds_0h_single,Ds_24h_single,xlim=c(0,80),ylim=c(0,80))
abline(a=0,b=1,col="red")
dev.off()

t.test(Ds_24h_single[!log2(reads24h+1)-log2(reads0h+1)>1&filters]-Ds_0h_single[!log2(reads24h+1)-log2(reads0h+1)>1&filters],
Ds_24h_single[log2(reads24h+1)-log2(reads0h+1)>1 & filters]-Ds_0h_single[log2(reads24h+1)-log2(reads0h+1)>1 & filters])

Drosaphila

smRNA data used in this part is showing as following.

Data

Accession number

Reference

Piwi IP

SRR1568759

[@han2015]

Aub IP

SRR1568760

[@han2015]

Ago3 IP

SRR1568761

[@han2015]

Piwi IP

SRR1746863

[@mohn2015]

Aub IP

SRR1746864

[@mohn2015]

Ago3 IP

SRR1746865

[@mohn2015]

Removing adapter

{bash}
NN="SRR1568759 SRR1568760 SRR1568761 SRR1746863 SRR1746864 SRR1746865"
 for N in  ${NN}; do
        echo $N
cutadapt -a NNNTGGAATTCTCGGGTGCCAAG --minimum-length 20 --maximum-length 30 -o ${N}_trim.fastq ${N}.fastq
fastq_to_fasta -Q33 -i ${N}_trim.fastq >  ${N}_trim.fasta
pigz ${N}_trim.fasta
rm ${N}_trim.fastq
done

Transposon length

{r}
library(Biostrings)
transposon <- readDNAStringSet("dmel-all-transposon-r6.27.fasta")
transposon2 <- transposon[width(transposon)> 999 & width(transposon) < 10001]
writeXStringSet(transposon2,"dmel-all-transposon-r6.27_2to3.fasta")


Mapping to TEs

{bash}
bowtie-build /home/Jie/Saturn/raw/dmel-all-transposon-r6.27_2to3.fasta Dro_transposon_2to3
L="Dro_transposon_2to3"
O="_trim.fasta"
NN="SRR1568759 SRR1568760 SRR1568761 SRR1746863 SRR1746864 SRR1746865"
for N in  ${NN}; do
        echo $N
(bowtie --offrate 3 -p 10 -a --best --strata -v 3 -t --sam ${L} -f /home/Jie/Saturn/raw/downloaded/${N}${O} > ${N}_${L}.sam )
samtools view -@ 10 -bS ${N}_${L}.sam> ${N}_${L}.bam
bamToBed -i ${N}_${L}.bam > ./bed/${N}_${L}_v3.bed
rm ${N}_${L}.sam
rm ${N}_${L}.bam
done

Defined 5ʹ end of piRNAs and calculation of Dscores

{bash}
NN="SRR1568759 SRR1568760 SRR1568761"
MM="SRR1746863 SRR1746864 SRR1746865"
i=1
for 1 in   do
for M in  ${MM}
	echo $N
	Rscript --slave --vanilla ./funs/Dro5endtab_v3.R ${N} ${M} &
done

{r}
#### Dro5endtab.R
#### Defined 5ʹ end of piRNAs
argv <- commandArgs(TRUE)

suppressMessages(library(Biostrings))
bed0 <- read.table(paste0("./bed/",argv[1],"_Dro_transposon_2to3_v3.bed"))
transposon <- readDNAStringSet("/home/Jie/Saturn/raw/dmel-all-transposon-r6.27_2to3.fasta")
tname <- sub("\\s.*", "", names(transposon))

bed <- bed0[bed0[, 3] - bed0[, 2] > 22 & bed0[, 3] - bed0[, 2] < 33, ]
bedp <- bed[bed[, 6] == "+", ]
bedm <- bed[bed[, 6] == "-", ]

tabp <- table(factor(bedp[, 1], levels = tname), factor(bedp[, 2] + 1, levels =
                                                          1:max(width(transposon))))
tabm <- table(factor(bedm[, 1], levels = tname), factor(bedm[, 3], levels =
                                                          1:max(width(transposon))))

mapped1 <- length(unique(bed$V4))
tab1 <- cbind(tabp,tabm)
tabsum1 <- apply(tab1,1,sum)/mapped1 *10^6

saveRDS(tab1,paste0("./rds/",argv[1],"_transposon_5tab.rds"))

bed0 <- read.table(paste0("./bed/",argv[2],"_Dro_transposon_2to3_v3.bed"))
transposon <- readDNAStringSet("/home/Jie/Saturn/raw/dmel-all-transposon-r6.27_2to3.fasta")
tname <- sub("\\s.*", "", names(transposon))

bed <- bed0[bed0[,3]-bed0[,2] > 22 & bed0[,3] - bed0[,2] < 33,]
bedp <- bed[bed[,6]=="+",]
bedm <- bed[bed[,6]=="-",]

tabp <- table(factor(bedp[,1],levels=tname),factor(bedp[,2]+1,levels=1:max(width(transposon))))
tabm <- table(factor(bedm[,1],levels=tname),factor(bedm[,3],levels=1:max(width(transposon))))

mapped2 <- length(unique(bed$V4))
tab2 <- cbind(tabp,tabm)
tabsum2 <- apply(tab2,1,sum)/mapped2 *10^6

saveRDS(tab2,paste0("./rds/",argv[2],"_transposon_5tab.rds"))
#### Calculation of Dscores

tabfil <- log2(tabsum1+tabsum2) > 1
percent1 <- sweep(tab1,1,apply(tab1,1,sum),"/")
percent2 <- sweep(tab2,1,apply(tab2,1,sum),"/")

Ds <- apply(abs(percent1-percent2),1,sum)/2*100
saveRDS(Ds,paste0("./rds/",argv[1],"VS",argv[2],"_transposon_dscores.rds"))

Consensus sequence

{bash}
NN="SRR1568759 SRR1568760 SRR1568761 SRR1746863 SRR1746864 SRR1746865"
i=1
for N in  ${NN}; do
	echo $N
	Rscript --slave --vanilla ./funs/vari_piRNA_v3.R ${N}
	i=$[$i+1]
	echo $i
done

{r}
##### Load Required Libraries #####
argv <- commandArgs(TRUE)
suppressMessages(library(pbapply))
suppressMessages(library(Biostrings))

##### Read Input Data #####
bed0 <- read.table(paste0("./bed/", argv[1], "_Dro_transposon_2to3_v3.bed"))
transposon <- readDNAStringSet("/home/Jie/Saturn/raw/dmel-all-transposon-r6.27_2to3.fasta")
fas <- readDNAStringSet(paste0("/home/Jie/Saturn/raw/downloaded/", argv[1], "_trim.fasta"))
tname <- sub("\\s.*", "", names(transposon))
fasnames <- sub("\\s.*", "", names(fas))

##### Filter Bed Data #####
bed <- bed0[bed0[, 3] - bed0[, 2] > 22 & bed0[, 3] - bed0[, 2] < 33, ]
bedp <- bed[bed[, 6] == "+", ]
bedm <- bed[bed[, 6] == "-", ]

mapped <- length(unique(bed$V4))

##### Compute RPM #####
bedpID <- paste(bedp[, 1], bedp[, 2], sep = "_")
bedmID <- paste(bedm[, 1], bedm[, 3], sep = "_")

RPMp <- table(bedpID) / mapped * 10^6
RPMm <- table(bedmID) / mapped * 10^6
tagp <- RPMp[RPMp > 3]
tagm <- RPMm[RPMm > 3]

bedp_sel <- bedp[is.element(bedpID, names(tagp)), ]
bedm_sel <- bedm[is.element(bedmID, names(tagm)), ]

##### Select Relevant Sequences #####
sel_names <- c(as.character(bedp_sel[, 4]), as.character(bedm_sel[, 4]))
fas_sel <- DNAStringSet(substr(fas[is.element(fasnames, sel_names)], 1, 22))
fasnames_sel <- fasnames[is.element(fasnames, sel_names)]
bedpID_sel <- bedpID[is.element(bedpID, names(tagp))]
bedmID_sel <- bedmID[is.element(bedmID, names(tagm))]

##### Initialize Result Tables #####
val_tabp <- matrix(0, nrow = length(tagp), ncol = 22)
val_tabm <- matrix(0, nrow = length(tagm), ncol = 22)
rownames(val_tabp) <- names(tagp)
rownames(val_tabm) <- names(tagm)

##### Calculate Variation Tables #####
pb <- txtProgressBar(min = 0, max = length(tagp), style = 3)
for (i in seq_len(length(tagp))) {
  fas_tmp <- fas_sel[is.element(fasnames_sel, bedp_sel[is.element(bedpID_sel, names(tagp)[i]), 4])]
  val_tabp[i, ] <- (length(fas_tmp) - apply(consensusMatrix(fas_tmp), 2, max)) / length(fas_tmp)
  setTxtProgressBar(pb, i)
}

pb <- txtProgressBar(min = 0, max = length(tagm), style = 3)
for (i in seq_len(length(tagm))) {
  fas_tmp <- fas_sel[is.element(fasnames_sel, bedm_sel[is.element(bedmID_sel, names(tagm)[i]), 4])]
  val_tabm[i, ] <- (length(fas_tmp) - apply(consensusMatrix(fas_tmp), 2, max)) / length(fas_tmp)
  setTxtProgressBar(pb, i)
}

##### Save Results #####
saveRDS(val_tabp, paste0("./rds/", argv[1], "_variated_piRNAs_posi.rds"))
saveRDS(val_tabm, paste0("./rds/", argv[1], "_variated_piRNAs_nega.rds"))

Mouse

smRNA data used in this part is showing as following.

Data

Accession number

Reference

Mili IP

SRR5304361

[@wenda2017]

Miwi2 IP

SRR5304364

[@wenda2017]

Removing adapter

{bash}
NN="SRR5304361 SRR5304364"
for N in  ${NN}; do
        echo $N
cutadapt -a AGATCGGAAGAG --minimum-length 20 --maximum-length 36 -o ${N}_trim.fastq ${N}.fastq
fastq_to_fasta -Q33 -i ${N}_trim.fastq >  ${N}_trim.fasta
pigz ${N}_trim.fasta
rm ${N}_trim.fastq
done

Mapping to TEs

{bash}
bowtie-build /home/Jie/Saturn/raw/mus_musculus_dfam.fasta mus_musculus_dfam
LL="mus_musculus_dfam"
NN="SRR5304361 SRR5304364"
O="_trim.fasta"
for L in ${LL}; do
for N in  ${NN}; do
echo $N
( bowtie --offrate 3 -p 20 -a --best --strata -v 3 -t --sam ${L} -f /home/Jie/Saturn/raw/downloaded/${N}${O}  > ${N}_${L}.sam )  
samtools view -@ 20 -bS ${N}_${L}.sam> ${N}_${L}.bam
bamToBed -i ${N}_${L}.bam > ./bed/${N}_${L}_v3.bed
rm ${N}_${L}.sam
rm ${N}_${L}.bam
awk '$3-$2>26 && $3-$2<40' ./bed/${N}_${L}_v3.bed > ./bed/${N}_${L}_piRNA_v3.bed
done
done

