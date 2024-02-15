---
title: "Autonomous Shaping of piRNA Repertoire by Adjacent Site Competition"
author: "Jie Yu, Natsuko Izumi, Yukihide Tomari, Keisuke Shoji"
date: "`r Sys.Date()`"
output:
  html_document:
  toc: true
---

# Global Configs on local computer

This section is the basic setup required to run this code locally. As mentioned in the code, it involves four packages:

ggplot2, Biostrings, tidyr, and DECIPHER.

```{r setup, include=FALSE}
# global settings
knitr::opts_chunk$set(show.error.messages = FALSE,
                      warning = -1, 
                      include = FALSE,
                      message = FALSE,
                      error = TRUE,
                      tidy = TRUE,
                      eval=FALSE)

#### Setting the working directory ####

setwd("~/Documents/Doctor/paper/v3")

#### load packages and files ####
pkgs <- c("ggplot2","Biostrings","tidyr")
sapply(pkgs,library,character.only=TRUE)

trans <- readDNAStringSet("./fasta/transposon.fasta")

load("data.RData")
#### figure size ####
Width = 3
Height = 2

#### figure theme ####
mytheme <- function() {
  theme_classic() +
    theme(
      text = element_text(size = 12,color = "black"),
      title = element_text(size = 10,color = "black"),
      axis.title = element_text(size = 12,color = "black"),
      axis.text = element_text(size = 12,color = "black"),
      legend.text = element_text(size = 12,color = "black"),
      panel.border = element_rect(fill = NA,color = "black",size = unit(1,"pt")),
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background = element_rect(fill = "transparent", colour = NA),
      legend.background = element_rect (fill='transparent')
    )
}
```

# Gloabl Configs on termial computer

In all titles, the code in italic font is a code snippet to be run on the terminal. The saved results are downloaded to the corresponding local folder using the scp command, for example, from ./results on the terminal to ./results locally.

```{r}
setwd("/home/Jie/Jupiter/paper/v3")
pkgs <- c("Biostrings","tidyr","DECIPHER")
sapply(pkgs,library,character.only=TRUE)
trans <-  readDNAStringSet("./fasta/transposon.fasta")
```

# Common Scripts

## smRNA seq

```{bash eval=FALSE, include=FALSE}
#### remove adapater ####
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

#### map to transposons ####
bowtie-build -f ./fasta/transposon.fasta transposon
MM="naive dsRluc"
NN="13 20"
LL="transposon"
for M in ${MM}; do
    for N in ${NN}; do
        for L in ${LL}; do
            echo "=========="$M$N"----->"$L"=========="   
            ( bowtie --offrate 3 -p 10 -a --best --strata -v 0 -t --sam ${L} -f ./fasta/${M}${N}_trim.fasta.gz  > ${M}${N}_${L}.sam ) 
            samtools view -@ 20  -bS ${M}${N}_${L}.sam > ${M}${N}_${L}.bam
            bamToBed -i ${M}${N}_${L}.bam > ./bed/${M}${N}_${L}.bed
            rm ${M}${N}_${L}.sam
            rm ${M}${N}_${L}.bam
        done
    done
done
rm *ebwt

#### map to fem & masc ####
#fem_mRNA: TTTCATTGTTACCTCTTTTTGTCAATTCATAAAGTCATTCAGTG
#Masc_mRNA:AAATGGCTTTGTGAATCGACAAAAAGAGGTAACAATTGAAGCTAATCAGAAGAAA
bowtie-build ./fasta/femmasc.fas femmasc
MM="naive dsRluc"
NN="13 20"
LL="femmasc"
for M in ${MM}; do
    for N in ${NN}; do
        for L in ${LL}; do
            echo "=========="$M$N"----->"$L"=========="
            ( bowtie --offrate 3 -p 10 -a --best --strata -v 0 -t --sam ${L} -f ./fasta/${M}${N}_trim.fasta.gz  > ${M}${N}_${L}.sam ) ##07.26 revised to v3m1
            samtools view -@ 20  -bS ${M}${N}_${L}.sam > ${M}${N}_${L}.bam
            bamToBed -i ${M}${N}_${L}.bam > ./bed/${M}${N}_${L}.bed
            rm ${M}${N}_${L}.sam
            rm ${M}${N}_${L}.bam
        done
    done
done
rm *ebwt

```

## RNA seq

```{bash}
#### mapping to transposon ####
hisat2-build -f /home/Jie/Jupiter/paper/fasta/transposon.fasta transposon
MM="close"
NN="13 20"
LL="transposon"
for M in ${MM}; do
    for N in ${NN}; do
        for L in ${LL}; do
            echo ${M}${N}
            hisat2 -x ${L} -1 ./fastq/${M}${N}_R1.fastq.gz  -2 ./fastq/${M}${N}_R2.fastq.gz -k 3 -p 10 -S ./sambam/${M}${N}_${L}_RNA_HS2.sam
            samtools view -@ 10 -bS ./sambam/${M}${N}_${L}_RNA_HS2.sam> ./sambam/${M}${N}_${L}_RNA_HS2.bam
            bamToBed -i ./sambam/${M}${N}_${L}_RNA_HS2.bam > ./bed/${M}${N}_${L}_RNA_HS2.bed
            grep '\+' <./bed/${M}${N}_${L}_RNA_HS2.bed | grep '/1' >tmp1
            grep '\-' <./bed/${M}${N}_${L}_RNA_HS2.bed | grep '/1' >tmp2
            grep '\+' <./bed/${M}${N}_${L}_RNA_HS2.bed | grep '/2' >tmp3
            grep '\-' <./bed/${M}${N}_${L}_RNA_HS2.bed | grep '/2' >tmp4
            cat tmp2 tmp4 > ./bed/${M}${N}_anti_${L}_RNA_HS2.bed
            cat tmp1 tmp3 > ./bed/${M}${N}_sense_${L}_RNA_HS2.bed
        done
    done
done
rm *.ht2
rm tmp*
```

## BEDtoRDS.R

change bed format to rds

```{r}
#### load files and args ####
args <- commandArgs(trailingOnly = TRUE)
bed_file <- args[1]

bed <- read.table(sprintf("./bed/%s.bed",bed_file))

#### save files in new format ####
if (!dir.exists("rds"))
{dir.create("rds")}
saveRDS(bed,sprintf("./rds/%s.rds",bed_file))
```

## piRNA.R

pre-processing piRNA bed files

```{r}
#### load args & packages & files ####
args <- commandArgs(trailingOnly = TRUE)
raw_file <- as.character(args[1])

suppressMessages(library("Biostrings"))
trans <- readDNAStringSet("./fasta/transposon.fasta")

#### transform bed to rds ####
message("Transforming ", raw_file, " to *.rds file")
rawdata <- read.table(sprintf("./bed/%s.bed", raw_file))
if (!dir.exists("rds"))
{
  dir.create("rds")
}
saveRDS(rawdata, sprintf("./rds/%s.rds", raw_file))
message(raw_file, ".rds file has saved")

#### length table ####
message("now handling ", raw_file)
rawdata$len <- rawdata$V3 - rawdata$V2
rawdata$V1 <- factor(rawdata$V1, levels = names(trans))
rawdata$V2 <- factor(rawdata$V2 + 1 , levels = 1:max(width(trans)))
len <- table(factor(rawdata$len, levels = 20:43))

#### save length files ####
if (!dir.exists("len"))
{
  dir.create("len")
}
saveRDS(len, sprintf("./len/%s_len.rds", raw_file))
message("saved ", raw_file, " length distribution")

#### defined piRNA ####
piRNA_plus <- subset(rawdata, 26 <= len & len <= 32 & V6 == "+")
piRNA_minus <- subset(rawdata, 26 <= len & len <= 32 & V6 == "-")

#### save piRNA files ####
if (!dir.exists("piRNA"))
{
  dir.create("piRNA")
}
saveRDS(piRNA_plus, sprintf("./piRNA/%s_piRNAplus.rds", raw_file))
saveRDS(piRNA_minus, sprintf("./piRNA/%s_piRNAminus.rds", raw_file))
message("saved ", raw_file, " piRNA")

#### count numbers of reads ####
sum_plus <- length(unique(piRNA_plus$V4))
sum_minus <- length(unique(piRNA_minus$V4))
sumreads <- sum_plus + sum_minus
message("Total reads in ", raw_file, "is ", sumreads)

#### table of reads distribution ####
message("Counting ", raw_file, " reads & rpm")
plus_reads <- table(piRNA_plus[, c(1, 2)])
piRNA_minus$V3 <- factor(piRNA_minus$V3, levels = 1:max(width(trans)))
minus_reads <- table(piRNA_minus[, c(1, 3)])

#### count rpm ####
plus_rpm <- plus_reads / sumreads * 10 ^ 6
minus_rpm <- minus_reads / sumreads * 10 ^ 6

#### save files ####
save(list = c("sumreads", "plus_reads", "minus_reads"),
     file = sprintf("./piRNA/%s_piRNA_reads.RData",raw_file))
saveRDS(plus_rpm, sprintf("./piRNA/%s_piRNA_plus_rpm.rds",raw_file))
saveRDS(minus_rpm, sprintf("./piRNA/%s_piRNA_minus_rpm.rds",raw_file))
message("saved ",raw_file," reads & rpm")
```

# Present in the Manuscript

## MAplot

### *code* MA

```{r eval=FALSE, message=TRUE, warning=TRUE, include=FALSE}
#### define args ####
sample_names <- c("naive13","naive20","dsRluc13","dsRluc20")
names(sample_names) <- sample_names

sample_set <- c(
  "naive13 vs dsRluc13",
  "naive20 vs dsRluc20",
  "naive13 vs naive20",
  "dsRluc13 vs dsRluc20"
)
names(sample_set) <- sample_set

#### load files ####
plus_rpm <-
  lapply(sample_names, function(i) {
    readRDS(sprintf("./piRNA/%s_transposon_piRNA_plus_rpm.rds", i))
  })
minus_rpm <-
  lapply(sample_names, function(i) {
    readRDS(sprintf("./piRNA/%s_transposon_piRNA_minus_rpm.rds", i))
  })

rpkm <- sapply(sample_names,function(i){
  rpm_tmp <- apply(plus_rpm[[i]],1,sum)+apply(minus_rpm[[i]],1,sum)
  (rpm_tmp / width(trans))*10^3
})

#### MA ####
# Apply a function to each element in the sample_set
MA <- lapply(sample_set,
             function(i) {
               
               # Split the string using " vs " as the delimiter and get the first and second parts
               aa <- unlist(strsplit(i, " vs "))[1]
               bb <- unlist(strsplit(i, " vs "))[2]
               
               # Create a data frame with the required columns
               df <- data.frame(
                 "trans" = factor(rownames(rpkm), levels = names(trans)),
                 "M" = log2(rpkm[, aa] + 0.01) - log2(rpkm[, bb] + 0.01),
                 "A" = log2(rpkm[, aa] + 0.01) / 2 + log2(rpkm[, bb] + 0.01) / 2,
                 "samples" = i
               )
               
               # Add a column 'range' based on conditions
               df$range <- ifelse(df$A >= 0 & -3 <= df$M & df$M <= 3,
                                  "center",
                                  ifelse(
                                    df$A < 0 & 3 < df$M ,
                                    "low_up",
                                    ifelse(
                                      df$A < 0 & -3 <= df$M & df$M <= 3,
                                      "low_center",
                                      ifelse(
                                        df$A < 0 & df$M < -3,
                                        "low_down",
                                        ifelse(df$A >= 0 &
                                                 df$M < -3, "down", "up")
                                      )
                                    )
                                  ))
               
               # Return the data frame
               df
             })

#### define color regions ####
MA1320 <- MA[["naive13 vs naive20"]]
targets <- as.character(MA1320[MA1320$range == "center", "trans"])

#### save files ####
saveRDS(rpkm, "./check/piRNA_rpkm.rds")
saveRDS(MA, "./results/MA_4samples.rds")
saveRDS(targets, "./results/MA_targets.rds")
```

### **code** dat.MA

```{R}
#### load files ####
MA <- readRDS("./results/MA_4samples.rds")
basicgroups <- readRDS("./results/basicgroups.rds")

#### defined bins in MA ####
MA_bin <- lapply(MA, function(i) {
  i$groups <- basicgroups[as.character(i$trans)]
  i$globlegroups <-
    ifelse(i$range == "center", as.character(i$groups), i$range)
  i$samples <-
    factor(
      i$samples,
      levels = c(
        "naive13 vs dsRluc13",
        "naive20 vs dsRluc20",
        "naive13 vs naive20",
        "dsRluc13 vs dsRluc20"
      ),
      labels = c(
        "Naive13 vs. dsRluc13",
        "Naive20 vs. dsRluc20",
        "Naive13 vs. Naive20",
        "dsRluc13 vs. dsRluc20"
      )
    )
  i
})

#### globelgroups ####
globlegroups <- MA_bin[["naive13 vs naive20"]]$globlegroups
names(globlegroups) <- MA_bin[["naive13 vs naive20"]]$trans

#### add to fig rdata ####
dat.MA <- MA_bin
save(list = ls(pattern = "dat"), file = "data.RData")
```

### **fig** MAplot

```{r}
lapply(dat.MA, function(i) {
  fig <- ggplot(i, aes(
    x = A,
    y = M,
    col = range,
    shape = range
  )) +
    geom_point(size = 1, shape = 21) +
    ggtitle(unique(i[, "samples"])) +
    xlab(expression(Log[2] ~ mean ~ expression)) +
    ylab(expression(Log[2] ~ fold ~ change)) +
    geom_hline(
      yintercept = c(-3, 3),
      linetype = "dashed",
      col = "black"
    ) +
    # geom_hline(yintercept = 0)+
    geom_vline(xintercept = -0,
               linetype = "dashed",
               col = "black") +
    scale_color_manual(values = c("black", rep("darkgrey", 6)),
                       guide = "none") +
    # scale_shape_manual(values = c(19, rep(21,6)),
    #                    guide = "none") +
    # facet_wrap(groups~.)+
    mytheme()
  print(fig)
  ggsave(
    filename = sprintf("./fig/MAplot_%s.pdf", unique(i[, "samples"])),
    width = Width,
    height = Height
  )
})
```

**fig**

```{r}
 ggplot(subset(dat.MA[["naive13 vs naive20"]],!is.na(groups)), aes(x = groups, y = M,fill = groups)) +
  geom_violin()+
    ggtitle(unique(dat.MA[["naive13 vs naive20"]][, "samples"])) +
    xlab("log2 mean expression") +
    ylab("Log2 fold change") +
    geom_hline(yintercept = c(-3, 3),
               linetype = "dashed",
               col = "#955251")+
  geom_hline(yintercept = 0)+
    geom_vline(xintercept = -0,
               linetype = "dashed",
               col = "#955251") +
    scale_fill_manual(values = c("blue", "#800080","red","#FF7A00","#FFD700")) +
    # facet_wrap(groups~.)+
    mytheme()
```

## TE RPM

### **code** example TEs

```{r}
#### define args ####
trans.names <-
  c("TE1_bm_228_LTR/Unknown", "TE1_bm_159_Unknown/Unknown")
samples.names <- c("naive13","naive20")

#### load files ####
rpkm <- readRDS( "./check/piRNA_rpkm.rds")

#### extract the two trans ####
rpkm_used_trans <- as.data.frame(rpkm[trans.names,samples.names])
rpkm_used_trans$trans <- row.names(rpkm_used_trans)

#### change to long frame ####
rpkm_used_trans_long <-
  pivot_longer(
    rpkm_used_trans,
    col = 1:2,
    names_to = "year",
    values_to = "RPM"
  )

#### add to rdata ####
dat.rpm <- rpkm_used_trans_long
save(list = ls(pattern = "dat"), file = "data.RData")
```

### **fig** rpm288159

```{r}
ggplot(dat.rpm) +
  geom_col(aes(x = trans, y = RPM, fill = year),
           position = "dodge",
           width = 0.6) +
  ggtitle("piRNA \nexpression levels") +
  xlab("") +
  ylab("RPKM") +
  scale_fill_manual(values = c("#22406a", "#962832"),
                    label = c("Naive13", "Naive20")) +
  scale_x_discrete(labels = c(),breaks = c()) +
  # scale_x_discrete(labels = c("Naive13", "Naive20")) +
  scale_y_log10(
    # limits = c(1,100000),
    breaks = c(1, 10, 100, 1000, 10000),
    labels = c(1, 10, expression(10 ^ 2), expression(10 ^ 3), expression(10 ^ 4))) +
  # facet_grid(.~trans)+
  mytheme() +
  # guides(fill=guide_legend(title=element_blank(),
  #                         nrow=1))+
  guides(fill="none")+
  theme(
    # legend.position = c (.5, .85),
    # legend.key.size = unit(0.5, "cm"),
    # legend.spacing.y = unit(0.001, "cm"),
    # legend.margin = margin(t = -18),
    # axis.text.x = element_text(angle=45,hjust=1),
    axis.text.y = element_text(hjust = 1),
    title = element_text(size = 10)
  )
ggsave(filename = "./fig/rpm_228159.pdf",dpi=300, width = 2, height = 1.5)

```

## TE distribution

### **code** dat.distri

```{r}
#### define args ####
sample_names <- c("naive13","naive20")
names(sample_names) <- sample_names

trans.names <-
  c("TE1_bm_228_LTR/Unknown", "TE1_bm_159_Unknown/Unknown")
names(trans.names) <- trans.names

#### load files ####
plus_rpm <-
  lapply(sample_names, function(i) {
    readRDS(sprintf("./piRNA/%s_transposon_piRNA_plus_rpm.rds", i))
  })
minus_rpm <-
  lapply(sample_names, function(i) {
    readRDS(sprintf("./piRNA/%s_transposon_piRNA_minus_rpm.rds", i))
  })

#### extract specific trans record ####
distribution <- lapply(trans.names,function(tn){
   trans_rpm <- lapply(sample_names,function(sn){
    trans_plus <- data.frame("rpm"=plus_rpm[[sn]][tn,],
                             "site" = 1:length(plus_rpm[[sn]][tn,]),
                             "direc" = "sense")
    trans_minus <- data.frame("rpm"=minus_rpm[[sn]][tn,]*(-1),
                              "site" = 1:length(minus_rpm[[sn]][tn,]),
                             "direc" = "anti") 
   tmp_rpm <- rbind(trans_plus,trans_minus)
    tmp_reads <- sum(trans_plus$rpm)+sum(abs(trans_minus$rpm))
    tmp_rpm$rpm <- tmp_rpm$rpm/tmp_reads *100
    tmp_rpm$trans <- tn
    tmp_rpm$year <- sn
    tmp_rpm
  })
   Reduce(rbind,trans_rpm)
})
distribution_long <-  Reduce(rbind,distribution)

#### add to rdata ####
dat.distri <- distribution_long
save(list = ls(pattern = "dat"), file = "data.RData")
```

### **fig** TE228

```{r}
ggplot(subset(dat.distri,trans=="TE1_bm_228_LTR/Unknown"))+
  geom_col(aes(x=site,y=rpm,fill=direc),col="black",size = 0.35)+
  # ylim(-2,6)+
  ggtitle("TE228")+
  ylab("Relative abundance (%)")+
  xlab("")+
  scale_x_continuous(limits = c(7200,7300),
                     labels = scales::comma)+
  # ylim(-20,80)+
  scale_fill_manual(values = c("white","black"))+
  guides(fill="none")+
  facet_grid(rows=vars(year),
             labeller = labeller(year=c("naive13"="Naive13","naive20"="Naive20"),
                                 TE=c("TE1_bm_228_LTR/Unknown"="TE228")))+
  mytheme()
ggsave(filename = "./fig/distribution_te228.pdf",width=Width,height=Height*1.8)
```

### **fig** TE159

```{r}
ggplot(subset(dat.distri,trans=="TE1_bm_159_Unknown/Unknown"))+
  geom_col(aes(x=as.numeric(site),y=rpm,fill=direc),col="black",size=0.35)+
  # ylim(-2,6)+
  ggtitle("TE159")+
  ylab("Relative abundance (%)")+
  xlab("")+
  scale_x_continuous(limits = c(1200,1300),
                     labels = scales::comma)+
  # ylim(-20,80)+
  scale_fill_manual(values = c("white","black"))+
  guides(fill="none")+
  facet_grid(rows=vars(year),
             labeller = labeller(year=c("naive13"="Naive13","naive20"="Naive20"),
                                 TE=c("TE1_bm_159_Unknown/Unknown"="TE159")))+
  mytheme()
ggsave(filename = "./fig/distribution_te159.pdf",width=Width,height=Height*1.8)
```

## Dscores

### *code* dscores

```{r eval=FALSE, include=FALSE}
#### define args ####
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

bin_names <- c("bin1", "bin2", "bin3", "bin4", "bin5")
names(bin_names) <- bin_names

#### load files ####
plus_rpm <-
  lapply(sample_names, function(i) {
    readRDS(sprintf("./piRNA/%s_transposon_piRNA_plus_rpm.rds", i))
  })
minus_rpm <-
  lapply(sample_names, function(i) {
    readRDS(sprintf("./piRNA/%s_transposon_piRNA_minus_rpm.rds", i))
  })
targets <- readRDS("./results/MA_targets.rds")

#### rpkm per TE ####
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
#### Dscores ####
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
         })

#### save checkpoint & reload ####
saveRDS(scores, "./check/Dscors_4samples.rds")
scores <- readRDS("./check/Dscors_4samples.rds")

#### select trans ####
trans_target <- trans[targets]
trans_used <-  trans_target[width(trans_target) >= 1000]

#### length ####
scores_used <-
  lapply(scores, function(i) {
    subset(i, trans %in% names(trans_used))
  })

#### grouping ####
basicgroups <-  cut(
  scores_used[["naive13 vs naive20"]]$scores,
  breaks = c(
    -Inf,
    quantile(scores_used[["naive13 vs naive20"]]$scores, 0.2),
    quantile(scores_used[["naive13 vs naive20"]]$scores, 0.4),
    quantile(scores_used[["naive13 vs naive20"]]$scores, 0.6),
    quantile(scores_used[["naive13 vs naive20"]]$scores, 0.8),
    Inf
  ),
  labels = bin_names
)
names(basicgroups) <-
  as.character(scores_used[["naive13 vs naive20"]]$trans)
scores_group <- lapply(scores_used,
                       function(i) {
                         i$groups <- basicgroups[as.character(i$trans)]
                         i
                       })
scores_group_long <- Reduce(rbind, scores_group, data.frame())
scores_group_long$len <-
  width(trans[as.character(scores_group_long$trans)])

#### save files ####
saveRDS(scores_group_long, "./results/scores_group_4samples.rds")
saveRDS(basicgroups, "./results/basicgroups.rds")
saveRDS(trans_used, "./check/transposon_uesed.rds")

```

### **code** dat.dscores

```{R}
#### define args ####
bin_names_2 <- c("Bin 1", "Bin 2", "Bin 3", "Bin 4", "Bin 5")

#### load files ####
scores_group_long <- readRDS("./results/scores_group_4samples.rds")

#### change to names shown in fig ####
scores_group_long$sample_set <- factor(scores_group_long$sample_set,
                                       levels = c("naive13 vs dsRluc13","naive20 vs dsRluc20",
                                                  "naive13 vs naive20", "dsRluc13 vs dsRluc20"),
                                       labels = c("Naive13 vs dsRluc13","Naive20 vs dsRluc20",
                                                  "Naive13 vs Naive20", "dsRluc13 vs dsRluc20"))

scores_group_long$groups <- factor(scores_group_long$groups, levels = bin_names,
                                   labels = bin_names_2)
#### add to rdata ####
dat.dscore <- scores_group_long
save(list = ls(pattern = "dat"), file = "data.RData")
```

### **fig** Dscore_naive13vsnaive20

```{r}
ggplot(subset(dat.dscore, sample_set == "Naive13 vs Naive20")) +
  geom_histogram(aes(x = as.numeric(scores), fill = groups),
                 binwidth = 0.3) +
  scale_x_continuous(limits = c(0, 100)) +
  scale_y_continuous(limits = c(0, 12), breaks = c(0, 2, 4, 6, 8, 10)) +
  ggtitle("") +
  xlab("Dscore") +
  ylab("Numbers of TE") +
  scale_fill_manual(values = c("#401962", "#466E8B", "#70BA73", "#A2D34F", "#dde136")) +
  # guides(fill = guide_legend(title=element_blank()))+
  guides(fill = "none") +
  mytheme() +
  theme(legend.position = c (.8, .75),
        legend.key.size = unit(0.5, "cm"))
  
ggsave(filename = "./fig/Dscore_naive13vsnaive20.pdf",
       width = Width,
       height = Height * 1.2)
```

### **fig** Dscore_dsRluc13vsdsRluc20

```{r}
ggplot(subset(dat.dscore, sample_set == "dsRluc13 vs dsRluc20")) +
  geom_histogram(aes(x = as.numeric(scores), fill = groups),
                 binwidth = 0.3) +
  scale_x_continuous(limits = c(0, 100)) +
  ggtitle("dsRluc13 vs. dsRluc20") +
  xlab("Dscore") +
  ylab("Numbers of TE") +
  scale_fill_manual(values = c("#401962", "#466E8B", "#70BA73", "#A2D34F", "#dde136")) +
  scale_y_continuous(breaks = c(0, 5, 10), limits = c(0, 12)) +
  guides(fill = "none") +
  facet_grid(groups ~ .) +
  mytheme() +
  # theme(strip.text.y = element_blank())+
  theme(legend.position = "right",
        legend.margin = margin(0, 0, 0, -10))

ggsave(filename = "./fig/Dscore_dsRluc13vsdsRluc20.pdf",
       width = Width,
       height = Height * 2)
```

### **fig** Dscore_nc

```{r}
ggplot(
  subset(
    dat.dscore,
    sample_set == "Naive13 vs dsRluc13" |
      sample_set == "Naive20 vs dsRluc20"
  )
) +
  geom_histogram(aes(x = as.numeric(scores)),
                 binwidth = 0.5,
                 fill = "black") +
  scale_x_continuous(limits = c(0, 100)) +
  xlab("Dscore") +
  ylab("Numbers of TE") +
  # scale_fill_manual(values ="black")+
  guides(fill = guide_legend(title = element_blank())) +
  facet_wrap(sample_set ~ ., nrow = 2) +
  mytheme()

ggsave(filename = "./fig/Dscore_nc.pdf",
       width = Width,
       height = Height * 1.7)
```

## piRNA & RNA vs Dscores

### *code* rpkms

```{r eval=FALSE, include=FALSE}
#### define args ####
sample_names <- c("naive13", "naive20")
names(sample_names) <- sample_names

# rna_names <-  c("close13", "close20")
# names(rna_names) <- rna_names
rna_names <-  c("close13", "close20_sampling_100nt")
names(rna_names) <- rna_names

#### load files ####
trans_used <- readRDS("./check/transposon_uesed.rds")
basicgroups <- readRDS("./results/basicgroups.rds")
scores_group_long <- readRDS("./results/scores_group_4samples.rds")
pirna_rpkm <- readRDS("./check/piRNA_rpkm.rds")

rna_reads <- lapply(rna_names, function(i) {
  read.table(sprintf("./bed/%s_transposon_RNA_HS2.bed", i))
})

# anti_reads <- lapply(rna_names, function(i) {
#   readRDS(sprintf("./rds/%s_anti_transposon_RNA_HS2.rds", i))
# })

#### rpkm ####
pirna_rpkm_used <- pirna_rpkm[names(trans_used), ]

# rna_rpkm <- sapply(rna_names, function(rn) {
#   sense <- length(unique(sense_reads[[rn]][, 4]))
#   anti <- length(unique(anti_reads[[rn]][, 4]))
#   sumreads <- sense + anti
#   reads_sense <-
#     table(factor(sense_reads[[rn]][, 1], levels = names(trans)))
#   reads_anti <-
#     table(factor(anti_reads[[rn]][, 1], levels = names(trans)))
#   rpm <- (reads_sense + reads_anti) / sumreads * 10 ^ 6
#   rpkm <- rpm / width(trans) * 10 ^ 3
#   rpkm
# })

rna_rpkm <- sapply(rna_names, function(rn) {
  sumreads <- length(unique(rna_reads[[rn]][, 4]))
  reads <- table(factor(rna_reads[[rn]][, 1], levels = names(trans)))
  rpm <- reads / sumreads * 10 ^ 6
  rpkm <- rpm / width(trans) * 10 ^ 3
  rpkm
})
rna_rpkm_used <- as.data.frame(rna_rpkm[names(trans_used), ])

#### reform to long frame ####
pirna_rpkm_used <- as.data.frame(pirna_rpkm_used)
pirna_rpkm_used$trans <- rownames(pirna_rpkm_used)
pirna_rpkm_used$type <- "piRNA"
pirna_rpkm_used$scores <-
  subset(scores_group_long, sample_set == "naive13 vs naive20")[pirna_rpkm_used$trans, "scores"]
pirna_rpkm_long <-
  pivot_longer(
    pirna_rpkm_used,
    col = 1:2,
    names_to = "year",
    values_to = "pirna_rpkm"
  )
pirna_rpkm_long$year <- rep(c(13, 20), nrow(pirna_rpkm_long) / 2)


rna_rpkm_long <-
  pivot_longer(
    rna_rpkm_used,
    col = 1:2,
    names_to = "year",
    values_to = "rna_rpkm"
  )

rpkm_long <-
  cbind(pirna_rpkm_long[, -c(1, 2)], "rna_rpkm" = rna_rpkm_long$rna_rpkm)
rpkm_long$group <- basicgroups[rpkm_long$trans]

#### save files ####
saveRDS(rpkm_long, "./results/piRNAandRNA_rpkm.rds")
saveRDS(rna_rpkm, "./check/rna_rpkm.rds")
```

### **code** dat.rpkms

```{R}
#### define args ####
bin_names_2 <- c("Bin 1", "Bin 2", "Bin 3", "Bin 4", "Bin 5")

#### load fiels ####
rpkm_long <- readRDS("./results/piRNAandRNA_rpkm.rds")

#### count the mean of rpkm ####
rpkm_long$year <- rep(c("naive13", "naive20"), nrow(rpkm_long) / 2)
rpkm_13 <- subset(rpkm_long, year == "naive13")
rpkm_20 <- subset(rpkm_long, year == "naive20")
rpkm_mean <- rpkm_13
rpkm_mean$pirna_rpkm <- (rpkm_13$pirna_rpkm + rpkm_20$pirna_rpkm) / 2
rpkm_mean$rna_rpkm <- (rpkm_13$rna_rpkm + rpkm_20$rna_rpkm) / 2

rpkm_mean$group <- factor(rpkm_mean$group, levels = bin_names,
                          labels = bin_names_2)

#### add to rdata ####
dat.rpkms <- rpkm_mean
save(list = ls(pattern = "dat"), file = "data.RData")
```

### **fig** data_gradient

```{r}
ggplot(dat.rpkms) +
  geom_point(aes(x = pirna_rpkm, y = rna_rpkm, color = scores)) +
  scale_x_log10(
    breaks = c(1, 10, 100, 1000, 10000),
    labels = c(1, 10, expression(10 ^ 2), expression(10 ^
                                                       3), expression(10 ^ 4))
  ) +
  scale_y_log10(
    breaks = c(1, 10, 100, 1000, 10000),
    labels = c(1, 10, expression(10 ^ 2), expression(10 ^
                                                       3), expression(10 ^ 4))
  ) +
  xlab(expression(RPKM[piRNA])) +
  ylab(expression(RPKM[transcripts])) +
  labs(color = "Dscore") +
  scale_color_gradientn(colors = c("#401962", "#466E8B", "#70BA73", "#A2D34F", "#dde136")) +
  mytheme() +
  theme(
    axis.text.x = element_text(hjust = 0.5, vjust = 0),
    # 调整x轴刻度标签水平居中，垂直底部对齐
    axis.text.y = element_text(hjust = 1, vjust = 0.5),
    axis.title.x = element_text(hjust = 0.5, vjust = 0),
  )

ggsave(filename = "./fig/data_gradient_mean.pdf",
       width = Width * 1.32 ,
       height = Height * 1.2)
```

### **fig** piRNA_scores

```{r}
ggplot(dat.rpkms) +
  geom_point(aes(x = pirna_rpkm, y = scores, col = group)) +
  scale_x_log10(
    # limits = c(0,500000),
    breaks = c(1, 10, 100, 1000, 10000),
    labels = c(1, 10, expression(10 ^ 2), expression(10 ^ 3), expression(10 ^
                                                                           4))
  ) +
  xlab(expression(RPKM[piRNA])) +
  ylab("Dscore") +
  # scale_color_gradientn(colors = c( "blue","red","#FFD700"),
  #                       # limits=c(0, 100),
  #                       guide = "none") +
  scale_color_manual(values = c("#401962", "#466E8B", "#70BA73", "#A2D34F", "#dde136")) +
  # guides(color = guide_legend(title=element_blank()))+
  guides (color = "none") +
  mytheme() +
  theme(
    # legend.position = c (.8, .75),
    #     legend.key.size = unit(0.5,"cm"),
    #     # legend.spacing.y = unit(0.01, "cm"),
    #     legend.margin = margin(0, 0, 0, -10),
    axis.text.x = element_text(hjust = 0.5, vjust = 0),
    axis.text.y = element_text(hjust = 1, vjust = 0.5),
    axis.title.x = element_text(hjust = 0.5, vjust = 0)
  )

ggsave(filename = "./fig/piRNA_scores_mean.pdf",
       width = Width,
       height = Height * 1.2)
```

## RNA MA

### preprocessing.sh

```{bash}
##close20_fastq_preparation.sh
NN="1 2"
for N in ${NN};do
echo ${N}
seqtk sample -s100 ./fastq/close20_R${N}.fastq.gz 0.1 > ./fastq/close20_sampling_R${N}.fastq
seqtk trimfq ./fastq/close20_sampling_R${N}.fastq -L 100 > ./fastq/close20_sampling_100nt_R${N}.fastq
gzip ./fastq/close20_sampling_100nt_R${N}.fastq > ./fastq/close20_sampling_100nt_R${N}.fastq.gz
done
##mapping
MM="close20_sampling_100nt close13"
LL="transposon"
hisat2-build -f ./fasta/transposon.fasta ${LL} &
for M in ${MM}; do
for L in ${LL}; do
#hisat2 -x ${L} -1 ./fastq/${M}_R1.fastq.gz  -2 ./fastq/${M}_R2.fastq.gz -k 3 -p 10 -S ./sambam/${M}_${L}_RNA_HS2.sam
#samtools view -@ 10 -bS ./sambam/${M}_${L}_RNA_HS2.sam > ./sambam/${M}_${L}_RNA_HS2.bam
bamToBed -i ./sambam/${M}_${L}_RNA_HS2.bam > ./bed/${M}_${L}_RNA_HS2.bed
#grep '\+' <./bed/${M}_${L}_RNA_HS2.bed | grep '/1' >./bed/${M}_${L}_RNA_R1_sense.bed
#grep '\-' <./bed/${M}_${L}_RNA_HS2.bed | grep '/1' >./bed/${M}_${L}_RNA_R1_anti.bed
#grep '\+' <./bed/${M}_${L}_RNA_HS2.bed | grep '/2' >./bed/${M}_${L}_RNA_R2_sense.bed
#grep '\-' <./bed/${M}_${L}_RNA_HS2.bed | grep '/2' >./bed/${M}_${L}_RNA_R2_anti.bed
rm ./sambam/${M}_${L}_RNA_HS2.sam
done
done 
```

### preprocessing_aftermapping.R

```{r}
#### extract piRNA ####
setwd("/home/Jie/Jupiter/paper/v3")
pkgs <- c("Biostrings","tidyr","DECIPHER")
sapply(pkgs, library, character.only = TRUE)

trans <- readDNAStringSet("./fasta/transposon.fasta")

sample_names <- c("close13", "close20_sampling_100nt")
read_number <- c("1", "2")
direc <- c("sense", "anti")
for (sn in sample_names) {
  for (rn in read_number) {
    for (d in direc) {
      bed <- read.table(sprintf("./bed/%s_transposon_RNA_R%s_%s.bed", sn, rn, d))
      message("readed bed file", sprintf("./bed/%s_transposon_RNA_R%s_%s.bed", sn, rn, d))
      bed$len <- bed[, 3] - bed[, 2]
      bed[, 1] <- factor(bed[, 1], levels = names(trans))
      bed[, 2] <- factor(bed[, 2] + 1 , levels = 1:max(width(trans)))
      bed[, 3] <- factor(bed[, 3] , levels = 1:max(width(trans)))
      mRNA <- subset(bed, len == 100)
      saveRDS(mRNA, sprintf("./rds/%s_R%s_%s_RNA_bed.rds", sn, rn, d))
    }
  }
}
```

### *code* rna_rpkm

```{r eval=FALSE, include=FALSE}
#### load files ####
rna_rpkm <- readRDS("./check/rna_rpkm.rds")

#### count the mean of rpkm ####
rna_rpkm <- as.data.frame(rna_rpkm)
rna_rpkm$trans <- rownames(rna_rpkm)
rna_rpkm$M <- log2(rna_rpkm[,"close13"] + 0.01) - log2(rna_rpkm[, "close20_sampling_100nt"] + 0.01)
rna_rpkm$A <- log2(rna_rpkm[, "close13"] + 0.01) / 2 + log2(rna_rpkm[, "close20_sampling_100nt"] + 0.01) / 2
rna_rpkm$group <- basicgroups[rna_rpkm$trans]
# rna_rpkm$globlegroups <- globlegroups[rna_rpkm$trans]

#### save files ####
saveRDS(rna_rpkm,"./results/rna_rpkm.rds")
```

### **code** dat.rnarpkm

```{r}
#### load files ####
rna_rpkm <- readRDS("./results/rna_rpkm.rds")

#### set color region to frame ####
rna_rpkm$range <- "total"
rna_rpkm[is.na(rna_rpkm$group),"range"] <- "out"
rna_rpkm <- dat.rnarpkm_100 %>% arrange(range)

#### add to rdata ####
dat.rnarpkm_100 <- rna_rpkm
save(list = ls(pattern = "dat"), file = "data.RData")
```

### **fig** MAplot_RNA

```{r}
ggplot(dat.rnarpkm_100) +
  geom_point(aes(x = A, y = M, color = range), size = 1, shape = 21) +
  # facet_wrap(.~group)+
  ggtitle("RNA expression levels") +
  xlab(expression(Log[2] ~ mean ~ expression)) +
  ylab(expression(Log[2] ~ fold ~ change)) +
  geom_hline(yintercept = c(-3, 3),
             linetype = "dashed",
             col = "black") +
  # geom_hline(yintercept = 0)+
  geom_vline(xintercept = -0,
             linetype = "dashed",
             col = "black") +
  scale_color_manual(values = c("darkgrey", "black"),
                     guide = "none") +
  mytheme()

ggsave(filename = "./fig/MAplot_RNA_100nt.pdf",
       width = Width,
       height = Height)
```

## ping-pong signature

### **code** dat.pingpong

```{r}
#### define args ####
sample_names <- c("naive13", "naive20")
names(sample_names) <- sample_names

bin_names_2 <- c("Bin 1", "Bin 2", "Bin 3", "Bin 4", "Bin 5")

#### load files ####
trans <-  readDNAStringSet("./fasta/transposon.fasta")
trans_used <- readRDS("./check/transposon_uesed.rds")
basicgroups <- readRDS("./results/basicgroups.rds")

plus_rpm <-
  lapply(sample_names, function(i) {
    readRDS(sprintf("./piRNA/%s_transposon_piRNA_plus_rpm.rds", i))
  })
minus_rpm <-
  lapply(sample_names, function(i) {
    readRDS(sprintf("./piRNA/%s_transposon_piRNA_minus_rpm.rds", i))
  })

#### count rpkm ####
plus_rpkm <- lapply(sample_names, function(sn) {
  rpkm <-   plus_rpm[[sn]]  / width(trans) * 10 ^ 3
  rpkm[names(trans_used),]
})
minus_rpkm <- lapply(sample_names, function(sn) {
  rpkm <-   minus_rpm[[sn]]  / width(trans) * 10 ^ 3
  rpkm[names(trans_used),]
})

#### ping-pong signature ####
maxlen <- max(width(trans_used))
pingpongpair <- lapply(sample_names, function(sn) {
  pp <- sapply(1:20, function(site) {
    signvalue <-
      plus_rpkm[[sn]][, 1:(maxlen - site + 1)] * minus_rpkm[[sn]][, site:maxlen]
    signvalue_sum <- apply(signvalue, 1, sum)
    signvalue_sum
  })
  pp.na <- subset(pp, apply(pp, 1, sum) != 0)
  
  pp_ratio <- sweep(pp.na, 1, apply(pp.na, 1, sum), "/")
  pp_mean <- as.data.frame(pp_ratio)
  
  pp_mean$year <- sn
  pp_mean$trans <- rownames(pp_mean)
  pp_mean$groups <- basicgroups[as.character(pp_mean$trans)]
  
  pp_group <- pp_mean[!is.na(pp_mean$groups), ]
  pp_range <-  arrange(pp_group, V10)
  
  pp_group$trans <- factor(pp_group$trans, levels = pp_range$trans)
  pivot_longer(
    pp_group,
    col = 1:20,
    values_to = "coverage",
    names_to = "site"
  )
})
#### change to long frame ####
pppair <- Reduce(rbind, pingpongpair, data.frame())
pppair$site <- rep(1:20, nrow(pppair) / 20)

#### change to names shown in fig ####
pppair$groups <- factor(pppair$groups,
                              levels = bin_names,
                              labels = bin_names_2)
#### add to rdata ####
dat.pingpong <- pppair
save(list = ls(pattern = "dat"), file = "data.RData")

```

### **fig** separat figuires pingpone_signature

```{r}
group_names <-  c("Bin 1", "Bin 2", "Bin 3", "Bin 4", "Bin 5")
names(group_names) <- group_names

lapply(group_names, function(gn) {
  fig <- ggplot(subset(dat.pingpong, groups == gn &
                         year == "naive13")) +
    geom_tile(aes(x = site,
                  y = trans,
                  fill = coverage)) +
    scale_fill_gradientn(
      colors = c("#EBF1F6", "#3D5E7E"),
      limits = c(0, 1),
      # guide = "colorbar",
      guide = "none"
    ) +
    # facet_grid(.~year)+
    ggtitle("") +
    xlab("Overlap (nt)") +
    ylab("") +
    scale_y_discrete() +
    scale_x_continuous(limits = c(0,20), breaks = seq(0, 20, by = 2)) +
    theme_void()
  # theme(axis.line.y = element_blank(),
  #       axis.text.y = element_blank(),
  #       axis.ticks.y = element_blank(),
  #       title = element_blank())
  
  # ggsave(
  #   filename = sprintf("./fig/figS2_B_%s.pdf", x),
  #   width = Width,
  #   height = Height * 0.7
  # )
  print(fig)
  ggsave(
    fig,
    filename = sprintf("./fig/pingpone_signature_%s.pdf", gn),
    width = Width * 0.8,
    height = Height * 0.6
  )
})
```

### **fig** single figure for pingpong_signature

```{r}
ggplot(subset(dat.pingpong, year == "naive13")) +
  geom_tile(aes(
    x = site,
    y = trans,
    fill = coverage,
    col = coverage
  )) +
  scale_fill_gradientn(
    colors = c("#EBF1F6", "#3D5E7E"),
    limits = c(0, 1),
    # guide = "colorbar",
    guide = "none"
  ) +
  scale_color_gradientn(
    colors = c("#EBF1F6", "#3D5E7E"),
    limits = c(0, 1),
    # guide = "colorbar",
    guide = "none"
  ) +
  facet_grid(groups ~ .) +
  ggtitle("") +
  xlab("Overlap (nt)") +
  ylab("") +
  scale_y_discrete() +
  scale_x_continuous(limits = c(0, 21), breaks = seq(1, 20, by =2)) +
  mytheme() +
  theme(
    axis.line.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

ggsave(
  filename = sprintf("./fig/pingpone_signatures.pdf"),
  width = Width * 1,
  height = Height * 3
)


```

## mismatches

### *code* seq library

trim fasta to piRNA length

```{r}
#### define args ####
sample_names <- c("naive13", "naive20", "dsRluc13", "dsRluc20")
names(sample_names) <- sample_names

#### trim seq to piRNA length ####
lapply(sample_names, function(sn) {
  # Print a message to the console indicating the current sample name
  message(sn)
  # Read DNA sequences from a gzipped FASTA file using Bioconductor's Biostrings package
  fas <- readDNAStringSet(sprintf("./fasta/%s_trim.fasta.gz", sn))
  # Select sequences with lengths between 26 and 32 (inclusive)
  fas_pi <- fas[26 <= width(fas) & width(fas) <= 32]
  # Extract the first 26 bases from each selected sequence
  fas_pi_26 <- subseq(fas_pi, 1, 26)
  # Write the extracted sequences to a new FASTA file
  writeXStringSet(fas_pi_26, sprintf("./fasta/%s_trim_pi26.fasta", sn))
})
```

### *code* zip fasta

```{bash}
gzip *_trim_pi26.fasta
```

### *code* split datas

divide bed files and fasta files into each bins

```{r}
#### define args ####
sample_names <- c("naive13", "naive20")
names(sample_names) <- sample_names

bin_names <- c("bin1", "bin2", "bin3", "bin4", "bin5")
names(bin_names) <- bin_names

#### load files ####
basicgroups <- readRDS("./results/basicgroups.rds")
basicgroups_list <- lapply(bin_names,
                           function(bn){
                             names(basicgroups[basicgroups == bn])
                           })

#### extract piRNA reads from bed file ####
lapply(sample_names, function(sn) {
  bed <- readRDS(sprintf("./rds/%s_transposon.rds", sn))
  message(sn, "readed bed file")
  bed$len <- bed[, 3] - bed[, 2]
  bed[, 1] <- factor(bed[, 1], levels = names(trans))
  bed[, 2] <- factor(bed[, 2] + 1 , levels = 1:max(width(trans)))
  bed[, 3] <- factor(bed[, 3] , levels = 1:max(width(trans)))
  piRNA <- subset(bed, 26 <= len & len <= 32)
  saveRDS(piRNA, sprintf("./piRNA/%s_piRNA_bed.rds", sn))
  rm("bed")
})

#### split fasta into each bins ####
for (years in c(13, 20)) {
  #### load ####
  piRNA <-
    readRDS(sprintf("./piRNA/naive%s_piRNA_bed.rds", years))
  message("naive", years, " readed bed file")
  fas <-
    readDNAStringSet(sprintf("./fasta/naive%s_trim_pi26.fasta.gz", years))
  message("naive", years, " readed fasta file")
  piRNA$V4 <- as.character(piRNA$V4)
  for (bn in bin_names) {
    piRNA_used <-
      unique(piRNA[piRNA$V1 %in% basicgroups_list[[bn]], "V4"])
    fas.piRNA <-  fas[piRNA_used]
    message("naive", years, " ", bn, " sequence extracted")
    writeXStringSet(fas.piRNA,
                    sprintf("./fasta/naive%s_pi26_%s.fasta", years, bn))
  }
}
```

### *code* consensus_sequence.R

```{r tmp }
#### basic config ####
setwd("/home/Jie/Jupiter/paper/v3")

#### load packages ####
pkgs <- c("Biostrings","tidyr","DECIPHER")
sapply(pkgs, library, character.only = TRUE)

#### define args ####
bin_names <- c("bin1", "bin2", "bin3", "bin4", "bin5")
names(bin_names) <- bin_names

args <- commandArgs(trailingOnly = TRUE)
# args <- c(13,"bin1","sense")
years <- args[1]
bn <- args[2]
direc <- args[3]

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

#### load files ####
trans <- readDNAStringSet("./fasta/transposon.fasta")

basicgroups <- readRDS("./results/basicgroups.rds")
basicgroups_list <- lapply(bin_names,
                           function(bn) {
                             namelist <- names(basicgroups[basicgroups == bn])
                             as.numeric(factor(namelist, levels = names(trans)))
                           })

piRNA_pm <- readRDS(sprintf("./piRNA/naive%s_piRNA_bed.rds", years))
message("naive", years, " ", bn, " ", direc, " readed bed file")

fas <-
  readDNAStringSet(sprintf("./fasta/naive%s_pi26_%s.fasta", years, bn))
message("naive", years, " ", bn, " ", direc, " readed fasta file")

#### find read peaks for TE in each bins ####

# Subset the 'piRNA_pm' data frame based on the condition V6 == c
piRNA <- subset(piRNA_pm, V6 == c)

# Convert 'V4' column to character and 'trans' column to numeric
piRNA$V4 <- as.character(piRNA$V4)
piRNA$trans <- as.numeric(piRNA$V1)

rm("piRNA_pm")
message("naive", years, " ", bn, " ", direc, " piRNA ", direc)

# Create a table of piRNA reads based on specified columns
piRNA.reads <- table(piRNA[, c(1, e)])

# Set the count for a specific row and column to 0
piRNA.reads["TE1_bm_1645_LINE/R4", ] <- 0
message("naive", years, " ", bn, " ", direc, " removed 1645")

# Find indices of elements in 'piRNA.reads' greater than 3
peaks <- which(piRNA.reads > 3, arr.ind = TRUE)

# Subset the 'peaks' matrix to include only rows corresponding to specific groups
peaks_bin <- peaks[peaks[, 1] %in% basicgroups_list[[bn]], ]
message("naive", years, " ", bn, " ", direc, " peaks : ", nrow(peaks_bin))

#### consensus sequence ####
conseq <- character(nrow(peaks_bin))

timestart <- Sys.time()
message("naive", years, " ", bn, " ", direc, " start time : ", timestart)

conseq <- apply(peaks_bin, 1, function(pb) {
  
  # Extract values from the current row of 'peaks_bin'
  tt <- pb[1]
  dd <- pb[2]
  
  # Create a logical mask based on the conditions
  mask <- piRNA$trans == tt & piRNA[, d] == dd
  
  # Get piRNA names based on the mask
  pinames <- piRNA$V4[mask]
  
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
    ifelse(length(max_seq) > 1,
           Reduce(function(x, y) {
             paste0(x, strrep("N", 24), y)
           }, max_seq),
           max_seq)
  } else {
    cat("-")
    conseq_tmp
  }
})


timeend <- Sys.time()
message("naive", years, " ", bn, " ", direc, "end time : ", timeend)

runningtime <- timeend - timestart
print(runningtime)

#### change the data.frame to Xstring ####
conseq <- DNAStringSet(conseq)
names(conseq) <- paste(peaks_bin[, 1], peaks_bin[, 2], f, sep = "_")

#### save files ####
writeXStringSet(conseq,sprintf("./check/naive%s_%s_%s_consensus_sequence.fasta",years,bn,direc))
message("naive",years," ",bn," ",direc," saved all results")
```

### *code* conbind sense and antisense

```{r}
#### define args ####
bin_names <- c("bin1", "bin2", "bin3", "bin4", "bin5")
names(bin_names) <- bin_names

sample_names <- c("naive13", "naive20")
names(sample_names) <- sample_names

#### bind ####
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
```

### *code* mismatches mapping

```{bash}
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

```

### *code* bash run

```{bash}
# set pwd
cd /home/Jie/Jupiter/paper/v3

# seq library 
nohup Rscript ./funs/seq_library.R > ./runlog/seq_library.log 2>&1 &

# split datas
nohup Rscript ./funs/split_datas.R > ./runlog/split_datas.log 2>&1 &

wait

# zip fasta
gzip ./fasta/*_trim_pi26.fasta
echo "raw fasta are ready"

# consesus_sequence.R
YY=(13 20)
BB=("bin1" "bin2" "bin3" "bin4" "bin5")
DD=("sense" "anti")
echo "Following commands starting..."
for Y in "${YY[@]}"; do
  # 启动多个任务
  for B in "${BB[@]}"; do
    for D in "${DD[@]}"; do
      echo ${Y}${B}${D}
      nohup Rscript ./funs/consensus_sequence.R ${Y} ${B} ${D} > ./runlog/0823_${Y}${B}${D}.out &
      pids+=($!)
    done
  done
  # 等待当前元素的所有任务完成
  for pid in "${pids[@]}"; do
    wait "$pid"
  done
done

# conbind sense and antisense
echo "bind datas"
nohup Rscript ./funs/bind.R > ./runlog/bind.log 2>&1 &

wait

# mismatches mapping
nohup bash ./funs/mapping.sh > ./runlog/mapping.log 2>&1 &
echo "missions done"

#### run 
nohup bash ./funs/0823.sh > ./runlog/script_output.log 2>&1 &
```

### code bed_con_all

```{r}
#### define args ####
bin_names <- c("bin1", "bin2", "bin3", "bin4", "bin5")
names(bin_names) <- bin_names

sample_names <- c("naive13", "naive20")
names(sample_names) <- sample_names

bed_con_all <- list()

#### con_seq ####
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

#### 159 ####
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
```

### <code> v123 total_reads

```{r}
#### misratio ####

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
```

### <code> postions_info

```{bash}
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
```

### <code> pos_results_long

```{r}
bin_names <- c("bin1", "bin2", "bin3", "bin4", "bin5")
names(bin_names) <- bin_names
sample_names <- c("naive13","naive20")
names(sample_names) <- sample_names
misnum <- c(1,2,3)

sn ="naive13"
results <- list()
  for (bn in bin_names){
    nn <- paste0(sn,"_",bn)
    results[[bn]] <- matrix(ncol = 3, nrow = 26,0)
      mis1 <- read.table(sprintf("./txt/%s_%s_v%s_%s.txt",sn,bn,1,1))
      mis1[,2] <-factor(mis1[,2],levels = 1:26)
      mis1 <- mis1[!as.character(mis1[,1]) %in% bed_con_all[[sn]], ]
      results[[bn]][,1] <- table(mis1[,2])
    
      mis2 <- read.table(sprintf("./txt/%s_%s_v%s_%s.txt",sn,bn,2,2))
      mis2[,2] <-factor(mis2[,2],levels = 1:26)
      mis2[,3] <-factor(mis2[,3],levels = 1:26)
      mis2 <- mis2[!as.character(mis2[,1]) %in% bed_con_all[[sn]], ]
      results[[bn]][,2] <- table(mis2[,2])+table(mis2[,3])
      
      mis3 <- read.table(sprintf("./txt/%s_%s_v%s_%s.txt",sn,bn,3,3))
      mis3[,2] <-factor(mis3[,2],levels = 1:26)
      mis3[,3] <-factor(mis3[,3],levels = 1:26)
      mis3[,4] <-factor(mis3[,4],levels = 1:26)
      mis3 <- mis3[!as.character(mis3[,1]) %in% bed_con_all[[sn]], ]
      results[[bn]][,3] <- table(mis3[,2])+table(mis3[,3])+table(mis3[,4])
      
      results[[bn]] <- as.data.frame(results[[bn]])
  }
  for (bn in bin_names){
      results[[bn]][,4] <- apply(results[[bn]],1,sum)
      results[[bn]] <- sweep( results[[bn]],2,apply( results[[bn]][c(2:9,11:26),],2,sum),"/") *100
  }

pos_results_13 <-  lapply(1:5,function(i){
  re <- as.data.frame(results[[i]])
    y <- data.frame(pos1 = re[1,4],
                    seed = sum(re[2:9,4]),
                  pos10 = re[10,4],
                 nonseed1 = sum(re[11:18,4]),
                 nonseed2 = sum(re[19:26,4]),
                 bins = names(bin_names[i]),
                 sample = sn)
  y
})
sn ="naive20"
for (bn in bin_names){
    nn <- paste0(sn,"_",bn)
    results[[bn]] <- matrix(ncol = 3, nrow = 26,0)
      mis1 <- read.table(sprintf("./txt/%s_%s_v%s_%s.txt",sn,bn,1,1))
      mis1[,2] <-factor(mis1[,2],levels = 1:26)
      mis1 <- mis1[!as.character(mis1[,1]) %in% bed_con_all[[sn]], ]
      results[[bn]][,1] <- table(mis1[,2])
    
      mis2 <- read.table(sprintf("./txt/%s_%s_v%s_%s.txt",sn,bn,2,2))
      mis2[,2] <-factor(mis2[,2],levels = 1:26)
      mis2[,3] <-factor(mis2[,3],levels = 1:26)
      mis2 <- mis2[!as.character(mis2[,1]) %in% bed_con_all[[sn]], ]
      results[[bn]][,2] <- table(mis2[,2])+table(mis2[,3])
      
      mis3 <- read.table(sprintf("./txt/%s_%s_v%s_%s.txt",sn,bn,3,3))
      mis3[,2] <-factor(mis3[,2],levels = 1:26)
      mis3[,3] <-factor(mis3[,3],levels = 1:26)
      mis3[,4] <-factor(mis3[,4],levels = 1:26)
      mis3 <- mis3[!as.character(mis3[,1]) %in% bed_con_all[[sn]], ]
      results[[bn]][,3] <- table(mis3[,2])+table(mis3[,3])+table(mis3[,4])
      
      results[[bn]] <- as.data.frame(results[[bn]])
  }
  for (bn in bin_names){
      results[[bn]][,4] <- apply(results[[bn]],1,sum)
      results[[bn]] <- sweep( results[[bn]],2,apply( results[[bn]][c(2:9,11:26),],2,sum),"/") *100
  }
pos_results_20 <-  lapply(1:5,function(i){
  re <- as.data.frame(results[[i]])
    y <- data.frame(pos1 = re[1,4],
                    seed = sum(re[2:9,4]),
                  pos10 = re[10,4],
                 nonseed1 = sum(re[11:18,4]),
                 nonseed2 = sum(re[19:26,4]),
                 bins = names(bin_names[i]),
                 sample = sn)
  y
})
pos_results <- rbind(Reduce(rbind,pos_results_13),Reduce(rbind,pos_results_20))
pos_results_long <- pivot_longer(pos_results,
                                 col=1:5,
                                 values_to = "percent",
                                 names_to = "regions")
pos_results_long$region <- rep(c(NA,"seed",NA,"nonseed","nonseed"),10)
saveRDS(pos_results_long, "./results/mismatches_positions_ratio.rds")

```

### <code> dat.pos_mis

```{r}
pos_mis <- readRDS("./results/mismatches_positions_ratio.rds")
pos_mis$region <- factor(pos_mis$region,levels=c("seed","nonseed"))
pos_mis$regions <- factor(pos_mis$regions,
                      levels=c("seed","nonseed1","nonseed2"),
                      labels=c("5'","Middle","3'"))
pos_mis$sample <- factor(pos_mis$sample,
                         levels = c("naive13","naive20"),
                         labels = c("Naive13","Naive20"))
pos_mis$bins <- factor(pos_mis$bins,
                       levels = bin_names,
                       labels = bin_names_2)
#### save ####
dat.pos_mis <- subset(pos_mis,region%in%c("seed","nonseed"))
save(list = ls(pattern = "dat"), file = "data.RData")
```

### \|fig\| snps_bins_regions

```{r}
ggplot(dat.pos_mis)+
  geom_bar(aes(x=bins,y=percent,fill = regions),
           stat="summary")+
  facet_grid(sample ~ . ,scales="free",space="free")+
  scale_fill_manual(values = c("black", "#666666","#cccccc")) +
  ggtitle("")+
  xlab("")+
  ylab(expression("Percentage of mismatches"))+
  labs(fill="Regions")+
  guides(fill="none")+
  mytheme()+
  # theme(axis.text.x=element_text(angle=45, hjust=1))
  theme(plot.margin = margin(l = 13))
ggsave(filename = "./fig/snps_bins_regions.pdf", width = Width*1.1, height = Height*1.8)
```

## Simulation

### **code** P0

```{r}
#### define args ####
Nmax <- 500
Gmax <- 10000 
Rep <- 100000 
Num <- 1

#### probability simulation ####
P0 <- sapply(1:Nmax,
             function(N) {
               results <- sapply(1:Rep,
                                 function(x) {
                                   Frequency <- (N - 1) / N
                                   i <- 0
                                   # Iterate until reaching the maximum number of generations (Gmax) or typeA frequency becomes 0 or 1
                                   while (i <= Gmax &
                                          Frequency != 0 & Frequency != 1) {
                                     Num <- rbinom(1, N, Frequency)
                                     Frequency <- Num / N
                                     i <- i + 1
                                   }
                                   Num
                                 })
               # Calculate the probability of complete elimination of typeA
               out <- sum(results == 0) / Rep
               cat(N)
               out
             })

#### save files ####
saveRDS(P0,"./results/simulation_P0.rds")

#### add to rdata ###
dat.pzero <- P0
save(list = ls(pattern = "dat"), file = "data.RData")
```

### **code** Preplace

```{r}
#### define args ####
Nmax <- 500
Gmax <- 10000 
Rep <- 100000 
Num <- 1
occurrences <- 100

#### load files ####
P0 <- readRDS("./results/simulation_P0.rds")

#### simulation Preplacement ####
Preplacement <- sapply(1:occurrences,
                           function(n){
                              1-(1-P0)^n
                             })

#### change to long frame ####
Preplacement <- as.data.frame(Preplacement)
colnames(Preplacement) <- 1:occurrences
Preplacement$N <- 1:Nmax
Preplace_long <- pivot_longer(Preplacement,cols = 1:occurrences,names_to="Oemergence",values_to = "P")

dat.simulate <- Preplace_long
```

### **fig** P0

```{r}
ggplot(dat.pzero)+
  geom_point(aes(x=N,y=pro),size=1)+
  xlab(expression(N[piRNA]))+
  ylab(expression(P[0]))+
  scale_y_log10(breaks = c(0,0.01,0.1,1),
                labels = c(0,0.01,0.1,1))+
  # scale_x_log10()+
  # labs(color = expression(P[0]))+
  guides(color="none")+
   mytheme()
ggsave(filename = "./fig/pzero.pdf",width=Width,height=Height)
```

### **fig** Preplacement

```{r}
ggplot(dat.simulate, ) +
  geom_tile(aes(
    x = N,
    y = as.numeric(Oemergence),
    fill = P
  )) +
  scale_fill_gradientn(
    colors = c("#401962", "#466E8B", "#70BA73", "#A2D34F", "#DDDD22"),
    limits = c(0, 1),
    labels = c(0, 0.25, 0.5, 0.75, 1),
    guide = "colorbar"
  ) +
  xlab(expression(N[piRNA])) +
  ylab(expression(O[emerging])) +
  labs(fill = expression(P[replacement])) +
  mytheme()

ggsave(filename = "./fig/simulation.pdf",
       width = Width * 1.34,
       height = Height * 1.2)
```

## Fem

### **code** dat.fem & dat.masc

```{r}
#### define args ####
sample_names <- c("naive13", "naive20", "dsRluc13", "dsRluc20")
names(sample_names) <- sample_names

#### load files ####
femmasc <- readDNAStringSet("./fasta/femmasc.fas")
samples <-
  lapply(sample_names, function(i) {
    read.table(sprintf("./bed/%s_femmasc.bed", i))
  })

#### define & extract piRNA ####
samples_pi <- 
  lapply(samples,function(i){
  i[,"len"] <- i[,3] - i[,2]
  i[,2] <- i[,2]+1
  subset(i,V6 == "+" & 26 <= len & len <= 32)
})

#### distribution of piRNA reads ####
table_fem <- lapply(samples_pi, function(sp) {
  s_fem  <- subset(sp, V1 == "fem")
  sum_fem <- length(unique(s_fem[, 4]))
  t_fem <-
    table(factor(s_fem[, 2], levels = 1:width(femmasc[1]))) / sum_fem * 100
  as.data.frame(t_fem)
})
table_masc <- lapply(samples_pi, function(sp) {
  s_masc <- subset(sp, V1 == "masc")
  sum_masc <- length(unique(s_masc[, 4]))
  t_masc <-
    table(factor(s_masc[, 2], levels = 1:width(femmasc[2]))) / sum_masc * 100
  as.data.frame(t_masc)
})

#### adjust data.frame to figure ####
table_fem[[1]][,"samples"] <- "Naive13"
table_fem[[2]][,"samples"] <- "Naive20"
table_masc[[1]][,"samples"] <- "Naive13"
table_masc[[2]][,"samples"] <- "Naive20"

fem <- rbind(table_fem[[1]],table_fem[[2]])
fem[,1] <- as.numeric(fem[,1])
fem_peak <- subset(fem, Freq >25)
masc <- rbind(table_masc[[1]],table_masc[[2]])
masc[,1] <- as.numeric(masc[,1])
masc_peak <- subset(masc,Freq > 90)

#### add to rdata ####
dat.fem <- fem
dat.masc <- masc
dat.fem_peak <- fem_peak
dat.masc_peak <- masc_peak
save(list = ls(pattern = "dat"), file = "data.RData")
```

### **fig** fem

```{r}
ggplot(dat.fem) +
  geom_col(aes(x = Var1, y = Freq), col = "black", size = 0.35) +
  # geom_col(data = dat.fem_peak, aes(x=Var1,y=Freq),col = "#3A88FC")+
  ggtitle(expression(italic(Fem))) +
  xlab("") +
  ylab("Relative abundance (%)") +
  facet_grid(samples ~ .) +
  mytheme()

ggsave(filename = "./fig/fem.pdf",
       width = Width,
       height = Height * 1.3)
```

### **fig** masc

```{r}
ggplot(dat.masc) +
  geom_col(aes(x = Var1, y = Freq), col = "black", size = 0.35) +
  # geom_col(data = dat.masc_peak, aes(x=Var1,y=Freq),col = "#17B506")+
  ggtitle(expression(italic(Masc))) +
  xlab("") +
  ylab("Relative abundance (%)") +
  ylim(0, 101) +
  scale_x_continuous(labels = scales::comma) +
  facet_grid(samples ~ .) +
  mytheme()
ggsave(filename = "./fig/masc.pdf",
       width = Width,
       height = Height * 1.3)
```

## Pelement

### *code* preprocessing of fasta

```{bash}
#SRR333509 #small RNAs from Har x w1 21 day
#SRR333507 #small RNAs from Har x w1 2-4 day
#SRR333506 #small RNAs from w1 x Har 21 day
#SRR333505 #small RNAs from w1 x Har 2-4 day

#SRR333512 #genomic DNA from w1 x Har 21 day
#SRR333513 #genomic DNA sequencing - (w1 x Har) x w1 21 day
#SRR333514 #genomic DNA sequencing - (w1 x Har) x w1 21 day

## find the files of fasta
NN="SRR333505 SRR333506 SRR333507 SRR333509"
for N in  ${NN}; do 
    echo $N
    # Output the value of N
    prefetch ${N} -O  /home/kshoji/Gawain/NGSdata/sra/
    # Download SRA files for the specified N and save them to the specified directory
    fastq-dump ${N}.sra
    # Convert SRA files to FASTQ format for the specified N
done

NN="SRR333512 SRR333513 SRR333514"
for N in  ${NN}; do 
        echo $N
prefetch ${N} -O  /home/kshoji/Gawain/NGSdata/sra/
fastq-dump  --split-files ${N}.sra
pigz ${N}_1.fastq
pigz ${N}_2.fastq
done

## cut adaptamer
NN="SRR333505 SRR333506 SRR333507 SRR333509"
for N in  ${NN}; do
    echo $N
    cutadapt -a TCGTATGCCGTCTTCTGCTTG --minimum-length 20 --maximum-length 30 -o ${N}_trim.fastq ${N}.fastq
    fastq_to_fasta -Q33 -i ${N}_trim.fastq >  ${N}_trim.fasta
    pigz ${N}_trim.fasta
    rm ${N}_trim.fastq
done

## mapping 
bowtie-build Pelement.fasta Pelement
PP="Pelement"
for P in ${PP};do
    NN="SRR333505 SRR333506 SRR333507 SRR333509"
    for N in ${NN};do
        echo $N
        ( bowtie --offrate 3 -p 20 -a --best --strata -v 3 -m 1 -t --sam ${P} -f /home/kshoji/Gawain/NGSdata/sra/${N}_trim.fasta.gz  > ${N}_${P}_bt_m1.sam )  
        samtools view -@ 20 -bS ${N}_${P}_bt_m1.sam > ${N}_${P}_bt_m1.bam
        bamToBed -i ${N}_${P}_bt_m1.bam > ${N}_${P}_bt_m1.bed
        rm ${N}_${P}_bt_m1.sam
        rm ${N}_${P}_bt_m1.bam
    done
done


for P in ${PP};do
    NN="SRR333512 SRR333513 SRR333514"
    for N in ${NN};do
        echo $N
        ( bowtie --offrate 3 -p 20 -a --best --strata -v 3 -m 1 -t --sam ${P} -1 /home/kshoji/Gawain/NGSdata/sra/${N}_1.fastq.gz -2 /home/kshoji/Gawain/NGSdata/sra/${N}_2.fastq.gz  > ${N}_${P}_bt_m1.sam )  
        samtools view -@ 20 -bS ${N}_${P}_bt_m1.sam > ${N}_${P}_bt_m1.bam
        bamToBed -i ${N}_${P}_bt_m1.bam > ${N}_${P}_bt_m1.bed
        rm ${N}_${P}_bt_m1.sam
    done
done
```

### *code* extract sequence information from bam file

```{r}
#### load packages ####
library(Biostrings)
library(Rsamtools)

#### define args ####
file_names <- c("SRR333512", "SRR333513", "SRR333514")


for (name in file_mames) {
  # load files 
  bam <- scanBam(paste0(name, "_Pelement_bt_m1.bam"))
  # Extract relevant information from BAM files
  chr3 <- bam[[1]]$rname == "Pelemenet"
  chr3[is.na(chr3)] <- FALSE
  pos <- bam[[1]]$pos[chr3]
  qname <- bam[[1]]$qname[chr3]
  qwidth <- bam[[1]]$qwidth[chr3]
  cigar <- bam[[1]]$cigar[chr3]
  seq <- bam[[1]]$seq[chr3]
  qual <- bam[[1]]$qual[chr3]
  strand <- bam[[1]]$strand[chr3]
  seq2 <- seq
  # Insert "N" characters at specific positions in the sequences
  for (i in 1:2907) {
    tag <- pos == i
    ins <- paste(rep("N", i - 1), collapse = "")
    seq2[tag] <- paste0(ins, seq[tag])
    cat(i)
  }
  # Write sense-strand variant table to file
  write.table(
    t(consensusMatrix(seq2[strand == "+"])[1:4, ]),
    paste0("231026_", name, "_Pelement_variant_table_sense.txt"),
    sep = "\t",
    quote = FALSE
  )
  # Write antisense-strand variant table to file
  write.table(
    t(consensusMatrix(seq2[strand == "-"])[1:4, ]),
    paste0("231026_", name, "_Pelement_variant_table_antisense.txt"),
    sep = "\t",
    quote = FALSE
  )
}

# Load smRNAseq data for different samples
data1 <- read.table("231026_SRR333512_Pelement_variant_table_sense.txt")
data2 <- read.table("231026_SRR333513_Pelement_variant_table_sense.txt")
data3 <- read.table("231026_SRR333514_Pelement_variant_table_sense.txt")

data4 <- read.table("231026_SRR333512_Pelement_variant_table_antisense.txt")
data5 <- read.table("231026_SRR333513_Pelement_variant_table_antisense.txt")
data6 <- read.table("231026_SRR333514_Pelement_variant_table_antisense.txt")

# Load RNAseq data for different samples
S5 <- read.table("SRR333505_Pelement_bt_m1.bed")
S6 <- read.table("SRR333506_Pelement_bt_m1.bed")
S7 <- read.table("SRR333507_Pelement_bt_m1.bed")
S9 <- read.table("SRR333509_Pelement_bt_m1.bed")

# Create a dummy matrix to handle missing data
dummy <- matrix(0, 2907, 4)
colnames(dummy) <- colnames(data1)
rownames(dummy) <- 1:2907
data1a <- rbind(data1, dummy[(nrow(data1) + 1):2907, ])
data2a <- rbind(data2, dummy[(nrow(data2) + 1):2907, ])
data3a <- rbind(data3, dummy[(nrow(data3) + 1):2907, ])

# Create binary matrices based on certain conditions
num <- cbind(
  !apply(data1a, 1, sum) == apply(data1a, 1, max),
  !apply(data2a, 1, sum) == apply(data2a, 1, max),
  !apply(data3a, 1, sum) == apply(data3a, 1, max),
  !apply(data4, 1, sum) == apply(data4, 1, max),
  !apply(data5, 1, sum) == apply(data5, 1, max),
  !apply(data6, 1, sum) == apply(data6, 1, max)
)

# Further calculations and analysis...
# (The code continues with various calculations, filtering, and plotting, 
# with comments indicating the main steps of the analysis)

# Filter and extract specific regions for smRNA data
S5pis <- S5[S5[,3] - S5[,2] > 22 & S5[,6] == "+",]
S6pis <- S6[S6[,3] - S6[,2] > 22 & S6[,6] == "+",]
S7pis <- S7[S7[,3] - S7[,2] > 22 & S7[,6] == "+",]
S9pis <- S9[S9[,3] - S9[,2] > 22 & S9[,6] == "+",]
S5pias <- S5[S5[,3] - S5[,2] > 22 & S5[,6] == "-",]
S6pias <- S6[S6[,3] - S6[,2] > 22 & S6[,6] == "-",]
S7pias <- S7[S7[,3] - S7[,2] > 22 & S7[,6] == "-",]
S9pias <- S9[S9[,3] - S9[,2] > 22 & S9[,6] == "-",]

# Organize data into lists
pis <- list(S5pis, S6pis, S7pis, S9pis)
pias <- list(S5pias, S6pias, S7pias, S9pias)

# Count reads for each sample
reads <- rep(0, 4)
for (i in 1:4) {
  reads[i] <- nrow(pis[[i]]) + nrow(pias[[i]])
}

# Perform additional calculations and analysis...
# (The code continues with various calculations and analyses, 
# with comments indicating the main steps of the process)

# Combine sum and max information for smRNA data
sumnum <- cbind(
  apply(data1a, 1, sum),
  apply(data2a, 1, sum),
  apply(data3a, 1, sum),
  apply(data4, 1, sum),
  apply(data5, 1, sum),
  apply(data6, 1, sum)
)
maxnum <- cbind(
  apply(data1a, 1, max),
  apply(data2a, 1, max),
  apply(data3a, 1, max),
  apply(data4, 1, max),
  apply(data5, 1, max),
  apply(data6, 1, max)
)

# Calculate percentage differences and identify sites of interest
diff_percentage <- ((sumnum - maxnum) / sumnum * 100)
sites_of_interest <- SNPsites1[apply((sumnum - maxnum) > 10, 1, sum) > 0,]

# Further processing of SNPs and regions...
# (The code continues with various calculations and analyses, 
# with comments indicating the main steps of the process)

# Assign regions based on SNPs for different samples
region5s <- matrix(0, length(SNPsites), 4)
centers <- matrix(0, length(SNPsites), 4)
region3s <- matrix(0, length(SNPsites), 4)
region5as <- matrix(0, length(SNPsites), 4)
centeras <- matrix(0, length(SNPsites), 4)
region3as <- matrix(0, length(SNPsites), 4)

for (j in 1:4) {
  for (i in 1:length(SNPsites)) {
    region5s[i, j] <- sum(SNPsites[i] - 9 <= pis[[j]][, 2] & pis[[j]][, 2] <= SNPsites[i] - 2)
    region5as[i, j] <- sum(SNPsites[i] + 1 <= pias[[j]][, 3] - 1 & pias[[j]][, 3] - 1 <= SNPsites[i] + 8)
    centers[i, j] <- sum(SNPsites[i] - 18 <= pis[[j]][, 2] & pis[[j]][, 2] <= SNPsites[i] - 11)
    centeras[i, j] <- sum(SNPsites[i] + 10 <= pias[[j]][, 3] - 1 & pias[[j]][, 3] - 1 <= SNPsites[i] + 17)
    region3s[i, j] <- sum(SNPsites[i] - 26 <= pis[[j]][, 2] & pis[[j]][, 2] <= SNPsites[i] - 19)
    region3as[i, j] <- sum(SNPsites[i] + 18 <= pias[[j]][, 3] - 1 & pias[[j]][, 3] - 1 <= SNPsites[i] + 25)
  }
}

# Further analysis and calculations...
# (The code continues with various calculations and analyses, 
# with comments indicating the main steps of the process)

# Perform additional calculations and create summary matrices
snpbox <- rbind(
  apply(region5s, 2, sum) + apply(region5as, 2, sum),
  apply(centers, 2, sum) + apply(centeras, 2, sum),
  apply(region3s, 2, sum) + apply(region3as, 2, sum)
)

# Further analysis and calculations...
# (The code continues with various calculations and analyses, 
# with comments indicating the main steps of the process)

# Plotting the results
pdf("231030_P-element_around652_plot.pdf")
par(mfrow = c(2, 2))
for (i in 1:4) {
  barplot(snpbox[, i], xlab = "", ylab = "", main = MAIN[i], ylim = YLIM)
  par(new = TRUE)
  barplot(-snpbox[, i], xlab = "P-element 622:712", ylab = "% in this region", main = MAIN[i], ylim = YLIM)
}
dev.off()

# Additional analysis using Biostrings library...
# (The code continues with the analysis using the Biostrings library, 
# with comments indicating the main steps of the process)

# Further analysis and calculations...
# (The code continues with various calculations and analyses, 
# with comments indicating the main steps of the process)

# Additional analysis using Biostrings library...
library(Biostrings)
S5 <- read.table("SRR333505_Pelement_bt_m1.bed")
S6 <- read.table("SRR333506_Pelement_bt_m1.bed")
S7 <- read.table("SRR333507_Pelement_bt_m1.bed")
S9 <- read.table("SRR333509_Pelement_bt_m1.bed")
S5fas <- readDNAStringSet("/home/kshoji/Gawain/NGSdata/sra/SRR333505_trim.fasta.gz")
S6fas <- readDNAStringSet("/home/kshoji/Gawain/NGSdata/sra/SRR333506_trim.fasta.gz")
S7fas <- readDNAStringSet("/home/kshoji/Gawain/NGSdata/sra/SRR333507_trim.fasta.gz")
S9fas <- readDNAStringSet("/home/kshoji/Gawain/NGSdata/sra/SRR333509_trim.fasta.gz")
S5fasnames <- sub(" .*", "", names(S5fas))
S6fasnames <- sub(" .*", "", names(S6fas))
S7fasnames <- sub(" .*", "", names(S7fas))
S9fasnames <- sub(" .*", "", names(S9fas))

# Further analysis and calculations...
# (The code continues with various calculations and analyses, 
# with comments indicating the main steps of the process)

# Analysis of RNAseq data...
sumnum <- c(1588335, 855750, 609023, 4214570)
S5fas[is.element(S5fasnames, S5[, 4])]
S6fas[is.element(S6fasnames, S6[, 4])]

# Create plots based on RNAseq data...
pdf("231025_Pelement_dis_recover_pattern_plot.pdf")
for (i in 1:4) {
  lib <- hlist[[i]]
  tabp <- table(factor(lib[lib[, 6] == "+", 2] + 1, levels = 1:len)) / sumnum[i] * 1000000
  tabm <- table(factor(lib[lib[, 6] == "-", 3], levels = 1:len)) / sumnum[i] * 1000000
  plot(tabp, ylab = "", ylim = YLIM, xlim = XLIM, xlab = "P-element", main = hs[i], type = "l")
  par(new = TRUE)
  plot(-tabm, ylab = "", ylim = YLIM, xlim = XLIM, xlab = "", type = "l")
}
dev.off()

# Additional analysis...
i <- 2
lib <- hlist[[i]]
tabp2 <- table(factor(lib[lib[, 6] == "+", 2] + 1, levels = 1:len))
tabm2 <- table(factor(lib[lib[, 6] == "-", 3], levels = 1:len))

# Further analysis and calculations...
# (The code continues with various calculations and analyses, 
# with comments indicating the main steps of the process)

# Analysis of consensus matrix...
consensusMatrix(S6fas2)[1:4,]
consensusMatrix(S7fas2)[1:4,]
consensusMatrix(S9fas2)[1:4,]

# Further analysis and calculations...
# (The code continues with various calculations and analyses, 
# with comments indicating the main steps of the process)

# More calculations...
num2 <- rep(0, 50)
num3 <- rep(0, 50)
num4 <- rep(0, 50)
for (i in 2:50) {
  num2[i] <- sum(tabp2[i:2907] * tabp2[1:(2907 - i + 1)]) + sum(tabm2[i:2907] * tabm2[1:(2907 - i + 1)])
  num3[i] <- sum(tabp3[i:2907] * tabp3[1:(2907 - i + 1)]) + sum(tabm3[i:2907] * tabm3[1:(2907 - i + 1)])
  num4[i] <- sum(tabp4[i:2907] * tabp4[1:(2907 - i + 1)]) + sum(tabm4[i:2907] * tabm4[1:(2907 - i + 1)])
}
mats <- cbind(num2, num3, num4)
sweep(mats, 2, apply(mats, 2, sum), "/") * 100
```

### **code** barplot for specific regions

```{r}
#### load files ####
dat.pelement <-
  read.table("./txt/Pelement.txt", sep = "\t", header = TRUE)

#### format the code ####
dat.pelement$Regions <- gsub ("region", "'", dat.pelement$Regions)
dat.pelement$Samples <- factor(dat.pelement$Samples,
                               levels = c("Har x w1", "w1 x Har"))
dat.pelement$Regions <-  factor(dat.pelement$Regions,
                                levels = c("5'", "Middle", "3'", "NR"))
dat.pelement[dat.pelement$direc == "antisense", "Percentage"] <-
  dat.pelement[dat.pelement$direc == "antisense", "Percentage"] * -1
dat.pelement[dat.pelement$direc == "antisense", "Sums"] <-
  dat.pelement[dat.pelement$direc == "antisense", "Sums"] * -1


#### add to rdata ####
save(list = ls(pattern = "dat"), file = "data.RData")
```

### **fig** barplot for specific regions

```{r}
p_day <- unique(dat.pelement$Days)

for (i in p_day) {
  fig <- ggplot(subset(dat.pelement, Days == i & !is.na(Regions))) +
    # geom_hline(yintercept= 0)+
    # geom_vline(xintercept = c(680,708),color="blue")+
    geom_vline(xintercept = c(660, 663, 666, 671, 679, 686, 687, 688),
               color = "pink") +
    geom_col(aes(x = Site, y = Percentage, fill =Regions),
             col = "black",
             size = 0.35) +
    # ggtitle("")+
    xlab("") +
    ylab("Relative abundance (%)") +
    ylim(-100, 100) +
    facet_grid(. ~ Samples) +
    scale_fill_manual(values = c("black", "#666666","#cdcdcd","white","white")) +
    scale_x_continuous(limits = c(657, 717), expand = c(0.1, 0.1)) +
    guides(fill = "none") +
    mytheme()
  print(fig)
  ggsave(
    fig,
    filename = sprintf("./fig/Pelement_%s.pdf", i),
    width = Width * 1.7,
    height = Height
  )
}
```

### **code** barplot for entire region

```{r}
#### load files ####
dat.pelement2 <-
  read.table("./txt/Pelement2.txt", sep = "\t", header = TRUE)

#### format the code ####
dat.pelement2$Regions <- gsub("region", "'", dat.pelement2$Regions)
dat.pelement2$Regions <-
  factor(dat.pelement2$Regions, levels = c("5'", "Middle", "3'"))
dat.pelement2$Samples <- factor(dat.pelement2$Samples,
                                levels = c("Har x w1", "w1 x Har"))
dat.pelement2$Percentage <- dat.pelement2$Percentage * 100
dat.pelement2$percent <- dat.pelement2$percent * 100

#### add to rdata ####
save(list = ls(pattern = "dat"), file = "data.RData")
```

### **fig** barplot for entire region

```{r}
p_day <- unique(dat.pelement2$Days)

for (i in p_day) {
  fig <- ggplot(subset(dat.pelement2, Days == i)) +
    geom_bar(aes(x = Samples, y = percent, fill = Regions),
             stat = "summary",
             wide = 0.5) +
    ggtitle("") +
    xlab("") +
    ylab("Percentage of\nmismatches") +
    # # ylim(-1,1)+
    scale_fill_manual(values = c("black", "#666666", "#cccccc")) +
    scale_x_discrete(expand = c(0.3, 0.3)) +
    guides(fill = "none") +
    mytheme()
  print(fig)
  ggsave(
    fig,
    filename = sprintf("./fig/Pelement_%s_summary.pdf", i),
    width = Width * 0.6,
    height = Height * 1.07
  )
}
```

# Backup

## RNA vs Dscores

### \<code\>

```{R}
#### save ####
rpkm_long <- readRDS("./results/piRNAandRNA_rpkm.rds")
rpkm_long$year <- rep(c("naive13", "naive20"),nrow(rpkm_long)/2)
dat.rpkms <- rpkm_long
save(list = ls(pattern = "dat"), file = "data.RData")
```

### \| fig \|

```{r}
lapply(c("naive13", "naive20"), function(i) {
  ggplot(subset(dat.rpkms, year == i)) +
    geom_point(aes(x = rna_rpkm, y = scores, color = scores)) +
    scale_x_log10() +
    ggtitle(i) +
    xlab(expression(RPKM[transcripts])) +
    ylab("Dscores") +
    scale_color_gradientn(colors = c( "blue","red","#FFD700"),
                          # limits=c(0, 100),
                          guide = "none") +
    mytheme()
})
```

## GC rich

### <code>

```{r eval=FALSE, include=FALSE}
#### read ####
trans_used <- readRDS("./check/transposon_uesed.rds")
basicgroups <- readRDS("./results/basicgroups.rds")

#### bases ####
bases <- c("AA","AT","TA","TT", #AT_AT
           "AC","AG","TC","TG", #AT_CG
           "CA","CT","GA","GT", #CG_AT
           "CC","CG","GC","GG") #CG_CG
bases_set <- rep(c("AT_AT", "AT_CG", "CG_AT ","CG_CG"),each = 4)
names(bases_set) <- bases

#### trans split ####
trans_used_split <- sapply(1:max(width(trans_used)),
                      function(x) {
                        message(x)
                        substr(trans_used, x, x + 1)
                      })

#### calculate ####
trans_used_table <- t(sapply(rownames(trans_used_split),
                      function(x){
                        table(factor(trans_used_split[x,],levels=bases))
                      }))

trans_used_bases <- sweep(trans_used_table,1,apply(trans_used_table,1,sum),"/")*100
trans_used_bases <-  as.data.frame(trans_used_bases)
trans_used_bases$trans <- rownames(trans_used_bases)

#### group ####
trans_used_bases$group <- basicgroups[trans_used_bases$trans]


trans_bases_long <- pivot_longer(trans_used_bases,cols=1:16,names_to = "bases", values_to = "percent")
trans_bases_long$bases_set <- bases_set[trans_bases_long$bases]

#### save ####
dat.transbases <-  trans_bases_long
save(list = ls(pattern = "dat"), file = "data.RData")

```

### \| fig \| transposon_group_bases_percent

```{r echo=FALSE}
ggplot(dat.transbases) +
  geom_bar(stat = "identity",
           position = "dodge",
           aes(x = bases_set, y = percent, fill = group)) +
  scale_fill_manual(values = c("blue", "#800080", "red", "#FF7A00", "#FFD700")) +
  ggtitle("transposon")+
  mytheme()
ggsave("./fig/transposon_group_bases_percent.pdf",
       width = Width *1.3 , height = Height * 1.3)
```

## piRNA seed GC

### \<code\> dat.pipercent

```{r}
#### samples ####
sn <- "naive13"

#### read ####
basicgroups <- readRDS("./results/basicgroups.rds")
piRNA_plus <- readRDS(sprintf("./piRNA/%s_transposon_piRNAplus.rds",sn))
piRNA_minus <- readRDS(sprintf("./piRNA/%s_transposon_piRNAminus.rds",sn))
fas <- readDNAStringSet(sprintf("./fasta/%s_trim.fasta.gz",sn))                    
#### group ####
bin_names <- c("bin1", "bin2", "bin3", "bin4", "bin5")
names(bin_names) <- bin_names

piRNA_total <- rbind(as.data.frame(piRNA_plus), as.data.frame(piRNA_minus))
piRNA_total$group <- basicgroups[piRNA_total$V1]
piRNA_group <- piRNA_total[!is.na(piRNA_total$group),]
                      
groups_names <- lapply(bin_names,function(bn){
  unique(as.character(piRNA_group[as.character(piRNA_group$group) == bn,"V4"]))
})
seq_group <- lapply(bin_names,function(bn){
  fas[groups_names[[bn]]]
})



#### bases  table ####
bases <- c("AA","AT","TA","TT", #AT_AT
           "AC","AG","TC","TG", #AT_CG
           "CA","CT","GA","GT", #CG_AT
           "CC","CG","GC","GG") #CG_CG
bases_set <- rep(c("AT_AT", "AT_CG", "CG_AT ","CG_CG"),each = 4)
names(bases_set) <- bases

piRNA_split <- lapply(bin_names,function(bn){
  message(bn)
  sapply(1:25,function(site) {
                        message(site)
                        bases_split <- substr(seq_group[[bn]], site, site + 1)
                        factor(bases_split, levels = bases)
                      })
  })
piRNA_table <- lapply(bin_names,function(bn){
  message(bn)
  bn_table <- t(sapply(c(1:25),function(nn){
    message(nn)
    nn_split <- factor(piRNA_split[[bn]][,nn],levels = bases)
    table(nn_split)
  }))
  sweep(bn_table,1,apply(bn_table[],1,sum),"/")*100
})

#### means ####
seed_mean <- as.data.frame(sapply(piRNA_table, function(i) {
  apply(i[2:8, ], 2, mean)
}))
seed_mean$bases <- rownames(seed_mean)
seed_mean$region <- factor("seed",levels = c("seed","nonseed"))

nonseed_mean <- as.data.frame(sapply(piRNA_table, function(i) {
  apply(i[11:25,], 2, mean)
}))
nonseed_mean$bases <- rownames(nonseed_mean)
nonseed_mean$region <- factor("nonseed",levels = c("seed","nonseed"))


seed_long <- pivot_longer(seed_mean,col=1:5, names_to = "groups",values_to = "percent")
nonseed_long <- pivot_longer(nonseed_mean,col=1:5, names_to = "groups",values_to = "percent")

piRNA_percent <- rbind(seed_long,nonseed_long)
piRNA_percent$sets <- factor(rep(c("AT_TA","AT_CG","CG_AT","CG_CG"),
                                 each = 20,2),
                             levels = c("AT_TA","AT_CG","CG_AT","CG_CG"))

#### save ####
saveRDS(piRNA_split,"./check/piRNA_split.rds")
saveRDS(piRNA_percent,"./results/piRNA_percent.rds")
```

```{R}
piRNA_percent <-  readRDS("./results/piRNA_percent.rds")
#### save ####
dat.pipercent <- piRNA_percent
save(list = ls(pattern = "dat"), file = "data.RData")
```

### \| fig \| piRNA_group_bases_percent\_%s

```{r}
sn <- "naive13"
ggplot(dat.pipercent)+
  geom_col(aes(x=bases,y=percent,fill=groups),
           position = "dodge")+
  ggtitle("piRNA")+
  xlab("")+
  facet_grid(region~sets,scales = "free")+
  scale_fill_manual(values = c("blue", "#800080","red","#FF7A00","#FFD700")) +
  mytheme()
ggsave(sprintf("./fig/piRNA_group_bases_percent_%s.pdf",sn),
       width = Width * 2 , height = Height * 1.7)
```

## mismatches

### mismatches.R

```{r}
#### setting ####
setwd("/home/Jie/Jupiter/paper/v3")
pkgs <- c("Biostrings","tidyr")
sapply(pkgs, library, character.only = TRUE)
trans <- readDNAStringSet("./fasta/transposon.fasta")

#### samples ####
args <- commandArgs(trailingOnly = TRUE)
years <- args[1]
direction <- args[2]
subfile <- args[3]
switch(direction,
       sense={
         b <- "plus"
         c <- "+"
         d <- "V2"
         e <- 2
       },
       anti={
         b <- "minus"
         c <- "-"
         d <- "V3"
         e <- 3
       })

#### load ####
bed <- readRDS(sprintf("./rds/naive%s_transposon.rds",years))
message("readed bed file")
fas <- readDNAStringSet(sprintf("./fasta/naive%s_trim.fasta.gz",years))
message("readed fasta file")

#### piRNA ####
bed$len <- bed[, 3] - bed[, 2]
bed[, 1] <- factor(bed[, 1], levels = names(trans))
bed[, 2] <- factor(bed[, 2] + 1 , levels = 1:max(width(trans)))
bed[, 3] <- factor(bed[, 3] , levels = 1:max(width(trans)))
piRNA <- subset(bed, 26 <= len & len <= 32 & V6 == c)
saveRDS(piRNA,sprintf("BmQ%s_%s.rds",years,b))
rm("bed")
saveRDS(piRNA,"./piRNA/naive%s_transposon_piRNAs.rds")
message("extracted piRNA")

#### sequences ####
piRNA.names <- as.character(unique(piRNA$V4))
fas.piRNA <- fas[piRNA.names]
rm("fas")
message("sequence")

#### distribution ####
piRNA.reads <- table(piRNA[,c(1,e)])
piRNA.reads["TE1_bm_1645_LINE/R4",] <- 0
message("removed 1645")

#### filter ####
peaks <- which(piRNA.reads>3,arr.ind=TRUE)
n <- ifelse(nrow(peaks)%% 2 == 0,nrow(peaks)/2 , (nrow(peaks)+1)/2)
m <- switch(subfile,"1"= c(1:n),"2"=c((n+1):nrow(peaks)))
message("selected 3 reads")

#### mismatches ratio ####
timestart<-Sys.time()
misratio <- matrix(0,nrow(peaks),27)
message(timestart)
for(i in m)
  {
  pinames <- piRNA[as.numeric(piRNA$V1) == peaks[i, 1] & piRNA[,d] == peaks[i, 2],4]
  misratio[i,1] <- piRNA.reads[peaks[i,1],peaks[i,2]]
  misratio[i,2:27] <- apply(consensusMatrix(fas.piRNA[pinames], as.prob = T)[1:4, 1:26], 2, max)
  message(i,"/",nrow(peaks))
  }
timeend<-Sys.time()
runningtime<-timeend-timestart
print(runningtime)

#### save ####
saveRDS(misratio,sprintf("./results/naive%s_%s_mis_%s.rds",years,direction,subfile))
saveRDS(peaks,sprintf("./results/naive%s_%s_peak.rds",years,direction))
message("finished ",years," ",direction)
```

```{bash}
nohup Rscript /home/Jie/Jupiter/paper/funs/mismatches.R  13 sense 1 > 0807_13sense1.out &
  nohup Rscript /home/Jie/Jupiter/paper/funs/mismatches.R  20 sense 1 > 0807_20sense1.out &
  nohup Rscript /home/Jie/Jupiter/paper/funs/mismatches.R  13 anti 1 > 0807_13anti1.out &
  nohup Rscript /home/Jie/Jupiter/paper/funs/mismatches.R  20 anti 1 > 0807_20anti1.out &
  nohup Rscript /home/Jie/Jupiter/paper/funs/mismatches.R  13 sense 2 > 0807_13sense2.out &
  nohup Rscript /home/Jie/Jupiter/paper/funs/mismatches.R  20 sense 2 > 0807_20sense2.out &
  nohup Rscript /home/Jie/Jupiter/paper/funs/mismatches.R  13 anti 2 > 0807_13anti2.out &
  nohup Rscript /home/Jie/Jupiter/paper/funs/mismatches.R  20 anti 2 > 0807_20anti2.out
```

### \<code\> dat.mismean

```{r}
#### read ####
basicgroups <- readRDS("./results/basicgroups.rds")
# basicgroups <- readRDS("/Users/yjjuly/Documents/Doctor/paper/v1/results/basicgroups.rds") #### v1's groups

#### assign ####
group_names <-  c("bin1", "bin2", "bin3", "bin4", "bin5")
names(group_names) <- group_names
readmin <- 100
readmax <- 10^6

#### mis reads ####
mis_mean <- lapply(c("naive13", "naive20"), function(years) {
  #### read
  reads <- lapply(c("sense", "anti"), function(direc) {
    reads_1 <-
      readRDS(sprintf("./results/%s_%s_mis_1.rds", years, direc))
    reads_2 <-
      readRDS(sprintf("./results/%s_%s_mis_2.rds", years, direc))
    reads_2[1:nrow(reads_1), ] <-  reads_1
  })
  mismatches <- rbind(as.data.frame(reads[[1]]), as.data.frame(reads[[2]]))
  colnames(mismatches) <- c("reads",1:26)
  #### peak
  plus <-
    readRDS(sprintf("./results/%s_sense_peak.rds", years))
  minus <-
    readRDS(sprintf("./results/%s_anti_peak.rds", years))
  peak <- names(trans[c(plus[, 1], minus[, 1])])
  #### group
  mismatches[, "groups"] <- as.character(basicgroups[peak])
  mis.mean <- sapply(group_names,
                     function(group_name) {
                       mis_group <- subset(mismatches, groups == group_name)
                       group_mean <-
                         apply(mis_group[mis_group[, "reads"] > readmin &
                                           mis_group[, "reads"] < readmax, 2:27], 2, mean)
                       group_mean
                     })
  mis.mean <- as.data.frame(mis.mean)
  #### info
  mis.mean$year <- years
  mis.mean$site <- 1:26
  mis.mean$region <-
    factor(c("1U", rep("seed", 8), "10A", rep("nonseed", 16)),
           levels = c("1U", "seed", "10A", "nonseed"))
  mis.mean$set <-
    c("1U", rep("2-9", 8), "10A", rep("11-18", 8), rep("19-26", 8))
  #### long
  mis.mean.long <-
    pivot_longer(
      mis.mean,
      cols = 1:5,
      names_to = "groups",
      values_to = "percent"
    )
  mis.mean.long
})

#### save ####
dat.mismean <-  Reduce(rbind,mis_mean,data.frame())
save(list = ls(pattern = "dat"), file = "data.RData")
```

### \| fig \| mismatches

```{r}
ggplot(subset(dat.mismean,region %in% c("seed","nonseed")),
       aes(x=region,y=1-percent,fill = groups))+
  geom_bar(stat = "summary",
           fun = "mean",
           position = "dodge",
           width=0.6)+
  facet_grid(year~.,scales ="free")+
  # geom_signif(test = "t.test",
  #             comparisons = list(c("insert", "nc")),
  #             map_signif_level = FALSE)+
  ggtitle(sprintf("%s < reads < %s",readmin,readmax))+
  ylab("percent of mismatches")+
  scale_fill_manual(values = c("blue", "#800080","red","#FF7A00","#FFD700")) +
  mytheme()
ggsave(filename = sprintf("./fig/mismatches_%s_%s.pdf",readmin,readmax), width = Width*1.7, height = Height*1.8)
```

### \<code\> updown

```{r}
#### load ####
updown <- readRDS("./check/updown_names.rds")

#### names ####
group_names <-  c("bin1", "bin2", "bin3", "bin4", "bin5")
names(group_names) <- group_names
updown_name <- c("up","down")
names(updown_name) <- updown_name

#### bins ####
readmin <- 100
readmax <- 10^6

mis_mean <- lapply(c("naive13", "naive20"),
                   function(years) {
                     #### read
                     reads <- lapply(c("sense", "anti"),
                                     function(direc) {
                                       reads_1 <-
                                         readRDS(sprintf("./results/%s_%s_mis_1.rds", years, direc))
                                       reads_2 <-
                                         readRDS(sprintf("./results/%s_%s_mis_2.rds", years, direc))
                                       reads_2[1:nrow(reads_1),] <-  reads_1
                                     })
                     mismatches <-
                       rbind(as.data.frame(reads[[1]]), as.data.frame(reads[[2]]))
                     colnames(mismatches) <- c("reads", 1:26)
                     #### peak
                     plus <-
                       readRDS(sprintf("./results/%s_sense_peak.rds", years))
                     minus <-
                       readRDS(sprintf("./results/%s_anti_peak.rds", years))
                     peak <- names(trans[c(plus[, 1], minus[, 1])])
                     mismatches[, "groups"] <- as.character(basicgroups[peak])
                     mismatches[, "updown"] <- as.character(updown[peak])
                     mis.mean <-  sapply(group_names,
                                         function(gn) {
                                           sapply(updown_name,
                                                  function(ud) {
                                                    mis_group <- subset(mismatches, updown == ud & groups == gn)
                                                    group_mean <-
                                                      apply(mis_group[mis_group[, "reads"] > readmin &
                                                                        mis_group[, "reads"] < readmax, 2:27], 2, mean)
                                                    group_mean
                                                  })
                                         })
                     mis.mean <- as.data.frame(mis.mean)
                     #### info
                     mis.mean$year <- years
                     mis.mean$site <- 1:26
                     mis.mean$region <-
                       factor(c("1U", rep("seed", 8), "10A", rep("nonseed", 16)),
                              levels = c("1U", "seed", "10A", "nonseed"))
                     mis.mean$set <-
                       c("1U", rep("2-9", 8), "10A", rep("11-18", 8), rep("19-26", 8))
                     mis.mean$updown <- rep(factor(updown_name,levels=updown_name), each = 26)
                     #### long
                     mis.mean.long <-
                       pivot_longer(
                         mis.mean,
                         cols = 1:5,
                         names_to = "groups",
                         values_to = "percent"
                       )
                     mis.mean.long$percent = 1 - mis.mean.long$percent
                     mis.mean.long
                   })

mismean_minus <- mis_mean[[1]]
mismean_minus$percent <- mis_mean[[2]]$percent - mis_mean[[1]]$percent

mismean_updown <-  Reduce(rbind,mis_mean,data.frame())

```

### \|fig\| bin5

```{r}
ggplot(subset(mismean_updown,region %in% c("seed","nonseed")),
       aes(x=region,y=percent,fill = groups))+
  geom_bar(stat = "summary",
           fun = "mean",
           position = "dodge",
           width=0.6)+
  facet_grid(updown~year,scales ="free")+
  # geom_signif(test = "t.test",
  #             comparisons = list(c("insert", "nc")),
  #             map_signif_level = FALSE)+
  ggtitle(sprintf("%s < reads < %s",readmin,readmax))+
  ylab("percent of mismatches")+
  # scale_fill_manual(values = c("blue", "#800080","red","#FF7A00","#FFD700")) +
  mytheme()

ggplot(subset(mismean_minus,region %in% c("seed","nonseed")),
       aes(x=region,y=percent,fill = updown))+
  geom_bar(stat = "summary",
           fun = "mean",
           position = "dodge",
           width=0.6)+
  geom_hline(yintercept = 0)+
  facet_wrap(groups~.)+
  # geom_signif(test = "t.test",
  #             comparisons = list(c("insert", "nc")),
  #             map_signif_level = FALSE)+
  ggtitle(sprintf("%s < reads < %s",readmin,readmax))+
  ylab("percent of mismatches")+
  # scale_fill_manual(values = c("blue", "#800080","red","#FF7A00","#FFD700")) +
  mytheme()


```

### <code> results_long

```{r}
bin_names <- c("bin1", "bin2", "bin3", "bin4", "bin5")
names(bin_names) <- bin_names
sample_names <- c("naive13","naive20")
names(sample_names) <- sample_names
misnum <- c(1,2,3)

sn ="naive13"
results <- list()
for (bn in bin_names){
    message(bn)
    nn <- paste0(sn,"_",bn)
    results[[bn]] <- matrix(ncol = 3, nrow = 26,0)
      mis1 <- read.table(sprintf("./txt/%s_%s_v%s_%s.txt",sn,bn,1,1))
      mis1 <- mis1[!as.character(mis1[,1]) %in% bed_con_all[[sn]], ]
      mis1[,2] <-factor(mis1[,2],levels = 1:26)
      results[[bn]][,1] <- table(mis1[,2])
    
      mis2 <- read.table(sprintf("./txt/%s_%s_v%s_%s.txt",sn,bn,2,2))
      mis2 <- mis2[!as.character(mis2[,1]) %in% bed_con_all[[sn]], ]
      mis2[,2] <-factor(mis2[,2],levels = 1:26)
      mis2[,3] <-factor(mis2[,3],levels = 1:26)
      results[[bn]][,2] <- table(mis2[,2])+table(mis2[,3])
      
      mis3 <- read.table(sprintf("./txt/%s_%s_v%s_%s.txt",sn,bn,3,3))
      mis3 <- mis3[!as.character(mis3[,1]) %in% bed_con_all[[sn]], ]
      mis3[,2] <-factor(mis3[,2],levels = 1:26)
      mis3[,3] <-factor(mis3[,3],levels = 1:26)
      mis3[,4] <-factor(mis3[,4],levels = 1:26)
      results[[bn]][,3] <- table(mis3[,2])+table(mis3[,3])+table(mis3[,4])
      
      results[[bn]] <- as.data.frame(results[[bn]])
  }
result13_percent <- list()
for(bn in bin_names){
      total.reads <-  readRDS(sprintf("./results/%s_%s_v3_total_reads.rds",sn,bn))
      # results[[bn]] <- results[[bn]]/total.reads
       result13_percent[[bn]] <- (apply(results[[bn]],2,sum)/total.reads)/c(1,2,3)*100

}
result13 <- Reduce(rbind,result13_percent,data.frame())
colnames(result13) <- c("1","2","3")
result13$bins <- bin_names
result13$samples <- "naive13"
result13_long <- pivot_longer(result13,col=1:3,values_to = "percent",names_to = "mismatches")


sn ="naive20"
for (bn in bin_names){
    message(bn)
    nn <- paste0(sn,"_",bn)
    results[[bn]] <- matrix(ncol = 3, nrow = 26,0)
      mis1 <- read.table(sprintf("./txt/%s_%s_v%s_%s.txt",sn,bn,1,1))
      mis1 <- mis1[!as.character(mis1[,1]) %in% bed_con_all[[sn]], ]
      mis1[,2] <-factor(mis1[,2],levels = 1:26)
      results[[bn]][,1] <- table(mis1[,2])
    
      mis2 <- read.table(sprintf("./txt/%s_%s_v%s_%s.txt",sn,bn,2,2))
      mis2 <- mis2[!as.character(mis2[,1]) %in% bed_con_all[[sn]], ]
      mis2[,2] <-factor(mis2[,2],levels = 1:26)
      mis2[,3] <-factor(mis2[,3],levels = 1:26)
      results[[bn]][,2] <- table(mis2[,2])+table(mis2[,3])
      
      mis3 <- read.table(sprintf("./txt/%s_%s_v%s_%s.txt",sn,bn,3,3))
      mis3 <- mis3[!as.character(mis3[,1]) %in% bed_con_all[[sn]], ]
      mis3[,2] <-factor(mis3[,2],levels = 1:26)
      mis3[,3] <-factor(mis3[,3],levels = 1:26)
      mis3[,4] <-factor(mis3[,4],levels = 1:26)
      results[[bn]][,3] <- table(mis3[,2])+table(mis3[,3])+table(mis3[,4])
      
      results[[bn]] <- as.data.frame(results[[bn]])
  }
result20_percent <- list()
for(bn in bin_names){
      total.reads <-  readRDS(sprintf("./results/%s_%s_v3_total_reads.rds",sn,bn))
      # results[[bn]] <- results[[bn]]/total.reads
       result20_percent[[bn]] <- (apply(results[[bn]],2,sum)/total.reads)/c(1,2,3)*100
}
result20 <- Reduce(rbind,result20_percent,data.frame())
colnames(result20) <- c("1","2","3")
result20$bins <- bin_names
result20$samples <- "naive20"
result20_long <- pivot_longer(result20,col=1:3,values_to = "percent",names_to = "mismatches")

result_long <- rbind(result13_long,result20_long)

saveRDS(result_long, "./results/snps_bins_percent.rds")
```

### <code> dat.snps_per

```{r}
snps_per <- readRDS("./results/snps_bins_percent.rds")
dat.snps_per <- snps_per 
save(list = ls(pattern = "dat"), file = "data.RData")
```

### [fig](#fig-1) snps_percent

```{r}
ggplot(snps_per)+
  geom_col(aes(x=bins,y=percent,fill = mismatches),col="black")+
  scale_fill_manual(values = c("white", "grey","black")) +
  facet_grid(samples ~ .)+
  mytheme()
ggsave(filename = "./fig/snps_percent.pdf", width = Width*1.3, height = Height*1.7)
```

## RNA snps

### consensus_sequence_mRNA.R

```{r tmp }
#### setting ####
setwd("/home/Jie/Jupiter/paper/v3")
pkgs <- c("Biostrings","tidyr","DECIPHER")
sapply(pkgs, library, character.only = TRUE)

trans <- readDNAStringSet("./fasta/transposon.fasta")
############## tmp for 
sample_names <- c("close13", "close20_sampling_100nt")
read_number <- c("1", "2")
directions <- c("sense", "anti")

for (sn in sample_names){
  for (rn in read_number){
    for ( direc in directions){
      if (sn == "close13"){
        
      }
      else{
 
# args <- commandArgs(trailingOnly = TRUE)
# args <- c(close20,1,"sense")
# sn <- args[1]
# rn <- args[2]
# direc <- args[3]

#### load ####
mRNA <- readRDS(sprintf("./rds/%s_R%s_%s_RNA_bed.rds", sn, rn, direc))
message("naive", sn, " ", rn, " ", direc, " readed bed file")

fas <- readDNAStringSet(sprintf("./fastq/%s_R%s.fastq.gz", sn, rn), format = "fastq")
message("naive", sn, " ", rn, " readed fasta file")

#### fas names ####
if (sn =="close20_sampling_100nt")  {
  seq_names <- names(fas)
  seq_names_processed <- gsub(" .*", "", seq_names)
  seq_names_final <- paste0(seq_names_processed,"/",rn)
  names(fas) <- seq_names_final
}

#### sense or anti ####
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

#### table ####
mRNA$V4 <- as.character(mRNA$V4)
mRNA$trans <- as.numeric(mRNA$V1)
message( sn, " ", rn, " ", direc," factor")

mRNAnames <- mRNA$V4
fas.mRNA <- fas[mRNAnames]
rm("fas")

mRNA.reads <- table(mRNA[, c(1, e)])
mRNA.reads["TE1_bm_1645_LINE/R4", ] <- 0
message("removed",sn, " ", rn, " ", direc)

#### reads peak ####
peaks <- which(mRNA.reads > 3, arr.ind = TRUE)
message("peaks", sn, " ", rn," ", direc, " peaks : ", nrow(peaks))

#### consensus sequence ####
conseq <- character(nrow(peaks))
timestart <- Sys.time()
message("started", sn, " ", rn," ", direc, " start time : ", timestart)
conseq <- apply(peaks, 1, function(pb) {
  tt <- pb[1]
  dd <- pb[2]
  mask <- mRNA$trans == tt & mRNA[, d] == dd
  mnames <- mRNA$V4[mask]
  fas.m <- fas.mRNA[mnames]
  conseq_tmp <-
    as.character(ConsensusSequence(fas.m, threshold = 0.5))
  if (grepl("[^ATCG]", conseq_tmp)) {
    cat("=")
    ast <-  table(fas.m)
    max_seq <- names(ast[ast == max(ast)])
    ifelse(length(max_seq) > 1,
           Reduce(function(x, y) {
             paste0(x, strrep("N", 24), y)
           }, max_seq),
           max_seq)
  } else {
    cat("-")
    conseq_tmp
  }
})
timeend <- Sys.time()
message("ended", sn, " ",rn, " ", direc, "end time : ", timeend)
runningtime <- timeend - timestart
print(runningtime)

conseq <- DNAStringSet(conseq)
names(conseq) <- paste(sn,rn,direc,peaks[, 1], peaks[, 2], sep = "_")

#### save ####
writeXStringSet(conseq, sprintf("./check/%s_R%s_%s_consensus_sequence.fasta", sn, rn, direc))
message( sn," ",rn, " ", direc, " saved all results")
   }}
  }
}
```

### run_0129.sh

```{bash}
echo "start preprocessing"
bash ./funs/preprocessing.sh  > ./runlog/preprocessing.out 

echo "start preprocessing_aftermapping"
Rscript ./funs/preprocessing_aftermapping.R  > ./runlog/preprocessing_aftermapping.out 

MM="close20_sampling_100nt close13"
NN="1 2"
DD="sense anti"
for M in ${MM};do
 for B in "${NN}"; do
    for D in "${DD}"; do
    # consensus seq
      echo ${M}${N}${D}
      nohup Rscript ./funs/consensus_sequence_mRNA.R ${M} ${N} ${D} > ./runlog/0129_${M}${N}${D}.out &
      pids+=($!)
    done
     for pid in "${pids[@]}"; do
    wait "$pid"
    done
  done
done

for M in ${MM}; do
  echo "bind datas"
  cat ./check/${M}*consensus_sequence.fasta  > ./fasta/${M}_consensus_sequence.fasta
  
  echo "build"
  hisat2-build -f ./fasta/${M}_consensus_sequence.fasta ${M}
  
  #mapping
  echo "mapping"
  hisat2 -x ${M} -1 ./fastq/${M}_R1.fastq.gz  -2 ./fastq/${M}_R2.fastq.gz -k 3 -p 10 -S ./sambam/${M}_consensus_RNA_HS2.sam
  samtools view -@ 10 -bS ./sambam/${M}_consensus_RNA_HS2.sam > ./sambam/${M}_consensus_RNA_HS2.bam
  for M in ${MM};do
  bamToBed -i ./sambam/${M}_consensus_RNA_HS2.bam > ./bed/${M}_consensus_RNA_HS2.bed
  ## position_info
#  time samtools view ./sambam/${M}_consensus_RNA_HS2.bam | awk 'BEGIN{
#    print "chromosome", "start", "mismatch", "type", "position"
#    }match($0, /XM:i:([0-9]+)/, m) && m[1] > 0 && $3 != "*"{ 
#    if (!and($2, 16) && $4 == 1) {                             # 检查第 5 位是否为 0 为正向链 ; 与运算：and & ，例$2 & 16  ## 64 Read1  ### 仅使用匹配到位置1的结果
#      for (i = 1; i <= 10; i++) {
#        match($0, /XM:i:([0-9]+)/, m) && m[1] >= i;
#        match($0, /MD:Z:([0-9A-Z^:]+)/, md);
#        split(md[1], pos, /[ATCG]/);
#        split($3,csid,"_");
#        sum = i - 1 ;
#        for (j = 1; j <= i; j++)
#            sum += pos[j];
#          ## samples_names, direc,  trans_num,  trans_pos,  Read, SNP_num, SNP_pos, SNP_site
#            print csid[1],  csid[3],  csid[4],  csid[5],  csid[2],  m[1], md[1],  sum,  csid[5]+sum -i
#        } 
#    }
#}' > ./results/${M}_consensus_RNA_HS2.txt 

## 2013
  time samtools view ./sambam/close13_consensus_RNA_HS2.bam |awk 'BEGIN{
             ## 1             2           3              4          5          6              7          8
    print "samples_names", "direc",  "trans_num",  "trans_pos",  "Read",  "reads_names",  "SNP_num", "SNP_pos"
    } $3 != "*" && !and($2, 16) && $4 == 1 { # 检查第 5 位是否为 0 为正向链 ; 与运算：and & ，例$2 & 16  ## 64 Read1  ### 仅使用匹配到位置1的结果
        match($0, /XM:i:([0-9]+)/, m);
        match($0, /MD:Z:([0-9A-Z^:]+)/, md);
        split(md[1], pos, /[ATCG]/);
        split($3,csid,"_");
        # processing mark
        if (m[1] > 0){
        printf("_") > "/dev/stderr";
        ## samples_names, direc,  trans_num,  trans_pos,  Read, reads_name, SNP_num, SNP_pos
        print csid[1],  csid[3],  csid[4],  csid[5],  csid[2],  $1, m[1], md[1]
        }else{
        printf("*") > "/dev/stderr";
        ## samples_names, direc,  trans_num,  trans_pos,  Read, reads_name, SNP_num, SNP_pos
        print csid[1],  csid[3],  csid[4],  csid[5],  csid[2],  $1, 0, "-"
        }
    }' > ./results/close13_consensus_RNA_HS2.txt 

## 2020
  time samtools view ./sambam/close20_sampling_100nt_consensus_RNA_HS2.bam |awk 'BEGIN{
             ## 1             2           3              4          5          6              7          8
    print "samples_names", "direc",  "trans_num",  "trans_pos",  "Read",  "reads_names",  "SNP_num", "SNP_pos"
    } $3 != "*" && !and($2, 16) && $4 == 1 { # 检查第 5 位是否为 0 为正向链 ; 与运算：and & ，例$2 & 16  ## 64 Read1  ### 仅使用匹配到位置1的结果
        match($0, /XM:i:([0-9]+)/, m);
        match($0, /MD:Z:([0-9A-Z^:]+)/, md);
        split(md[1], pos, /[ATCG]/);
        split($3,csid,"_");
        # processing mark
        if (m[1] > 0){
        printf("_") > "/dev/stderr";
        ## samples_names, direc,  trans_num,  trans_pos,  Read, reads_name, SNP_num, SNP_pos
        print csid[1],  csid[5],  csid[6],  csid[7],  csid[4],  $1, m[1], md[1]
        }else{
        printf("*") > "/dev/stderr";
        ## samples_names, direc,  trans_num,  trans_pos,  Read, reads_name, SNP_num, SNP_pos
        print csid[1],  csid[5],  csid[6],  csid[7],  csid[4],  $1, 0, "-"
        }
    }' > ./results/close20_sampling_100nt_consensus_RNA_HS2.txt 

```

### <code> run

```{bash}
####
cd /home/Jie/Jupiter/paper/v3
nohup bash ./funs/run_0129.sh > ./runlog/script_output_0129.log 2>&1 &
tailf ./runlog/script_output_0129.log
```

### <code> positions_analysis

```{r}
#### define args ####
sample_names <- c("close13", "close20_sampling_100nt")
names(sample_names) <- sample_names
read_number <- c("1", "2")
names(read_number) <- sample_names
direc <- c("sense", "anti")
names(direc) <- direc

#### load files ####
samples <-
  lapply(sample_names, function(i) {
    read.table(sprintf("./results/%s_consensus_RNA_HS2.txt", i),header = TRUE)
  })
# bed <-
#   lapply(sample_names,function(i){
#     read.table(sprintf("./bed/%s_consensus_RNA_HS2.bed",i))
#   })
basicgroups <- readRDS("./results/basicgroups.rds")

#### add info of names and bins ####
samples_trans <- lapply(sample_names,function(sn){
  samples[[sn]]$trans <- factor(names(trans)[samples[[sn]]$trans_num],
                                levels = names(trans))
  samples[[sn]]
})

####
samples_ratio <- lapply(sample_names,function(sn){
  sums <- table(samples_trans[[sn]]$trans)
  mis <- tapply(samples_trans[[sn]]$SNP_num,samples_trans[[sn]]$trans,sum)
  ratio <- as.data.frame(mis/sums)
  ratio$samples <- sn
  ratio$bins <- basicgroups[ratio[,1]]
  ratio
})
#### reform ####
mis_ratio <- Reduce(rbind,samples_ratio)
mis_13vs20 <- Reduce(cbind,samples_ratio)
colnames(mis_13vs20) <- c("trans","num13","samples","bins","trans2","num20","samples2","bins2")
mis_13vs20$dif <- mis_13vs20$num20 - mis_13vs20$num13
mis_13vs20$aver <- (mis_13vs20$num20 + mis_13vs20$num13)/2

```

### **fig** {#fig-1}

```{r}
ggplot(subset(mis_13vs20,!is.na(bins)))+
  geom_boxplot(aes(x=bins,y=aver))+
  xlab("")+
  ylab("Numbers of mismatch \nper read")+
  # scale_y_log10()+
  # facet_grid(years~.)
  mytheme()
ggsave(
    filename = sprintf("./fig/rna_snps_boxplot1.pdf"),
    width = Width ,
    height = Height)
ggplot(subset(mis_13vs20,!is.na(bins)))+
  geom_boxplot(aes(x=bins,y=dif))+
  xlab("")+
  ylab("Difference of numbers of \nmismatch per read")+
  # scale_y_log10()+
  # facet_grid(years~.)
  mytheme()
ggsave(
    filename = sprintf("./fig/rna_snps_boxplot2.pdf", i),
    width = Width*1.4,
    height = Height*1.4)
ggplot(subset(mis_13vs20,!is.na(bins)))+
  geom_point(aes(x=num13,y=num20,col=bins))+
  xlab("Numbers of mismatch \nper read in 2013")+
  ylab("Numbers of mismatch \nper read in 2020")+
  scale_color_manual(values = c("#401962", "#466E8B", "#70BA73", "#A2D34F", "#dde136")) +
  xlim(0,1.2)+
  ylim(0,1.2)+
  geom_abline(intercept = 0, slope = 1) +
  mytheme()
ggsave(
    filename = sprintf("./fig/rna_snps_scatter.pdf"),
    width = Width*1.4,
    height = Height*1.4)
```
