## Overview

This file contains the all code needed to reproduce the result of the following manuscript:

**Autonomous Shaping of the piRNA Sequence Repertoire by Competition between Adjacent Ping-Pong Amplification Sites**

*Jie Yu, Fumiko Kawasaki, Natsuko Izumi, Takashi Kiuchi, Susumu Katsuma, Yukihide Tomari, Keisuke Shoji*

## Requirement

To perform the code, `cutadapt`, `Bowtie`, `bedtools`, `samtools` and `R` are required.

`pkgs` contains all packages required in R.

```{r}
pkgs <- c("ggplot2","BiocManager","Biostrings","tidyr","ggsignif","DECIPHER")
suppressMessages(sapply(pkgs,library,character.only=TRUE))
```

## Note

Each code block is labeled with the running environment (bash or R) in bold font before it.

All results start with `fig_` were saved in the same .RData file, `figdata.RData`, for picturing figures.

```{r}
save(list = ls(pattern = "fig_"), file = "figdata.RData")
```
