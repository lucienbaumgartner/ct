library(dplyr)
library(ggplot2)
library(ggrepel)
library(data.table)
library(RANN)
library(Rtsne)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
rm(list=ls())

load("../../output/loco/data/subsample_loco.RDS")

df <- fread("../../output/loco/embeddings/embeddings_sentences.tsv", sep = "\t")
meta <- readLines("../../output/loco/embeddings/metadata_sentences.tsv")

df$lemma <- meta
eos_indices <- which(meta == "@@@") + 1
rep_lengths <- lengths(split(meta, f = cumsum(seq_along(meta) %in% eos_indices)))
rep_lengths <- unname(rep_lengths)
subsample <- subsample[rep(1:nrow(subsample), rep_lengths),]
df <- cbind(df, subsample)

ggplot(df, aes(x = date)) +
  geom_histogram()

rm(subsample)
rm(meta)

df <- df %>% filter(!lemma == "@@@")
ct.index <- which(df$lemma == "conspiracy" & lead(df$lemma) == "theory")
ct.indices <- c()
for(i in ct.index){
  ct.indices <- append(ct.indices, c(i, i+1))
}
df.sub <- df[ct.indices,1:768]
df.sub$g <- rep(1:length(ct.index), 2)
df.sub <- df.sub %>% group_by(g) %>% summarise(across(starts_with('V'), mean))
df[ct.index,1:768] <- df.sub[2:769]
df$lemma[ct.index] <- "conspiracy theory"
df <- df[-(ct.index+1),]

table(df$subcorpus[df$lemma == "conspiracy theory"])

sem <- df %>% filter(!grepl("\\#", lemma))
#sem <- df %>% filter(!grepl("[[:punct:]]", lemma))
sem$lemma
sem <- sem %>% group_by(lemma) %>% summarise(across(starts_with('V'), mean))
sem.mat <- as.matrix(sem[2:769])
rownames(sem.mat) <- sem$lemma

set.seed(12345)
set.seed(123)
rtsne.sem <- Rtsne(sem.mat)
rtsne.sem <- as.data.frame(rtsne.sem$Y)
rtsne.sem$lemma <- sem$lemma
rownames(rtsne.sem) <- sem$lemma

nn.indices <- nn2(rtsne.sem[1:2], as.matrix(rtsne.sem["conspiracy theory",1:2]), k = 30)
rownames(sem.mat)[nn.indices$nn.idx]
rtsne.sem$group <- ifelse(1:nrow(rtsne.sem) %in% nn.indices$nn.idx, "Top 30", "Rest") 

rtsne.sem <- filter(rtsne.sem, !lemma %in% c("prosperity", "prosperous", "vera"))

ggplot(rtsne.sem, aes(V1, V2, label = lemma, color = group)) +
  geom_point(alpha = .1) +
  geom_text_repel() +
  scale_y_continuous(limits = c(-13, -12)) +
  scale_x_continuous(limits = c(-21, -14))

ggplot(rtsne.sem[!rtsne.sem$lemma == "conspiracy theory",], aes(V1, V2, label = lemma, color = group)) +
  geom_point(alpha = .8) +
  geom_text_repel(seed = 2637) +
  geom_point(data = rtsne.sem[rtsne.sem$lemma == "conspiracy theory",], color = "black") +
  geom_text(data = rtsne.sem[rtsne.sem$lemma == "conspiracy theory",], color = "black", hjust = -0.1) +
  scale_y_continuous(limits = c(-13.5, -9.5)) +
  scale_x_continuous(limits = c(-23, -17)) +
  theme_classic()

ggplot(rtsne.sem, aes(V1, V2, label = lemma, color = group)) +
  geom_point(alpha = .1) +
  geom_text_repel() +
  scale_y_continuous(limits = c(-20, -13)) +
  scale_x_continuous(limits = c(7, 12))
