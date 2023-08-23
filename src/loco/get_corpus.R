library(dplyr)
library(jsonlite)
library(tokenizers)
library(ggplot2)
library(tidytext)
library(textstem)
library(foreach)
library(doParallel)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
rm(list=ls())

df <- fromJSON("~/phd/data/LOCO.json", flatten = T)
names(df)
table(is.na(df$date))
df <- mutate(df, date = as.Date(date))
head(df$date)
df <- filter(df, !is.na(date))
range(df$date)

ggplot(df, aes(x = date)) +
  geom_histogram()

sentences <- tokenize_sentences(df$txt)

#setup parallel backend to use many processors
cores = detectCores()
cl <- makeCluster(cores[1]-4) #not to overload your computer
registerDoParallel(cl)

clusterEvalQ(cl, library("textstem"))
sentences <- parSapply(cl, sentences, lemmatize_strings)
sentences[[1]]

sentence.ids <- rep(1:length(sentences), lengths(sentences))

sentences <- tibble(txt_lemma = unlist(sentences), id = sentence.ids)
sentences <- cbind(sentences, df[sentence.ids,])
sentences <- sentences %>% mutate(ct_dummy = grepl("conspiracy theory", txt_lemma))

table(sentences$ct_dummy, sentences$subcorpus)

set.seed(46659)
ct_sentences <- filter(sentences, ct_dummy)
ct_sentences <- ct_sentences[sample(1:nrow(ct_sentences), 5000),]
subsample <- sentences[!sentences$ct_dummy,]
set.seed(46659)
subsample <- subsample %>% group_by(subcorpus) %>% sample_n(10000)
table(subsample$subcorpus)
subsample <- rbind(subsample, ct_sentences)
subsample <- mutate(subsample, txt_lemma = tolower(txt_lemma))

save(subsample, file = "../../output/loco/data/subsample_loco.RDS")

write.table(subsample$txt_lemma, file = "../../output/loco/data/sentences.txt", row.names = F, quote = F, col.names = F)

