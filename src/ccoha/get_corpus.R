library(dplyr)
library(stringr)
library(ggplot2)
library(foreach)
library(doParallel)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
rm(list=ls())

options(dplyr.summarise.inform = FALSE)

folders <- list.files("~/phd/data/ccoha", pattern = "CCOHA", full.names = T)
files <- list.files(paste0(folders, "/tagged"), full.names = T, pattern = "\\.txt", recursive = T)
genyears <- gsub(".*\\/|_[0-9]+\\.txt", "", files)
genyears <- strsplit(genyears, "_")
genyears <- do.call(rbind, genyears)
genyears <- as.data.frame(genyears)
table(genyears$V1, genyears$V2)

files <- files[grepl("(mag|news)_(19|20)", files)]

searchwords <- c("conspiracy", "theory")

#setup parallel backend to use many processors
cores = detectCores()
cl <- makeCluster(cores[1]-4) #not to overload your computer
registerDoParallel(cl)

clusterExport(cl, c("searchwords"))
clusterEvalQ(cl, library("stringr"))
clusterEvalQ(cl, library("dplyr"))
out.par <- parSapply(cl, files, function(file){
  #file = files[1]
  temp <- readLines(file)
  temp <- strsplit(temp, "\t")
  temp <- as.data.frame(do.call(rbind, temp))
  colnames(temp) <- c("txt", "token", "fo")
  id <- temp$txt[1]
  temp <- temp[-1,]
  
  # Find indices where <eos> occurs
  eos_indices <- which(temp$fo == "<eos>")
  
  # Split the dataframe based on <eos> indices
  temp_list <- split(temp, cumsum(seq_along(temp$fo) %in% eos_indices))
  
  # Remove the <eos> rows from each dataframe in the list
  temp_list <- lapply(temp_list, function(x) x[x$fo != "<eos>", ])
  
  # Log the sentences containing the search words
  logs <- sapply(temp_list, function(x) any(x$token %in% searchwords))
  logs.terms <- sapply(temp_list, function(x) searchwords[searchwords %in% x$token])
  
  temp <- do.call(rbind, temp_list[logs])
  temp$sentence_id <- as.numeric(gsub("\\..*", "", rownames(temp)))
  temp$doc_id <- id
  temp$category <- str_extract(file, "(mag|news|fic|nf)(?=_)")
  temp$year <- as.numeric(str_extract(file, "(?<=(mag|news|fic|nf)_)[0-9]+"))
  
  if(!is.data.frame(temp)) return(NULL)
  
  #container[[file]] <- as_tibble(temp)
  temp.sparse <- temp %>% 
    group_by(category, year, doc_id, sentence_id) %>% 
    summarise(s.plain = paste0(token, collapse = " "),
              s.pos = paste0(token, "__", fo, collapse = " "),
              pos = paste0(fo, collapse = " "))
  temp.sparse$target <- logs.terms[logs]
  #container.sparse[[file]] <- temp.sparse
  list(long = as_tibble(temp), wide = temp.sparse)
})

container <- sapply(out.par, "[[", 1)
container.sparse <- sapply(out.par, "[[", 2)
rm(out.par)
container <- do.call(rbind, container)
container.sparse <- do.call(rbind, container.sparse)

saveRDS(container, file = "../../output/ccoha/data/ct_long.RDS")
saveRDS(container.sparse, file = "../../output/ccoha/data/ct_wide.RDS")

container.sparse <- readRDS("../output/data/ct_wide.RDS")
container.sparse <- container.sparse %>% filter(!grepl("<nul>", s.plain))
write.table(container.sparse$s.plain, file = "../../output/ccoha/data/sentences_ct.txt", row.names = F, quote = F, col.names = F)

stopCluster(cl)

ggplot(container.sparse, aes(x = year)) +
  geom_histogram()

table(grepl("conspiracy theory", container.sparse$s.plain))

library(umap)
library(plotly)
library(factoextra)
library(cluster)
library(Rtsne)
l <- read.table("../output/embeddings/embeddings_ct.tsv", sep = "\t")
meta <- readLines("../output/embeddings/metadata_ct.tsv")

l$lemma <- meta
eos_indices <- which(meta == "@@@") + 1
rep_lengths <- lengths(split(meta, f = cumsum(seq_along(meta) %in% eos_indices)))
rep_lengths <- unname(rep_lengths)
sup <- container.sparse[1:1000,]
l <- cbind(l, sup[rep(1:nrow(sup), rep_lengths),])
grep(l$lemma)
school <- filter(l, lemma == "school")
school.umap <- umap(school[1:768])

cl <- hclust(dist(scale(school[1:768])), method = "ward.D")
fviz_nbclust(school[1:768], FUN = hcut, method = "wss")
fviz_nbclust(school[1:768], FUN = hcut, method = "silhouette")
gap_stat <- clusGap(school[1:768], FUN = hcut, nstart = 25, K.max = 10, B = 50)
fviz_gap_stat(gap_stat)
sub_grp <- cutree(cl, k = 5)
school$cluster <- sub_grp

p <- ggplot(as.data.frame(school.umap$layout), aes(x = V1, y = V2, label = school$s.plain, color = as.factor(school$cluster))) +
  geom_point()
ggplotly(p)

rtsne <- Rtsne(school[!duplicated(school[1:768]), 1:768])
q <- ggplot(as.data.frame(rtsne$Y), aes(x = V1, y = V2, label = school$s.plain[!duplicated(school[1:768])], color = as.factor(school$cluster[!duplicated(school[1:768])]))) +
  geom_point()
ggplotly(q)
