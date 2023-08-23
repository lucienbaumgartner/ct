library(dplyr)
library(ggplot2)
library(ggrepel)
library(data.table)
library(RANN)
library(Rtsne)
library(factoextra)
library(stopwords)
library(quanteda)
library(plotly)
library(udpipe)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
rm(list=ls())

udmodel <- udpipe::udpipe_load_model("/Users/lucienbaumgartner/phd/data/udpipe/english-gum-ud-2.5-191206.udpipe")
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
## combine "conspiracy" & "theory"
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

## combine ".*" & "##.*"
ct.index <- which(grepl("\\#\\#", df$lemma))
#ct.index+1 == lead(ct.index)
#ct.index[ct.index+1 == lead(ct.index)]
#df$lemma[ct.index+1 == lead(ct.index)]
ct.indices <- c()
for(i in ct.index){
  ct.indices <- append(ct.indices, c(i-1, i))
}
ct.indices <- ct.indices[!duplicated(ct.indices)]
df.sub <- df[ct.indices,1:768]
eos_indices <- c(1:length(df$lemma[ct.indices]))[!grepl("\\#\\#", df$lemma[ct.indices])]
temp_list <- split(df.sub, cumsum(seq_along(c(1:length(df$lemma[ct.indices]))) %in% eos_indices))
g_rep <- rep(1:length(eos_indices), unname(sapply(temp_list, nrow)))
df.sub$g <- g_rep
df.sub <- df.sub %>% group_by(g) %>% summarise(across(starts_with('V'), mean))
length(eos_indices) == nrow(df.sub)
df[eos_indices,1:768] <- df.sub[2:769]
temp_list <- split(df$lemma[ct.indices], g_rep)
temp_list <- lapply(temp_list, function(x) paste0(gsub("\\#", "", x), collapse = ""))
df$lemma[eos_indices] <- unlist(temp_list)
df <- df[-ct.index,]

rm(df.sub)
rm(temp_list)



as.data.frame(udpipe::udpipe_annotate(udmodel, df$lemma[1:100]))

tibble(doc_id = df$id[1], text = df$txt[1]) %>% 
  udpipe(., object = udmodel) %>% as_tibble() %>% 
  mutate(phrase_tag = as_phrasemachine(upos, type = 'upos')) %>% 
  group_by(paragraph_id, sentence_id) %>% 
  do(keywords_phrases(.$phrase_tag, .$token, pattern = "N+N", 
                      is_regex = TRUE, 
                      ngram_max = 4, 
                      detailed = TRUE))

# mono-semic approach NEW
sem <- df %>% group_by(subcorpus, lemma) %>% summarise(across(starts_with('V'), mean))
sem <- sem %>% filter(!grepl("\\#", lemma))
sem <- sem %>% filter(!grepl("[[:punct:]]", lemma))
sem <- sem %>% filter(!lemma %in% stopwords::stopwords())
sem <- sem %>% filter(!grepl("[0-9]+", lemma))

sem.sub <- split(sem, sem$subcorpus)
write.table(sem.sub$conspiracy[-c(1:2)], file = "../../output/loco/data/data_conspiracy.txt", row.names = F, quote = F, col.names = F, sep = ",")
write.table(sem.sub$mainstream[-c(1:2)], file = "../../output/loco/data/data_mainstream.txt", row.names = F, quote = F, col.names = F, sep = ",")

tsne <- lapply(sem.sub, function(x){
  # PCA
  x = sem.sub$conspiracy
  pca.res <- prcomp(x[grepl("V", colnames(x))], center = T, scale. = T)
  expl_var <- pca.res$sdev^2/sum(pca.res$sdev^2)
  N_perm <- 10
  expl_var_perm <- matrix(NA, ncol = length(pca.res$sdev), nrow = N_perm)
  set.seed(12345)
  for(k in 1:N_perm){
    expr_perm <- apply(x[grepl("V", colnames(x))],2,sample)
    PC_perm <- prcomp(expr_perm, center=TRUE, scale=TRUE)
    expl_var_perm[k,] <- PC_perm$sdev^2/sum(PC_perm$sdev^2)
  }
  tmp_df1 <- data.frame(explained_variance = expl_var[1:500], PCs = seq(1:500),
                        var_perm = colMeans(expl_var_perm)[1:500])
  
  pval <- apply(t(expl_var_perm) >= expl_var,1,sum) / N_perm
  optPC <- head(which(pval>=0.05),1)-1
  tmp_df2 <- data.frame(p_value = pval[1:500], PCs = seq(1:500))
  
  # t-SNE
  init <- pca.res$x[,1:2]/sd(pca.res$x[,1])*0.0001
  tsne <- Rtsne(pca.res$x[,1:optPC], is_distance = FALSE, 
                verbose = F, num_threads=3, 
                perplexity = nrow(pca.res$x)**0.5,
                #perplexity = 30,
                Y_init = init,
                #eta = 200,
                pca = F,
                eta = nrow(pca.res$x)/12
  )
  tsne <- as.data.frame(tsne$Y)
  tsne$lemma <- x$lemma
  rownames(tsne) <- x$lemma
  return(tsne)
})

### Nearest Neighbours (k-NN)
nn.indices <- nn2(tsne$conspiracy[1:2], as.matrix(tsne$conspiracy["conspiracy theory",1:2]), k = 50)
rownames(tsne$conspiracy)[nn.indices$nn.idx]
tsne$conspiracy$group <- ifelse(1:nrow(tsne$conspiracy) %in% nn.indices$nn.idx, "Top 50", "Rest")

nn.indices <- nn2(tsne$mainstream[1:2], as.matrix(tsne$mainstream["conspiracy theory",1:2]), k = 50)
rownames(tsne$mainstream)[nn.indices$nn.idx]
tsne$mainstream$group <- ifelse(1:nrow(tsne$mainstream) %in% nn.indices$nn.idx, "Top 50", "Rest") 

sz <- 3
xdot <- tsne$conspiracy[tsne$conspiracy$lemma == "conspiracy theory",]$V1
xlim <- sort(c(xdot - sz, xdot + sz))
ydot <- tsne$conspiracy[tsne$conspiracy$lemma == "conspiracy theory",]$V2
ylim <- sort(c(ydot - sz, ydot + sz))
g <- ggplot(tsne$conspiracy, aes(V1, V2, label = lemma, color = group)) +
  geom_point(alpha = .8) +
  geom_text_repel(seed = 2637) +
  scale_y_continuous(limits = ylim) +
  scale_x_continuous(limits = xlim)
ggsave(g, file = "../../output/loco/plots/tsne_conspiracy_local.png", width = 10, height = 10)
p <- ggplot(tsne$conspiracy, aes(V1, V2, label = lemma, color = group)) +
  geom_point(alpha = .8)
ggsave(p, file = "../../output/loco/plots/tsne_conspiracy_global.png", width = 10, height = 10)




sz <- 3
xdot <- tsne$mainstream[tsne$mainstream$lemma == "conspiracy theory",]$V1
xlim <- sort(c(xdot - sz, xdot + sz))
ydot <- tsne$mainstream[tsne$mainstream$lemma == "conspiracy theory",]$V2
ylim <- sort(c(ydot - sz, ydot + sz))
g <- ggplot(tsne$mainstream, aes(V1, V2, label = lemma, color = group)) +
  geom_point(alpha = .8) +
  geom_text_repel(seed = 2637) +
  scale_y_continuous(limits = ylim) +
  scale_x_continuous(limits = xlim)
g
ggsave(g, file = "../../output/loco/plots/tsne_mainstream_local.png", width = 10, height = 10)
p <- ggplot(tsne$mainstream, aes(V1, V2, label = lemma, color = group)) +
  geom_point(alpha = .8)
ggsave(p, file = "../../output/loco/plots/tsne_mainstream_global.png", width = 10, height = 10)



####### PACMAP

consp_pacmap <- read.csv("../../output/loco/data/pacmap_data_conspiracy.txt")
consp_pacmap <- cbind(consp_pacmap, sem.sub$conspiracy[1:2])
rownames(consp_pacmap) <- consp_pacmap$lemma

nn.indices <- nn2(consp_pacmap[1:2], as.matrix(consp_pacmap["conspiracy theory",1:2]), k = 50)
rownames(consp_pacmap)[nn.indices$nn.idx]
consp_pacmap$group <- ifelse(1:nrow(consp_pacmap) %in% nn.indices$nn.idx, "Top 50", "Rest")

sz <- 1
xdot <- consp_pacmap[consp_pacmap$lemma == "conspiracy theory",]$Dim1
xlim <- sort(c(xdot - sz, xdot + sz))
ydot <- consp_pacmap[consp_pacmap$lemma == "conspiracy theory",]$Dim2
ylim <- sort(c(ydot - sz, ydot + sz))
g <- ggplot(consp_pacmap, aes(Dim1, Dim2, label = lemma, color = group)) +
  geom_point(alpha = .8) +
  geom_text_repel(seed = 2637) +
  scale_y_continuous(limits = ylim) +
  scale_x_continuous(limits = xlim)
g
ggsave(g, file = "../../output/loco/plots/pacmap_conspiracy_local.png", width = 10, height = 10)
p <- ggplot(consp_pacmap, aes(Dim1, Dim2, label = lemma, color = group)) +
  geom_point(alpha = .8)
p
ggsave(p, file = "../../output/loco/plots/pacmap_conspiracy_global.png", width = 10, height = 10)



mainstr_pacmap <- read.csv("../../output/loco/data/pacmap_data_mainstream.txt")
mainstr_pacmap <- cbind(mainstr_pacmap, sem.sub$mainstream[1:2])
rownames(mainstr_pacmap) <- mainstr_pacmap$lemma

nn.indices <- nn2(mainstr_pacmap[1:2], as.matrix(mainstr_pacmap["conspiracy theory",1:2]), k = 50)
rownames(mainstr_pacmap)[nn.indices$nn.idx]
mainstr_pacmap$group <- ifelse(1:nrow(mainstr_pacmap) %in% nn.indices$nn.idx, "Top 50", "Rest")

sz <- 1
xdot <- mainstr_pacmap[mainstr_pacmap$lemma == "conspiracy theory",]$Dim1
xlim <- sort(c(xdot - sz, xdot + sz))
ydot <- mainstr_pacmap[mainstr_pacmap$lemma == "conspiracy theory",]$Dim2
ylim <- sort(c(ydot - sz, ydot + sz))
g <- ggplot(mainstr_pacmap, aes(Dim1, Dim2, label = lemma, color = group)) +
  geom_point(alpha = .8) +
  geom_text_repel(seed = 2637) +
  scale_y_continuous(limits = ylim) +
  scale_x_continuous(limits = xlim)
g
ggsave(g, file = "../../output/loco/plots/pacmap_mainstream_local.png", width = 10, height = 10)
p <- ggplot(mainstr_pacmap, aes(Dim1, Dim2, label = lemma, color = group)) +
  geom_point(alpha = .8)
p
ggsave(p, file = "../../output/loco/plots/pacmap_mainstream_global.png", width = 10, height = 10)










r <- ggplot(tsne, aes(V1, V2, label = lemma, color = group)) + 
  geom_point(alpha = .8)
ggplotly(r)

rm(df)
rm(df.sub)
sem <- sem.proto %>% filter(!grepl("\\#", lemma))
sem <- sem.proto %>% filter(!grepl("[[:punct:]]", lemma))
rm(sem.proto)
sem$lemma
#sem <- sem %>% group_by(lemma) %>% summarise(across(starts_with('V'), mean))
sem.mat <- as.matrix(sem[2:769])
rownames(sem.mat) <- sem$lemma

nn.indices <- nn2(sem.mat[,1:768], t(as.matrix(sem.mat["conspiracy theory",1:768])), k = 500)
rownames(sem.mat)[nn.indices$nn.idx]

### t-SNE
seeds <- c(678,7324,998, 87797)
for(seed in seeds){
  set.seed(seed)
  rtsne.sem <- Rtsne(sem.mat, num_threads = 3, pca = F)
  rtsne.sem <- as.data.frame(rtsne.sem$Y)
  rtsne.sem$lemma <- sem$lemma
  rownames(rtsne.sem) <- sem$lemma
  
  ### Nearest Neighbours (k-NN)
  nn.indices <- nn2(rtsne.sem[1:2], as.matrix(rtsne.sem["conspiracy theory",1:2]), k = 50)
  rownames(sem.mat)[nn.indices$nn.idx]
  rtsne.sem$group <- ifelse(1:nrow(rtsne.sem) %in% nn.indices$nn.idx, "Top 50", "Rest") 
  
  #rtsne.sem <- filter(rtsne.sem, !lemma %in% c("prosperity", "prosperous", "vera"))
  
  q <- ggplot(rtsne.sem[!rtsne.sem$lemma == "conspiracy theory",], aes(V1, V2, label = lemma, color = group)) +
    geom_point(alpha = .8)
  
  ggsave(q, file = paste0("../../output/loco/plots/tsne_", seed,".png"), width = 6, height = 5)
  write.table(rownames(sem.mat)[nn.indices$nn.idx], file = paste0("../../output/loco/recovery/tsne_", seed,".txt"))
}

seed <- 4712874
set.seed(seed)
rtsne.sem <- Rtsne(sem.mat)
rtsne.sem <- as.data.frame(rtsne.sem$Y)
rtsne.sem$lemma <- sem$lemma
rownames(rtsne.sem) <- sem$lemma

### Nearest Neighbours (k-NN)
nn.indices <- nn2(rtsne.sem[1:2], as.matrix(rtsne.sem["conspiracy theory",1:2]), k = 17)
rownames(sem.mat)[nn.indices$nn.idx]
rtsne.sem$group <- ifelse(1:nrow(rtsne.sem) %in% nn.indices$nn.idx, "Top 50", "Rest") 

q <- ggplot(rtsne.sem[!rtsne.sem$lemma == "conspiracy theory",], aes(V1, V2, label = lemma, color = group)) +
  geom_point(alpha = .8) +
  geom_text_repel(seed = 2637) +
  geom_point(data = rtsne.sem[rtsne.sem$lemma == "conspiracy theory",], color = "firebrick1") +
  geom_text(data = rtsne.sem[rtsne.sem$lemma == "conspiracy theory",], color = "firebrick1", hjust = -0.1) +
  scale_y_continuous(limits = c(-25, -19)) +
  scale_x_continuous(limits = c(-7, 0)) +
  #scale_y_continuous(limits = c(14, 18)) +
  #scale_x_continuous(limits = c(10, 15)) +
  scale_color_manual(values = c("bisque4", "darkblue")) +
  theme_bw() + 
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "Dimension 2", x = "Dimension 1")
q
ggsave(q, file = "../../output/loco/plots/tsne___4712874.png", width = 6, height = 5)


q <- ggplot(rtsne.sem[!rtsne.sem$lemma == "conspiracy theory",], aes(V1, V2, label = lemma, color = group)) +
  geom_point(alpha = .8) +
  geom_text_repel(seed = 2637) +
  geom_point(data = rtsne.sem[rtsne.sem$lemma == "conspiracy theory",], color = "firebrick1") +
  geom_text(data = rtsne.sem[rtsne.sem$lemma == "conspiracy theory",], color = "firebrick1", hjust = -0.1) +
  scale_y_continuous(limits = c(-13.5, -9.5)) +
  scale_x_continuous(limits = c(-23, -17)) +
  scale_color_manual(values = c("bisque4", "darkblue")) +
  theme_bw() + 
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "Dimension 2", x = "Dimension 1")

ggsave(q, file = "../../output/loco/plots/tsne.png", width = 6, height = 5)


# mono-semic approach OLD
sem.proto <- df %>% group_by(lemma) %>% summarise(across(starts_with('V'), mean))
sem <- sem.proto %>% filter(!grepl("\\#", lemma))
#sem <- sem.proto
#sem <- sem.proto %>% filter(!grepl("[[:punct:]]", lemma))
#sem <- sem %>% filter(!lemma %in% stop_words$word)
sem$lemma
sem.mat <- as.matrix(sem[2:769])
rownames(sem.mat) <- sem$lemma

#nn.indices <- nn2(sem.mat[,1:768], t(as.matrix(sem.mat["conspiracy theory",1:768])), k = 30)
#rownames(sem.mat)[nn.indices$nn.idx]

### t-SNE
set.seed(12345)
rtsne.sem <- Rtsne(sem.mat)
rtsne.sem <- as.data.frame(rtsne.sem$Y)
rtsne.sem$lemma <- sem$lemma
rownames(rtsne.sem) <- sem$lemma

### Nearest Neighbours (k-NN)
nn.indices <- nn2(rtsne.sem[1:2], as.matrix(rtsne.sem["conspiracy theory",1:2]), k = 17)
rownames(sem.mat)[nn.indices$nn.idx]
rtsne.sem$group <- ifelse(1:nrow(rtsne.sem) %in% nn.indices$nn.idx, "Top 50", "Rest") 

rtsne.sem <- filter(rtsne.sem, !lemma %in% c("prosperity", "prosperous", "vera"))

q <- ggplot(rtsne.sem[!rtsne.sem$lemma == "conspiracy theory",], aes(V1, V2, label = lemma, color = group)) +
  geom_point(alpha = .8) +
  geom_text_repel(seed = 2637) +
  geom_point(data = rtsne.sem[rtsne.sem$lemma == "conspiracy theory",], color = "firebrick1") +
  geom_text(data = rtsne.sem[rtsne.sem$lemma == "conspiracy theory",], color = "firebrick1", hjust = -0.1) +
  scale_y_continuous(limits = c(-25, -19)) +
  scale_x_continuous(limits = c(-7, 0)) +
  #scale_y_continuous(limits = c(14, 18)) +
  #scale_x_continuous(limits = c(10, 15)) +
  scale_color_manual(values = c("bisque4", "darkblue")) +
  theme_bw() + 
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "Dimension 2", x = "Dimension 1")
ggsave()

q <- ggplot(rtsne.sem[!rtsne.sem$lemma == "conspiracy theory",], aes(V1, V2, label = lemma, color = group)) +
  geom_point(alpha = .8) +
  geom_text_repel(seed = 2637) +
  geom_point(data = rtsne.sem[rtsne.sem$lemma == "conspiracy theory",], color = "firebrick1") +
  geom_text(data = rtsne.sem[rtsne.sem$lemma == "conspiracy theory",], color = "firebrick1", hjust = -0.1) +
  scale_y_continuous(limits = c(-13.5, -9.5)) +
  scale_x_continuous(limits = c(-23, -17)) +
  #scale_y_continuous(limits = c(14, 18)) +
  #scale_x_continuous(limits = c(10, 15)) +
  scale_color_manual(values = c("bisque4", "darkblue")) +
  theme_bw() + 
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "Dimension 2", x = "Dimension 1")
q


ggsave(q, file = "../../output/loco/plots/tsne.png", width = 6, height = 5)





### poly-semic approach
ct.sem <- df %>% filter(lemma == "conspiracy theory")
set.seed(12345)
dist_mat <- dist(ct.sem[,1:768])
hclust <- hclust(dist_mat, method = 'ward')
#fviz_nbclust(ct.sem[,1:768], FUN = hcut, method = "wss")
sub_grp <- cutree(hclust, k = 4)
ct.sem$lemma <- paste0("conspiracy theory_", sub_grp)
df[df$lemma == "conspiracy theory",] <- ct.sem

sem <- df %>% filter(!grepl("\\#", lemma))
sem <- sem %>% group_by(lemma) %>% summarise(across(starts_with('V'), mean))
sem.mat <- as.matrix(sem[2:769])
rownames(sem.mat) <- sem$lemma
nn.indices <- nn2(sem.mat[,1:768], as.matrix(sem.mat[paste0("conspiracy theory_", 1:4),1:768]), k = 50)
rownames(sem.mat)[nn.indices$nn.idx]

set.seed(12345)
rtsne.sem <- Rtsne(sem.mat)
rtsne.sem <- as.data.frame(rtsne.sem$Y)
rtsne.sem$lemma <- sem$lemma
rownames(rtsne.sem) <- sem$lemma
nn.indices <- nn2(rtsne.sem[1:2], as.matrix(rtsne.sem[paste0("conspiracy theory_", 1:4),1:2]), k = 50)
rownames(sem.mat)[nn.indices$nn.idx]
rtsne.sem$group <- ifelse(1:nrow(rtsne.sem) %in% unique(nn.indices$nn.idx), "Top 50", "Rest") 

ggplot(rtsne.sem[!rtsne.sem$lemma %in% paste0("conspiracy theory_", 1:4),], aes(V1, V2, label = lemma, color = group)) +
  geom_point(alpha = .8) +
  geom_point(data = rtsne.sem[rtsne.sem$lemma %in% paste0("conspiracy theory_", 1:4),], color = "firebrick1") +
  geom_text(data = rtsne.sem[rtsne.sem$lemma %in% paste0("conspiracy theory_", 1:4),], color = "firebrick1", hjust = -0.1) +
  geom_text_repel(seed = 2637) +
  scale_y_continuous(limits = c(-22, -15)) +
  scale_x_continuous(limits = c(-18, -13))

