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
sem <- df %>% group_by(lemma) %>% summarise(across(starts_with('V'), mean))
sem <- sem %>% filter(!grepl("\\#", lemma))
sem <- sem %>% filter(!grepl("[[:punct:]]", lemma))
sem <- sem %>% filter(!lemma %in% stopwords::stopwords())
sem <- sem %>% filter(!grepl("[0-9]+", lemma))
nn <- nn2(as.matrix(sem[,2:769]), as.matrix(sem[sem$lemma == "conspiracy theory",2:769]), k = round(nrow(sem)/10))
sem <- sem[as.vector(nn$nn.idx),]

# PCA
pca.res <- prcomp(sem[2:769], center = T, scale. = T)
expl_var <- pca.res$sdev^2/sum(pca.res$sdev^2)
N_perm <- 10
expl_var_perm <- matrix(NA, ncol = length(pca.res$sdev), nrow = N_perm)
set.seed(12345)
for(k in 1:N_perm){
  expr_perm <- apply(sem[2:769],2,sample)
  PC_perm <- prcomp(expr_perm, center=TRUE, scale=TRUE)
  expl_var_perm[k,] <- PC_perm$sdev^2/sum(PC_perm$sdev^2)
}
tmp_df1 <- data.frame(explained_variance = expl_var[1:50], PCs = seq(1:50),
                      var_perm = colMeans(expl_var_perm)[1:50])
gp1 <- ggplot(tmp_df1) +
  geom_point(aes(PCs, explained_variance, shape = "." )) +
  geom_line(aes(PCs, explained_variance)) +
  geom_point(aes(PCs, var_perm, colour = "red", shape = ".")) +
  geom_line(aes(PCs, var_perm, colour = "red")) +
  theme(legend.position = "none")

pval <- apply(t(expl_var_perm) >= expl_var,1,sum) / N_perm
optPC <- head(which(pval>=0.05),1)-1
tmp_df2 <- data.frame(p_value = pval[1:500], PCs = seq(1:500))
gp2 <- ggplot(tmp_df2, aes(PCs, p_value, shape = ".")) +
  geom_point() +
  geom_line() +
  ggtitle(paste("optimal number of PCs = ", optPC)) +
  theme(legend.position = "none")

gp1
gp2

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
tsne$lemma <- sem$lemma
rownames(tsne) <- sem$lemma

### Nearest Neighbours (k-NN)
nn.indices <- nn2(tsne[1:2], as.matrix(tsne["conspiracy theory",1:2]), k = 100)
rownames(tsne)[nn.indices$nn.idx]
tsne$group <- ifelse(1:nrow(tsne) %in% nn.indices$nn.idx, "Top 50", "Rest") 

xdot <- tsne[tsne$lemma == "conspiracy theory",]$V1
xlim <- sort(c(xdot - 10, xdot + 10))
ydot <- tsne[tsne$lemma == "conspiracy theory",]$V2
ylim <- sort(c(ydot - 10, ydot + 10))
ggplot(tsne, aes(V1, V2, label = lemma, color = group)) +
  geom_point(alpha = .8) +
  geom_text_repel(seed = 2637) +
  scale_y_continuous(limits = ylim) +
  scale_x_continuous(limits = xlim)

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

