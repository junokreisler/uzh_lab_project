# Read count data (normalized) from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80744
# RNA-seq profiling of 728 Arabidopsis thaliana accessions
data <- read.table('./data/GSE80744_ath1001_tx_norm_2016-04-21-UQ_gNorm_normCounts_k4.tsv/GSE80744_ath1001_tx_norm_2016-04-21-UQ_gNorm_normCounts_k4.tsv',
                   header = TRUE)

rownames(data) <- data$gene_id
data<- data[,-1]
# the numbers are actually characters. so gotta change them to numbers
for (i in c(1:length(colnames(data)))) {
  data[,i] <- as.numeric(data[,i])
}
class(data[,32])
data <- t(data)
# testing on one gene pair
plot(data[,'AT1G01010'],data[,'AT1G01020'])
cor(data[,'AT1G01010'],data[,'AT1G01020'])

# generating empty matrix
gene_correlation <- matrix(nrow = length(colnames(data)), 
                           ncol = length(colnames(data)))

colnames(gene_correlation) <- colnames(data)
rownames(gene_correlation) <- colnames(data)

##### calculate PCC value between two genes and insert into the correlation matrix

# test - 100x100

test_matrix <- matrix(nrow = 100, ncol = 100)
colnames(test_matrix) <- colnames(data)[1:100]
rownames(test_matrix) <- colnames(data)[1:100]

for (n in c(1:100)) {
  for (m in c(1:100)) {
    test_matrix[colnames(data)[n],colnames(data)[m]] <-
      cor(data[,n],data[,m])
  }
}

### make the real matrix

for (n in c(1:length(colnames(data))-1)) {
  for (m in c(2:length(colnames(data)))) {
    gene_correlation[colnames(data)[n],colnames(data)[m]] <-
      cor(data[,n],data[,m])
  }
}

