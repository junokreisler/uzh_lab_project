#########################################
# size-dependent reproduction in plants #
#########################################

# load data
d = read.csv("./test/plant_data.csv", header=TRUE)

# standard linear model
res_lm = lm(dry_weight~height_cm, data=d)
summary(res_lm)

# Poisson GLM
res_glm = glm(seeds~height_cm, data=d, family="poisson")
summary(res_glm)

# make .pdf figure and save it in ./figure
pdf("./figure/plant_data_analysis.pdf", width=8, height=4)
par(mfcol=c(1,2))
plot(d$height_cm, d$dry_weight, pch=20, las=1, ylab="dry weight (g)", xlab="plant height (cm)", bty="l")
abline(a=res_lm$coefficients[1], b=res_lm$coefficients[2], lty=2)

plot(d$height_cm, d$seeds, pch=20, las=1, ylab="No. of seeds", xlab="plant height (cm)", bty="l")
curve(exp(res_glm$coefficients[1]+x*res_glm$coefficients[2]), lty=2, add=TRUE)
dev.off()

########################################

# Read count data (normalized) from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80744
# RNA-seq profiling of 728 Arabidopsis thaliana accessions
data <- read.table('./data/GSE80744_ath1001_tx_norm_2016-04-21-UQ_gNorm_normCounts_k4.tsv/GSE80744_ath1001_tx_norm_2016-04-21-UQ_gNorm_normCounts_k4.tsv',
                   header = TRUE)

data <- as.matrix(data)

data_transposed <- t(data)
data_transposed <- as.data.frame(data_transposed)
colnames(data_transposed) <- data_transposed[1,]
data_transposed <- data_transposed[-1,]
# the numbers are actually characters. so gotta change them to numbers
for (i in c(1:length(colnames(data_transposed)))) {
  data_transposed[,i] <- as.numeric(data_transposed[,i])
}
class(data_transposed[,3513])
# testing on one gene pair
plot(data_transposed$AT1G01010,data_transposed$AT1G01020, pch=3)
cor(data_transposed$AT1G01010,data_transposed$AT1G01020)


# generating empty matrix
gene_correlation <- matrix(nrow = length(colnames(data_transposed)), 
                           ncol = length(colnames(data_transposed)))

colnames(gene_correlation) <- colnames(data_transposed)
rownames(gene_correlation) <- colnames(data_transposed)

##### calculate PCC value between two genes and insert into the correlation matrix

# test - 100x100

test_matrix <- matrix(nrow = 100, ncol = 100)
colnames(test_matrix) <- colnames(data_transposed)[1:100]
rownames(test_matrix) <- colnames(data_transposed)[1:100]

for (n in c(1:100)) {
  for (m in c(1:100)) {
    test_matrix[colnames(data_transposed)[n],colnames(data_transposed)[m]] <-
      cor(data_transposed[,n],data_transposed[,m])
  }
}

### make the real matrix

for (n in c(1:length(colnames(data_transposed))-1)) {
  for (m in c(2:length(colnames(data_transposed)))) {
    gene_correlation[colnames(data_transposed)[n],colnames(data_transposed)[m]] <-
      cor(data_transposed[,n],data_transposed[,m])
  }
}

