acc_names <- t(read.csv('acc_excerpt.csv')) 
files <- list.files(path="/srv/kenlab/ejansone/SNP_Table", # location of folder with all chunks to be added
                    pattern="*.csv.gz", full.names=TRUE, recursive=FALSE)
cat(length(files), 'files found')

print('Opening initial table...')

SNP_table <- read.table('SNP_table_init.csv.gz', header = FALSE, sep= ',') # chunk no. 000 renamed
SNP_table <- cbind(acc_names, SNP_table)
cat('Initial table with dimensions', dim(SNP_table), 'opened, saving as RDS...')
saveRDS(SNP_table, 'SNP_Table_init.rds')
print('Initial SNP table saved as SNP_Table_init.rds')

SNP_table <- readRDS('SNP_Table_init.rds')

for (file in files) {
  cat('Now adding', file, '\n')
  temp <- read.table(file, header=FALSE, sep= ',') # load file
  cat('adding dimensions', dim(temp), 'to the SNP table..')
  SNP_table <- cbind(SNP_table, temp)
  cat('current table dimensions', dim(SNP_table))
}

print('Saving SNP table into as SNP_Table_full_AraGWAS.rds ...')
saveRDS(SNP_table, 'SNP_Table_full_AraGWAS.rds')

print('Table saved successfully.')
