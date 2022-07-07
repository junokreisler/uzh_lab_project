import h5py
import pandas as pd # only for exports and descriptive statistics
import numpy as np
import time
import wget
from matplotlib import pyplot as plt # only for histograms

# import GSEA data from url
GSEA_url = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE80nnn/GSE80744/suppl/GSE80744_ath1001_tx_norm_2016-04-21-UQ_gNorm_normCounts_k4.tsv.gz'
GSEA_filename = wget.download(GSEA_url) # downloads the file into the project directory, var contains filename

GSEA_data = pd.read_csv(GSEA_filename, delimiter='\t')

GSEA_accs = list(GSEA_data.columns)[1:]

############################ import AraGWAS data
data = h5py.File('4.hdf5', 'r') # hdf5 file obtained, unzipped from https://aragwas.1001genomes.org/api/genotypes/download
# data.keys() => ['accessions', 'positions', 'snps'];
# snps (10709466, 2029)

############################
# finding common accessions between AraGWAS and GSEA
############################

accessions = [0 for _ in range(0,len(data['accessions']))]

for i in range(0, len(accessions)):
    accessions[i] = int(data['accessions'][i])
    accessions[i] = 'X' + str(accessions[i])

common_accessions = np.array(list(set(accessions).intersection(GSEA_accs)))

acc_index = []

for i in range(0, len(accessions)):
    if accessions[i] in common_accessions:
        acc_index.append(i)

### export data for GAPIT (SNP indices and accessions + their indices)
acc_excerpt = pd.DataFrame([acc_index,common_accessions])
acc_excerpt.to_csv('acc_excerpt.csv', header = False, index = False)
############################
# trimming minor alleles at two-sided 5%
############################

# splitting all SNPs into chunks of 10% because my laptop can't handle calculating the full thing at once
chunks = [i*1070947 for i in range(0,11)]
chunks[-1] = len(data['positions'])

# initializing the SNP count array
SNP_counts = np.empty((0,1))

# filling the SNP count array with values, takes 20-25min
i_start = time.perf_counter()
for i in range(0,len(chunks)-1):
    print('Calculating sums for SNPs', chunks[i], 'to', chunks[i+1],'...')
    SNP_counts = np.append(SNP_counts, np.sum(data['snps'][chunks[i]:chunks[i+1],acc_index], axis = 1))
    i_stop = time.perf_counter()
    print('Sums for SNPs', chunks[i], 'to', chunks[i+1],'were calculated in', round((i_stop-i_start)/60,4), 'minutes')
    i_start = time.perf_counter()

# since sums of 0 and 1 mean that there is very little variance for a locus within all accessions,
# these loci can be excluded.

pd.DataFrame(SNP_counts).hist(bins=100)
plt.savefig('Hist_raw_counts_all_SNP.png')

SNP_diverse_ind = np.where(SNP_counts > 1)
SNP_diverse_ind_int = np.array(SNP_diverse_ind).astype(int)
SNP_counts_wo_zeros_ones = SNP_counts[SNP_diverse_ind]
SNP_counts_wo_zeros_ones_int = SNP_counts_wo_zeros_ones.astype(int)

pd.DataFrame(SNP_counts_wo_zeros_ones).hist(bins=100)
plt.savefig('Hist_raw_counts_no_01_SNP.png')

cutoff_bottom, cutoff_top = 0.05, 0.95

# check which SNPs to cut and which to keep
delete_SNP = []
keep_SNP = []

j = 0 # for performance timing

for i in list(SNP_diverse_ind)[0]:
    j += 1
    i_start = time.perf_counter()
    # obtaining the original index of SNP positions to keep
    if SNP_counts[i]/681 <= cutoff_bottom:
        delete_SNP.append(i)
    if SNP_counts[i]/681 >= cutoff_top:
        delete_SNP.append(i)

    if (SNP_counts[i]/681 > cutoff_bottom) and (SNP_counts[i]/681 < cutoff_top):
        keep_SNP.append(i)

    i_stop = time.perf_counter()
    if j % len(SNP_diverse_ind_int) == 10:
        print(round(j / i,2),'% done in', round(i_stop-i_start, 2), 'seconds.')

#pd.DataFrame(SNP_counts[keep_SNP]).hist(bins=100)
#plt.savefig('Hist_raw_counts_remaining_SNP.png')

SNP_counts_wo_zeros_ones_tails = SNP_counts[keep_SNP]

##################################
# make marker information table
##################################
# ID = indices of SNPs (keep_SNP) /// Chr = Chromosome number /// Pos = Index value (data['positions'][keep_SNP][i]

chr = 5
pos_all = np.array(data['positions'])
pos_excerpt = pos_all[keep_SNP]
chr_of_pos = pos_excerpt

# find chromosome starting indices, based on original / MAF-untrimmed data
new_chr_ind = []
for i in range(0, len(data['positions'])):
    if data['positions'][i]<data['positions'][i-1]:
        print('Chromosome change at index', i)
        new_chr_ind.append(i)

# Chr 1: 0    Chr 2: 2597735   Chr 3: 4466530   Chr 4: 6660782   Chr 5: 8427786

new_chr_ind = [0, 2597735, 4466530, 6660782, 8427786, len(pos_all)]


for i in range(len(new_chr_ind)-1,-1,-1):
    print(i)
    chr_of_pos[np.where(np.array(keep_SNP) < new_chr_ind[i])] = chr
    chr -= 1

np.savetxt('marker_information_ID.csv.gz', np.array(keep_SNP), delimiter=',')
np.savetxt('marker_information_CHR.csv.gz', chr_of_pos, delimiter=',')
np.savetxt('marker_information_POS.csv.gz', pos_excerpt.astype(int), delimiter=',')

###################################
# exporting the SNP table in 100 chunks (R might crash trying to load larger chunks)
###################################

def proper_number(number): # to ensure ordered looping through SNPs
    if (number // 10 == 0):
        return '00' + str(number)
    elif (number // 10 <= 9):
        return '0' + str(number)

# 1) creating the chunks and filling them with required indices
chunks_SNP = [i for i in range(0, len(keep_SNP), int(len(keep_SNP) / 100))]
chunks_SNP[-1] = len(keep_SNP)
chunks_SNP_keep = [keep_SNP[chunks_SNP[i]:chunks_SNP[i+1]] for i in range(0, len(chunks_SNP)-1)]

# 2) extracting the chunks, takes 2-3h !!!
for i in range(0, len(chunks_SNP_keep)):
    print('Accessing chunk', i+1,'...')
    i_start = time.perf_counter()

    # selecting the selected SNPs of a chunk for all accessions
    snp_all_accs = np.array(data['snps'][chunks_SNP_keep[i]])
    i_stop = time.perf_counter()
    print('Obtained chunk', i+1, 'in', round((i_stop-i_start)/60,2), 'minutes; trimming accessions now...')

    # trimming accessions
    snp_accs = snp_all_accs[:,acc_index]

    # saving new chunk
    print('Saving chunks', i+1, 'as SNP_table_'+proper_number(i)+'.csv.gz')
    np.savetxt('SNP_table_100_'+proper_number(i)+'.csv.gz', snp_accs.T, delimiter=',') # saving for GAPIT format

    # snp_excerpt = np.append(snp_excerpt, snp_all_accs[:,acc_index])
    i_stop = time.perf_counter()
    print('SNP chunk', i+1, 'out of', len(chunks_SNP_keep),'was extracted in a total of',
          round((i_stop-i_start)/60,2), 'minutes')
