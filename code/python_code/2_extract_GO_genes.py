import pandas as pd

go_genes = pd.read_csv('GO_Groups.csv') 
# obtained from ShinyGO as GO_Gropus.csv by inputing the list of all genes from GSEA data
# http://bioinformatics.sdstate.edu/go/ ver. 0.76, Arabidopsis thaliana (not STRINGdb)

# select categories of interest and extract them into a list
go_genes_shortlist = go_genes.iloc[[47,49,67,71,72,79,87,90,92,98]] # GO categories of interest
go_genes_shortlist.to_csv('GO_Groups_shortlist.csv') # export shortlist of GO categories and represented genes

go_genes_shortlist_genelist = go_genes_shortlist['Genes'] # extract only the gene strings, separated by ' '

gene_lists = [[] for i in range(len(go_genes_shortlist))]
i = 0
gene = ''
for list in go_genes_shortlist_genelist:
    print(list)
    for c in list:
        if c != ' ':
            gene += c
        else:
            gene_lists[i].append(gene)
            gene = ''

    gene_lists[i].append(gene)
    gene = ''

    print(gene_lists[i])
    i += 1

redundant_list_index = 0
final_gene_list = []

for i in range(len(gene_lists)-1):
    a = set(gene_lists[i])
    for j in range(1, len(gene_lists)):
        b = set(gene_lists[j])
        temp = a.intersection(b)
        if i != j:
            print('Common between lists', i, 'and', j, ':', len(temp))
            if len(temp) != 0:
                print(temp)
            if len(temp) == len(gene_lists[j]):
                print(go_genes_shortlist['High level GO category'].iloc[i], 'includes',
                      go_genes_shortlist['High level GO category'].iloc[j], '!')
                redundant_list_index = j


# list 3 (activation of immune response) is redundant
# remove the redundant category from shortlisted datasets
gene_lists.pop(redundant_list_index)
go_genes_shortlist = go_genes_shortlist.reset_index()
go_genes_shortlist_genelist = go_genes_shortlist_genelist.reset_index()
go_genes_shortlist = go_genes_shortlist.drop(go_genes_shortlist.index[redundant_list_index])
go_genes_shortlist_genelist = go_genes_shortlist_genelist.drop(go_genes_shortlist_genelist.index[redundant_list_index])

# extract the shortlist of genes, avoiding repeated entries
final_gene_list = []

for list in gene_lists:
    for gene in list:
        if gene not in final_gene_list:
            final_gene_list.append(gene)
        else: print(gene, 'redundant')

all_genes_list = []

with open('gene_output.txt', 'r') as all_genes:
    for line in all_genes:
        all_genes_list.append(line[0:9])


# adding selected housekeeping genes for variability reference
housekeeping_genes = ['AT1G58050', 'AT1G59830', 'AT5G60390', 'AT1G07920', 'AT1G07930', 'AT1G07940']
# helicase, PP2A, 4x EF-1alpha. suggestions based on paper https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1203353/

for gene in housekeeping_genes:
    if gene in all_genes_list:
        print(gene, 'is in the full gene list')


# making the final list of genes, including housekeeping genes
final_gene_list = housekeeping_genes

for list in gene_lists:
    for gene in list:
        if gene not in final_gene_list:
            final_gene_list.append(gene)
        else: print(gene, 'redundant')

for i in range(len(gene_lists)):
    temp = pd.DataFrame(gene_lists[i])
    temp.to_csv(go_genes_shortlist['High level GO category'].iloc[i]+'GENE_LIST.csv')

# exporting file for GAPIT phenotype selection
with open('gene_output_for_gapit.csv', 'w') as gene_list:
    for gene in final_gene_list:
        gene_list.write(gene+',')
