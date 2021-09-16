import pandas as pd
import subprocess
import statsmodels.api as sm
import numpy as np
import math

'''
This function prcesses the gene file
Output is a one-row file for a gene
Each individual is in a column

Input file must have rowname
gene: gene ENSG ID of interest
start_col: column number which the gene exp value starts
gene_col: column name for the gene column
gene_start_col: column name for the gene start position
chr_col: column name for the gene chromosome
'''
def process_input(gene_file, vcf_file, cov_file, gene, start_col, gene_col, chr_col, gene_start_col):
    all_gene = pd.read_csv(gene_file, sep='[\t,]', header=0)  #sep='[\t,]' allows read in both , and tab delimited files'''
    gene=all_gene.loc[all_gene[gene_col]==gene,]
    
    gene_start=int(gene.loc[:,gene_start_col])
    chrom=int(gene.loc[:,chr_col])
    gene=gene.iloc[:,start_col:gene.shape[1]]
    
    start=int(gene_start-1e6)
    if start < 0:start = 0
    end=int(start+1e6)

    cmd='tabix '+ vcf_file + ' ' + str(chrom) + ':' + str(start) + '-' + str(end)
    s = subprocess.check_output(cmd, shell=True)
    s = s.decode().strip()
    s = s.split('\n')
    gt=[]
    for i in s:
        gt.append(i.split('\t')) 
    s1=pd.DataFrame(gt)
    info=s1.iloc[:,0:9]
    s1=s1.drop([0,1,2,3,4,5,6,7,8],axis=1)
    s1.index=info.iloc[:,2]

    s2= pd.DataFrame()
    for i in s1.columns:
        s2[i] = s1[i].apply(lambda x: x.split(':')[1])

    sample_ids = subprocess.check_output('/usr/local/bin/bcftools query -l {}'.format(vcf_file),    shell=True).decode().strip().split()
    s2.columns=sample_ids
    s3=s2[gene.columns]

    cov = pd.read_csv(cov_file, sep='\t', index_col=0, header=0)  
    cov=cov[gene.columns]

    return gene, s3, cov

'''This function takes the input from the previous function
   Fit linear model 
   Return beta and pvalues for the SNPs
'''
def lm_res(snps,gene,cov):
    res = pd.DataFrame(np.zeros([snps.shape[0],2], dtype=np.float32))
    res.index=snps.index
    res.columns=['beta','pval']                                          

    for i in range(snps.shape[0]):
        X=pd.concat([snps.iloc[i,].T, cov.T], axis=1)
        X = X.apply(pd.to_numeric)
        X = sm.add_constant(X)
        est = sm.OLS(pd.to_numeric(gene.T.iloc[:,0]), X).fit()
        res.iloc[i,0]=est.params[1]
        res.iloc[i,1]=est.pvalues[1]        
    return res
