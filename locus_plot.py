import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import math
import pandas as pd

class Locus:   
    def __int__(self, chrom, locus_start, locus_end):
        self.chr=chrom
        self.locus_start=locus_start
        self.locus_end=locus_end
    
    def get_genes(self, gene_info_file):
        '''return a dataframe where each row is a gene with gene_start, gene_end and strand'''
        gene_info = pd.read_csv(gene_info_file, sep='[\t,]', header=0)
        gene_info = gene_info.loc[gene_info['chr']==self.chr,]
        gene_info = gene_info.loc[gene_info['start']>= self.locus_start,]
        gene_info = gene_info.loc[gene_info['end']<= self.locus_end,]        
        if (gene_info.shape[0] < 1):
            print('No genes in this locus')
        else:
            return gene_info

def locus_plot(res, pval_colname, snp_pos, pos_colname, snps, snp, tick_spacing=100):
    res=pd.merge(res,snp_pos,left_on=res.index, right_on=snp_pos.index)
    res.index=res.key_0
    res=res.drop(columns=['key_0'])
    res['log10p']=res[pval_colname].apply(lambda x: -math.log(x,10))
    
    #calculate LD R2
    ld=pd.DataFrame(snps.corrwith(snps.loc[snp], axis=1, method='pearson')**2)
    ld.columns=['ld_R2']
    
    res=pd.merge(res, ld, left_on=res.index, right_on=ld.index)
    res.index=res.key_0
    res=res.drop(columns=['key_0'])

    res['ld_cat'] = pd.cut(res.ld_R2, [0,0.2,0.4,0.6,0.8,1], include_lowest=True)
    
    fig, ax = plt.subplots(1, 1) 
    ax.scatter(res[pos_colname], res['log10p'], alpha=0.5, c=res.ld_cat.cat.codes) 
    ax.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
    plt.xlabel('Genomic positions')
    plt.ylabel('pvalue')
    return plt 
