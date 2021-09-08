import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import math
import pandas as pd

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
