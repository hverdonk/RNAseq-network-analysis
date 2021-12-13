import pandas as pd
import os

class DEgenes():
    '''
    A special type of Pandas dataframe containing differential gene expression data.
    '''
    def __init__(self, DE_genes):
        '''
        DE_genes: a results table output by DEseq2 containing log2 fold change data for each gene and associated 
                  p-value for measuring differential expression.
        
        '''
        assert os.path.isfile(DE_genes)
        self.full_data = __extractfulldata__(DE_genes)
        self.DEGs = __extractdata__(DE_genes)

    def get_DEGs(self):
        '''
        Returns an ndarray of the gene names in self.DEGs
        '''
        return self.DEGs['geneNames']

    def get_DEGs_string(self):
        '''
        Returns a string of comma-separated gene names in self.DEGs
        Intended for easy search using the STRING database app in cytoscape
        '''
        return ','.join(self.genes())


    def get_overexpressed(self):
        '''
        Returns a new DEgenes object where object.DEGs is a copy of self.DEGs containing only overexpressed genes
        '''
        return DEgenes(self.DEGs[self.DEGs['expressionChange']=='overexpressed'])

    def get_underexpressed(self):
        '''
        Returns a new DEgenes object where object.DEGs is a copy of self.DEGs containing only underexpressed genes
        '''
        return DEgenes(self.DEGs[self.DEGs['expressionChange']=='underexpressed'])

    def reduce_data(self, data, cols):
        '''
        Returns a subset of either self.full_data or self.DEGs containing only the columns in COLS.
        DATA may be 'full' or 'deg'
        COLS is a list of column names to keep
        '''
        if data=='deg':
            return self.DEGs[cols]
        elif data=='full':
            return self.full_data[cols]
        else:
            print("Please specify which data to reduce: 'full' or 'deg'")
            
def __extractfulldata__(file):
    '''
    Given a file of differentially expressed genes (yielded by DESeq2), extract all genes and their respective data.
    '''
    df = pd.read_csv(file, sep=' ', header=0, index_col=0)
    df.reset_index(inplace=True)
    df = df.rename(columns = {'index':'geneNames'})
    
    return df
        
def __extractdata__(file, pthresh=0.05, logthresh=1.0):
    '''
    Given a file of differentially expressed genes (yielded by DESeq2), extract all genes and their respective data. Filter genes 
    by adjusted p-value (padj) and log2 fold change, keeping only those above the appropriate threshold (pthresh and logthresh, respectively)
    '''
    DEgenes = pd.read_csv(file, sep=' ', header=0, index_col=0)
    # Genes with low count numbers and outlier samples are prefiltered by DESeq2 and have a padj of "NA"
    DEgenes.dropna(subset=['padj'], inplace=True)
    DEgenes.reset_index(inplace=True)
    DEgenes = DEgenes.rename(columns = {'index':'geneNames'})

    # find differentially expressed genes for each comparison
    # sort by lfc, keep only stuff below lower -logthresh and above upper logthresh
    DEgenes = DEgenes[(DEgenes['log2FoldChange'] < -logthresh) | (DEgenes['log2FoldChange'] > logthresh)]
    # sort by padj (lowest to highest), keep only genes with padj values lower than pthresh
    DEgenes = DEgenes[DEgenes['padj'] < pthresh]
    DEgenes.reset_index(inplace=True)
    # note which genes are over- and underexpressed
    numrows = DEgenes.shape[0]
    emptycol = ['blank']*numrows
    DEgenes['expressionChange'] = emptycol
    for row in range(numrows):
        if DEgenes.at[row, 'log2FoldChange'] < 0:
            DEgenes.at[row, 'expressionChange'] = 'underexpressed'
        elif DEgenes.at[row, 'log2FoldChange'] > 0:
            DEgenes.at[row, 'expressionChange'] = 'overexpressed'

    return DEgenes

class ArachneData():
    '''
    Instances of ArachneData store gene expression data and transcription factors passed to the ARACHNe network inference algorithm.
    '''
    def __init__(self, data, factors):
        '''
        DATA: A standard gene expression matrix provided as a csv file where rows are genes and columns are samples.
              Each row contains a gene's transcript counts per sample.
        FACTORS: A csv file containing a list of transcription factors.
        '''
        assert os.path.isfile(data)
        assert os.path.isfile(factors)
        self.expression = pd.read_csv(data)
        self.TFs = pd.read_csv(factors).iloc[:, 0]
