import pandas as pd

# pd.DataDEgenes
class CleanData():
    '''
    A special type of Pandas dataframe containing differential gene expression data.
    '''
    def __init__(self, filename='', dataframe=None):
        '''
        FILENAME: a file containing the genes we want to analyze.
        TODO: improve description
        '''
        if dataframe:
            self.df=dataframe
        else:
            self.df = __extractdata__(filename)


    def genes(self):
        '''
        Returns an ndarray of the gene names in self.df
        '''
        return self.df['geneNames']

    def genes_string(self):
        '''
        Returns a string of comma-separated gene names in self.df
        Intended for easy search using the STRING database app in cytoscape
        '''
        return ','.join(self.genes())


    def get_overexpressed(self):
        '''
        Returns a new CleanData object where object.df is a copy of self.df containing only overexpressed genes
        '''
        return CleanData(self.df[self.df['expressionChange']=='overexpressed'])

    def get_underexpressed(self):
        '''
        Returns a new CleanData object where object.df is a copy of self.df containing only underexpressed genes
        '''
        return CleanData(self.df[self.df['expressionChange']=='underexpressed'])

    def reduce_data(self, *args):
        '''
        Returns a subset of self.dataframe containing only the columns in ARGS.
        '''
        return self.df[list(args)]
        
def __extractdata__(file, pthresh=0.05, logthresh=1.0):
    '''
    Given a file of differentially expressed genes (yielded by DESeq2), extract all genes and their respective data. Filter genes 
    by adjusted p-value (padj) and log2 fold change, keeping only those above the appropriate threshold (pthresh and logthresh, respectively)
    '''
    DEgenes = pd.read_csv(file, sep=' ', header=0, index_col=0)
    # Genes with low count numbers and outlier samples are prefiltered by DESeq2 and have a padj of "NA"
    DEgenes.reset_index(inplace=True)
    DEgenes = DEgenes.rename(columns = {'index':'geneNames'})

    # TODO: look up the default p_adj threshold for DESeq2

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
