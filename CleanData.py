import pandas as pd

class CleanData(pd.DataFrame):
    '''
    A special type of Pandas dataframe.

    NOTE: for timepoint-based analyses, make separate CleanData objects for each timepoint's comparison
    '''
    def __init__(self, filenames):
        '''
        FILENAMES: a list of files containing the genes we want to analyze.
        TODO: improve description
        '''
        pd.DataFrame.__init__(self)
        self.dataframe = __extractdata__(filenames)


    def get_overexpressed(self):
        '''
        Returns a copy of self.df containing only overexpressed genes
        '''
        return self.df[self.df['expressionChange']=='overexpressed']

    def get_underexpressed(self):
        '''
        Returns a copy of self.df containing only underexpressed genes
        '''
        return self.df[self.df['expressionChange']=='underexpressed']
        
def __extractdata__(files):
    '''
    Given a list of files, extract all genes and their respective data. Filter genes 
    '''
    DEgenes = pd.DataFrame(columns=['geneNames', "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"])
    for file in files:
        try:
            frame = pd.read_csv(file, sep=' ', header=0, index_col=0)
            # Genes with low count numbers and outlier samples have a padj of "NA"

            # The genes with NA are the ones DESeq2 has filtered out. From DESeq2 manual: “The results function of the DESeq2 package 
            # performs independent filtering by default using the mean of normalized counts as a filter statistic. A threshold on the 
            # filter statistic is found which optimizes the number of adjusted p values lower than a [specified] significance level”.

            # I used the default False Discovery Rate (FDR) cutoff of 0.1 as the specified significance level for thresholding
            # adjusted p-values in the independent filtering step.
            frame.dropna(subset=['padj'], inplace=True)
            frame.reset_index(inplace=True)
            frame = frame.rename(columns = {'index':'geneNames'})

            # TODO: look up the default p_adj threshold for DESeq2

            # find overexpressed and underexpressed genes for each comparison
            # sort by lfc, keep only stuff above upper lfc_thr=-1.0 and below lower lfc_thr=1.0
            frame = frame[(frame['log2FoldChange'] < -1.0) | (frame['log2FoldChange'] > 1.0)]
            # sort by padj (lowest to highest), keep only stuff above pv_thr=0.05
            frame = frame[frame['padj'] < 0.05]
            frame.reset_index(inplace=True)
            # note which genes are over- and underexpressed
            numrows = frame.shape[0]
            emptycol = ['blank']*numrows
            frame['expressionChange'] = emptycol
            for row in range(numrows):
                if frame.at[row, 'log2FoldChange'] < 0:
                    frame.at[row, 'expressionChange'] = 'underexpressed'
                elif frame.at[row, 'log2FoldChange'] > 0:
                    frame.at[row, 'expressionChange'] = 'overexpressed'

            DEgenes = DEgenes.append(frame)

        except Exception as e:
            print(Exception, e)
            print(file)

    # if a gene is an outlier at multiple time points, only keep the first entry in the dataframe
    # TODO: should I keep duplicate with most significant padj value regardless of timepoint? Or just 
    # stick to the first timepoint for consistency?

    print(DEgenes.head())
    # check if duplicates are overexpressed at one timepoint and underexpressed at another
    # print(DEgenes[DEgenes.duplicated(subset='geneNames', keep=False)])

    DEgenes.drop_duplicates(subset='geneNames', inplace=True)
    return DEgenes




filepaths = ['/home/hannah/Dropbox/Temple/kulathinal-rotation/DE_results/FCGXXM-vs-FCGXXF-105-diffexp']

XXMvsXXF_DEGs = CleanData(filepaths)
print(XXMvsXXF_DEGs.head())