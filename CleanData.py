import pandas as pd

# pd.DataFrame
class CleanData():
    '''
    A special type of Pandas dataframe.

    NOTE: for timepoint-based analyses, make separate CleanData objects for each timepoint's comparison
    '''
    def __init__(self, filenames):
        '''
        FILENAMES: a list of files containing the genes we want to analyze.
        TODO: improve description
        '''
        self.df = __extractdata__(filenames)


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
            # Genes with low count numbers and outlier samples are prefiltered by DESeq2 and have a padj of "NA"
            frame.reset_index(inplace=True)
            frame = frame.rename(columns = {'index':'geneNames'})

            # TODO: look up the default p_adj threshold for DESeq2

            # find differentially expressed genes for each comparison
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

    # TODO: should I keep duplicate with most significant padj value regardless of timepoint? Or just 
    # stick to the first timepoint for consistency? Will I even need to drop duplicates in the final data, 
    # or will each timepoint be its own thing?

    # check if duplicates are overexpressed at one timepoint and underexpressed at another
    # print("Duplicates:")
    # temp = DEgenes[DEgenes.duplicated(subset='geneNames', keep=False)]
    # grouped = temp.groupby('geneNames')

    # for name,group in grouped:
    #     print(name)
    #     print(group)

    # TODO: at least one gene is overexpressed at one timepoint and underexpressed at another. Even for genes with similar trends,
    # how should I deal with different l2fc or different padj values? Which gene should I keep, if any?
    # PERHAPS: keep duplicate with most significant padj value 
        # remove all dups from one df, keep all dups in another frame for processing. Then at the end, merge frames together
        # write helper function to handle dups df - keep only gene with highest padj value
    '''
    Gm22711
        geneNames   baseMean  log2FoldChange     lfcSE      stat    pvalue      padj    index expressionChange
    122   Gm22711  17.088073        1.729882  0.597521  2.829907  0.004656  0.032593  13156.0    overexpressed
    37    Gm22711  17.088073       -1.949715  0.550194 -3.456516  0.000547  0.035195  13156.0   underexpressed
    '''

    # if a gene is an outlier at multiple time points, only keep the first entry in the dataframe
    DEgenes.drop_duplicates(subset='geneNames', inplace=True)
    DEgenes.reset_index(inplace=True)
    return DEgenes
