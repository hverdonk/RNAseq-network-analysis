import pandas as pd

class DataCleaner(pd.DataFrame):
    '''
    A special type of Pandas dataframe.
    '''
    def __init__(filenames):
        '''
        FILENAMES: a list of files containing the genes of
        '''
        
    def __extractdata__(files):
        for file in files:
            with open(file, 'r') as f:
                pass
