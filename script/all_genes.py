import pandas as pd


cdrs_sep = ';'
gap_character = '.'

class TCR_Gene:
    def __init__( self, l ):
        # l comes from pandas dataframe.itertuples
        self.id = l.id
        self.organism = l.organism
        self.chain = l.chain
        self.region = l.region
        self.nucseq = l.nucseq
        self.alseq = l.aligned_protseq
        if pd.isna(l.cdrs):
            self.cdrs = []
            self.cdr_columns = []
        else:
            self.cdrs = l.cdrs.split(cdrs_sep)
            ## these are still 1-indexed !!!!!!!!!!!!!!
            self.cdr_columns = [ list(map( int,x.split('-'))) for x in l.cdr_columns.split(cdrs_sep) ]
        frame = l.frame
        #assert frame in ['+1','+2','+3','1','2','3']
        assert frame in [1,2,3] # now parsed by pandas, so string converted to int
        #self.nucseq_offset = int( frame[-1] )-1 ## 0, 1 or 2 (0-indexed for python)
        self.nucseq_offset = frame-1 ## 0, 1 or 2 (0-indexed for python)
        self.protseq = translation.get_translation( self.nucseq, f'+{frame}' )
        assert self.protseq == self.alseq.replace(gap_character,'')
        # sanity check
        if self.cdrs:
            assert self.cdrs == [ self.alseq[ x[0]-1 : x[1] ] for x in self.cdr_columns ]

