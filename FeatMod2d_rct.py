"""
2D Feature Model Reaction Module.

Pandas data structure is used.
"""

from FeatMod2d_readin import sp_name, mat_name

import pandas as pd

class React():
    """Define reaction structure."""
    
    def __init__(self, sp_name=list(), mat_name=list()):
        """Init the data structure."""
        self.sp_name = sp_name
        self.mat_name = mat_name
        df = list()
        col_name=['sp', 'mat', 'prod1', 'prod2', 'type', 'prob']
        for sp in sp_name:
            for mat in mat_name:
                df.append([sp, mat, sp, mat, 'rflct', 1.00])
        self.df = pd.DataFrame(df, columns=col_name)
    
    def readin(self, csv_file=None):
        """Read in reactions from csv file."""
        df_readin = pd.read_csv(csv_file)  
        self.df = self.df.append(df_readin)
        self.df = self.df.sort_values(by=['sp','mat'])
        self.df = self.df.reset_index(drop=True)
    
    # Normalization is not necessary for random pick
    def _norm_prob(self):
        """Normalize the prob for same reactants."""
        for sp in self.sp_name:
            for mat in self.mat_name:
                temp_df = self.df[
                    (self.df['sp'] == sp) & (self.df['mat'] == mat)
                    ]
                temp_sum = temp_df['prob'].sum()
                print(self.df.loc[temp_df.index]['prob'])
                # self.df.loc[temp_df.index]['prob'] = \
                #     self.df.loc[temp_df.index]['prob'].apply(lambda x:x/temp_sum)
                print(temp_sum)
                self.df.loc[temp_df.index]['prob'].apply(lambda x:x/temp_sum)
                print(self.df.loc[temp_df.index]['prob'])
                

if __name__ == '__main__':
    react = React(sp_name=sp_name, mat_name=mat_name)
    react.readin('Chem.csv')
    print(react.df)
    
    
