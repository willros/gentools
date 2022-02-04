import subprocess
from gentools.config import Configfile, CreateFolders
import pandas as pd
from pathlib import Path

#WORKS
class DifferentialExpression:
    '''Class that handle differential expression with gentools. Relies on DESeq2.'''
    
    def __init__(self, config_file):
        self.config = Configfile(config_file)
        self.folders = CreateFolders(config_file)
        self.program = 'DESeq2'
    
    def _create_df(self):
        
        counts = self.folders.feature_counts / 'count_matrix.csv'
        df = pd.read_csv(counts)
        samples = df.columns.to_list()[4:]
        coldata = pd.DataFrame({'samples':samples, 'condition': 'NA'})
        return coldata

    def _create_coldata(self):
        coldata = self._create_df()
        
        for params in self.config['DESeq2']:
            for key, value in params.items():
                if key == 'coldata':
                    for cond in value:
                        for key, value in cond.items():
                            start, _, end = value.split(' ')
                            start = int(start) - 1
                            end = int(end)
                            coldata['condition'].iloc[start:end] = key
        
        self.coldata = self.folders.results / 'coldata.csv'
        coldata.to_csv(self.coldata, index=False)
                            
    def run_command(self):
        #create the coldata file
        self._create_coldata()
        
        counts = self.folders.feature_counts / 'count_matrix.csv'
        de_genes = self.folders.results / 'de_genes.csv'
       
        #descript = os.path.join(path, 'DEanalysis.R')
        path = Path(__file__, '..').resolve().parent
        R_script = path / 'scripts/DESeq2.R'
        
        subprocess.call(['Rscript', R_script, counts, self.coldata, de_genes])
        
