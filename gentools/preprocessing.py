import re
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from .config import Configfile, CreateFolders
from pathlib import Path
 
class ProgramGraphs:
    '''Super class for making graphs and so on'''
    def __init__(self, config_file):
        self.config = Configfile(config_file)
        self.folders = CreateFolders(config_file)
        self.program = 'program'
        self.pattern_pre = ''
        self.pattern_post = ''
        
    def list_of_files(self):
        folders_dict = {'umi_tools_extract': self.folders.umi_tools_log, 
                        'umi_tools_dedup': self.folders.umi_tools_log_dedup, 
                        'cutadapt': self.folders.cutadapt_log, 
                        'bowtie2': self.folders.bowtie2_log}
        folder = folders_dict[self.program]
        return sorted([file for file in folder.iterdir() if file.is_file()])
        
    def extract_info(self):
        files = self.list_of_files()
        pre_processing = 0
        post_processing = 0
        dataframe_dict = {'name': [], 'pre_filtering': [], 'post_filtering': []}
        
        for file in files:
            with open(file, 'r') as f:
                content = f.read()
                pre = re.findall(self.pattern_pre, content)
                pre = int(pre[0].replace(',', ''))
                post = re.findall(self.pattern_post, content)
                post = int(post[0].replace(',', ''))
                pre_processing += pre
                post_processing += post
                dataframe_dict['name'].append(file.stem)
                dataframe_dict['pre_filtering'].append(pre)
                dataframe_dict['post_filtering'].append(post)
        
        name = f'{self.program}.csv'
        csv_file = Path(self.folders.results, name)
        pd.DataFrame(dataframe_dict).to_csv(csv_file, index=False)
        return pre_processing, post_processing
    
    def make_pie_graph(self):
        explode = (0.1,0)
        pre, post = self.extract_info()
        fig_name = self.folders.results / f'{self.program}.pdf'
        
        plt.pie((pre - post, post), explode=(0.1, 0), shadow=False, startangle=120, labels=('', 'Survived filtering'), autopct='%1.1f%%')
        plt.title(f'Procent of reads filtered by {self.program}')
        plt.savefig(fig_name)
        plt.close()
        
        
class UmitoolsExtractGraphs(ProgramGraphs):
    '''Returns information about umi_tools'''
    def __init__(self, config_file):
        super().__init__(config_file)
        self.program = 'umi_tools_extract'
        self.pattern_pre = re.compile(r'INFO Input Reads: (\d+)')
        self.pattern_post = re.compile(r'INFO Reads output: (\d+)')

        
class CutadaptGraphs(ProgramGraphs):
    '''Returns information about cutadapt'''
    def __init__(self, config_file):
        super().__init__(config_file)
        self.program = 'cutadapt'
        self.pattern_pre = re.compile(r'Total reads processed: +(\d+,\d+)')
        self.pattern_post = re.compile(r'Reads written \(passing filters\): +(\d+,\d+)')
        
        
class Bowtie2Graphs(ProgramGraphs):
    '''Returns information about bowtie2'''
    def __init__(self, config_file):
        super().__init__(config_file)
        self.program = 'bowtie2'
        self.pattern_pre = re.compile(r'(\d+) reads; of these:')
        self.pattern_post = re.compile(r'(\d*) \(\d*\.\d{2}%\) aligned exact')

# This works    
class UmitoolsDedupGraphs(ProgramGraphs):
    '''Returns information about umi_tools'''
    def __init__(self, config_file):
        super().__init__(config_file)
        self.program = 'umi_tools_dedup'
        self.pattern_pre = re.compile(r'INFO Reads: Input Reads: (\d+)')
        self.pattern_post = re.compile(r'INFO Number of reads out: (\d+)')
        

class FeatureCountsGraphs(ProgramGraphs):
    '''Returns information about featureCounts'''
    def __init__(self, config_file):
        super().__init__(config_file)
        self.program = 'featureCounts'
        self.pattern_pre = re.compile(r'Assigned\t(\d+.*)') 
        self.pattern_post = re.compile(r'Unassigned_NoFeatures\t(\d+.*)')
    
    def extract_info(self):
        names = sorted([file.stem for file in self.folders.raw_reads.iterdir() if file.is_file()])
        summary = self.folders.feature_counts / 'count_matrix.txt.summary'
        
        with open(summary, 'r') as f:
            content = f.read()
            assigned = re.findall(self.pattern_pre, content)
            assigned = "".join(assigned).split('\t')
            assigned = [int(number) for number in assigned]
            not_assigned = re.findall(self.pattern_post, content)
            not_assigned = "".join(not_assigned).split('\t')
            not_assigned = [int(number) for number in not_assigned]
            pre = np.array(assigned) + np.array(not_assigned)

        pre_processing = sum(pre)
        post_processing = sum(assigned)   
        pd.DataFrame({'sample':names, 'mapped_reads': pre, 'assigned_reads': assigned})

        name = f'{self.program}.csv'
        csv_file = Path(self.folders.results, name)
        pd.DataFrame({'sample':names, 'mapped_reads': pre, 'assigned_reads': assigned}).to_csv(csv_file, index=False)
        
        return pre_processing, post_processing
