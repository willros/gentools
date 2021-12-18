import yaml
from pathlib import Path

class Configfile:
    '''Class that hold configfile'''
    def __init__(self, configfile):
        self._yaml = self.__read_yaml(configfile)
    
    def __read_yaml(self, file):
        with open(file, 'r') as f:
            return yaml.full_load(f)
        
    def __getitem__(self, item):
        return self._yaml[item]
    
# Include all folders here!
class CreateFolders:
    '''Class that creates and holds correct folder structure'''
    def __init__(self, config_file):
        #set up
        self.config = Configfile(config_file)
        self.working_directory = Path('.').absolute() if self.config['working_directory'] == 'this_folder' else Path(self.config['working_directory']).absolute()
        self.raw_reads = Path(self.config['raw_reads']).absolute()
        #umi_tools
        self.umi_tools_processed = Path(self.working_directory, 'umi_tools', 'processed')
        self.umi_tools_processed_dedup = Path(self.working_directory, 'umi_tools', 'processed', 'dedup')
        self.umi_tools_log = Path(self.working_directory, 'umi_tools', 'log')
        self.umi_tools_log_dedup = Path(self.working_directory, 'umi_tools', 'log', 'dedup')
        #cutadapt
        self.cutadapt_processed = Path(self.working_directory, 'cutadapt', 'processed')
        self.cutadapt_log = Path(self.working_directory, 'cutadapt', 'log')
        #bowtie2
        self.bowtie2_processed_aligned = Path(self.working_directory, 'bowtie2', 'processed', 'aligned')
        self.bowtie2_processed_unaligned = Path(self.working_directory, 'bowtie2', 'processed', 'unaligned')
        self.bowtie2_log = Path(self.working_directory, 'bowtie2', 'log')
        #featureCounts
        self.feature_counts = Path(self.working_directory, 'results', 'featureCounts')
        #results
        self.results = Path(self.working_directory, 'results')
        #list of all directories
        self.__directory_list = [self.umi_tools_processed, self.umi_tools_log, self.umi_tools_processed_dedup, 
                                 self.umi_tools_log_dedup, self.cutadapt_processed, self.cutadapt_log, 
                                self.bowtie2_processed_aligned, self.bowtie2_processed_unaligned,
                                self.bowtie2_log, self.results, self.feature_counts]
                
    def make_directories(self):
        for i in self.__directory_list:
            i.mkdir(exist_ok=True, parents=True)
      
