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
        self.working_directory = Path('.') if self.config['working_directory'] == 'this_folder' else Path(self.config['working_directory'])
        self.raw_reads = Path(self.config['raw_reads'])
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
      
    
# function that creates default config file. FIX THIS 
def create_default_config(work_dir):
    path = Path(work_dir, 'gentool_config.yaml')
    content = '### test comment\nraw_reads: test_data\nworking_directory: testar_lite\ngenomes: ../scripts/genomes\n\numi_tools_extract:\n- mode: extract\n- extract-method: regex\n- bc-pattern: (?P<discard_1>.*)(?P<umi_1>[ACT]{8}CA)\n\ncutadapt:\n- adapter: TGGAATTCTCGGGTGCCAAGG\n- minimum-length: 18\n- maximum-length: 41\n- error-rate: 0.1\n- overlap: 1\n\nbowtie2:\n- k: 100\n- local: very-sensitive-local\n- x: indices/bowtie2/rna_central/rna_central\n\numi_tools_dedup:\n- mode: dedup\n- method: unique\n\n\norder:\n- cutadapt\n- umi_tools_extract\n- bowtie2'
    
    with open(path, 'w') as file:
        file.write(content)
