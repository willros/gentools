import yaml
from pathlib import Path
import subprocess
from .config import Configfile, CreateFolders
import pandas as pd
import os
import re
from datar.all import select, contains

#parallellization for umi_tools
from joblib import Parallel, delayed

class ProgramCommando:
    ''' Superclass for different programs that can be run in gentools '''
    def __init__(self, config_file):
        self.config = Configfile(config_file)
        self.folders = CreateFolders(config_file)
        self.threads = str(self.config['threads'])
    
    def new_out_file(self, read, new_dir, suffix=None):
        suffix = read.suffix if not suffix else suffix
        return new_dir / (read.stem + '_' + self.program + suffix)
    
    def reads_in(self, folder, suffix=None):
        folder = self._reads_in_dict(folder)
        if suffix:
            return sorted([read for read in folder.iterdir() if read.suffix == suffix])
        return sorted([read for read in folder.iterdir() if read.is_file()])
    
    def _reads_in_dict(self, folder):
        input_folder_dict = {'raw': self.folders.raw_reads,
                            'umi_tools_extract': self.folders.umi_tools_processed,
                            'umi_tools_dedup': self.folders.umi_tools_processed_dedup,
                            'cutadapt': self.folders.cutadapt_processed,
                            'bowtie2_aligned': self.folders.bowtie2_processed_aligned,
                            'bowtie2_unaligned': self.folders.bowtie2_processed_unaligned,
                            'featureCounts': self.folders.feature_counts,
                            'fastp': self.folders.fastp_processed}
        return input_folder_dict[folder]
    
    def extract_parameters(self):
        pass
    
    def create_command(self):
        pass
    
    def run_command(self):
        pass

# WORKS LIKE A CHARM! :D 
class CutadaptCommando(ProgramCommando):
    def __init__(self, config_file):
        super().__init__(config_file)
        self.program = 'cutadapt'
       
    def extract_parameters(self):
        command = [self.program, '-j', self.threads]
        for params in self.config[self.program]:
            for key, value in params.items():
                if key == 'input':
                    self.input_files = self.reads_in(value)
                elif key == 'trimmed_only?':
                    if value == 'Y':
                        command.append('--trimmed-only')
                elif len(key) == 1:
                    command.append(f'-{key}')
                    command.append(str(value))
                else:
                    command.append(f'--{key}')
                    command.append(str(value))             
        return command
        
    def create_command(self):
        # to get the attribute self.input_files, self.extract_parameters needs to be run
        self.extract_parameters()
        list_of_commands = []
        log_files = []
        for read in self.input_files:
            command = self.extract_parameters()
            read_out = self.new_out_file(read, self.folders.cutadapt_processed)
            log_file = self.folders.cutadapt_log / f'{read.stem}.log'
            log_files.append(log_file)
            command.append('-o')
            command.append(read_out)
            command.append(read)
            list_of_commands.append(command)
        return zip(list_of_commands, log_files)
             
    def run_command(self):
        for command, logfile in self.create_command():
            with open(logfile, 'w') as file:
                subprocess.call(command, stdout=file)
 

 #WORKS
class UmitoolsExtractCommando(ProgramCommando):
    def __init__(self, config_file):
        super().__init__(config_file)
        self.program = 'umi_tools_extract'
    
    def extract_parameters(self):
        command = ['umi_tools']
        for params in self.config[self.program]:
            for key, value in params.items():
                if key == 'input':
                    self.input_files = self.reads_in(value)
                elif len(key) == 1:
                    command.append(f'-{key}')
                    command.append(str(value))
                elif key == 'mode':
                    command.append(value)
                elif key == 'extract-method':
                    command.append(f'--extract-method={value}')  
                elif key == 'bc-pattern':
                    command.append(f'--bc-pattern={value}') 
                else:
                    command.append(f'--{key}')
                    command.append(str(value))             
        return command
        
    def create_command(self):
        # to get the attribute self.input_files, self.extract_parameters needs to be run
        self.extract_parameters()
        list_of_commands = []
        for read in self.input_files:
            command = self.extract_parameters()
            read_out = self.new_out_file(read, self.folders.umi_tools_processed)
            log_file = self.folders.umi_tools_log / f'{read.stem}.log'
            command.append('-I')
            command.append(read)
            command.append('-S')
            command.append(read_out)
            command.append('-L')
            command.append(log_file)
            list_of_commands.append(command)
        return list_of_commands
            
    def _run_command(self, command):
        subprocess.call(command)

    def run_command(self):
        #parallellizes umi_tools
        commands = self.create_command()
        Parallel(n_jobs=int(self.threads))(delayed(self._run_command)(command) for command in commands)

        
#WORKS
class UmitoolsDedupCommando(ProgramCommando):
    def __init__(self, config_file):
        super().__init__(config_file)
        self.program = 'umi_tools_dedup'
    
    def extract_parameters(self):
        command = ['umi_tools']
        for params in self.config[self.program]:
            for key, value in params.items():
                if key == 'input':
                    self.input_files = self.reads_in(value, suffix='.bam')
                elif len(key) == 1:
                    command.append(f'-{key}')
                    command.append(str(value))
                elif key == 'mode':
                    command.append(value)
                elif key == 'method':
                    command.append(f'--method={value}')  
                else:
                    command.append(f'--{key}')
                    command.append(str(value))             
        return command
        
    def create_command(self):
        # to get the attribute self.input_files, self.extract_parameters needs to be run
        self.extract_parameters()
        list_of_commands = []
        for read in self.input_files:
            command = self.extract_parameters()
            read_out = self.new_out_file(read, self.folders.umi_tools_processed_dedup)
            log_file = self.folders.umi_tools_log_dedup / f'{read.stem}.log'
            command.append('-I')
            command.append(read)
            command.append('-S')
            command.append(read_out)
            command.append('-L')
            command.append(log_file)
            list_of_commands.append(command)
        return list_of_commands
    
    def _run_command(self, command):
        subprocess.call(command)

    def run_command(self):
        #parallellizes umi_tools
        commands = self.create_command()
        Parallel(n_jobs=int(self.threads))(delayed(self._run_command)(command) for command in commands)

        #index with samtools index
        files = [file for file in self.folders.umi_tools_processed_dedup.iterdir() if file.suffix == '.bam']
        for file in files:
            subprocess.call(['samtools', 'index', '-@', self.threads, file])

            
# WORKS LIKE A CHARM! :D 
class Bowtie2Commando(ProgramCommando):
    def __init__(self, config_file):
        super().__init__(config_file)
        self.program = 'bowtie2'
    
    def extract_parameters(self):
        command = [self.program, '-p', self.threads]
        for params in self.config[self.program]:
            for key, value in params.items():
                if key == 'input':
                    self.input_files = self.reads_in(value)
                elif len(key) == 1:
                    command.append(f'-{key}')
                    command.append(str(value))
                elif key == 'local':
                    command.append(f'--{key}')
                    command.append(f'--{value}') 
                elif key == 'end-to-end':
                    command.append(f'--{key}')
                    command.append(f'--{value}')   
                else:
                    command.append(f'--{key}')
                    command.append(str(value))             
        return command
        
    def create_command(self):
        # to get the attribute self.input_files, self.extract_parameters needs to be run
        self.extract_parameters()
        list_of_commands = []
        log_files = []
        list_of_out = []
        for read in self.input_files:
            command = self.extract_parameters()
            read_aligned = self.new_out_file(read, self.folders.bowtie2_processed_aligned, suffix='.bam')
            read_unaligned = self.new_out_file(read, self.folders.bowtie2_processed_unaligned)
            log_file = self.folders.bowtie2_log / f'{read.stem}.log'
            command.append('-U')
            command.append(read)
            command.append('--un')
            command.append(read_unaligned)
            log_files.append(log_file)
            list_of_commands.append(command)
            list_of_out.append(read_aligned)
        return zip(list_of_commands, log_files, list_of_out)
            
    def run_command(self):
        for command, logfile, reads_out in self.create_command():
            with open(logfile, 'w') as file:
                #convert to bam file
                bowtie2 = subprocess.Popen(command, stderr=file, stdout=subprocess.PIPE)
                sam_view = subprocess.Popen(['samtools', 'view', '-b'], stdin=bowtie2.stdout, stdout=subprocess.PIPE)
                subprocess.call(['samtools', 'sort', '-o', reads_out], stdin=sam_view.stdout)
                #index
                subprocess.call(['samtools', 'index', '-@', self.threads, reads_out])
      
    
#this works
class FeatureCountsCommando(ProgramCommando):
    def __init__(self, config_file):
        super().__init__(config_file)
        self.program = 'featureCounts'
        
    def extract_parameters(self):
        command = [self.program, '-T', self.threads]
        for params in self.config[self.program]:
            for key, value in params.items():
                if key == 'input':
                    self.input_files = self.reads_in(value, suffix='.bam')
                elif key in ['O', 'f']:
                    command.append(f'-{key}')
                elif len(key) == 1:
                    command.append(f'-{key}')
                    command.append(str(value)) 
                else:
                    command.append(f'--{key}')
                    command.append(str(value))             
        return command
    
    def create_command(self):
        # to get the attribute self.input_files, self.extract_parameters needs to be run
        self.extract_parameters()
        command = self.extract_parameters()
        self.count_matrix = Path(self.folders.feature_counts, 'count_matrix.txt')
        command.append('-o')
        command.append(self.count_matrix)
        command.extend(self.input_files)
        
        return command
    
    def clean_matrix(self):
        # read in data 
        counts_table = pd.read_csv(self.count_matrix, delimiter='\t', skiprows=1)
        # filter for correct columns
        counts_table = counts_table >> select(contains('Geneid') | contains('bam'))
        # rename the columns
        raw_names = ['gene_id']
        samples = sorted([file.stem for file in self.folders.raw_reads.iterdir() if file.is_file()])
        raw_names.extend(samples) 
        counts_table.columns = raw_names
        # save and remove original file
        out_file = self.folders.feature_counts / 'count_matrix.csv'
        counts_table.to_csv(out_file, index=False)
        os.remove(self.count_matrix)
    
    def run_command(self):
        command = self.create_command()
        subprocess.call(command)
        self.clean_matrix()
 

class FastpCommando(ProgramCommando):
    def __init__(self, config_file):
        super().__init__(config_file)
        self.program = 'fastp'
       
    def extract_parameters(self):
        command = [self.program, '--thread', self.threads]
        for params in self.config[self.program]:
            for key, value in params.items():
                if key == 'input':
                    self.input_files = self.reads_in(value)
                elif len(key) == 1:
                    command.append(f'-{key}')
                    command.append(str(value))
                else:
                    command.append(f'--{key}')
                    command.append(str(value))             
        return command
      
    def create_command(self):
        # to get the attribute self.input_files, self.extract_parameters needs to be run
        self.extract_parameters()
        list_of_commands = []
        for read in self.input_files:
            command = self.extract_parameters()
            read_out = self.new_out_file(read, self.folders.fastp_processed)
            log_file = self.folders.fastp_log / f'{read.stem}.json'
            command.append('-i')
            command.append(read)
            command.append('-o')
            command.append(read_out)
            command.append('--json')
            command.append(log_file)
            list_of_commands.append(command)
        return list_of_commands
            
    def run_command(self):
        commands = self.create_command()
        for command in commands:
            subprocess.call(command)