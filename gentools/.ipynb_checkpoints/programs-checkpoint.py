import yaml
from pathlib import Path
import subprocess
from .utils import Configfile


def create_default_config(work_dir):
    default = {'raw_reads': 'test_data',
     'working_directory': 'new_test',
     'bowtie': [{'v': 10}, {'n': '--norc'}],
     'umi_tools': [{'method': 'extract'},
      {'extract-method': 'regex'},
      {'bc-pattern': '(?P<discard_1>.*)(?P<umi_1>[ACT]{8}CA)'}],
     'cutadapt': [{'adapter': 'TGGAATTCTCGGGTGCCAAGG'},
      {'minimum-length': 18},
      {'error-rate': 0.1},
      {'overlap': 1}]}
    
    path = Path(work_dir, 'gentool_config.yaml')
    
    with open(path, 'w') as file:
        yaml.dump(default, file)
        
class CreateFolders:
    '''Class that creates correct folder structure'''
    def __init__(self, config_file):
        self.config = Configfile(config_file)
        self.working_directory = Path('.') if self.config['working_directory'] == 'this_folder' else Path(self.config['working_directory'])
               
    def __path_to_parent_folder(self, program):
        base = self.working_directory / program
        log = base / 'log'
        processed = base / 'processed'
        return log, processed
        
    def __make_directories(self, *directory):
        for i in directory:
            Path(i).mkdir(exist_ok=True, parents=True)
      
    def create_directories(self):
        for program in self.config._yaml:
            if program not in ['raw_reads', 'working_directory']:
                log, processed = self.__path_to_parent_folder(program)
                self.__make_directories(log, processed)
        Path(self.working_directory, 'results').mkdir(exist_ok=True)
        
class ProgramCommando:
    ''' Superclass for different programs that can be run in gentools '''
    def __init__(self, config_file):
        self.config = Configfile(config_file)
        self.working_dir = '.' if self.config['working_directory'] == 'this_folder' else self.config['working_directory']
        self.program = ''
    
    def new_out_file(self, read, new_dir, suffix=None):
        suffix = read.suffix if not suffix else suffix
        return new_dir / (read.stem + '_' + self.program + suffix)
    
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
        self.reads_in = sorted([file for file in Path(self.working_dir, 'umi_tools', 'processed').iterdir() if file.is_file()])
        self.reads_out = Path(self.working_dir, 'cutadapt', 'processed')
        self.reads_log = Path(self.working_dir, 'cutadapt', 'log')

    def extract_parameters(self):
        command = [self.program]
        for params in self.config[self.program]:
            for key, value in params.items():
                if len(key) == 1:
                    command.append(f'-{key}')
                    command.append(str(value))
                else:
                    command.append(f'--{key}')
                    command.append(str(value))             
        return command
        
    def create_command(self):
        list_of_commands = []
        log_files = []
        for read in self.reads_in:
            command = self.extract_parameters()
            read_out = self.new_out_file(read, self.reads_out)
            log_file = self.reads_log / f'{read.stem}.log'
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
        
# WORKS LIKE A CHARM! :D       
class UmitoolsCommando(ProgramCommando):
    def __init__(self, config_file):
        super().__init__(config_file)
        self.program = 'umi_tools'
        self.reads_in = sorted([file for file in Path(self.config['raw_reads']).iterdir() if file.is_file()])
        self.reads_out = Path(self.working_dir, 'umi_tools', 'processed')
        self.reads_log = Path(self.working_dir, 'umi_tools', 'log')
    
    def extract_parameters(self):
        command = [self.program]
        for params in self.config[self.program]:
            for key, value in params.items():
                if len(key) == 1:
                    command.append(f'-{key}')
                    command.append(str(value))
                elif key == 'method':
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
        list_of_commands = []
        for read in self.reads_in:
            command = self.extract_parameters()
            read_out = self.new_out_file(read, self.reads_out)
            log_file = self.reads_log / f'{read.stem}.log'
            command.append('-I')
            command.append(read)
            command.append('-S')
            command.append(read_out)
            command.append('-L')
            command.append(log_file)
            list_of_commands.append(command)
        return list_of_commands
            
    def run_command(self):
        for command in self.create_command():
            subprocess.call(command)
