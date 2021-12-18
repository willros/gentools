from .utils import Configfile, CreateFolders
from .programs import Bowtie2Commando, CutadaptCommando, UmitoolsDedupCommando, UmitoolsExtractCommando, FeatureCountsCommando
from .preproccessing import Bowtie2Graphs, CutadaptGraphs, UmitoolsDedupGraphs, UmitoolsExtractGraphs, FeatureCountsGraphs


class GenTools:
    def __init__(self, config_file):
        self.config_file = config_file
        self.config = Configfile(config_file)
        self.folders = CreateFolders(config_file)
        self.order = self.config['order']
        
    def _run_program(self, program):
        program.run_command()
        
    def _run_analysis(self, program):
        program.make_pie_graph()
    
    # add featureCounts
    def _choose_program(self, program):
        program_dict = {'cutadapt': CutadaptCommando(self.config_file), 
                        'bowtie2': Bowtie2Commando(self.config_file),
                       'umi_tools_dedup': UmitoolsDedupCommando(self.config_file), 
                        'umi_tools_extract': UmitoolsExtractCommando(self.config_file),
                       'featureCounts': FeatureCountsCommando(self.config_file)}   
        return program_dict[program]
    
    def _choose_analysis(self, program):
        analysis_dict = {'cutadapt': CutadaptGraphs(self.config_file), 
                        'bowtie2': Bowtie2Graphs(self.config_file),
                       'umi_tools_dedup': UmitoolsDedupGraphs(self.config_file), 
                        'umi_tools_extract': UmitoolsExtractGraphs(self.config_file),
                        'featureCounts': FeatureCountsGraphs(self.config_file)}   
        return analysis_dict[program]

    def run_pipeline(self):
        #create folders first
        self.folders.make_directories()
        
        # create object only when needed and run it
        for program in self.order:
            print(f'Running {program}')
            tool = self._choose_program(program)
            self._run_program(tool)
            print(f'{program} done!')
            
    def run_analysis(self):
        for program in self.order:
            tool = self._choose_analysis(program)
            self._run_analysis(tool)
            print(f'Look in {self.folders.results} for analysis of {program}')


        
