import click
from .default_config import default_config
from .config import CreateFolders
from .preprocessing import UmitoolsExtractGraphs, CutadaptGraphs, UmitoolsDedupGraphs, Bowtie2Graphs, FeatureCountsGraphs
from .programs import UmitoolsExtractCommando, UmitoolsDedupCommando, CutadaptCommando, Bowtie2Commando, FeatureCountsCommando, FastpCommando
from .main import GenTools
from .differential_expression import DifferentialExpression

@click.group()
def cli():
    '''WELCOME TO GENTOOLS'''
    
    pass


### SET UP COMMAND TO CREATE CONFIG FILE AND FOLDER STRUCTURE ###

#WORKS!               
@cli.command()
@click.argument('directory')
def config(directory):
    '''Creates default config file in the directory specified in the directory argument'''
    
    default_config(directory)
    click.echo(f'default config file created in {directory}')

    
#WORKS!    
@cli.command()
@click.argument('config_file')
def folders(config_file):
    ''' Creates folder structure for the pipeline. Enter the work directory for the folder tree in the config file, 
    and pass the config file as an argument to this function.'''
    
    folders = CreateFolders(config_file)
    folders.make_directories()
    click.echo(f'Folder structure created in {folders.working_directory}')
    

### RUN WHOLE PIPELINE AND ANALYSIS ###


#WORKS
@cli.command()
@click.argument('config_file')
def run_pipe(config_file):
    '''Runs the pipeline in the order specified in the config file.'''
    pipe = GenTools(config_file)
    click.echo(f'Starting to run the pipeline')
    pipe.run_pipeline()
    click.echo(f'The pipeline is done!')
    
#WORKS
@cli.command()
@click.argument('config_file')
def run_analysis(config_file):
    '''Runs analysis of the step runned in the pipeline and saves images and csv files in
    in the results folder'''
    
    pipe = GenTools(config_file)
    pipe.run_analysis()
    click.echo(f'The analysis is done!')

    
### RUNNING INVIDIVUAL PROGRAMS GRAPHS###

#WORKS
@cli.command()
@click.argument('config_file')
def umi_extract_info(config_file):
    '''Collects information about the umi_tools extraction and creates. 
    Generates a graph and a csv file in the results folder.'''
    
    umi_graph = UmitoolsExtractGraphs(config_file)
    umi_graph.make_pie_graph()
    click.echo('Look for result files in the results folder!')


#WORKS
@cli.command()
@click.argument('config_file')
def umi_dedup_info(config_file):
    '''Collects information about the umi_tools deduplication and creates. 
    Generates a graph and a csv file in the results folder.'''
    
    umi_graph = UmitoolsDedupGraphs(config_file)
    umi_graph.make_pie_graph()
    click.echo('Look for result files in the results folder!')


#WORKS
@cli.command()
@click.argument('config_file')
def cutadapt_info(config_file):
    '''Collects information about the cutadapt and creates. 
    Generates a graph and a csv file in the results folder.'''
    
    cutadapt_graph = CutadaptGraphs(config_file)
    cutadapt_graph.make_pie_graph()
    click.echo('Look for result files in the results folder!')
    
#WORKS  
@cli.command()
@click.argument('config_file')
def bowtie2_info(config_file):
    '''Collects information about the bowtie2 and creates. 
    Generates a graph and a csv file in the results folder.'''
    
    bowtie2_graph = Bowtie2Graphs(config_file)
    bowtie2_graph.make_pie_graph()
    click.echo('Look for result files in the results folder!')

#WORKS      
@cli.command()
@click.argument('config_file')
def featurecounts_info(config_file):
    '''Collects information about the featureCounts and creates. 
    Generates a graph and a csv file in the results folder.'''
    
    feat_graph = FeatureCountsGraphs(config_file)
    feat_graph.make_pie_graph()
    click.echo('Look for result files in the results folder!')
    
    
### RUN INVIDIDUAL PROGRAMS ###

#WORKS
@cli.command()
@click.argument('config_file')
def umi_extract_run(config_file):
    '''Runs umi_tools extract only with the input files specified in the config file.'''
    
    tool = UmitoolsExtractCommando(config_file)
    
    click.echo(f'Starting to run {tool.program}')
    tool.run_command()
    click.echo(f'{tool.program} done!')

#WORKS
@cli.command()
@click.argument('config_file')
def umi_dedup_run(config_file):
    '''Runs umi_tools dedup only with the input files specified in the config file.'''
    
    tool = UmitoolsDedupCommando(config_file)
    
    click.echo(f'Starting to run {tool.program}')
    tool.run_command()
    click.echo(f'{tool.program} done!')
    
#WORKS   
@cli.command()
@click.argument('config_file')
def cutadapt_run(config_file):
    '''Runs cutadapt only with the input files specified in the config file.'''
    
    tool = CutadaptCommando(config_file)
    
    click.echo(f'Starting to run {tool.program}')
    tool.run_command()
    click.echo(f'{tool.program} done!')

#WORKS
@cli.command()
@click.argument('config_file')
def bowtie2_run(config_file):
    '''Runs bowtie2 only with the input files specified in the config file.'''
    
    tool = Bowtie2Commando(config_file)
    
    click.echo(f'Starting to run {tool.program}')
    tool.run_command()
    click.echo(f'{tool.program} done!')
    

#WORKS
@cli.command()
@click.argument('config_file')
def featurecounts_run(config_file):
    '''Runs featureCounts only with the input files specified in the config file.'''
    
    tool = FeatureCountsCommando(config_file)
    
    click.echo(f'Starting to run {tool.program}')
    tool.run_command()
    click.echo(f'{tool.program} done!')
    
    
### fastp
@cli.command()
@click.argument('config_file')
def cutadapt_run(config_file):
    '''Runs fastp only with the input files specified in the config file.'''
    
    tool = FastpCommando(config_file)
    
    click.echo(f'Starting to run {tool.program}')
    tool.run_command()
    click.echo(f'{tool.program} done!')

### DESeq2

@cli.command()
@click.argument('config_file')
def deseq2_run(config_file):
    '''Runs DESeq2 only with the input files specified in the config file.'''
    
    tool = DifferentialExpression(config_file)
    
    click.echo(f'Starting to run {tool.program}')
    tool.run_command()
    click.echo(f'{tool.program} done!')
    
               
if __name__ == '__main__':
    cli()
