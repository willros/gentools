from pathlib import Path

def default_config(directory: str):
    '''creates default config file in the path specified'''
    
    path = Path(directory, 'config.yaml').absolute()
    content = '### Setup ###\n\n# Enter the directory containing the raw reads and a directory name that you want to use as a working directory\n\nraw_reads: data_smallseq_10_files\nworking_directory: gentools_smallseq_10_files\nthreads: 6\n\n### Preprocessing ###  \n\n# The default configuration is optimized to work with miRNA data. Especially libraries prepared with the \n# Small-seq protocol. \n\numi_tools_extract:\n- input: raw\n- mode: extract\n- extract-method: regex\n- bc-pattern: (?P<discard_1>.*)(?P<umi_1>[ACT]{8}CA)\n\n# For trimmed only enter Y for yes and N for no if you want cutadapt to only keep reads that have been filtered.\n\ncutadapt:\n- input: umi_tools_extract\n- adapter: TGGAATTCTCGGGTGCCAAGG\n- minimum-length: 18\n- maximum-length: 41\n- error-rate: 0.1\n- overlap: 1\n- trimmed_only?: N\n\n# enter folder and prefix of index name, e.g. index/human_genome/hg38\n\nbowtie2:\n- input: cutadapt\n- k: 100\n- local: very-sensitive-local\n- x: /bowtie2_index/genome/GRCh38\n\n\numi_tools_dedup:\n- input: bowtie2_aligned\n- mode: dedup\n- method: unique\n\n\n# enter the path to the annotation file (a), GTF or GFF\n\nfeatureCounts:\n- input: bowtie2_aligned\n- M:\n- O:\n- a: gencode.v39.annotation.gtf.gz\n\n# Enter the coldata for the two conditions of the dataset. For example, cancer: 1 - 10, normal: 11 - 20\n\nDESeq2:\n- input: featureCounts\n- coldata:\n    - cancer: 1 - 10\n    - normal: 11 - 20\n\n\n# Write the programs in the order that you want to run them in. This order will be run by the run-pipe command.\n\norder:\n- umi_tools_extract\n- cutadapt\n- bowtie2\n- umi_tools_dedup\n- featureCounts\n- DESeq2\n'
    with open(path, 'w+') as file:
        file.write(content)
   
    
    
