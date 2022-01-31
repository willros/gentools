from pathlib import Path

def default_config(directory: str):
    '''creates default directory in the path specified'''
    
    path = Path(directory, 'config.yaml').absolute()
    content = '#SETUP \nraw_reads: test_data\nworking_directory: testar_lite\nthreads: 6\n\n# choose input between: raw, cutadapt, umi_tools_extract, umi_tools_dedup,\n# bowtie2_aligned, bowtie2_unaligned, star_aligned, star_unaligned...\numi_tools_extract:\n- input: raw\n- mode: extract\n- extract-method: regex\n- bc-pattern: (?P<discard_1>.*)(?P<umi_1>[ACT]{8}CA)\n\n# Y for yes N for no\ncutadapt:\n- input: raw\n- adapter: TGGAATTCTCGGGTGCCAAGG\n- minimum-length: 18\n- maximum-length: 41\n- error-rate: 0.1\n- overlap: 1\n- trimmed_only?: N\n\nfastp:\n- input: raw\n- adapter_sequence: TGGAATTCTCGGGTGCCAAGG\n\n# enter folder and prefix of index name, e.g. index/human_genome/hg38\nbowtie2:\n- input: cutadapt\n- k: 100\n- local: very-sensitive-local\n- x: indices/bowtie2/GRCh38/GRCh38_noalt_as\n\numi_tools_dedup:\n- input: bowtie2_aligned\n- mode: dedup\n- method: unique\n\n\n# enter the path to the annotation file (a), GTF or GFF\nfeatureCounts:\n- input: umi_tools_dedup\n- t: miRNA\n- g: ID\n- f:\n- O:\n- a: annotation/mirbase_hsa.gff\n\n# create csv with only samples and empty condition column. User should be able to enter wich row number should be what:\n# e.g. 1-10 = cancer. 11-20 = normal, and the program should then insert that information\nDESeq2:\n- input: featureCounts\n- coldata:\n    - sick: 1 - 1\n    - normal: 2 - 3\n\norder:\n- umi_tools_extract\n- cutadapt\n- bowtie2\n- umi_tools_dedup\n- featureCounts\n#- DESeq2\n'
    with open(path, 'w+') as file:
        file.write(content)
   
    
    
