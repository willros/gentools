from pathlib import Path

def default_config(directory: str):
    '''creates default directory in the path specified'''
    
    path = Path(directory, 'config.yaml').absolute()
    content = 'raw_reads: test_data\nworking_directory: testar_lite\n#genomes: ../scripts/genomes\n\n## choose input between: raw, cutadapt, umi_tools_extract, umi_tools_dedup bowtie2_aligned, bowtie2_unaligned, star_aligned, star_unaligned\numi_tools_extract:\n- input: raw\n- mode: extract\n- extract-method: regex\n- bc-pattern: (?P<discard_1>.*)(?P<umi_1>[ACT]{8}CA)\n\ncutadapt:\n- input: umi_tools_extract\n- adapter: TGGAATTCTCGGGTGCCAAGG\n- minimum-length: 18\n- maximum-length: 41\n- error-rate: 0.1\n- overlap: 1\n\n# enter folder and prefix of index name, e.g. index/human_genome/hg38\nbowtie2:\n- input: cutadapt\n- k: 100\n- end-to-end: very-sensitive\n- x: indices/bowtie2/GRCh38/GRCh38_noalt_as\n\n#remove mode:dedup?? un necessary? \numi_tools_dedup:\n- input: bowtie2_aligned\n- mode: dedup\n- method: unique\n\n\n# enter the path to the annotation file (a) (GTF or GFF)\nfeatureCounts:\n- input: umi_tools_dedup\n- t: miRNA\n- g: ID\n- f:\n- O:\n- a: annotation/mirbase_hsa.gff\n\n\norder:\n- umi_tools_extract\n- cutadapt\n- bowtie2\n- umi_tools_dedup\n- featureCounts\n'
    with open(path, 'w+') as file:
        file.write(content)
   
    
    
