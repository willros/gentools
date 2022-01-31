#### information and script to build STAR index for smallseq analysis. 
### will use --sjdbOverhang to 49, because apparently it is supposed to be len(read) - 1
#### will use parameters from this site  https://www.encodeproject.org/documents/b4ec4567-ac4e-4812-b2bd-e1d2df746966/@@download/attachment/ENCODE_miRNA-seq_STAR_parameters_v2.pdf
# beacuse the parameters used in smallseq article does not align at all. 


# STAR index command:
STAR --runThreadN 30 --runMode genomeGenerate --genomeDir indexes/star_index_ensembl --sjdbGTFfile genomes/Homo_sapiens.GRCh38.104_ensembl.gtf --sjdbOverhang 49 --genomeFastaFiles genomes/Homo_sapiens.GRCh38.dna_sm.primary_assembly_ensembl.fa


