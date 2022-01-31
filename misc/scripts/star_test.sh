
params=' --runThreadN 16
--sjdbGTFfile genomes/gencode_miRNA_subset.gtf
--alignEndsType EndToEnd
--outFilterMismatchNmax 1
--outFilterMultimapScoreRange 0
--quantMode TranscriptomeSAM GeneCounts
--outReadsUnmapped Fastx
--outSAMtype BAM SortedByCoordinate
--outFilterMultimapNmax 10
--outSAMunmapped Within
--outFilterScoreMinOverLread 0
--outFilterMatchNminOverLread 0
--outFilterMatchNmin 16
--alignSJDBoverhangMin 1000
--alignIntronMax 1
--outWigType wiggle
--outWigStrand Stranded
--outWigNorm RPM
'

STAR --genomeDir indexes/star_index_gencode/ --readFilesIn snakemake_cutadapt/SRR3495727_umi_cutadapt.fastq --outFileNamePrefix STAR_TEST/XXX_ $params
