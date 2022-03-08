# gentools - a CLI Python wrapper for analyzis of miRNA sequencing data. 

gentools is built with Python and can be runned entirely from the command line.
geentools uses well known bioinformatic programs under the hood. These programs needs to be installed separately. The programs are: 

* bowtie2
* cutadapt
* umi_tools
* featureCounts

Recommended way to install with conda. See the `conda_env.yaml` file for dependencies and for an easy installation through conda:
```conda create -n <name_of_new_environment> -f conda_env.yaml```


### Installation

Clone this repository:
`git clone https://github.com/willros/gentools.git`

cd into the folder and install via pip:
```pip install .```

### Download of indicies and GFF files

For gentools to work it needs to have an index and GFF file available. The recommended files are listed below, but you can use the index file most suitable for your analysis. 

* Index for bowtie2
Download the following prebuilt index from bowtie2:
https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip

* GFF file for featureCounts
Download the following GFF file from gencode: 
https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.annotation.gff3.gz

### Usage:
gentools relies on a config file. A default config file can be generated through the following command:
```gentools config <path_of_folder_to_put_file>```

1. This results in a config file which has optimzed parameters for miRNA analysis. 

2. Edit the file with approppriate path to indexes and number of threads. Specify input files to each program.

3. Specify the order of execution. 

4. Run gentools through the command line, with the command: gentools run-pipe <config.yaml>

4. To obtains a summary report from the preprocessing, run: gentools run-analysis <config.yaml>

5. To further analyze the data, go to https://share.streamlit.io/willros/gentools_streamlit/main/app.py. Take the files generated from gentools (coldata, de_genes and countmatrix) and drag them into the correct place. This will generate plots and tables about the sequencing data. 



### Python dependencies

* numpy
* pandas
* plotnine
* click

```All of the above is present in the requirements.txt and is automatically installed when pip installed```



### Author
William Rosenbaum 


