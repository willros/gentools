# gentools - a package to help you with every day sequencing analysis processes!

gentools is built with python and can be runned entirely from the command line.
geentools relies on many different and highly regarded programs, such as bowtie2, cutadapt etc..

* bowtie2
* cutadapt
* umi_tools
* featureCounts
* STAR

These programs need to be installed separatly. Recommended way to install with conda. See conda_config.yaml file for easily installing everything through conda:
```conda create -n <name_of_new_environment> -f conda_config.yaml```

### Default settings

The default settings is set to be optimal for miRNA analysis.

### Installation

Install gentools via pip:
```pip install gentools```

### Usage:
gentools relies on a config file. A default config file can be generated through:
```gentools config <path_of_folder_to_put_file>```

Edit the file with approppriate path to indexes and so on. Specify input files to each program. 

### Python dependencies

* numpy
* pandas
* datar
* click


```All of the above is present in the requirements.txt and is automatically installed when pip installed```


## TODO:
should I do so that every program creates its own folders instead of the class folders doing that? 

### Author
William Rosenbaum 


