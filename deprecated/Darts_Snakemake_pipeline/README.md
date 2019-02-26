# Readme for Darts Snakemake pipeline
### Author: "Zijun Zhang"
### Date: "2.15.2018"

### Table of Contents
- [Installation](#installation)
- [Testing](#testing)
- [Using Darts_Snakemake_pipeline](#usage)

### Installation
The Darts Snakemake pipeline runs in python3 environment and depends on the python 
package [snakemake](#).

You will also need to put the following software in your environment variables. Click the link
will navigate you to the software page if you don't have it yet:
  - [STAR](#)
  - [kallisto](#)
  - [rMATS](#)

### Testing
The `Darts_Snakemake_pipeline` contains a dummy test folder called `test_snakemake`.
If you have successfully installed snakemake, type the following
```
PROJECT='test_snakemake' snakemake -n
```
.. and you should see a list of jobs to run. Since this is a dry-run, this indicates
the pipeline itself is runnable.

### Usage
The input for the pipeline is a configuration file and a set of fastq files (either gzipped or plain text). And the pipeline
will generate the splicing analysis results by `Darts_BHT` and `Darts_DNN`, which utilized our trained DNN model (as of Feburary 2018).

To set up a new project, say named 'new_project', first make the file directories:
```
mkdir new_project
mkdir new_project/config
mkdir new_project/fq
```
Then make a configuration file, similar to the test file in 'test_snakemake/config/config.yaml'.
Basically, you will need to specify the sample names, and all the comparisons.

Finally, to run the whole pipeline,
```
PROJECT='new_project' snakemake
```
We also provide the option for submitting all jobs to a SGE cluster
```
qsub -N new ./submit.sh new_project
```
