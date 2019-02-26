# Getting Started

## 1. Installation

Installation of ``Darts_BHT`` is made easy through [Anaconda](https://anaconda.org/Darts-comp-bio/darts_bht).
It's recommended to start by creating a new environment:

```bash
conda create -n darts python=2.7 # optional
source activate darts
conda install -c darts-comp-bio darts_bht
```

Upon finish, type in the following command in shell to show the help messages:
```bash
> Darts_BHT -h
usage: Darts_BHT [-h] [--version] {rmats_count,bayes_infer} ...

Darts_BHT -- DARTS - Deep-learning Augmented RNA-seq analysis of Transcript
Splicing

positional arguments:
  {rmats_count,bayes_infer}
    rmats_count         Darts_BHT rmats_count: run rMATS-turbo to count
                        junction reads in BAM files
    bayes_infer         Dart_BHT bayes_infer: perform Bayesian hypothesis
                        testing inference

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit

For command line options of each sub-command, type: Darts_BHT COMMAND -h
```

## 2. Using Darts BHT on Empirical RNA-seq data

There are two steps in running ``Darts_BHT`` on empirical RNA-seq data, i.e. 1) Counting, and 2) Statistical inference.

The counting part of the ``Darts_BHT`` is currently powered by the fast [rMATS-turbo](http://rnaseq-mats.sourceforge.net/) counting pipeline. As such, you can run the counting
step of Darts_BHT the same way as in rMATS-turbo. In the simplest case, use the following command
on two bam files as input:

```bash
Darts_BHT rmats_count --b1 /path/to/b1.bam --b2 /path/to/b2.bam \
--gtf /path/to/gtfFile.gtf --od outDir -t paired \
--nthread 8 --readLength 150
```
**Note**: it's important to have the correct `--readLength` specified for your own data; for details, see rMATS [user guide](http://rnaseq-mats.sourceforge.net/user_guide.htm). In this example run, the data has 150bp paired-end reads.

This will generate two types of outputs in the folder {outDir}:
 - 5 types annotation files named "fromGTF.{eventType}.txt", where eventType={SE,A5SS,A3SS,RI,MXE} ;
 - the counted RNA-seq reads named "JC.raw.input.{eventType}.txt".

Next to run the ``Darts_BHT bayes_infer`` on the generated rMATS count files:
```bash
Darts_BHT bayes_infer --rmats-count JC.raw.input.A5SS.txt --od ./ --annot fromGTF.A5SS.txt -t A5SS --cutoff 0.05 --nthread 8
```

This command will prepare the rMATS count in Darts format in "A5SS.input.txt", and then generate the output inference reults in "A5SS.darts_bht.flat.txt". If you have ``openpyxl``, another spreadsheet will be appended to the Excel file "Darts_BHT.results.xlsx". You can install it by `pip install openpyxl`.

The output in the Excel "Darts_BHT.results.xlsx" should be fairly human readable. The primary quantity of interest is the column "Posterior", which depicts the how confident ``Darts_BHT`` is that the difference in the observed data exceeds the user-specified cutoff.


## 3. Using Darts BHT with Darts DNN

The above procedure uses `Darts_BHT` as a conventional computational tool for splicing analysis based on RNA-seq. Given the interface of `bayes_infer`, it is fairly easy to incorporate any informative prior into the differential splicing detection, e.g., use our ``Darts_DNN`` deep learning model to predict the probability of differential splicing based on exon-specific cis-features and sample-specific trans features. You can follow the documentations for Darts_DNN [here](#) to generate an informative prior file named `pred.txt`.

Suppose we have already made the pred.txt file, simple run the same inference command again, but with the `--prior` option specified the predicted prior file:

```bash
Darts_BHT bayes_infer --rmats-count JC.raw.input.A5SS.txt --od ./ --annot fromGTF.A5SS.txt --prior pred.txt -t A5SS --cutoff 0.05 --nthread 8
```

When the prior file is specified, the output file name will become "A5SS.darts_bht.info.txt" instead of "A5SS.darts_bht.flat.txt". Also a new spreadsheet named "A5SS-info" will be appended to "Darts_BHT.results.xlsx" in the output directory, in order to facilitate downstream analyses.


## 4. A Note on Perform Deep-learning Augmented analysis

While it's generally useful to incorporate the informative prior to perform the augmented analysis, there are cases where the prediction is not very accurate and could in turn harm the inference if incorporated as an informative prior (e.g. AQR shRNA knockdown in ENCODE). Hence users should keep an eye on the evaluation results (i.e. AUROC and AUPR) when generating the pred.txt file, and choose accordingly between using 

 1. Augmented results by combining empirical RNA-seq junction counts with prediction (**Informative Posterior**)
 2. Ordinal results using just empirical RNA-seq junction counts (**Flat Posterior**)
 3. Darts_DNN predicted scores (**Informative Prior**)

In general, an AUROC>0.8 indicates good performance when incorporating the informative prior.

## 5. FAQ

Waiting to be updated.