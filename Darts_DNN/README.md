# Readme for Darts_DNN
### Author: "Zijun Zhang"
### Date: "3.17.2018"

### Table of Contents
- [Installation](#installation)
- [Using Darts DNN](#using-darts-dnn)
- [Training from scratch](#training-from-scratch)

### Installation
To install `Darts_DNN` python package, navigate to this folder, then type
```
> cd Darts_DNN
> make install
```

The `Makefile` will also connect to the Internet to download the latest exon cis- feature set as well as
the trained parameter data files. Currently, these files will be released upon paper publication or by request.

There are a few Deep-learning packages that `Darts_DNN` requires, including
the popular interfaces [Keras](#), [Theano](#). 

To test whether you have successfully installed `Darts_DNN`, type the following command in your shell:

```
> Darts_DNN -h
usage: Darts_DNN [-h] [--version] {train,predict,build_feature} ...

Darts_DNN -- DARTS - Deep-learning Augmented RNA-seq analysis of Transcript
Splicing

positional arguments:
  {train,predict,build_feature}
    train               Darts_DNN train: train a DNN model using Darts
                        Framework from scratch
    predict             Darts_DNN predict: make predictions on a built feature
                        sets in h5 format
    build_feature       Darts_DNN build_feature: build feature file given
                        required information

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit

For command line options of each sub-command, type: Darts_DNN COMMAND -h
```

`Darts_DNN` was trained using an older version of `Keras`, and we are actively working on
adapting to the new version.


### Using Darts DNN

Assume you have already run `Darts_BHT` and get the Darts-flat inference output file, say `"darts_flat.out.txt"`. 
There are two simple steps to run `Darts_DNN` prediction on it:

#### Darts_DNN build_feature
You will need to build the feature file for your target Darts-flat output. The input is `"darts_flat.out.txt"`, and the
output is a feature set in hdf5 data store. Below is an example:

```
Darts_DNN build_feature -i darts_flat.out.txt \
-c /path/to/ENCODE_sequenceFeature_absmax_normalized.h5 \
-e /path/to/condition1/kallisto/ /path/to/condition2/kallisto/ \
-o data.h5
```

#### Darts_DNN predict
Now run the `predict`. This will also estimate the prediction accuracy by
the significant Darts-flat events, so that the users can decide whether to
proceed running `Darts_BHT` with the informative prior.


```
Darts_DNN predict -i data.h5 -o pred.txt \
-m /path/to/model_param.h5
```


#### Darts_Snakemake_pipeline
`Darts_DNN` is best used in `Darts_Snakemake_pipeline`, given the complicated nature of the
whole procedure, which includes construction of feature sets for any given specific data,
calling deep-learning predictions, evaluation of the prediction accuracy, and subsequently
incorporation of prediction to run `Darts_BHT` with informative prior.


### Training from scratch

You can train your own `Darts_DNN` model from scratch by first building a set a training data files in hdf5 format, then use
the `train` subcommand in `Darts_DNN` package:

```
> Darts_DNN train -i /some/path/to/train/folder -o /some/path/to/output/folder
```

The training data of ENCODE and Roadmap will be released upon paper publication.
