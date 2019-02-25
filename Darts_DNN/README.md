# README

[![Anaconda-Server Badge](https://anaconda.org/darts-comp-bio/darts_dnn/badges/version.svg)](https://anaconda.org/darts-comp-bio/darts_dnn)
[![Anaconda-Server Badge](https://anaconda.org/darts-comp-bio/darts_dnn/badges/installer/conda.svg)](https://conda.anaconda.org/darts-comp-bio)
[![Anaconda-Server Badge](https://anaconda.org/darts-comp-bio/darts_dnn/badges/latest_release_date.svg)](https://anaconda.org/darts-comp-bio/darts_dnn)
[![Anaconda-Server Badge](https://anaconda.org/darts-comp-bio/darts_dnn/badges/platforms.svg)](https://anaconda.org/darts-comp-bio/darts_dnn)
[![Anaconda-Server Badge](https://anaconda.org/darts-comp-bio/darts_dnn/badges/downloads.svg)](https://anaconda.org/darts-comp-bio/darts_dnn)
[![Anaconda-Server Badge](https://anaconda.org/darts-comp-bio/darts_dnn/badges/license.svg)](https://anaconda.org/darts-comp-bio/darts_dnn)

Version: "v0.1.0"

Author: "Zijun Zhang"

Date: "2.17.2019"

### Table of Contents
- [Installation](#installation)
- [Testing](#testing)
- [Using Darts DNN](#using-darts-dnn)
- [Training from scratch](#training-from-scratch)

### Installation

#### Install via Anaconda
The recommended way to install `Darts_DNN` is through [Anaconda](https://anaconda.org/darts-comp-bio).
You can also create a new environment for Darts, because currently DARTS works in Python 2.7.

```bash
conda create -n darts python=2.7 # optional
source activate darts
conda install -c darts-comp-bio darts_dnn
```

This will allow conda to do all the heavy-lifting and most often the easiest way to get things spinning.


#### From GitHub
Alternatively, to install `Darts_DNN` python package from Github, navigate to this folder, then type
```sh
> cd Darts_DNN
> make install
```

There are a few Deep-learning packages that `Darts_DNN` requires, including
the popular high-level interface [Keras](#). 

To test whether you have successfully installed `Darts_DNN`, type the following command in your shell:

```sh
> Darts_DNN -h
usage: Darts_DNN [-h] [--version] {train,predict,build_feature,get_data} ...

Darts_DNN -- DARTS - Deep-learning Augmented RNA-seq analysis of Transcript
Splicing

positional arguments:
  {train,predict,build_feature,get_data}
    train               Darts_DNN train: train a DNN model using Darts
                        Framework from scratch
    predict             Darts_DNN predict: make predictions on a built feature
                        sets in h5 format
    build_feature       Darts_DNN build_feature: build feature file given
                        required information
    get_data            Darts_DNN get_data: connects online to get Darts_DNN
                        data for the current version.

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit

For command line options of each sub-command, type: Darts_DNN COMMAND -h
```


### Using Darts DNN

Assume you have already run `Darts_BHT` and get the Darts-flat inference output file, say `"darts_bht.flat.txt"`. 
There are two simple steps to run `Darts_DNN` prediction on it:

#### Darts_DNN build_feature
You will need to build the feature file for your target Darts-flat output. The input is `"darts_flat.out.txt"`, and the
output is a feature set in hdf5 data store. Below is an example:

```sh
Darts_DNN build_feature -i darts_flat.out.txt \
-c /path/to/ENCODE_sequenceFeature_absmax_normalized.h5 \
-e /path/to/condition1/kallisto/ /path/to/condition2/kallisto/ \
-o data.h5
```

#### Darts_DNN predict
Now run the `predict`. This will also estimate the prediction accuracy by
the significant Darts-flat events, so that the users can decide whether to
proceed running `Darts_BHT` with the informative prior.


```sh
Darts_DNN predict -i data.h5 -o pred.txt \
-m /path/to/model_param.h5
```


### Training from scratch

You can train your own `Darts_DNN` model from scratch by first building a set a training data files in hdf5 format, then use
the `train` subcommand in `Darts_DNN` package:

```
> Darts_DNN train -i /some/path/to/trainListFile -o /some/path/to/output/folder
```
