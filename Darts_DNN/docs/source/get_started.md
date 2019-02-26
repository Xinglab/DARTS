# Getting Started

## 1. Installation

Installation of ``Darts_DNN`` is made easy through [Anaconda](https://anaconda.org/Darts-comp-bio/darts_dnn).
It's recommended to start by creating a new environment:

```bash
conda create -n darts python=2.7 # optional
source activate darts
conda install -c darts-comp-bio darts_dnn
```

Upon finish, type in the following command in shell to show the help messages:
```bash
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

## 2. Using Predict
``Darts_DNN predict`` is probably the most used utility. Please note, **using prediction DOES NOT require a GPU machine**, and you can totally do it on your laptop!

In the simplest case, the ``predict`` function can be invoked by providing a labelled input file (generated from ``Darts_BHT bayes_infer``) and a trans gene expression file. 

If you have not installed the Darts_DNN previously, you will need to download the cis-Features and trained 
model parameters, etc. through ``Darts_DNN get_data``. get_data function will automatically resume previous run and check md5sum - so don't worry about doubled storage space.

For the purpose of this walk-through tutorial, since our test data is A5SS, we only need to download the files for A5SS splicing events.

```bash
Darts_DNN get_data -d transFeature cisFeature trainedParam -t A5SS
```

Next as an example, download the [test_data](#) from GitHub then run:

```bash
wget https://github.com/zj-zhang/DARTS-BleedingEdge/raw/master/Darts_DNN/test_data/A5SS.thymus_adipose.tgz
tar -xvzf A5SS.thymus_adipose.tgz
Darts_DNN predict -i darts_bht.flat.txt -e RBP_tpm.txt -o pred.txt -t A5SS
```

In the screen log output, you should see something like:
```
2019-02-25 15:02:32,659 - Darts_DNN.predict - INFO -
 AUROC=0.8686118716025868
2019-02-25 15:02:32,659 - Darts_DNN.predict - INFO -
 AUPR=0.5410178835754661
```

The output of the predictions is in the user-specified filename, in this case "pred.txt". The output file is a three-column text file, with

ID   Y_true   Y_pred

The `ID` is a unique identifier for an alternative splicing event. The `Y_true` is the observed posterior probability for differential splicing, and the `Y_pred` is the predicted probability differential splicing. In computing the AUROC and AUPR, only the high-confidence events (i.e. Y_true>0.9 as positive, Y_true<0.1 as negative) are used.

This prediction output file `pred.txt` can be further utilized to perform deep-learning augmented analysis in ``Darts_BHT``. See the user guide page for ``Darts_BHT`` [here](#).


## 3. Using Train to train a model from scratch

If you want to train a new model from scratch, you can download our pre-processed training data by using 
``Darts_DNN get_data`` utilies like below; in this case, we download the training and held-out data for A5SS:
```bash
mkdir A5SS_train
Darts_DNN get_data -d trainingDataSet -t A5SS -o A5SS_train/
cd A5SS_train/
tar -xvzf A5SS.trainSet.tgz
```

The ``Darts_DNN train`` function takes in a training summary file that lists all associated files, like the example below:

```bash
Darts_DNN train -i Darts_DNN-train_data/trainSet/A5SS/A5SS_Roadmap_trainList.txt Darts_DNN-train_data/trainSet/A5SS/A5SS_ENCODE_trainList.txt -o ./ -t A5SS
```

Training usually takes a few hours for A5SS/A3SS/RI, and a few days for SE, depending on the amount of training data and the processing speed for your local machine. For training it is recommended to get a GPU server, especially for training SE events.

## 4. Working with other human genome assemblies

Currently the Darts DNN cis-features are compiled on hg19 genome assembly, so the exon coordinates are based on hg19. If you are using other human genome assembly, say hg38, please use the UCSC [liftOver](https://genome.ucsc.edu/cgi-bin/hgLiftOver) to convert the exon coordinates so that DNN can work correctly.

## 5. FAQ

Wait to be updated.