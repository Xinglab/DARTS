# Readme for Darts_DNN
### Author: "Zijun Zhang"
### Date: "2.15.2018"

### Table of Contents
- [Installation](#installation)
- [Training from scratch](#training-from-scratch)
- [Using Darts DNN](#using-darts-dnn)

### Installation
The `Darts_DNN` does need installation; once you download/clone the repo into your
local folder, you can directly call it for training and prediction. 

However, there are a few Deep-learning packages that `Darts_DNN` requires, including
the popular interfaces [Keras](#), [Theano](#). 

`Darts_DNN` was trained using an older version of `Keras`, and we are actively working on
adapting to the new version.

### Training from scratch

To train your own `Darts_DNN` on your local machine, first run the
`configure.py`(TODO) to configure the paths to your local training data, then 
type following command in your shell:

```
python deep_neural_net_gpu.py
```

The training data will be released upon publication.

### Using Darts DNN

`Darts_DNN` is best used in `Darts_Snakemake_pipeline`, given the complicated nature of the
whole procedure, which includes construction of feature sets for any given specific data,
calling deep-learning predictions, evaluation of the prediction accuracy, and subsequently
incorporation of prediction to run `Darts_BHT` with informative prior.

