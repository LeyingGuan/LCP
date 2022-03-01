# LCP

## Installation

install_github("LeyingGuan/LCP")

Make sure that you installed the pre-requisites. The neuralnet work related packages are only used to perform empirical comparisons as described in the paper.

## Performing LCP construction
You can find two examples using LCP under man/examples. Example 1 performs LCP with user-specified h and example 2 performs LCP with auto-tuned LCP. If you want to rerun all numerical examples [1], please also install LCPexperiment, and in which case, please make sure that you have properly installed  torch, torchvision, and XRPython.

1. example_intro.r: illustration of LCP with simple example (Fig 1 in [1])



[1] Localized Conformal Prediction: A Generalized Inference Framework for Conformal Prediction https://arxiv.org/abs/2106.08460.


