################################################################################
# Script Name:        fit_hmsc_hpc_models.R
# Description:        This script fits hmsc models. The models are defined in the script build_hmsc_hpc_models.R
#                     The fitting is done using the HPC implementation. See github: https://github.com/hmsc-r/hmsc-hpc
#                     and paper https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011914
#
# Author:             Jonathan Jupke
# Date Created:       2025-09-26
# Last Modified:      2025-09-26
#
# R Version:          R 4.3.2
# Required Packages: 
#
# Notes:              Build to be called from within a singularity container. 
################################################################################

#define python 
python = "C:/HYAPP/Anaconda3-2023.07-2/python.exe"
system2(python, "--version")
Sys.setenv(TF_CPP_MIN_LOG_LEVEL=3)  # reduce debug output from tensorflow
system2(python, "-c \"import tensorflow as tf; print(tf.constant(1))\"")
system2(python, "-c \"import hmsc\"")

post_file_path = file.path(getwd(), "post_file.rds")
python_cmd_args = paste("-m hmsc.run_gibbs_sampler",
                        "--input", shQuote(init_file_path),
                        "--output", shQuote(post_file_path),
                        "--samples", nSamples,
                        "--transient", transient,
                        "--thin", thin,
                        "--verbose", verbose)
cat(paste(shQuote(python), python_cmd_args), "\n")
system2(python, python_cmd_args)
importFromHPC = from_json(readRDS(file = post_file_path)[[1]])
postList = importFromHPC[1:nChains]
fitTF = importPosteriorFromHPC(m, postList, nSamples, thin, transient)
cat(sprintf("fitting time %.1f sec\n", importFromHPC[[nChains+1]]))
savename <- paste0()
saveRDS(fitTF)