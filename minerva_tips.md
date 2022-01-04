# The purpose of this to give a basic overview / introduction on the various ways to use minerva:

## R Studio:
- R-Studio can be started with the following: `sh minerva-rstudio-web-r4.sh `
- R-Studio CANNOT download data because it is behind the singularity firewall. If you want to download data, such as using a BioLinks package, then please run your script via the command line and not piecemeal in R-Studio
- If you want to download R packages that are NOT included in the 
## Jupyter Notebooks:
- Jupyter notebooks can be run with the following command
  - minerva-jupyter-web.sh -h 
  - the bash script accepts a --image flag that allows you to specify a docker image. that way you can run docker pull and specify a custom image
- Jupyter uses singularity, which is a more secure version of docker, but you can still access images via dockerhub. You need to pull the image that you wouild like to use. 
  - pulling base python image 3.7: `singularity pull --name "jupyter-base-376.img" docker://jupyter/base-notebook:python-3.7.6 # referencing image here`
  - pulling jupyter tensorflow if you want to access GPUs: `singularity pull --name "jupyter-tf.img" docker://jupyter/tensorflow-notebook`
- You can then use these singularity images when starting the jupyter notebook. 
  - minerva-jupyter-web.sh --image jupyter-base-376.img
  - minerva-jupyter-web.sh --image jupyter-tf.img
## Scripts
- Scripts must be submitted to queues via `bsub`
  - `bjobs`: status of jobs
  - `bjobs -u rantid01`: status of MY jobs
  - `bqueues`: queues    
  - `bpeek`: See the output of a job submission so far
  - `bkill`: kill jobs in the queue
- Submitting to a queue requires the following:
  - `bsub -P acc_mnibc_bcg -n 5 -W 05:00 -R rusage[mem=50000] "python spatial_script.py" -o spatial_010421.txt -e spatial_errors_010421.txt`
  - `-P`: the project to submit to. The spelling of our project is messed up, so it is `acc_mnibc_bcg`
  - `-W`: wall time. Each queue has a max that you can find online
  - `-n`: number of CPUs to use
  - `-R rusage[mem=50000]`: RAM usage in megabytes
  - `-o`: output file
  - `-e`: error file

## Random notes on issues I have encountered
### For the annoy installation
- You cannot install annoy via pip. You must create a new conda environment and install it via conda-forge. 
  - Python just released Python 3.10. It doesnt work with most things. Stick with 3.7 - 3.9 until compatibility is fixed. 
    - `conda activate uroconda381` (this is 3.8.1 environment)
    - `module load gcc `
    - `conda install -c conda-forge python-annoy`
    - `pip install -r requirements.txt` for the rest of the requirements
    - `bsub -P acc_mnibc_bcg -n 5 -W 05:00 -R rusage[mem=50000] "python spatial_script.py" -o spatial_010421.txt -e - spatial_errors_010421.txt`
- THE PROJECT NAME IS MISSPELLED!! ITS acc_mnibc_bcg --> NOT NMIBC
