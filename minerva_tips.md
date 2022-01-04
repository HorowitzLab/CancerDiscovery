# The purpose of this to give a basic overview / introduction on the various ways to use minerva:

## R Studio:
- R-Studio can be started with the following: `sh minerva-rstudio-web-r4.sh `
- In order to run packages that utilize data downloads (such as BioLinks) you need to set up a proxy to get around the singularity firewall. Run the following commands before launching R-Studio
  - `export http_proxy=http://172.28.7.1:3128`
  - `export https_proxy=http://172.28.7.1:3128`
  - `export all_proxy=http://172.28.7.1:3128`
  - `export no_proxy=localhost,*.hpc.mssm.edu,*.chimera.hpc.mssm.edu,172.28.0.0/16`
- If you would like to install packages that are not default, you need to do the following
  - Open R-Studio via the link that the bash command outputs 
  - Select the shell terminal console WITHIN R-STUDIO
  - Do the following: 
    `[INFO] $ export http_proxy=http://172.28.7.1:3128`
    `[INFO] $ export https_proxy=http://172.28.7.1:3128`
    `[INFO] $ export all_proxy=http://172.28.7.1:3128`
    `[INFO] $ export no_proxy=localhost,*.hpc.mssm.edu,*.chimera.hpc.mssm.edu,172.28.0.0/16`
    `[INFO] $ R`
    `[INFO] >>> install.packages(name_of_package)`
- Troubleshooting: 
  - If you try to run the job, it successfully submits to a queue, but it then dies immediately, this is likely due to your home directory not having enough space. Home directories are limited to 20 GB. 
  - Delete any large files in your home directory ***not the project directory***
  - Run `du -hs .[^.]*` to see the hidden files and their size if any of the non-hidden files are not big. 

## Jupyter Notebooks:
- Jupyter notebooks can be run with the following command
  - minerva-jupyter-web.sh -h 
  - the bash script accepts a --image flag that allows you to specify a docker image. that way you can run docker pull and specify a custom image
- Jupyter uses singularity, which is a more secure version of docker, but you can still access images via dockerhub. You need to pull the image that you wouild like to use. 
  - pulling base python image 3.7: `singularity pull --name "jupyter-base-376.img" docker://jupyter/base-notebook:python-3.7.6 # referencing image here`
  - pulling jupyter tensorflow if you want to access GPUs: `singularity pull --name "jupyter-tf.img" docker://jupyter/tensorflow-notebook`
- You can then use these singularity images when starting the jupyter notebook. 
  - `minerva-jupyter-web.sh --image jupyter-base-376.img`
  - `minerva-jupyter-web.sh --image jupyter-tf.img`

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
- In order to submit a script to a queue successfully, your username must be successfully associated with a project
  - To check which projects you're associated with run `mybalance`. Your projects will be in the second column in the form `acc_xxx`

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
