# Data analysis related to the Ultra-Deep Sequencing

## Usage

### Step 1:

Setup Bioconda as shown [here](https://bioconda.github.io).

### Step 2:

Install Snakemake into an isolated environment

    conda create -n snakemake snakemake

### Step 3:

Checkout this git repository

    git clone https://github.com/bioconda/bioconda-paper.git
    cd bioconda-paper

### Step 4:

Execute the analysis workflow with Snakemake

    snakemake --use-conda -j

This tells Snakemake to run at most as many jobs in parallel as you have CPU cores.

Results can be found in the folder `plots/`.
