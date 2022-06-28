# Prerequisites

## Java
To check if you have Java available:
```
which java
```

If is not available, you can either:
1. [install](https://java.com/en/download/) it
2. load with:
    ```
    module load java
    ```

## Nextflow

[Install](https://www.nextflow.io/docs/latest/getstarted.html) nextflow:
1. Navigate to where you would download nextflow (in this example, `/home/users/kaitlin3`)
2. Download the executable with either:
    1. `wget -qO- https://get.nextflow.io | bash`
    2. `curl -s https://get.nextflow.io | bash`
3. Make the binary executable
    1. `chmod +x nextflow`
4. Optional: make the nextflow file accessible, by adding either of these to your `~/.bashrc` or `~/.bash_profile`, depending on setup. Be sure to source the file so that changes take effect, with `source ~/.bashrc`.
    1. Add file to $PATH
        1. `export PATH="$PATH:/home/users/kaitlin3/nextflow"`
    2. Set up alias
        1. `alias nextflow='~/nextflow'`

Confirm your nextflow installation with:
```
which nextflow
```
## Environments

### Singularity

Singularity is a container management system, that is the recommended method of enviornment management with this pipeline. It is preinstalled on Sherlock, so no further action is needed.

### Docker

Docker is another container management system, but it cannot run on Sherlock (or most HPCs) because it requires root access. However, Docker can be used for runs on local or cloud environments.

### Conda

For those who choose to use conda over Singularity/Docker, conda must be [installed](https://docs.conda.io/en/latest/miniconda.html) and available (miniconda is fine, and is less bloated than anaconda).

# Setting Up Your First Run

Generally, the pipeline is launched from a head node, which will launch child processes for each step/branch of the pipeline. Thus, the pipeline can be launched in one of several ways:
1. Interactively
    1. (Recommended but optional) Start a [tmux](https://www.hamvocke.com/blog/a-quick-and-easy-guide-to-tmux/) session, so that your head node is not interrupted by timeouts or signal iterruptions
    2. Start a [interative/compute](https://www.sherlock.stanford.edu/docs/user-guide/running-jobs/#interactive-jobs) node on Sherlock
        1. Even though the head node requires very little memory, it is not recommended to run computation on login nodes
    3. Launch your pipeline in the command line
2. Launch the pipeline in a sbatch script (described below)

## Pipeline sbatch scripts

Here is an example of a sbatch script that would launch a nextflow pipeline. Note that the time is set to max allowed on Sherlock, with a relatively small memory request.

**IMPORTANT:**
You cannot have more than one nextflow pipeline running in the same directory! This is because there are hidden files being created which a pipeline is launched, so they cannot be overwritten. Thus, if you are running more than one pipeline at once, you **MUST** create separate work directories into which you can launch each pipeline.

This sbatch script is doing the following:
1. Set up a work directory
2. Navigate into the work directory
3. Launch a pipeline
    1. We can run pipelines by providing the account and repo name of pipeline.
        1. So for example, `nextflow run kaitlinchaung/nomad` refers to the pipeline that is stored at https://github.com/kaitlinchaung/nomad
    2. Use the following profiles (which are defined in `nextflow.config`):
        1. `test`, which runs the test dataset that is stored on github
        2. `horence`, which calls the parameters required to run this pipeline on sherlock on the horence, quake, and owners partitions
        3. `singularity`, which calls the containter required to use Singularity to run this pipeline. Note that this is interchangeable with `docker` or `conda`
    3. Running the `main` branch of the repo (`-r main`)
    4. Running the latest revision/commit of the repo (`-latest`)
    5. Turn on the resume flag as a safe default(`-resume`)
        1. If this is the first you are running the pipeline, it will provide a warning message that can be overlooked
        2. If this is NOT the first time you are running the pipeline, it will [resume](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html) at the step of the last code change or input change.


```
#!/bin/bash
#
#SBATCH --job-name=launch_nf
#SBATCH --output=launch_nf.%j.out
#SBATCH --error=launch_nf.%j.err
#SBATCH --time=48:00:00
#SBATCH -p owners
#SBATCH --nodes=1
#SBATCH --mem=4Gb
#SBATCH --requeue

work_dir=${PWD}/test_run
mkdir -p ${work_dir}
cd ${work_dir}

nextflow run kaitlinchaung/nomad \
    -profile test,horence,singularity \
    -r main \
    -latest \
    -resume
```

# Setting Up Other Runs

Here is an example of an sbatch script that runs a real dataset, with all default parameters (see: README):

**NOTE** Pipeline parameters (detailed in README) always get a double hyphen, while pipeline run settings get a single hyphen. If no pipeline parameters are passed, default parameters are used (detailed in README and `nextflow.config`).

```
#!/bin/bash
#
#SBATCH --job-name=launch_nf
#SBATCH --output=launch_nf.%j.out
#SBATCH --error=launch_nf.%j.err
#SBATCH --time=48:00:00
#SBATCH -p owners
#SBATCH --nodes=1
#SBATCH --mem=4Gb
#SBATCH --requeue

samplesheet=/home/users/kaitlin3/data/samplesheet_viral.csv
outdir=/home/users/kaitlin3/results/viral

work_dir=${PWD}/viral
mkdir -p ${work_dir}
cd ${work_dir}


nextflow run kaitlinchaung/nomad \
    -profile horence,singularity \
    -r main \
    -latest \
    -resume \
    --input ${samplesheet} \
    --output ${outdir}
```


Here is an sbatch script, that modifies a couple of pipeline parameters. Please note that Groovy booleans are lower case.
```
#!/bin/bash
#
#SBATCH --job-name=launch_nf
#SBATCH --output=launch_nf.%j.out
#SBATCH --error=launch_nf.%j.err
#SBATCH --time=48:00:00
#SBATCH -p owners
#SBATCH --nodes=1
#SBATCH --mem=4Gb
#SBATCH --requeue

samplesheet=/home/users/kaitlin3/data/samplesheet_viral.csv
outdir=/home/users/kaitlin3/results/viral

work_dir=${PWD}/viral
mkdir -p ${work_dir}
cd ${work_dir}


nextflow run kaitlinchaung/nomad \
    -profile horence,singularity \
    -r main \
    -latest \
    -resume \
    --input ${samplesheet} \
    --output ${outdir} \
    --K_num_hashes 200 \
    --L_num_random_Cj 100 \
    --run_pvals_only true \
    --gene_bed genes.bed
```

# Troubleshooting

## Repo issues
According to nextflow docs:
>When you launch a script execution with Nextflow, it will look for a file with the pipeline name you’ve specified. If that file does not exist, it will look for a public repository with the same name on GitHub (unless otherwise specified). If it is found, the repository is automatically downloaded to your computer and executed. This repository is stored in the Nextflow home directory, that is by default the $HOME/.nextflow path, and thus will be reused for any further executions.

Sometimes, there can be issues with commits not being up to date in the local repo. This can happen when you are getting errors that are not in sync with the version of the pipeline you are running (ie your version seems to be a commit behind the remote version), or if you get an error like this:
```
Unknown error accessing project `kaitlinchaung/nomad` -- Repository may be corrupted: /home/users/kaitlin3/.nextflow/assets/kaitlincha
ung/nomad
```

In that case, you can always refresh your local version with:
```
rm /home/users/kaitlin3/.nextflow/assets/kaitlinchaung/nomad
nextflow pull kaitlinchaung/nomad
```

Note that the above may be different for you, depending on where nextflow is installed.
