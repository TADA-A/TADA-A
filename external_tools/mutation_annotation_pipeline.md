
# Annotation of allele specific base level mutation rate

This [SoS workflow](https://vatlab.github.io/sos-docs/) is used to generate allele-specific base level mutation rate.

## Software tools

We have developed a
[Docker container](https://hub.docker.com/r/yuwenliu/tadaa-tools) that
includes all software components necessary to run the analyses.
If you do not have [Docker](https://www.docker.com/community-edition) 
please download and install it to your system, following the instructions
provided on the Docker website. Once you have installed Docker, check
that Docker is working correctly by following [this tutorial](https://docs.docker.com/get-started).

If your docker works, you need to setup your terminal with this `alias` command (run it in terminal):

```
alias tadaa-tools='docker run --rm --security-opt label:disable -t -P -h "TADA-A Tools" '\
'-w $PWD -v $HOME:/home/$USER -v /tmp:/tmp -v $PWD:$PWD '\
'-u $UID:${GROUPS[0]} -e HOME=/home/$USER -e USER=$USER yuwenliu/tadaa-tools'
```

Then run:

```
tadaa-tools uname -sn
```

This command will download the Docker image if it is the first time you run it.

*Note:* If you get error "Cannot connect to the Docker daemon. Is the
docker daemon running on this host?" in Linux or macOS, see
[here for Linux](https://askubuntu.com/questions/477551/how-can-i-use-docker-without-sudo)
or [here for Mac](https://github.com/wodby/docker4drupal/issues/15) for
suggestions on how to resolve this issue.

## Run the mutation rate pipeline

We provide workflows to download hg19 reference genome that will be triggered by the default analysis pipeline.
To run the analysis, under [the TADA-A github root repo](https://github.com/TADA-A/TADA-A) for example after you download it,

```
tadaa-tools sos run external_tools/mutation_annotation_pipeline.ipynb \
    --window_file test_data/test_windows.txt \
```

This will reproduce the mutation files we have provided for a default TADA-A run. 

To check on what command were executed exactly:

```
tadaa-tools sos dryrun external_tools/mutation_annotation_pipeline.ipynb \
    --window_file test_data/test_windows.txt
```

Optionally you can also provide mutation file reference etc, for example:

```
    --mutation_ref test_data/fordist_1KG_mutation_rate_table.txt
    --hg_summary test_data/hg19.genome
```

For more advance usage please checkout [SoS documentation](https://vatlab.github.io/sos-docs/) if you have troubles customizing this pipeline.

## Pipeline in detail


```sos
[global]
parameter: wd = path('./mutation_pipeline_workdir')
parameter: resource_dir = f"{wd:a}/hg19"
parameter: window_file = file_target()
parameter: mutation_ref = file_target('/data/fordist_1KG_mutation_rate_table.txt')
parameter: hg_summary = file_target('/data/hg19.genome.txt')
ref_fa = "hg19.fasta"
```

### Step 0: prepare human reference genome data


```sos
[hg19: provides = file_target(f"{resource_dir}/{ref_fa}")]
ucsc_url = "http://hgdownload.cse.ucsc.edu"
output: f"{resource_dir}/{ref_fa}"
download: dest_dir = resource_dir, expand = True
    {ucsc_url}/goldenPath/hg19/bigZips/hg19.2bit
bash: expand = True
    twoBitToFa {resource_dir}/hg19.2bit {_output}
```


```sos
Export workflow to Markdown file:
```


```sos
[0]
input: [item for item in paths(sys.argv) if item.suffix == '.ipynb'], group_by = 1
output: [f'{wd:a}/{item:bn}.md' for item in paths(sys.argv) if item.suffix == '.ipynb'], group_by = 1
bash: expand = True, stderr = False
  sos convert {_input} {_output}
```

### Step 1: prepare extended genomic windows file


```sos
[1]
fail_if(not window_file.is_file(), msg = f'Please provide valid window file via ``--window_file``')
fail_if(not mutation_ref.is_file(), msg = f'Please provide valid mutation reference file via ``--mutation_ref``')
depends: file_target(f"{resource_dir}/{ref_fa}")
input: f'{window_file:a}'
output: f'{wd:a}/{_input:bn}.mutrate.bed'
bash: expand = "${ }"
    sed '1d' ${_input} | awk {'print $1"\t"$2-1"\t"$3+1"\t"$4'} > ${_output}
```

### Step 2: Get the nucleotide sequence of each interval in tab format


```sos
[2]
depends: executable("bedtools")
output: f'{_input:n}.fasta'
bash: expand = True
    bedtools getfasta -fi {resource_dir}/{ref_fa} -bed {_input} -fo {_output} -tab
```

### Step 3: Use the output file to extract tri-nuleotide sequence of each base within the window intervals


```sos
[3]
depends: executable('tri_extract_for_TADA-A.py'), executable('MutRateBase_for_TADA-A_v2.py')
output: f'{_input}.tri'
bash: expand = True
    tri_extract_for_TADA-A.py {_input} > {_output}
_input.zap()
```

### Step 4:Assesing mutation rate
Use .tri file as an input file to get the allele-specific mutation rate. 


```sos
[4]
depends: executable('MutRateBase_for_TADA-A_v2.py'), mutation_ref
output: f'{_input:n}.tri.mr'
bash: expand = True
    MutRateBase_for_TADA-A_v2.py {mutation_ref:a} {_input} > {_output}
_input.zap()
```

### Step 5: generte mutrate file base on the alternative nucleotide


```sos
[5]
alt = list('GCTA')
input: for_each = 'alt', group_by = 1, concurrent = True
output: [f'{_input[0]:n}.alt_{x}.mr' for x in alt], group_by = 1
bash: expand = True
    awk '$5 == "{_alt}"' {_input} > {_output}
```

### Step 6: Generate alternative-allele-specific Wig file


```sos
[6]
depends: executable('base_mutarate_to_wiggle_file.sh')
input: group_by = 1, concurrent = True
output: f'{_input}.wiggle'
bash: expand = True
    base_mutarate_to_wiggle_file.sh {_input}
```

### Step 7: tranform Wig file to bigwig file


```sos
[7]
depends: executable('wigToBigWig'), hg_summary
input: group_by = 1, concurrent = True
output: f'{_input:n}.bw'
bash: expand = True
    wigToBigWig {_input} {hg_summary:a} {_output}
```
