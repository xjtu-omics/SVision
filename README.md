## SVision

SVision is a deep learning-based structural variants caller that takes aligned reads or contigs as input. Especially, SVision implements a targeted multi-objects recognition framework, detecting and characterizing both simple and complex structural variants from three-channel similarity images.

<img src="https://github.com/xjtu-omics/SVision/blob/master/supports/workflow.png" alt="SVision workflow" align=center/>


Please check the wiki page for more details. 


## Install and run

### Install from source
Step1: Create a python environment with conda

```
conda create -n svision-env python=3.6
```

Step2: Install required packages:

```
conda install pysam
conda install tensorflow
conda install opencv-python
```
Step3: Install from source code

```
git clone https://github.com/songbowang125/SVision.git
cd SVision
python setup.py install
```

### Install from Anaconda

To be continue ...


### Usage

```
SVision [parameters] -o <output path> -b <input bam path> -g <reference> -m <model path>
```


## Change Logs

2021.1.3 (v1.2.1):

1. Fixing insertion length for detailed breakpoints.

2020.11.16 (v1.2)

1. Adding function for calling from minimap2 aligned BAM, where CIGAR operator is different from NGMRL.
2. Adding Graph representation for detected complex structural variants.
3. Adding a prameter for detecting from contig aligned BAMs.

2020.10.18 (v1.1.6)

1. Making changes to the formation mechanism inference module.
2. Adding GT, DV and DF to the standard VCF output.

2020.9.17 (v1.1.5)

1. Fixed bug: The function process_cigars() in collect_signatures.py affect the breakpoints' precision of short (about 50bp) DEL and INS.

2020.8.26 (v1.1.4)

1. Adding a SV formation mechanism inference module.

2020.8.26 (v1.1.3)

1. Adding internal breakpoints refine module.

2020.8.12 (v1.1.2)

1. Fixing bug while processing alternative contigs, such as chrUn_JTFH01001938v1_decoy
2. Adding breakpoint left shift operation
3. Fixing bug while distinguish major and minor segments at src/analyze_reads.py line 20

2020.7.20 (v1.0.3)

1. Adding log file in the output

## Contact
If you have any questions, please feel free to contact: jiadong324@gmail.com, songbowang125@163.com
