<img src="https://github.com/xjtu-omics/SVision/blob/master/supports/svision-logo.png" alt="svision_logo" width="30%" height="30%" align=center/>


SVision is a deep learning-based structural variants caller that takes aligned reads or contigs as input. Especially, SVision implements a targeted multi-objects recognition framework, detecting and characterizing both simple and complex structural variants from three-channel similarity images.

<img src="https://github.com/xjtu-omics/SVision/blob/master/supports/workflow.png" alt="SVision workflow" width="60%" height="60%" align=center/>


Please check the [wiki](https://github.com/xjtu-omics/SVision/wiki) page for more details. 


## License

SVision is free for non-commercial use by academic, government, and non-profit/not-for-profit institutions. A commercial version of the software is available and licensed through Xi’an Jiaotong University. 
For more information, please contact with Jiadong Lin (jiadong324@stu.xjtu.edu.cn) or Kai Ye (kaiye@xjtu.edu.cn).

## Install and run

### Install from PyPI

Step1: Create a python environment with conda

```
conda create -n svision-env python=3.6
```
Step2: Install deep-learning related packages

```
conda install -c conda-forge opencv==4.5.1
conda install -c conda-forge tensorflow==1.14.0
```


### Install from source
Step1: Create a python environment with conda

```
conda create -n svision-env python=3.6
```
step2: Install basic packages
```
conda install -c anaconda scipy, pysam, numpy, beautifulsoup4
```

Please install numpy=1.16.4 to avoid feature warnings raised by tensorflow

Step3: Install deep-learning related packages

```
conda install -c conda-forge opencv==4.5.1
conda install -c conda-forge tensorflow==1.14.0
```
Step4: Install from source code

```
git clone https://github.com/xjtu-omics/SVision.git
cd SVision
python setup.py install
```

### Usage

```
SVision [parameters] -o <output path> -b <input bam path> -g <reference> -m <model path>
```

## Change Logs

**V1.2.1**

Fixing insertion length for detailed breakpoints.

**V1.2**

1. Adding function for calling from minimap2 aligned BAM, where CIGAR operator is different from NGMRL.
2. Adding Graph representation for detected complex structural variants.
3. Adding a prameter for detecting from contig aligned BAMs.

**V1.1.6**

1. Making changes to the formation mechanism inference module.
2. Adding GT, DV and DF to the standard VCF output.

**V1.1.5**

Fixed bug: The function process_cigars() in collect_signatures.py affect the breakpoints' precision of short (about 50bp) DEL and INS.

**V1.1.4**

Adding a SV formation mechanism inference module.

**V1.1.3**

Adding internal breakpoints refine module.

**V1.1.2**

1. Fixing bug while processing alternative contigs, such as chrUn_JTFH01001938v1_decoy
2. Adding breakpoint left shift operation
3. Fixing bug while distinguish major and minor segments at src/analyze_reads.py line 20


## Contact
If you have any questions, please feel free to contact: jiadonglin324@163.com, songbowang125@163.com
