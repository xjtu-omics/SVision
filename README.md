<img src="https://github.com/xjtu-omics/SVision/blob/master/supports/svision-logo.png" alt="svision_logo" width="30%" height="30%" align=center/>


SVision is a deep learning-based structural variants caller that takes aligned reads or contigs as input. Especially, SVision implements a targeted multi-objects recognition framework, detecting and characterizing both simple and complex structural variants from three-channel similarity images.

<img src="https://github.com/xjtu-omics/SVision/blob/master/supports/workflow.png" alt="SVision workflow" width="60%" height="60%" align=center/>


## License

SVision is free for non-commercial use by academic, government, and non-profit/not-for-profit institutions. A commercial version of the software is available and licensed through Xi’an Jiaotong University. 
For more information, please contact with Jiadong Lin (jiadong324@stu.xjtu.edu.cn) or Kai Ye (kaiye@xjtu.edu.cn).

## Install

Step1: Create a python environment with conda

```
conda create -n svision-env python=3.6
```
Step2: Install required packages of specific versions

```
conda install -c anaconda pysam==0.16.0
conda install -c conda-forge opencv==4.5.1
conda install -c conda-forge tensorflow==1.14.0
```
Step3: Install SVision from PyPI

```
pip install SVision
```

(Optional) Install from source code

```
git clone https://github.com/xjtu-omics/SVision.git
cd SVision
python setup.py install
```

## Usage

```
SVision [parameters] -o <output path> -b <input bam path> -g <reference> -m <model path>
```

Please check the [wiki](https://github.com/xjtu-omics/SVision/wiki) page for more usage details. 

#### Input/output parameters

```
-o OUT_PATH           Absolute path to output
-b BAM_PATH           Absolute path to bam file
-m MODEL_PATH         Absolute path to CNN predict model
-g GENOME             Absolute path to your reference genome (.fai required in the directory)
-n SAMPLE             Name of the BAM sample name
```

```-g``` path to the reference genome, the index file should under the same directory.

```-m``` path to the pre-trained deep learning model, which is available at https://drive.google.com/drive/folders/1j74IN6kPKEx9hy3aENx3zHYPUnyYWGvj?usp=sharing.

#### General parameters
```
-t THREAD_NUM         Thread numbers [1]
-s MIN_SUPPORT        Min support read number for an SV [1]
-c CHROM              Specific region to detect, format: chr1:xxx-xxx or 1:xxx-xxx
--hash_table          Activate hash table to align unmapped sequences
--cluster_callset     Cluster original callset to merge uncovered event
--report_mechanism    Report mechanisms for DEL event
--report_graph        Report graph for events
--contig              Activate contig mode
```

```--hash_table``` enables the image subtraction process, which is activated by default. 

```--report_graph``` enables the program to create the CSV graph in GFA format, which is not activated by default. 

```--report_mechanism``` is used to infer the formation mechansim according to the breakpoint sequence features. 
This is still underdevelopment, which is not recommended to use for current version.

```--contig``` is used for calling from assemblies, which currently uses minimap2 aligned BAM file as input.

#### Other parameters

```--partition_max_distsance``` maximum distance allowed of a group of feature sequences.

```--cluster_max_distance``` maximum distance for feature sequence clustering. This is implemented via Scipy hierarchical clustering.

```--k_size``` size of kmer used in hash-table realignment, only used when ```--hash_table``` is activated.

```--min_accept``` minimum matched segment length, default is 50bp.

## Contact
If you have any questions, please feel free to contact: jiadonglin324@163.com, songbowang125@163.com
