<img src="https://github.com/xjtu-omics/SVision/tree/master/supports/svision-logo.png" alt="svision_logo" width="30%" height="30%" align=center/>


SVision is a deep learning-based structural variants caller that takes aligned reads or contigs as input. Especially, SVision implements a targeted multi-objects recognition framework, detecting and characterizing both simple and complex structural variants from three-channel similarity images.

<img src="https://github.com/xjtu-omics/SVision/tree/master/supports/workflow.png" alt="SVision workflow" width="60%" height="60%" align=center/>


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
conda install -c bioconda pysam
conda install -c conda-forge opencv==4.5.1
conda install -c conda-forge tensorflow==1.14.0
```
Step3: Install SVision from source

```
python setup.py install
```

The Pip and Conda install would be available after the **-beta* version.

## Usage

```
SVision [parameters] -o <output path> -b <input bam path> -g <reference> -m <model path>
```
Check all parameters with

```
SVision -h
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

```-m``` path to the pre-trained deep learning model *svision-cnn-model.ckpt* ([download link](https://drive.google.com/drive/folders/1j74IN6kPKEx9hy3aENx3zHYPUnyYWGvj?usp=sharing)). Please download all files and save them to your local directory.


#### General parameters
```
-t THREAD_NUM         Thread numbers [1]
-s MIN_SUPPORT        Min support read number for an SV [1]
-c CHROM              Specific region to detect, format: chr1:xxx-xxx or 1:xxx-xxx
--hash          Activate hash table to align unmapped sequences
--qname    Report support names for each events
--cluster     Cluster calls that might occur together
--graph        Report graph for events
--contig              Activate contig mode
```

```--hash``` enables kmer based alignment for unmapped sequences. 

```--graph``` enables the program to create the CSV graph in GFA format. This function requires ```--qname``` enabled.

```--contig``` is used for calling from assemblies, which currently uses minimap2 aligned BAM file as input.

## SVision output

### VCF

The SV ```ID``` column is given in the format of ```a_b```, where ```b``` indicates site ```a``` contains other type of SVs. 

Filters used in the output.

```Covered```: The entire SV is spanned by long-reads, producing the most confident calls.

```Uncovered```: SV is partially spanned by long-reads, i.e. reads spanning one of the breakpoints.

```Clustered```: SV is partially spanned by long-reads, but can be spanned through reads clusters.

We add extra attributes in the ```INFO``` column of VCF format for SVision detected structural variants.

```BRPKS```: The CNN recognized breakpoint junctions through tMOR.

```GraphID```: The graph index used to indicate the graph structure, which requires ```--report_graph``` and is obtained by calculating isomorphic graphs. 
The ID for simple SVs is -1. 

```VAF```: The estimated variant allele fraction, which is calculated by DV/DR. Note that SVision does not calculate the genotypes in the current version.

### CSV graph 

#### CSV graph compare

SVision classify the graph of each CSV instances by comparing their graph topologies. 
It will create two text (.txt) file along with the VCF output.

1. ```sample.graph_exactly_match.txt```: Unique graphs for all CSV instances.
2. ```sample.graph_symmetry_match.txt```: Topological similar graph found from unique graphs.

#### Graph format

The graph output requires ```--graph``` activated. 
The below example is an CSV in rGFA format (node sequence is omitted for display purpose), which is detected by SVision at chr11:99,819,283-99,820,576 in HG00733. 
The graph output is saved in separated files for each CSV events.

```
S	S1	SN:Z:chr11	SO:i:99819338	SR:i:0	LN:i:2990
S	I0	SN:Z:m54329U_190827_173812/140708091/ccs	SO:i:15813	SR:i:0	LN:i:1113
S	I1	SN:Z:m54329U_190827_173812/140708091/ccs	SO:i:16927	SR:i:0	LN:i:466
S	I2	SN:Z:m54329U_190827_173812/140708091/ccs	SO:i:17400	SR:i:0	LN:i:377	DP:S:S1:99820198
S	I3	SN:Z:m54329U_190827_173812/140708091/ccs	SO:i:17778	SR:i:0	LN:i:838
S	I4	SN:Z:m54329U_190827_173812/140708091/ccs	SO:i:18617	SR:i:0	LN:i:61	DP:S:S0:99819276
L	S0	+	I0	+	0M	SR:i:0
L	I0	+	I1	+	0M	SR:i:0
L	I1	+	I2	-	0M	SR:i:0
L	I2	-	I3	+	0M	SR:i:0
L	I3	+	I4	+	0M	SR:i:0
L	I4	+	S1	+	0M	SR:i:0
```

Besides the information included in standard [rGFA](https://github.com/lh3/gfatools/blob/master/doc/rGFA.md) format, 
we add another ```DP:S``` column to indicate sequence with detected origins via local realignment, 
such as node ```I2``` is duplicated from node ```S1```. 

#### Graph genotyping

**Note**: This is a post-processing step that tries to validate the detected CSVs. 

**Step1: Extract HiFi raw reads**

```
samtools view -b HG00733.ngmlr.sorted.bam chr11:99810000-99830000 > tmp.bam
samtools fasta tmp.bam > tmp.fasta
```

**Step2: Align with GraphAligner**

Please check [GraphAligner](https://github.com/maickrau/GraphAligner) for the detailed usage.

```
GraphAligner -g chr11-99819283-99820576.gfa -f tmp.fasta -a aln.gaf -x vg
```

Example of CSV path supporting reads

```
m54329U_190827_173812/140708091/ccs     21668   0       21668   +       >S0>I0>I1<I2>I3>I4>S1
m54329U_190617_231905/88145984/ccs      13612   0       13612   +       >S0>I0>I1<I2>I3>I4>S1
m54329U_190617_231905/88145984/ccs      13612   0       13612   +       >S0>I0>I1<I2>I3>I4>S1
```


## Contact
If you have any questions, please feel free to contact: jiadong66@stu.xjtu.edu.cn, songbowang125@163.com
