![header](supports/svision-logo.png)  


SVision is a deep learning-based structural variants caller that takes aligned reads or contigs as input. 
Especially, SVision implements a targeted multi-objects recognition framework, detecting and characterizing both simple and complex structural variants from three-channel similarity images.

<img src="https://github.com/xjtu-omics/SVision/tree/master/supports/workflow.png" alt="SVision workflow" width="60%" height="60%" align=center/>


## License

SVision is free for non-commercial use by academic, government, and non-profit/not-for-profit institutions. A commercial version of the software is available and licensed through Xi’an Jiaotong University. 
For more information, please contact with Jiadong Lin (jiadong324@stu.xjtu.edu.cn) or Kai Ye (kaiye@xjtu.edu.cn).

## Installation

### Install from source

```
## Get the source code
git clone https://github.com/xjtu-omics/SVision.git
cd SVision
## Create a conda environment for SVision
conda env create -f ./environment.yml 
## Install from source
conda activate svisionenv
python setup.py install
```

The Pip and Conda install would be available later.

### Docker

#### Pull docker image

```
docker pull jiadongxjtu/svision:1.3.6
```

#### Run docker image

```
docker run jiadongxjtu/svision:1.3.6 SVision -h
```

**Note:** Please ensure you have the permission to write into docker.


## Usage

Please visit our wiki page for [performance evaluation](https://github.com/xjtu-omics/SVision/wiki/Performance-evaluation).

### Short usage

```
SVision [parameters] -o <output path> -b <input bam path> -g <reference> -m <model path>
```

#### Check all parameters with

```
SVision -h
```

Required Input/Ouput parameters

```
-o OUT_PATH           Absolute path to output
-b BAM_PATH           Absolute path to bam file
-m MODEL_PATH         Absolute path to CNN predict model
-g GENOME             Absolute path to your reference genome (.fai required in the directory)
-n SAMPLE             Name of the BAM sample name
```

```-g``` path to the reference genome, the index file should under the same directory.

```-m``` path to the pre-trained deep learning model *svision-cnn-model.ckpt* ([external download link](https://drive.google.com/drive/folders/1j74IN6kPKEx9hy3aENx3zHYPUnyYWGvj?usp=sharing)).


Please check the [wiki](https://github.com/xjtu-omics/SVision/wiki) page for more usage and parameter details.


### Run demo data

The demo data is ./supports/HG00733.svision.demo.bam. 
The HiFi whole genome sequencing data of [HG00733](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/working/20190925_PUR_PacBio_HiFi/) is published on [Science](https://www.science.org/doi/10.1126/science.abf7117?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed).

#### Run with graph option

1. Download the reference genome GRCh38

2. Run SVision with your reference

```
SVision -o ./output_dir -b ./supports/HG00733.svision.demo.bam -m /path/to/svision-cnn-model.ckpt -g ./reference.fa -n HG00733 -s 5 --graph --qname
```

Please use the same parameter settings if you use the docker image.

3. Output files

``` *.graph.vcf ``` The standard VCF output with CSV graph info columns.

```*.graph_exactly_match.txt ``` CSV graphs of exactly identical structure.

```*.graph_symmetry_match.txt``` Identified isomorphic graphs from all CSV graphs.

```graphs``` The directory for CSV graph in rGFA format.

Please check the [wiki](https://github.com/xjtu-omics/SVision/wiki) page for more details of output format. 

## Contact
If you have any questions, please feel free to contact: jiadong66@stu.xjtu.edu.cn, songbowang125@163.com
