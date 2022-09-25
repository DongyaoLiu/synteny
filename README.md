# synteny

#### This pipeline is designed to visualization of small range of mulit-species whole genome alignment. Base on the file MAF. 

## Usage 

### make index of maf
if you need to creat index file of your own .maf
```shell
  python3 01.maf_make_index.py <your.maf.index>
```
### truncate maf base on coordinates and scaffold
Base on the index, we need to extract the lines in maf which contain expect alignment information. This script only could do truncate base on ref genome. If you wannt to truncate base on other genomes. DO hal2maf with parameter of your wanted genome. And follow this from begining. <start> <end> should be offered as gtf(1 based and fullly closed coordinates). The out is stand out. And all following script will take <species_Chr_start_end.maf> as input.
```shell
  python3 02.truncated_maf.py <your.maf.index> <Chr> <start> <end> > <species_Chr_start_end.maf>
```
### extract all verse all (ava) alignment relationship base on maf
Becase I add gap detection in the alignment of a single block. Thus, it will offer the nogap alignment coodinates ot .nogap.ava file
```shell
  python3 04.syntny_data_prepare.py <species_Chr_start_end.maf> <species_Chr_start_end.ava> <species_Chr_start_end.nogap.ava>
 ```
### extract annotation from gtf base on the ava
```shell
  python3 <species_Chr_start_end.ava> <species_Chr_start_end.ana>
```
### visualization
```shell
  Rscript synteny.R
```
