
# Index_investigator
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Index switching is potentially a large problem in Illumina HiSeq data. These scripts are designed to test if index switching is occuring in your dataset. 
#### IMPORTANT POINT: 
This will give a false positive signal if genetic grouping corresponds to lane grouping. For example, if populations are intentionally sequenced together in the same lane. Additionally, this requires samples sequenced on multiple different lanes since it relies on greater similarity for samples sequenced on the same lane, versus other lanes. 

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Installing

These scripts need no installation. Just download the perl scripts using git.

```
git clone https://github.com/owensgl/index_investigator.git
```

### Input files

Both scripts require two input files:
1. A vcf file containing SNPs. This script was written based on vcfs from freebayes v1.1.0 and may not correctly parse files from other programs.
2. An tab separated info file containing 4 columns: Samplename, Lane1, Lane2, Machine. There are two lanes in case a single sample was sequenced on multiple lanes. If it was only sequenced on one lane, please repeat the same identifier in both columns.

***

## vcf2indexswitcher.pl
This script takes a vcf file (formatted from Freebayes), an info file that tells it what technology (i.e. sequencing machine) and lane each sample was sequenced on, and a decimal fraction of reads to switch (i.e. 0.05 for 5%) 
It bioinformatically switches n percent of reads to different samples of the same lane. Genotypes are recalled with the new read depths.

### Options:
* info=FILENAME; The name of your info file.
* max_sites=INTEGER; The number of sites to process before stopping.
* switch_rate=SCALAR; The minimum read depth to consider an unbalanced heterozygote (0 < N < 1).
### Example:
```
zcat < example_data.vcf.gz | perl ./vcf2indexswitcher.pl --info example_infofile.txt --max_sites 1000 --switch_rate 0.01 > example_data.switched.vcf
```
### Output:
A vcf file fit for using on vcf2indexinvestigator.pl. Note: Much metadata has been stripped from this vcf file.

***

## vcf2indexinvestigator.pl
This script takes a vcf file (formatted from Freebayes) and an info file that tells it what technology (i.e. sequencing machine) and lane each sample was sequenced on. 
It will only use di-allelic SNPs from the vcf file. It looks for unbalanced heterozygotes where one allele has one read and the other has multiple. By default genotypes need to have >=5 reads to be considered unbalanced.

Options:
* info=FILENAME; The name of your info file.
* max_sites=INTEGER; The number of sites to process before stopping. More is better, but around 500,000 you're maxing out.
* min_dp=INTEGER; The minimum read depth to consider an unbalanced heterozygote.

```
zcat < example_data.vcf.gz | perl ./vcf2indexinvestigator.pl --info example_infofile.txt --max_sites 1000 --min_dp 10 > out.txt
```
### Output:
A text file with nine columns and two rows per unbalanced heterozygote.
1. site: The position in the genome.
2. sample: The sample with the unbalanced heterozygote.
3. technology: The machine used to sequenced the sample.
4. lane1: The lane used to sequence the sample.
5. lane2: The other lane used to sequence the sample.
6. depth: The read depth of the unbalanced heterozygote.
7. percent: The predicted percent of allele sharing within the lane (p_hat).
8. type: Whether it's testing samples within the lane or control samples outside the lane.
9. value: Whether there is allele sharing (1) or not (0).

Each site is represented by two rows, one testing within the lane and another testing a control set of samples outside the lane.

***

## plot_indexinvestigator.R
A simple R script is provided to plot the output of vcf2indexinvestigator.pl. It takes two arguments, the name of the output file vcf2indexinvestigator.pl produces and the name of your final resulting pdf.

*NOTE*: This script requires the tidyverse package, so it will attempt to install it if it is not already installed. This may take a few minutes.

```
Rscript plot_indexinvestigator.R out.txt out.pdf 
```


***

## Authors

* **Gregory L. Owens** 


## License

This project is licensed under the MIT License.
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)



