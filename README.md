# index_switching
Scripts and dataset for index switching paper by Owens et al., 2017

There are two main programs in this repository. 
-vcf2laneallelesharing.pl: This script takes a vcf file (formatted from Freebayes) and an info file that tells it what technology (i.e. sequencing machine) and lane each sample was sequenced on. 
It will only use di-allelic SNPs from the vcf file. It looks for unbalanced heterozygotes where one allele has one read and the other has multiple. By default genotypes need to have >=5 reads to be considered unbalanced.
-vcf2barcodeslippage.pl: This script takes a vcf file (formatted from Freebayes), an info file that tells it what technology (i.e. sequencing machine) and lane each sample was sequenced on, and a decimal fraction of reads to switch (i.e. 0.05 for 5%) 
It bioinformatically switches n percent of reads to different samples of the same lane. Genotypes are recalled with the new read depths.
