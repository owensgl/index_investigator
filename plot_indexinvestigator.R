#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if(!require(tidyverse)){
    install.packages("tidyverse",repos = "http://cran.us.r-project.org")
    library("tidyverse")
}
input_file = args[1]
output_file = args[2]

data <- read.delim(input_file,header=T,row.names=NULL,sep="\t")

pdf(output_file)
data %>% 
  ggplot(.,aes(x=percent,y=value)) +
  geom_smooth(aes(color=type),method="auto") + 
  coord_cartesian(ylim=c(0,1)) + geom_abline() +
  ylab(expression(paste( italic(p),": observed allele sharing" ))) + 
  xlab(expression(paste( italic(hat(p)),": predicted allele sharing" )))
dev.off()
