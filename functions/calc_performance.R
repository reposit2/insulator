# This code was created based on codes in 'Tasaki, S., Gaiteri, C., Mostafavi, S. & Wang, Y. Deep learning decodes the principles of differential gene expression. Nature Machine Intelligence (2020)'

outloc = commandArgs(trailingOnly=TRUE)[1] 
options(stringsAsFactors = FALSE)
library(tidyverse)
library(ggpubr)
library(broom)

#message(outloc)
#message(dirname(outloc))

#sample_names = read_delim(paste0(dirname(outloc),"/feature_norm_stats.txt"),
sample_names = read_delim(paste0(outloc,"/feature_norm_stats.txt"),
                          delim = ",",col_names = T)%>%
  filter(feature_type=="deg_stat")%>%.$feature_name
geneids = read_delim(paste0(outloc,"/test_data/geneid.txt.gz"),
                     delim = "\t",col_names = F)%>%.$X1

actual = read_delim(paste0(outloc,"/test_data/actual.txt.gz"),
                    delim = "\t",col_names = F)
colnames(actual)=sample_names

prediction = read_delim(paste0(outloc,"/test_data/prediction.txt.gz"),
                        delim = "\t",col_names = F)
colnames(prediction)=sample_names

actual%>%
  as.data.frame()%>%
  mutate(geneid=geneids)%>%
  gather(sample,actual,-geneid)%>%
  inner_join(.,
             prediction%>%
               as.data.frame()%>%
               mutate(geneid=geneids)%>%
               gather(sample,prediction,-geneid),
             by=c("geneid","sample"))%>%
  group_by(sample)%>%
#  tryCatch({
      do(cor.test(.$actual,.$prediction,method = "spearman")%>%broom::tidy()) -> cor_tbl
#  do(cor.test(.$actual,.$prediction,method = "pearson")%>%broom::tidy()) -> cor_tbl
#    },
#    error = function(e) {
#      message(e)
#    },
#    warning = function(e) {
#      message(e)
#    },
#    silent = TRUE
#  )

#  tryCatch({
      write.table(cor_tbl,file=paste0(outloc,"/test_data/cor_tbl.txt"),
            append = F,quote = F,sep = "\t",row.names = F,col.names = T)  
#    },
#    error = function(e) {
#      message(e)
#    },
#    warning = function(e) {
#      message(e)
#    },
#    silent = TRUE
#  )
