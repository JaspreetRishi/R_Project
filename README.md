## Part I:  Data inspection

Loading (tidyverse) package:
  

library(tidyverse)


## Downloading files


# downloaded files from github into Assignment2 folder
genotypes_data<-read_tsv("fang_et_al_genotypes.txt")

snp_data<-read_tsv("snp_position.txt")


## Data Inspection


# Load data:
str(genotypes_data)
str(snp_data)
unique(genotypes_data$Group)
unique(snp_data$Chromosome)
nrow(genotypes_data)
nrow(snp_data)
ncol(dgenotypes_data)
ncol(snp_data)
dim(genotypes_data)
dim(snp_data)

### Data Processing

#Rearranged the snp_data file so that the Chromosome column is followed by position, followed by transposition 

snp_data<-snp_data[c(1,3,4,2,5:15)]

#Subset the genotypes_data dataframe for the maize and teosinte group seperately

maize_genotypes_data<-genotypes_data %>% filter(Group=="ZMMIL"|Group=="ZMMLR"|Group=="ZMMMR") 

teosinte_genotypes_data<-genotypes_data %>%filter(Group=="ZMPJA"|Group=="ZMPIL"|Group=="ZMPBA")

### Transposed the subset of genotypes_data to merge with the snp_data data frame


maize_genotypes_data<-as.data.frame(t(maize_genotypes_data), stringsAsFactors = F)
teosinte_genotypes_data<-as.data.frame(t(teosinte_genotypes_data), stringsAsFactors = F)


#Column names transposed to row names. Created another row from the row names called SNP_ID. First row converted to row names and removed the first three columns. 
SNP_ID <- rownames(maize_genotypes_data)
rownames(maize_genotypes_data) <- NULL
maize_genotypes_data<- cbind(SNP_ID,maize_genotypes_data,stringsAsFactors = FALSE)
names(maize_genotypes_data)<-  c("SNP_ID",maize_genotypes_data[1,-1])
maize_genotypes_data <- maize_genotypes_data[-c(1,2,3),]

SNP_ID <- rownames(teosinte_genotypes_data)
rownames(teosinte_genotypes_data) <- NULL
teosinte_genotypes_data<- cbind(SNP_ID,teosinte_genotypes_data,stringsAsFactors = FALSE)
names(teosinte_genotypes_data)<-  c("SNP_ID",teosinte_genotypes_data[1,-1])
teosinte_genotypes_data <- teosinte_genotypes_data[-c(1,2,3),]

#merged the transposed genotype file with snp file

merged_maize<-merge(snp_data, maize_genotypes_data, by="SNP_ID")
merged_teosinte<-merge(snp_data, teosinte_genotypes_data, by="SNP_ID")
data<-merged_maize %>% 
  filter(Chromosome %in% c(1))

#subset the merged file according to Chromosome number(numerical values) then  sorted it based on position(ascending values). Exported as csv file based on Chrom number, order and grouping(maize/teosinte)
for(i in c(1:10)){
  data<-merged_maize %>% filter(Chromosome==i)%>% mutate(Pos=as.numeric(Position))%>%arrange(Pos)
  data$Pos<-NULL
  write.csv( data,paste0("Maize_Chromosome_",i,"_ascending.csv"), row.names = F)
  data<-merged_maize%>% filter(Chromosome==i)%>% mutate(Pos=as.numeric(Position))%>%
    arrange(-Pos)
  data$Pos<-NULL
  data[data=="?/?"]<-"-/-"
  write.csv( data,paste0("Maize_Chromosome_",i,"_descending.csv"), row.names = F)
  
  data<-merged_teosinte %>% filter(Chromosome==i)%>% mutate(Pos=as.numeric(Position))%>%arrange(Pos)
  data$Pos<-NULL
  write.csv( data,paste0("Teosinte_Chromosome_",i,"_ascending.csv"), row.names = F)
  data<-merged_teosinte %>% filter(Chromosome==i)%>% mutate(Pos=as.numeric(Position))%>%
    arrange(-Pos)
  data$Pos<-NULL
  data[data=="?/?"]<-"-/-"
  write.csv( data,paste0("Teosinte_Chromosome_",i,"_descending.csv"), row.names = F)
  
}

##Part 2

#Data Visualization

#loaded ggplot
library(ggplot2)

#visualized the number of polymorphism positions in each chromosome

ggplot(data = snp_data[!is.na(as.numeric(snp_data$Chromosome)),]) +   geom_bar(mapping = aes(x = as.numeric(Chromosome), fill=Chromosome)) + scale_x_discrete(limit=c(1:10))+ labs(x = "Chromosome", y="No. of polymorphism position") 

#tidying up the data using pivot_longer

#loaded tidyr

library(tidyr)

tidy_file<-genotypes_data %>% pivot_longer(!c("Sample_ID", "JG_OTU", "Group"),names_to = "SNP_ID",values_to = "base")
    
tidy_file<-merge(tidy_file, snp_data, by="SNP_ID")
    
ggplot(data = tidy_file[!is.na(as.numeric(tidy_file$Chromosome)),]) +   geom_bar(mapping = aes(x = as.numeric(Chromosome), fill=Group)) + scale_x_discrete(limit=c(1:10))+ labs(x = "Chromosome"
                                                                                                                                                                      , y="No. of polymorphism position")

#Representation of polymorphisms contributed by each group

graph <- tidy_file
  
#Mutate Chromosone number into numeric
mutate(Chromosome=as.numeric(Chromosome)) %>% 
  
#Select required fields
select(Group, SNP_ID, Chromosome, base) %>%  
  
#Filter chromosome 1 to 10
filter(Chromosome %in% c(1:10)) %>%
  
#remove all unknown bases  
filter(base!="?/?")%>%
  
  #group by Group and SNP_ID  
group_by(Group,SNP_ID) %>% 
  
  #remove all SNPS in a group that has only one base i.e no polymorphism in the group 
filter(length(unique(base))>1) %>%  
  
  #select requred fields
select(Group, SNP_ID, Chromosome) 

graph<-graph[!duplicated(graph),] #remove duplicated records

#plotted the filtered data based on polymorphism contributed by each group
ggplot(data=graph) +
geom_bar(mapping=aes(x=Chromosome, fill=Group)) + 
scale_x_discrete(limit=c(1:10), label=c(1:10))

#created a new column based on homozygousity of the polymorphism

tidy_file$homozygous<-TRUE
tidy_file$homozygous[tidy_file$base=="?/?"]<-NA

tidy_file$homozygous[substr(tidy_file$Base,1,1)!=substr(tidy_file$base,3,3)]<-FALSE

# Graphed the SNPs by Group, filling by Homozygosity:
graphs<-tidy_file
ggplot(data=tidy_file) +
  geom_bar(mapping=aes(x=Group, fill=homozygous), position="fill")

#Graphs based on polymorphism in every chromosome
fig.width=10
fig.height=11


SNP_graph<-tidy_file[!is.na(as.numeric(tidy_file$Chromosome)),]
SNP_graph$Chromosome<-as.numeric(SNP_graph$Chromosome)
SNP_graph<-SNP_graph%>% select(SNP_ID, Group, base, Position, Chromosome, homozygous)%>% filter(base!="?/?")%>% unique()%>%filter(!is.na(as.numeric(Position)))

final_graph<-ggplot(data = SNP_graph) + geom_point(mapping=aes(x=as.numeric(Position), y=Group, color=homozygous)) +labs(y = "Group" , x="Chromosome position")

final_graph + facet_wrap(~ Chromosome,ncol=2, scales = "free") +labs(color="Genotype")


 
  



