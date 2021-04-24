####################################################################################################################
### 1. MUTATION DATA COLLECTION & PRE-PROCESSING
####################################################################################################################

#Load the dataset
df <- read.csv("INPUT_TCGA_data.txt", header = TRUE, sep='\t', stringsAsFactors=FALSE)

#remove duplicated rows
df <- unique(df)

#convert the gene names from hugo to entrez ID
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("mygene")
library(mygene)


unique_genes <- as.data.frame(unique(df$Hugo_Symbol))
unique_genes_vector <- as.vector(unique_genes$`unique(df$Hugo_Symbol)`)

df_hugo_entrez <- queryMany(unique_genes_vector, scopes="symbol", fields="entrezgene", species="human")
n <- as.data.frame(df_hugo_entrez)
nn <- n[,c(1,4)]
nnn <- nn[complete.cases(nn), ]

df_hugo_uniprot <- queryMany(unique_genes_vector, scopes="symbol", fields=c("entrezgene", "uniprot"), species="human")
s <- as.data.frame(df_hugo_uniprot)
ss <- s[,c(4,6)]
ss$uniprot.Swiss.Prot <- as.character(ss$uniprot.Swiss.Prot)
sss <- ss[complete.cases(ss), ]

#create a new dataframe with the entrez gene IDs
df_2 <- merge(df, nnn, by.x="Hugo_Symbol", by.y="query")

#Remove p.
df_2$HGVSp_Short <- gsub("^.{0,2}", "", df_2$HGVSp_Short)

#write the output to a file
write.table(df_2, "TCGA_only_missense.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote = FALSE)





#########################################################################################################
###  GENE-LEVEL FEATURE EXTRACTION (max no of PPI for HQ-interactions from Interactome INSIDER)
#########################################################################################################

HQ_interf <- read.csv("GENE_level_HQinterfaces.txt", header = TRUE, sep='\t', stringsAsFactors=FALSE)

#remove the rows which show protein self-self interactions
HQ_interf_diff <- HQ_interf[-which(HQ_interf$P1 == HQ_interf$P2),]

#keep the columns of HQ interactions
PPI_data <- unique(HQ_interf_diff[,1:2])

#convert the SWISSPROT to Entrez ID for PPI data
#extract mapping IDs from UniProt
mapping_table <-  read.csv("mapping_table_uniprot_entrez.txt", header = TRUE, sep='\t')

#convert uniprot IDs to entrez IDs
mer_1 <- merge(PPI_data, mapping_table, by.x="P1", by.y="From")
mer_2 <- merge(mer_1, mapping_table, by.x="P2", by.y="From")
#calculate the max number of interactions for each protein
PPI_pairs_by_entrez <- unique(mer_2[,c(4,3)])

#change the column names
colnames(PPI_pairs_by_entrez) <- c("P1_entrez", "P2_entrez")


#count how many numbers of interactions do each protein have with dplyr package
library(dplyr)
#summarise the maximum no of interactions per protein and write the other interactors as a vector
PPI_pairs_by_entrez <- PPI_pairs_by_entrez %>% group_by(P1_entrez)
PPI_max_no_entrez <- PPI_pairs_by_entrez %>% summarise(P2_entrez = paste(P2_entrez, collapse=","), times=count.fields(textConnection(P2_entrez), sep=",")) %>% arrange(desc(times), P2_entrez)
#write the max interactions output to a file
write.table(PPI_max_no_entrez, "max_no_inter.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote = FALSE)

#merge the total number of interactions with the mutation data & remove duplicated data
df_only_missense_and_PPI_entrez <- merge(PPI_max_no_entrez, df_2, by.x="P1_entrez", by.y="entrezgene")
#check if there are any duplicated rows -> NO
dim(unique(df_only_missense_and_PPI_entrez))

#write the output to a file 
write.table(df_only_missense_and_PPI_entrez, "mutations_and_PPI.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote = FALSE)





##########################################################################################################
### 2.2.2 GENE-LEVEL FEATURE EXTRACTION (enriched pathway for groups of genes which interact together)
##########################################################################################################

#create a new column with a vector including sets of proteins required for pathway searches
df_only_missense_and_PPI_entrez$interact_genes_entrez <- paste(df_only_missense_and_PPI_entrez$P1_entrez, df_only_missense_and_PPI_entrez$P2_entrez, sep=",")

#extract the vector to a new dataframe & replace name of the column
genes_as_vector <- unique(as.data.frame(df_only_missense_and_PPI_entrez$interact_genes_entrez))
colnames(genes_as_vector) <- c("int_ent")

#write the data to a new file
write.table(genes_as_vector, "groups_of_genes_as_vector.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote = FALSE)




#install and load package reactomePA (based on Reactome database)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ReactomePA")
library(ReactomePA)
#install other packages required for reactomePA
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)

#read the file with genes as vector
gns_as_vec <- read.csv("groups_of_genes_as_vector.txt", header = TRUE, sep='\t', stringsAsFactors = FALSE)

#load the package and split genes as columns
library(splitstackshape)
split_entrez_genes <- cSplit(gns_as_vec, "int_ent", sep=",")



#install package msingdbr for molecular signatures & dependencies
install.packages("msigdbr")
library(msigdbr)
library(clusterProfiler)
#extract molecular signatures for human by category H (hallmarks)
human_df <- msigdbr(species="Homo sapiens", category="H")
#check what types of hallmarks are present
types_hallmarks <- as.data.frame(table(human_df$gs_name))

#create a dataframe with molecular signatures for each gene
T_2_G <- human_df[,c(2,7)]




#for each row of the split_entrez_gene
for (i in 1:dim(split_entrez_genes)[1]) {
  #vectorise each row with entrez gene names so that it can inputed into ReactomePA
  vec_i <- dput(as.character(na.omit(as.vector(unlist(split_entrez_genes[i,])))))
  #if the no of genes are more than 0
  if (length(vec_i)>0) {
    #compute the enrichment pathway analysis
    p_i <- enricher(vec_i, pvalueCutoff = 0.05, pAdjustMethod = "BH",
                    universe,minGSSize = 20, maxGSSize = 500,
                    qvalueCutoff = 0.05, TERM2GENE = T_2_G)
    
    #if any results are NULL put the result as NA
    if (!is.null(p_i)) {
      #add the name of the pathway to the genes_as_vector dataframe
      genes_as_vector[i,2] <-  p_i@result$ID[1]
    }
    else {
      genes_as_vector[i,2] <- NA
    }
  }
  #if the no of genes are <0 put the result as NA
  else {
    genes_as_vector[i,2] <- NA
  }
}


#get a summary of the hallmarks
unique(genes_as_vector$V2)
length(genes_as_vector$V2)

types_hallmarks_by_freq <- as.data.frame(table(genes_as_vector$V2))
sum(types_hallmarks_by_freq$Freq)





#load the hallmarks with categories data
hallmarks_and_categ <- read.table("hallmarks_process_category.txt", sep="\t")
#add the categories to our groups of genes
genes_as_vector_4.2 <- genes_as_vector[,1:2]
genes_as_vector_4.2$V2 <- gsub("^.{0,9}", "", genes_as_vector_4.2$V2)
genes_as_vector_4_categ <- merge(genes_as_vector_4.2, hallmarks_and_categ, by.x="V2", by.y="Hallmark_name", all.x=TRUE)
#change the column names and order
colnames(genes_as_vector_4_categ)[which(names(genes_as_vector_4_categ) == "V2")] <- "Hallmark_name"
colnames(genes_as_vector_4_categ)[which(names(genes_as_vector_4_categ) == "int_ent")] <- "interact_genes_entrez"
genes_as_vector_4_categ <- genes_as_vector_4_categ[, c(
  "interact_genes_entrez", "Process_category", "Hallmark_name", 
  "Description", "Number_of_founder_sets", "Number_of_genes")]
colnames(genes_as_vector_4_categ)[which(names(genes_as_vector_4_categ) == "interct_genes_entrez")] <- "interact_ent"

genes_as_vector_5_categ <- as.data.frame(genes_as_vector_4_categ)
colnames(genes_as_vector_5_categ)[which(names(genes_as_vector_5_categ) == "interct_genes_entrez")] <- "interact_genes_Entrez"
names(genes_as_vector_5_categ)[1] <- "interact_genes_Entrez"

#merge the full datasets with the mutations
df_PPI_path <- merge(df_only_missense_and_PPI_entrez, genes_as_vector_5_categ, 
                     by.x="interact_genes_entrez", 
                     by.y="interact_genes_Entrez")


#reorder the columns of the dataframe to match the initial data structure
df_PPI_path_2 <- df_PPI_path[,c(2, 5:12,1,13, 4,15,14)]



#write the output with max no of PPI and hallmarks+ pathways to a file
write.table(df_PPI_path_2, "mutations_and_PPI_and_proc+hallm.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote = FALSE)

df_PPI_path_2 <-  read.csv("mutations_and_PPI_and_proc+hallm.txt", header = TRUE, sep='\t', stringsAsFactors=FALSE)


#make pathways as columns (dichotomous variable)
genes_and_process_unique <- unique(df_PPI_path_2[,c(10,13)])
library(data.table)
genes_and_process_as_cols <- dcast(genes_and_process_unique, interact_genes_entrez~Process_category, fill=0, fun.aggregate=length)

#create a new dataframe with pathways_as_cols
df_PPI_path_3 <- merge(df_PPI_path_2, genes_and_process_as_cols, by.x="interact_genes_entrez", by.y="interact_genes_entrez")

#keep only needed columns
df_PPI_path_4 <- df_PPI_path_3[,c(2:12,15:23,14)]
colnames(df_PPI_path_4)[which(names(df_PPI_path_4) == "NA")] <- "NA_hallmarks"
colnames(df_PPI_path_4)[which(names(df_PPI_path_4) == "times")] <- "max_no_interactions"
colnames(df_PPI_path_4)[which(names(df_PPI_path_4) == "class")] <- "labels_driver_mutation"
colnames(df_PPI_path_4)[which(names(df_PPI_path_4) == "P1_entrez")] <- "EntrezID"




#save the output to a file
write.table(df_PPI_path_4, "mut_and_PPI_and_proc.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)




df_PPI_path_5 <- unique(merge(df_PPI_path_4, sss, by.x="EntrezID", by.y="entrezgene"))





#########################################################################################################
### 2.2.1 GENE-LEVEL FEATURE EXTRACTION (PTMs from PhosphoSitePlus)
#########################################################################################################

#make a list of proteins for which PTMs will be extracted
proteins_list <- as.data.frame(as.character(unique(df_PPI_path_5$uniprot.Swiss.Prot)))
proteins_list$`as.character(unique(df_PPI_path_5$uniprot.Swiss.Prot))` <- as.character(proteins_list$`as.character(unique(df_PPI_path_5$uniprot.Swiss.Prot))`)
#remove rows with non-uniprot IDs or multiple IDs
proteins_list_2 <- as.data.frame(proteins_list[(which(nchar(proteins_list$`as.character(unique(df_PPI_path_5$uniprot.Swiss.Prot))`) == 6)),])
#rename column name
colnames(proteins_list_2)[1] <- "proteins"

#for each 300 proteins, create a file which will be used an input for PhosphoSitePlus
#split the protein list
proteins_list_split <- split(proteins_list_2, rep(1:ceiling(nrow(proteins_list_2)/300), each=300)[1:nrow(proteins_list_2)])
#for each set write a new file
for (i in seq_along(proteins_list_split)) {
  filename <- paste(i, ".txt", sep="")
  write.table(proteins_list_split[[i]], filename, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
}



#DELETE THE ACCESSION NUMBER FROM EACH FILE BEFORE LOADING THE XLS files

#extract PTM for each set of proteins
install.packages("readxl")
library("readxl")

#put all the xls files into a list of files to be looped
files <- list.files(pattern=".xls")
out.file <- ""

#for each file extract the PTMs and bind them together 
for (i in files) {
  PTM_i <- read_excel(i)
  ptm_i <- data.frame(sapply(PTM_i,c))
  p_i <- ptm_i[ptm_i$ORGANISM=="human",]
  P_i <- p_i[,c(4,7)]
  out.file <- rbind(out.file, P_i)
}

total_PTM <- unique(as.data.frame(out.file))
#remove rows with NA
total_PTM2 <- total_PTM[complete.cases(total_PTM), ]
  

#load the key for the PTM description
key_PTM <- read_excel("key_for_PTM.xlsx")

#split the PTMs by comma
install.packages("splitstackshape")
library(splitstackshape)
total_PTM_split <- cSplit(total_PTM2, "MODIFICATION.SUMMARY", sep=",")

#reshape the structure of the data
library(reshape2)
#make each corresponding PTM to one protein
ptm_reshape_1 <- melt(total_PTM_split, id=c("ACC."))
#delete the second column which is not needed
ptm_reshape_2 <- ptm_reshape_1[,c(1,3)]
#add the description of the PTM and amino acid residue
ptm_reshape_3 <- unique(merge(ptm_reshape_2, key_PTM, by.x="value", by.y="PTM"))
#match each protein to the general description of the PTM (some proteins may undergo multiple PTMs)
ptm_reshape_4 <- unique(ptm_reshape_3[,c(2:3)])

#make PTMs as columns (dichotomous variable)
library(data.table)
#general PTM
proteins_and_general_PTM_as_cols <- dcast(ptm_reshape_4, ACC.~Description, fill=0, fun.aggregate=length)
#amino acid specific PTM
proteins_and_specific_PTM_as_cols <- dcast(ptm_reshape_3, ACC.~Description + Amino_acid_by_code, fill=0, fun.aggregate=length)




#copy the swissprot column for merging the data 
df_PPI_path_5$SWISSPROT <- df_PPI_path_5$uniprot.Swiss.Prot

#merge the general PTMs to the final dataframe & remove column used for merging
df_PPI_pathways_PTM_gen <- merge(df_PPI_path_5, proteins_and_general_PTM_as_cols, by.x="SWISSPROT", by.y="ACC.")
df_PPI_pathways_PTM_gen$SWISSPROT <- NULL
#write an output file for general PTM
write.table(df_PPI_pathways_PTM_gen, "mutations_PPI_proc_PTM_gen.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)


#nerge the specific PTMs to the final dataframe & remove column used for merging
df_PPI_pathways_PTM_spec <- merge(df_PPI_path_5, proteins_and_specific_PTM_as_cols, by.x="SWISSPROT", by.y="ACC.")
df_PPI_pathways_PTM_spec$SWISSPROT <- NULL
#write an output file for specific PTM
write.table(df_PPI_pathways_PTM_spec, "mutations_PPI_proc_PTM_spec.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)



#create a duplicate dataframe for merging all features with the initial mutation dataframe
df_novel_features <- df_PPI_pathways_PTM_gen

#create a new column for merging VEP features
df_novel_features$VEP_merge_col <- paste(df_novel_features$Chromosome, sep="_", df_novel_features$Start_Position, df_novel_features$Reference_Allele)
df_novel_features$VEP_merge_col2 <- paste(df_novel_features$VEP_merge_col, sep="/", df_novel_features$Tumor_Seq_Allele2)

#extract only needed rows
mutations_and_novel_features <- df_novel_features[,c(1,2,4,11,38:40,46,47:55,56:68,70,75)]


df_novel_features_ratio$genes <- NULL
df_novel_features_ratio[,c(2:10, 35)] <- NULL
df_novel_features_ratio <- merge(df_novel_features_ratio, df_only_missense_and_PPI_entrez, by.x="EntrezID", by.y="P1_entrez")




###############################################################################################################################
### GENE-LEVEL FEATURES (ratio-metric features from 20/20+)
###############################################################################################################################

#load the dataset
feat_2020 <- read.csv("GENE_level_ratio_metric_2020plus.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)

#extract needed columns
ratio_metric_features <- feat_2020[,c(1:16,21:29)]

#create a new column for merging purposes
df_novel_features$genes <- df_novel_features$Hugo_Symbol

#merge the data
df_novel_features_ratio <- merge(df_novel_features, feat_2020, by.x="Hugo_Symbol", by.y="gene")
write.table(df_novel_features_ratio, "gene_based_features.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)




####################################################################################################################################
### MUTATION-LEVEL FEATURES (using VEP) & preparation for merging with the gene-level features
######################################################################################################################

#VEP INPUT
VEP_input <- df_novel_features_ratio[,c(4,5,6,8,9,7)]
VEP_input$allele <- paste(VEP_input$Reference_Allele, sep="/", VEP_input$Tumor_Seq_Allele2)
VEP_input2 <- VEP_input[,c(1:3,6:7)]
VEP_input2 <- VEP_input2[,c(1:3,5,4)]
#sort the dataframe
library(dplyr)
VEP_input2 <- VEP_input2 %>% arrange(Chromosome, Start_Position, End_Position)
#write input to a file
write.table(VEP_input2, "MUTATION_level_TCGA_VEP_input.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

#VEP OUTPUT
VEP_output <- read.csv("MUTATION_level_TCGA_VEP_output.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)

#remove unwanted columns
VEP_output2 <- VEP_output[,c(1,29,30,40:76)]

#extract just what is between paratheses for PolyPhen
VEP_output2$PolyPhen <- gsub(".*\\((.*)\\).*", "\\1", VEP_output2$PolyPhen)

#transform columns as numeric
VEP_output2[,c(2:40)] <- sapply(VEP_output2[,c(2:40)], as.numeric)

#remove column 16 which has only NA
VEP_output2$LINSIGHT_rankscore <- NULL

#for each VEP mutation, keep only max score columns wise for each unique mutation
library(dplyr)
VEP_output3 <- 
    VEP_output2 %>% 
    group_by(X.Uploaded_variation) %>% #group by the VEP identifier column
    summarise_each(funs(max(., na.rm=TRUE))) #summarise each group so that it contains the max value col-wise for that group, also skipping any NA value

#due to some columns having only NA, the max function will return -Inf (which is the max from NA values), therefore replace -Inf with NA
VEP_output3[VEP_output3 == -Inf] <- NA

#compute the average from all columns of ranks
VEP_output3$average_rank <- rowMeans(VEP_output3[,c(5:39)], na.rm=TRUE)

#extract columns only with SIFT, PolyPhen, Condel and average_rank
VEP_output4 <- VEP_output3[,c(1:4, 40)]





######################################################################################################################
### MERGE GENE-BASED FEATURES, MUTATION-BASED FEATURES & LABELS and write the output to a file
######################################################################################################################

#merge gene-based features + mutation-based features
gene_based_features <- unique(df_novel_features_ratio)
  
mutation_based_features <- VEP_output4

all_features <- merge(df_novel_features_ratio, VEP_output4, by.x="VEP_merge_col2", by.y="X.Uploaded_variation")

#extract genes and mutations
genes_and_mut_TCGA <- all_features[,c(1,2)]
write.table(genes_and_mut_TCGA, "genes_and_mutations_TCGA.txt", row.names=FALSE, col.names=TRUE, quote=FALSE)

summary(as.factor(all_features2$labels_driver_mutation))

#remove unwanted columns, change order of the columns
all_features2 <- all_features[,c(1, 12:21, 24:36, 38:52, 57:65, 75:78, 22)]


#write data to a file for ML script
write.table(all_features2, "data_ready_for_ML_TCGA.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)