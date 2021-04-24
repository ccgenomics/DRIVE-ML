#######################################################################################################################################
###  BENCHMARK DATA AQUISITION & PRE-PROCESSING
#######################################################################################################################################

#Load the benchmark dataset
df_BENCHMARK <- read.csv("INPUT_benchmark_data.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)

#load the dataset from transvar (containing position)
transvar_result <- read.csv("transvar_results.csv", header=TRUE, sep=",")
#extract the chromosome no
transvar_result$chromosome <- gsub(":.*$", "", transvar_result$coordinates.gDNA.cDNA.protein.)
#extract the start/end position with qdap package
library(qdap)
transvar_result$Start_Position <- genXtract(transvar_result$coordinates.gDNA.cDNA.protein., ":g.", "/c.")
#extract reference allele
library(stringi)
transvar_result$Ref_allele <- stri_sub(transvar_result$Start_Position, from=-3, to=-3)
#extract alternative allele
transvar_result$Alt_allele <- stri_sub(transvar_result$Start_Position, from=-1, to=-1)
transvar_result$Start_Position <- gsub('.{3}$', '', transvar_result$Start_Position) #delete last 3 characters from position
transvar_result$chromosome <- substring(transvar_result$chromosome, 4) #delete first 3 characters from chromosome (i.e. delete chr and start from the 4th character)

#merge the benchmark dataset with transvar
df_BENCHMARK$merge_col <- paste(":p.", df_BENCHMARK$mutation, sep="")#add :p.
df_BENCHMARK$merge_col <- paste(df_BENCHMARK$gene, df_BENCHMARK$merge_col, sep="")

df_BENCHMARK_position <- merge(df_BENCHMARK, transvar_result, by.x="merge_col", by.y="input")
#remove columns which are not needed & rename columns
df_BENCHMARK_position2 <- unique(df_BENCHMARK_position[,c(2,3,10,14,15,15,16,17,7)])
#remove rows which have deletions and other types of mutations than missense
df_BENCHMARK_position3 <- subset(df_BENCHMARK_position2, nchar(as.character(df_BENCHMARK_position2$chromosome)) <= 2)
df_BENCHMARK_position4 <- subset(df_BENCHMARK_position3, nchar(as.character(df_BENCHMARK_position3$Start_Position)) < 12)
df_BENCHMARK_position5 <- subset(df_BENCHMARK_position4, nchar(as.character(df_BENCHMARK_position4$Start_Position.1)) < 12)
#change column names
colnames(df_BENCHMARK_position5) <- c("Hugo_Symbol", "Mutation", "Strand", "Chromosome", "Start_Position", 
                                      "End_Position", "Reference_Allele", "Alternative_Allele", "Class_label_MutaGene")


#extract Hugo gene symbols for conversion to UniProt_ID and Entrez_ID
write.table(unique(df_BENCHMARK_position2$gene.x), "genes_hugo_BENCHMARK.txt", sep="\t", row.names = FALSE, col.names = TRUE, quote=FALSE)
#map hugo to uniprot and entrez
mapping_hugo_uniprot <- read.csv("mapping_hugo_uniprot_BENCHMARK.tab", sep="\t", header=TRUE)
mapping_hugo_uniprot2 <- (mapping_hugo_uniprot[-c(13,35,54),c(1,2)])
mapping_uniprot_entrez <-  read.csv("mapping_table_uniprot_entrez.txt", header = TRUE, sep='\t')
mapping_hugo_uniprot3 <- merge(mapping_hugo_uniprot2, mapping_uniprot_entrez, by.x="Entry", by.y="From", all.x=TRUE)
#change column names
colnames(mapping_hugo_uniprot3) <- c("UniProt_ID", "Hugo_Symbol", "Entrez_Gene_ID")

#merge the original dataframe adding the UniProt_IDs and Entrez_IDs (needed for feature extraction)
df_BENCHMARK_position6 <- merge(df_BENCHMARK_position5, mapping_hugo_uniprot3, all.x=TRUE)


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
df_only_missense_and_PPI_entrez <- merge(PPI_max_no_entrez, df_BENCHMARK_position6, by.x="P1_entrez", by.y="Entrez_Gene_ID")
#check if there are any duplicated rows -> NO
dim(unique(df_only_missense_and_PPI_entrez))

#write the output to a file 
write.table(df_only_missense_and_PPI_entrez, "mutations_and_PPI.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote = FALSE)






##########################################################################################################
###  GENE-LEVEL FEATURE EXTRACTION (enriched pathway for groups of genes which interact together)
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
split_entrez_genes <- cSplit(genes_as_vector, "int_ent", sep=",")



#install package msingdbr & clusterProfiler for molecular signatures & dependencies
install.packages("msigdbr")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")

library(msigdbr)
library(clusterProfiler)

#extract molecular signatures for human by category H (hallmarks)
human_df <- msigdbr(species="Homo sapiens", category="H")
#check what types of hallmarks are present
types_hallmarks <- as.data.frame(table(human_df$gs_name))

#create a dataframe with molecular signatures for each gene
T_2_G <- human_df[,c(3,6)]




#for each row of the split_entrez_gene
for (i in 1:dim(split_entrez_genes)[1]) {
  #vectorise each row with entrez gene names so that it can inputed into ReactomePA
  vec_i <- dput(as.character(na.omit(as.vector(unlist(split_entrez_genes[i,])))))
  #if the no of genes are more than 0
  if (length(vec_i)>0) {
    #compute the enrichment pathway analysis
    p_i <- enricher(vec_i, pvalueCutoff = 0.05, pAdjustMethod = "BH",
                    minGSSize = 20, maxGSSize = 500,
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
colnames(genes_as_vector_4_categ)[which(names(genes_as_vector_4_categ) == "interact_genes_entrez")] <- "interact_ent"

genes_as_vector_5_categ <- as.data.frame(genes_as_vector_4_categ)
colnames(genes_as_vector_5_categ)[which(names(genes_as_vector_5_categ) == "interact_genes_entrez")] <- "interact_genes_Entrez"
names(genes_as_vector_5_categ)[1] <- "interact_genes_Entrez"

#merge the full datasets with the mutations
df_PPI_path <- merge(df_only_missense_and_PPI_entrez, genes_as_vector_5_categ, 
                     by.x="interact_genes_entrez", 
                     by.y="interact_genes_Entrez")





df_PPI_path_2 <-  df_PPI_path


#make pathways as columns (dichotomous variable)
genes_and_process_unique <- unique(df_PPI_path_2[,c(1,15)])
library(data.table)
genes_and_process_as_cols <- dcast(genes_and_process_unique, interact_genes_entrez~Process_category, fill=0, fun.aggregate=length)

#create a new dataframe with pathways_as_cols
df_PPI_path_3 <- merge(df_PPI_path_2, genes_and_process_as_cols, by.x="interact_genes_entrez", by.y="interact_genes_entrez")

df_PPI_path_4 <- df_PPI_path_3



#save the output to a file
write.table(df_PPI_path_4, "mut_and_PPI_and_proc.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)






#########################################################################################################
###  GENE-LEVEL FEATURE EXTRACTION (PTMs from PhosphoSitePlus)
#########################################################################################################

#make a list of proteins for which PTMs will be extracted
proteins_list <- as.data.frame(unique(df_PPI_path_4$UniProt_ID))

#write the list to a file with max 300 proteins/file/search
write.table(proteins_list[1:57,], "unique_proteins_BENCHMARK.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

#DELETE THE ACCESSION NUMBER FROM EACH FILE BEFORE LOADING THE XLS files

#extract PTM for each set of proteins
install.packages("readxl")
library("readxl")       
PTM_1 <- read_excel("protein_BENCHMARK_PTMs.xls")
ptm_1 <- data.frame(sapply(PTM_1,c))
p_1 <- ptm_1[ptm_1$ORGANISM == "human",]
P_1 <- p_1[,c(4,7)]


#load the key for the PTM description
key_PTM <- read_excel("key_for_PTM.xlsx")

#split the PTMs by comma
install.packages("splitstackshape")
library(splitstackshape)
total_PTM_split <- cSplit(P_1, "MODIFICATION.SUMMARY", sep=",")

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
df_PPI_path_4$UniProt_ID2 <- df_PPI_path_4$UniProt_ID

#merge the general PTMs to the final dataframe & remove column used for merging
df_PPI_pathways_PTM_gen <- merge(df_PPI_path_4, proteins_and_general_PTM_as_cols, by.x="UniProt_ID2", by.y="ACC.")
df_PPI_pathways_PTM_gen$UniProt_ID2 <- NULL
#write an output file for general PTM
write.table(df_PPI_pathways_PTM_gen, "mutations_PPI_proc_PTM_gen.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)


#nerge the specific PTMs to the final dataframe & remove column used for merging
df_PPI_pathways_PTM_spec <- merge(df_PPI_path_4, proteins_and_specific_PTM_as_cols, by.x="UniProt_ID2", by.y="ACC.")
df_PPI_pathways_PTM_spec$SWISSPROT_2 <- NULL
#write an output file for specific PTM
write.table(df_PPI_pathways_PTM, "mutations_PPI_proc_PTM_spec.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)



#create a duplicate dataframe for merging all features with the initial mutation dataframe
df_novel_features <- df_PPI_pathways_PTM_gen

#create a new column for merging VEP features
df_novel_features$VEP_merge_col <- paste(df_novel_features$Chromosome, sep="_", df_novel_features$Start_Position, df_novel_features$Reference_Allele)
df_novel_features$VEP_merge_col2 <- paste(df_novel_features$VEP_merge_col, sep="/", df_novel_features$Alternative_Allele)

#extract only needed rows
mutations_and_novel_features <- df_novel_features[,c(1,2,4,38:40,46,47:55,56:68,70,75)]






###############################################################################################################################
### GENE-LEVEL FEATURES (ratio-metric features from 20/20+)
###############################################################################################################################

#load the dataset
feat_2020 <- read.csv("output_2020plus_pan_cancer.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)

#extract needed columns
ratio_metric_features <- feat_2020[,c(1:16,21:29)]

#create a new column for merging purposes
df_novel_features$genes <- df_novel_features$Hugo_Symbol

#merge the data
df_novel_features_ratio <- merge(df_novel_features, feat_2020, by.x="genes", by.y="gene")
write.table(df_novel_features_ratio, "gene_based_features.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)





####################################################################################################################################
### MUTATION-LEVEL FEATURES (using VEP) & preparation for merging with the gene-level features
######################################################################################################################

#VEP INPUT
VEP_input <- df_PPI_pathways_PTM_gen[,c(8:12,7)]
VEP_input$allele <- paste(VEP_input$Reference_Allele, sep="/", VEP_input$Alternative_Allele)
VEP_input2 <- VEP_input[,c(1,2,3,7,6)]

#sort the dataframe
library(dplyr)
VEP_input2 <- VEP_input2 %>% arrange(Chromosome, Start_Position, End_Position)
#write output to a file
write.table(unique(VEP_input2), "MUTATION_level_BENCHMARK_VEP_input.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

#VEP OUTPUT
VEP_output <- read.csv("MUTATION_level_BENCHMARK_VEP_output.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
#remove unwanted columns
VEP_output2 <- VEP_output[,c(1,29,30,77, 41:76)]

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
VEP_output3$average_rank <- rowMeans(VEP_output3[,c(4, 6:39)], na.rm=TRUE)

#extract columns only with SIFT, PolyPhen, Condel and average_rank
VEP_output4 <- VEP_output3[,c(1,2,3,5, 40)]
VEP_output5 <- VEP_output4[complete.cases(VEP_output4),]





######################################################################################################################
### MERGE GENE-BASED FEATURES, MUTATION-BASED FEATURES & LABELS and write the output to a file
######################################################################################################################

#merge gene-based features + mutation-based features
  
all_features <- merge(df_novel_features_ratio, VEP_output5, by.x="VEP_merge_col2", by.y="X.Uploaded_variation", all.x=TRUE)

all_features <- all_features[complete.cases(all_features),]





#####################################################################################################################
### EXTRACT INPUT & OUTPUT FOR CHASMplus
#####################################################################################################################

#extract file for CHASMplus input
CHASMplus_input <- all_features[,c(10,11,9,13,14)]
#modify the Chromosome column to contain "chr"
CHASMplus_input$Chromosome <- paste("chr", CHASMplus_input$Chromosome, sep="")
write.table(CHASMplus_input, "CHASMplus_input.txt", quote=FALSE, row.names = FALSE, col.names=FALSE)

#read the output file from CHASMplus
CHASMplus_output <- read.csv("CHASMplus_output.tsv", sep="\t")
#change name of the column with CHASMplus scores (i.e. probabilities)
colnames(CHASMplus_output)[13] <- "CHASMplus_score"
CHASMplus_output$Chrom <- substring(CHASMplus_output$Chrom, 4) #delete first 3 characters from chromosome (i.e. delete chr and start from the 4th character)


#add CHASMplus scores
all_features2 <- merge(all_features, CHASMplus_output, 
                       by.x=c("Chromosome", "Start_Position", "Reference_Allele", "Alternative_Allele"), 
                       by.y=c("Chrom", "Position", "Ref_Base", "Alt_Base"))

#order features for machine learning and keep only needed columns
all_features3 <- all_features2[,c(5, 10, 22:43, 45:59, 64:72, 82:85, 15, 94)]
#keep only columns with complete data (i.e. no missing data)
all_features4 <- all_features3[complete.cases(all_features3),]




##########################################################################################################################
### DATA PREPARATION FOR ML
##########################################################################################################################

#change column names to match the other datasets of GENIE and TCGA
colnames(all_features4)[2] <- "max_no_interactions"
colnames(all_features4)[3] <- "cellular.component"
colnames(all_features4)[5] <- "DNA.damage"
colnames(all_features4)[11] <- "NA_hallmarks"
colnames(all_features4)[13] <- "caspase.cleavage"
colnames(all_features4)[14] <- "di.methylation"
colnames(all_features4)[16] <- "mono.methylation"
colnames(all_features4)[17] <- "N.Glycosylation"
colnames(all_features4)[18] <- "O.GalNAc"
colnames(all_features4)[19] <- "O.GlcNAc"
colnames(all_features4)[23] <- "tri.methylation"

#replace 1 and 0 with yes and no for hallmarks and PTMs features
all_features4[,c(3:24)] <- factor(ifelse(all_features4[,c(3:24)] == 1, "yes", "no"))

#replace non-neutral with driver and neutral with passenger
all_features4$Class_label_MutaGene <- ifelse(all_features4$Class_label_MutaGene == "Non-neutral", "driver", "passenger")

#write data to a file for ML script
write.table(all_features4, "data_ready_for_ML_BENCHMARK.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
