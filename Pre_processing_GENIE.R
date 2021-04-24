####################################################################################################################
### 1. MUTATION DATA COLLECTION & PRE-PROCESSING
####################################################################################################################

#Load the dataset
df <- read.csv("INPUT_GENIE_data.txt", header = TRUE, sep='\t', stringsAsFactors=FALSE)

#remove duplicated rows
df <- unique(df)

#Keep only missense mutations
df_only_missense <- unique(df[grep("Missense_Mutation", df$Variant_Classification), ])

#Remove p.
df_only_missense$HGVSp_Short <- gsub("^.{0,2}", "", df_only_missense$HGVSp_Short)

#Duplicate first column for 3 times
df_only_missense$Position <- df_only_missense$HGVSp_Short
df_only_missense$First_amino_acid <- df_only_missense$HGVSp_Short
df_only_missense$Second_amino_acid <- df_only_missense$HGVSp_Short

#Keep only the first letter for the First amino acid
df_only_missense$First_amino_acid <- strtrim(df_only_missense$First_amino_acid, 1)

#Keep only the number for Position
df_only_missense$Position <- gsub('^.|.$', '', df_only_missense$Position)


#Keep only the last letter for the Second amino acid
install.packages("stringr")
library(stringr)
df_only_missense$Second_amino_acid <- str_sub(df_only_missense$Second_amino_acid, -1)

#Remove rows with more than one SwissProt identifier (isoforms)
df_only_missense <- df_only_missense[(which(nchar(as.character(df_only_missense$SWISSPROT)) == 6)),]

#Keep only rows with numeric mutation position
df_only_missense <- subset(df_only_missense, grepl('^\\d+$', df_only_missense$Position))

#check how many types of variants and their numbers we have
types_variant <- as.data.frame(table(df_only_missense$Variant_Type))

#check how many affected proteins/genes do we have
length(unique(df_only_missense$SWISSPROT))

#remove any other missense mutations than SNPs (e.g. DNP, DEL, ONP, TNP)
df_only_missense <- df_only_missense[which(df_only_missense$Variant_Type != "DNP"),]
df_only_missense <- df_only_missense[which(df_only_missense$Variant_Type != "DEL"),]
df_only_missense <- df_only_missense[which(df_only_missense$Variant_Type != "ONP"),]
df_only_missense <- df_only_missense[which(df_only_missense$Variant_Type != "TNP"),]

#write the output to a file
write.table(df_only_missense, "GENIE_only_missense.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote = FALSE)





#########################################################################################################
###  2.1.1.1 GENE-LEVEL FEATURE EXTRACTION (max no of PPI for HQ-interactions from Interactome INSIDER)
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
df_only_missense_and_PPI_entrez <- merge(PPI_max_no_entrez, df_only_missense, by.x="P1_entrez", by.y="Entrez_Gene_Id")
#check if there are any duplicated rows -> NO
dim(unique(df_only_missense_and_PPI_entrez))

#write the output to a file 
write.table(df_only_missense_and_PPI_entrez, "mutations_and_PPI.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote = FALSE)





##########################################################################################################
### 2.1.1.2 GENE-LEVEL FEATURE EXTRACTION (enriched pathway for groups of genes which interact together)
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


#reorder the columns of the dataframe to match the initial GENIE data structure
colnames(df_only_missense)
df_PPI_path$Entrez_Gene_Id <- df_PPI_path$P1_entrez
colnames(df_PPI_path)
df_PPI_path_2 <- df_PPI_path[,c("Hugo_Symbol",                "Entrez_Gene_Id",                  "Center",                       
                                "NCBI_Build",                    "Chromosome",                    "Start_Position"   ,            
                                "End_Position",                  "Strand"  ,                      "Variant_Classification" ,      
                                "Variant_Type" ,                 "Reference_Allele" ,             "Tumor_Seq_Allele1"  ,          
                                "Tumor_Seq_Allele2" ,            "dbSNP_RS"    ,                  "dbSNP_Val_Status" ,            
                                "Tumor_Sample_Barcode"  ,        "Matched_Norm_Sample_Barcode" ,  "Match_Norm_Seq_Allele1",       
                                "Match_Norm_Seq_Allele2",        "Tumor_Validation_Allele1"  ,    "Tumor_Validation_Allele2"  ,   
                                "Match_Norm_Validation_Allele1", "Match_Norm_Validation_Allele2", "Verification_Status"  ,        
                                "Validation_Status",             "Mutation_Status"     ,          "Sequencing_Phase"  ,           
                                "Sequence_Source",               "Validation_Method" ,            "Score"  ,                      
                                "BAM_File",                      "Sequencer"  ,                   "HGVSp_Short" ,                 
                                "t_ref_count",                   "t_alt_count" ,                  "n_ref_count" ,                 
                                "n_alt_count",                   "Protein_position" ,             "Codons"     ,                  
                                "SWISSPROT",                     "RefSeq"  ,                      "t_depth"   ,                   
                                "n_depth",                       "FILTER",                        "mutationInCis_Flag" ,          
                                "Position",                      "First_amino_acid" ,             "Second_amino_acid",
                                "P1_entrez", "P2_entrez", 
                                "times", "interact_genes_entrez", "Process_category", 
                                "Hallmark_name", "Description", "Number_of_founder_sets", 
                                "Number_of_genes")]



#write the output with max no of PPI and hallmarks+ pathways to a file
write.table(df_PPI_path_2, "mutations_and_PPI_and_proc+hallm.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote = FALSE)

df_PPI_path_2 <-  read.csv("mutations_and_PPI_and_proc+hallm.txt", header = TRUE, sep='\t', stringsAsFactors=FALSE)


#make pathways as columns (dichotomous variable)
genes_and_process_unique <- unique(df_PPI_path_2[,52:53])
library(data.table)
genes_and_process_as_cols <- dcast(genes_and_process_unique, interact_genes_entrez~Process_category, fill=0, fun.aggregate=length)

#create a new dataframe with pathways_as_cols
df_PPI_path_3 <- merge(df_PPI_path_2, genes_and_process_as_cols, by.x="interact_genes_entrez", by.y="interact_genes_entrez")

#keep only needed columns
df_PPI_path_4 <- df_PPI_path_3[,c(2:46,52,58:66)]
colnames(df_PPI_path_4)[which(names(df_PPI_path_4) == "NA")] <- "NA_hallmarks"
colnames(df_PPI_path_4)[which(names(df_PPI_path_4) == "times")] <- "max_no_interactions"



#save the output to a file
write.table(df_PPI_path_4, "mut_and_PPI_and_proc.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)





#########################################################################################################
### 2.1.1.3 GENE-LEVEL FEATURE EXTRACTION (PTMs from PhosphoSitePlus)
#########################################################################################################

#make a list of proteins for which PTMs will be extracted
proteins_list <- as.data.frame(unique(df_PPI_path_4$SWISSPROT))

#write the list to a file with max 300 proteins/file/search
write.table(proteins_list[1:300,], "unique_proteins_1_300.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(proteins_list[301:600,], "unique_proteins_301_600.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(proteins_list[601:900,], "unique_proteins_601_900.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(proteins_list[901:1200,], "unique_proteins_901_1200.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(proteins_list[1201:1396,], "unique_proteins_1201_1396.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

#DELETE THE ACCESSION NUMBER FROM EACH FILE BEFORE LOADING THE XLS files

#extract PTM for each set of proteins
install.packages("readxl")
library("readxl")       
PTM_1 <- read_excel("1-300.xls")
ptm_1 <- data.frame(sapply(PTM_1,c))
p_1 <- ptm_1[ptm_1$ORGANISM == "human",]
P_1 <- p_1[,c(4,7)]


PTM_2 <- read_excel("301-600.xls")
ptm_2 <- data.frame(sapply(PTM_2,c))
p_2 <- ptm_2[ptm_2$ORGANISM == "human",]
P_2 <- p_2[,c(4,7)]


PTM_3 <- read_excel("601-900.xls")
ptm_3 <- data.frame(sapply(PTM_3,c))
p_3 <- ptm_3[ptm_3$ORGANISM == "human",]
P_3 <- p_3[,c(4,7)]


PTM_4 <- read_excel("901-1200.xls")
ptm_4 <- data.frame(sapply(PTM_4,c))
p_4 <- ptm_4[ptm_4$ORGANISM == "human",]
P_4 <- p_4[,c(4,7)]


PTM_5 <- read_excel("1201-1396.xls")
ptm_5 <- data.frame(sapply(PTM_5,c))
p_5 <- ptm_5[ptm_5$ORGANISM == "human",]
P_5 <- p_5[,c(4,7)]


#bind all 4 files together
total_PTM <- rbind(P_1, P_2, P_3, P_4, P_5)

#load the key for the PTM description
key_PTM <- read_excel("key_for_PTM.xlsx")

#split the PTMs by comma
install.packages("splitstackshape")
library(splitstackshape)
total_PTM_split <- cSplit(total_PTM, "MODIFICATION.SUMMARY", sep=",")

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
df_PPI_path_4$SWISSPROT_2 <- df_PPI_path_4$SWISSPROT

#merge the general PTMs to the final dataframe & remove column used for merging
df_PPI_pathways_PTM_gen <- merge(df_PPI_path_4, proteins_and_general_PTM_as_cols, by.x="SWISSPROT_2", by.y="ACC.")
df_PPI_pathways_PTM_gen$SWISSPROT_2 <- NULL
#write an output file for general PTM
write.table(df_PPI_pathways_PTM_gen, "mutations_PPI_proc_PTM_gen.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)


#nerge the specific PTMs to the final dataframe & remove column used for merging
df_PPI_pathways_PTM_spec <- merge(df_PPI_path_4, proteins_and_specific_PTM_as_cols, by.x="SWISSPROT_2", by.y="ACC.")
df_PPI_pathways_PTM_spec$SWISSPROT_2 <- NULL
#write an output file for specific PTM
write.table(df_PPI_pathways_PTM, "mutations_PPI_proc_PTM_spec.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)



#create a duplicate dataframe for merging all features with the initial mutation dataframe
df_novel_features <- df_PPI_pathways_PTM_gen

#create a new column for merging VEP features
df_novel_features$VEP_merge_col <- paste(df_novel_features$Chromosome, sep="_", df_novel_features$Start_Position, df_novel_features$Tumor_Seq_Allele1)
df_novel_features$VEP_merge_col2 <- paste(df_novel_features$VEP_merge_col, sep="/", df_novel_features$Tumor_Seq_Allele2)

#extract only needed columns
mutations_and_novel_features <- df_novel_features[,c(1,2,4,38:40,46,47:55,56:68,70,75)]

write.csv(df_novel_features, "GENIE_novel_feat.csv", row.names = FALSE)





###############################################################################################################################
### 2.1.2 GENE-LEVEL FEATURES (ratio-metric features from 20/20+)
###############################################################################################################################

#load the dataset
feat_2020 <- read.csv("GENE_level_ratio_metric_2020plus.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)

#extract needed columns
ratio_metric_features <- feat_2020[,c(1:16,21:29)]

#create a new column for merging purposes
df_novel_features2$genes <- df_novel_features2$Hugo_Symbol

#merge the data
df_novel_features_ratio <- merge(df_novel_features2, feat_2020, by.x="genes", by.y="gene")
write.table(df_novel_features_ratio, "gene_based_features.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)




####################################################################################################################################
### 2.2 MUTATION-LEVEL FEATURES (using VEP) & preparation for merging with the gene-level features
######################################################################################################################

#VEP INPUT
VEP_input <- df_PPI_pathways_PTM_gen[,c(5,6,7,12,13,8)]
VEP_input$allele <- paste(VEP_input$Tumor_Seq_Allele1, sep="/", VEP_input$Tumor_Seq_Allele2)
VEP_input2 <- VEP_input[,c(1:3,6:7)]
VEP_input2 <- VEP_input2[,c(1:3,5,4)]
#sort the dataframe
library(dplyr)
VEP_input2 <- VEP_input2 %>% arrange(Chromosome, Start_Position, End_Position)
#write output to a file
write.table(VEP_input2, "MUTATION_level_GENIE_VEP_input.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

#VEP OUTPUT
VEP_output <- read.csv("MUTATION_level_GENIE_VEP_output.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
#remove unwanted columns
VEP_output2 <- VEP_output[,c(1,29,30,40:80)]

#extract just what is between parentheses for PolyPhen
VEP_output2$PolyPhen <- gsub(".*\\((.*)\\).*", "\\1", VEP_output2$PolyPhen)

#transform columns as numeric
VEP_output2[,c(2:44)] <- sapply(VEP_output2[,c(2:44)], as.numeric)

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
VEP_output3$average_rank <- rowMeans(VEP_output3[,c(5:43)], na.rm=TRUE)

#extract columns only with SIFT, PolyPhen, Condel and average_rank
VEP_output4 <- VEP_output3[,c(1:4, 44)]





#####################################################################################################################
### 3. LABELS (based on statistical models according to frequency of mutations)
#####################################################################################################################
hotspots <- read.csv("hotspots.tsv", header = TRUE, sep='\t', stringsAsFactors=FALSE)

#re-write the amino acids ad positions for merging with the labels
#Duplicate first column for 3 times
df_novel_features$Position <- df_novel_features$HGVSp_Short
df_novel_features$First_amino_acid <- df_novel_features$HGVSp_Short
df_novel_features$Second_amino_acid <- df_novel_features$HGVSp_Short

#Keep only the first letter for the First amino acid
df_novel_features$First_amino_acid <- strtrim(df_novel_features$First_amino_acid, 1)

#Keep only the number for Position
df_novel_features$Position <- gsub('^.|.$', '', df_novel_features$Position)


#Kepp only the last letter for the Second amino acid
install.packages("stringr")
library(stringr)
df_novel_features$Second_amino_acid <- str_sub(df_novel_features$Second_amino_acid, -1)


#create a new column in the labels dataframe
hotspots$codon_id2 <- hotspots$codon_id
#split the gene name and mutation and re-name columns
library(stringr)
hotspots[,13:14] <- str_split_fixed(hotspots$codon_id2, "_", 2)
colnames(hotspots)[which(names(hotspots) == "codon_id2")] <- "driver_mutation"
colnames(hotspots)[which(names(hotspots) == "V13")] <- "driver_gene"


#create a new column for merging labels/'ground truth'
df_novel_features$labels_merge <- paste(df_novel_features$Hugo_Symbol, sep="_", df_novel_features$First_amino_acid)
df_novel_features$labels_merge2 <- paste(df_novel_features$labels_merge, sep="", df_novel_features$Position)


#delete unwanted columns
df_novel_features2 <- df_novel_features[,c(1:68,70,75)]

write.table(hotspots, "hotspots_2.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(all_features_and_labels, "all_features_and_labels.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

#extract needed columns
labels <- hotspots[,c(12:13)]






######################################################################################################################
### MERGE GENE-BASED FEATURES, MUTATION-BASED FEATURES & LABELS and write the output to a file
######################################################################################################################

#merge gene-based features + mutation-based features
gene_based_features <- merge(mutations_and_novel_features, ratio_metric_features, by.x="Hugo_Symbol", by.y="gene", all.x=TRUE)
  
mutation_based_features <- merge(mutations_and_novel_features, VEP_output4, by.x="VEP_merge_col2", by.y="X.Uploaded_variation", all.x=TRUE)

all_features <- merge(unique(gene_based_features), unique(mutation_based_features))

#add labels to the mutations+features
labels$driver_mutation2 <- labels$driver_mutation
all <- merge(all_features, labels, by.x="labels_merge2", by.y="driver_mutation", all.x=TRUE)

colnames(all)[which(names(all) == "labels_merge2")] <- "Mutation"
colnames(all)[which(names(all) == "driver_mutation2")] <- "labels_driver_mutation"


all$driver_gene <- NULL

#write object "all" to a file for merging with benchmark dataset
write.csv(all, "all_features_GENIE.csv", row.names = FALSE, quote=FALSE)


#extract gene names and mutations
gene_and_mut <- all2[,c(2,31)]
write.table(gene_and_mut, "genes_and_mutations_GENIE.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)


summary (as.factor(all2$labels_driver_mutation))


#write data to a file for ML script
write.table(all2, "data_ready_for_ML_GENIE.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)




##########################################################################################################################
### DATA PREPARATION FOR ML (if needed)
##########################################################################################################################

#replace all rows with a driver mutation with 1
all$labels_driver_mutation[!is.na(all$labels_driver_mutation)] <- 1
#replace all rows with NA with 0
all$labels_driver_mutation[is.na(all$labels_driver_mutation)] <- 0


#remove rows with NA
all2 <- na.omit(all)
#write the final data frame to a file
write.table(all2, "all2.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

#replace all rows with a driver gene with 1
all$driver_gene[!is.na(all$driver_gene)] <- 1
#replace all rows with NA with 0
all$driver_gene[is.na(all$driver_gene)] <- 0

colnames(all)[which(names(all) == "codon_id")] <- "labels_driver_mutation"
colnames(all)[which(names(all) == "driver_gene")] <- "labels_driver_gene"