![CCG](/other/ccg.png)



# DRIVE-ML


**A Feature-Based Machine Learning Model for Pan-Cancer Assessment of Somatic Missense Mutations**


**Release Version:** 1.0

**Author(s):** Ionut Dragomir, Adnan Akbar, John W Cassidy, Nirmesh Patel, Harry 
W Clifford, Gianmarco Contino


## Description

Project completed by Ionut Dragomir during his internship at Cambridge Cancer Genomics, Cambridge, United Kingdom, 2020-2021.

Simple Summary: Genes dictate the grounds of life by comprising molecular bases which encode proteins. A mutation represents a gene modification that may influence the protein function. Cancer occurs when the mutation triggers uncontrolled cellular growth. Judging by the cancer expansion, mutations labelled as drivers confer a growth advantage, while passengers do not contribute to this augmentation. The aim of this study is methodological, which assesses the usefulness of a classification method for distinguishing between driver and passenger mutations. Based on 51 molecular characteristics of mutations and genes, including 3 novel features, multiple machine learning algorithms were used to determine whether these characteristics biologically represent the driver mutations and how they impact the classification procedure. To test the ability of the present methodology, the same steps were applied to an independent dataset. The results showed that both gene and mutation level characteristics are representative of the driver mutations, and the proposed approach achieved more than 80% accuracy in finding the true type of mutation. The evidence suggests that machine learning methods can be used to gain knowledge from mutational data seeking to deliver more targeted cancer treatment.

Abstract: Sporadic cancer occurs due to gene damage which is translated at a molecular level as a somatic mutation. Out of all small-scale somatic aberrations, 95% are base substitutions, with 90% being missense mutations. While multiple studies focused on the importance of this mutation type, a solid characterisation based on cancer growth has not been thoroughly adopted. To address these drawbacks, this study aims to develop an improved computational method for driver identification, validation and evaluation (DRIVE), which is compared to other methods for assessing its performance. DRIVE aims at distinguishing between driver and passenger mutations using a feature-based learning approach comprising two levels of biological classification for a pan-cancer assessment of somatic mutations. Gene-level features include the maximum number of protein–protein interactions (PPIs), the biological process and the type of post-translational modifications (PTMs) while mutation-level features are based on pathogenicity scores. Multiple supervised classification algorithms were trained on Genomics Evidence Neoplasia Information Exchange (GENIE) project data and then tested on an independent dataset from The Cancer Genome Atlas (TCGA) study. finally, the most powerful classi- fier using DRIVE was evaluated on a benchmark dataset, which showed a better overall performance compared to other state-of-the-art methodologies, however, considerable care must be taken due to the reduced size of the dataset. DRIVE outlines the outstanding potential that multiple levels of a feature-based learning model will play in the future of oncology-based precision medicine.


## Data pre-processing 

* uses ***stringr***, ***dplyr***, ***splitstackshape***, ***data.table***, ***readxl*** and ***reshape2*** R packages. 
1. **FILTERING** (keeping only the missense mutations)
2. **ADDRESS THE PROBLEM OF MISSING DATA**
3. **EXPLORATORY DATA ANALYSIS**


## Feature extraction 

* uses ***mygene***, ***ReactomePA***, ***msigdbr***, ***org.Hs.eg.db*** and ***clusterProfiler*** R packages. 
1. **GENE-LEVEL FEATURES EXTRACTION** 
   * NOVEL STRUCTURAL FEATURE (max no of PPIs from Interactome INSIDER)
   * NOVEL STRUCTURAL FEATURE (enriched pathways from MSigDB)
   * NOVEL STRUCTURAL FEATURE (PTMs from PhosphoSitePlus)
   * EXISTENT RATIOMETRIC FEATURES (from 20/20+)
3. **MUTATION-LEVEL FEATURES EXTRACTION**
   * SIFT, PolyPhen, Condel and average score of multiple rank scores (from VEP)
4. **TRUE CLASS LABELS COMPILATION** (based on statistical models)


## Machine learning 

* uses ***caret*** + dependencies, ***varImp***, ***pROC***, ***ggplot2*** and ***ggpubr*** R packages.
1. **DATA PRE-PROCESSING FOR A ML FORMAT** 
   * pre-process training data (GENIE) - one-hot encoding necessary for *mlp* models
   * pre-process independent prediction dataset (TCGA)
   * pre-process benchmark data (MutaGene) - one-hot encoding necessary for *mlp* models
3. **IMBALANCED CLASSIFICATION PROBLEM** - addressed by undersampling using ***UBL*** R package
4. **K-FOLD CROSS-VALIDATION** (random forest, decision tree, EGB, logistic regression, SVM, KNN, MLP)
5. **PREDICTION ON INDEPENDENT DATA** (random forest, decision tree, EGB, logistic regression, SVM, KNN, MLP)
6. **BENCHMARKING** (random forest vs. state-of-the-art-model)
