Project: Skin UVB SKH1 mouse model treated with UA/SFN
Scientist: Yuqing (Anne) Yang
Data Analysis: Ran Yin, Renyi Wu, Davit Sargsyan

---

## Daily Logs
### 04/05/2018
* Changed hitmaps to donut plots; added more comparisons
* ToDo: tumor vs. normal tissue at week 25

### 03/29/2018
* Pairwise comaprisons   
* Added MA plots, heatmaps, and counts for Venn diagrams

### 03/22/2018
* Added DSS script

### 03/15/2015
* Analysis plan (Ran)    

1. DNA methylation of SFN, UVB, Control
2. DNA methylation of UA, UVB, Control
3. RNA seq of SFN, UVB, Control
4. RNA seq of UA, UVB, Control
5. DNA methyl seq with RNA seq UVB induced carcinogenesis.
6. Aging and methylation

Figures and Tables in Paper 1 (SFN), 2 (UA) DNA papers.
1. Table 1: Up and down regulated top 10 genes in all time points and tumors.
2. Figure 1: Top 20 pathways with z-score indication for all time points and tumor
3. Figure 2: Heatmap and cluster for all samples (4 samples from Control, 4 samples from UVB,
4 samples from SFN with UVB).
4. Figure 3: PCA or column chart of all samples (samples are same in 3)
5. Figure 4: MA ploting (DEseq2)
6. Figure 5: Gene heatmap and Venn diagram.
6. Figure 6: Time vs CpG region with ??? (still working on this figure)

Figures and Tables in Paper 3 (SFN), 4 (UA) RNA papers.
1. Table 1: Up and down regulated top 10 genes in all time points and tumors.
2. Figure 1: Top 20 pathways with z-score indication for all time points and tumor
3. Figure 2: Heatmap and cluster for all samples (4 samples from Control, 4 samples from UVB,
4 samples from SFN with UVB).
4. Figure 3: PCA or column chart of all samples (samples are same in 3)
5. Figure 4: MA ploting (DEseq2)
6. Figure 5: Genes heatmap and Venn diagram.
Figures and Tables in paper

* Analysis plan (Davit)    

1. For Skin UVB:
   a. We are splitting data into two sets: Ctrl, UVB and UA; Ctrl, UVB and SFN, and rerunning DESeq2 analysis of RNA-seq
   b. Same for Methyl-seq data and analyzing using John's DMRfinder to  cluster CpG followed by DSS two-factor analysis
   c. DNA methylation correlation with RNA expression
   d. Separately, examine aging effect in control samples only over time 

### 03/17/2018
* Added Methyl-seq analysis with package DSS (*skin_uvb_methylseq_dss_v1.R*)

### 03/10/2018
* Added RNA-seq data analysis with DESeq2

### 01/06/2017
* Added analysis for weeks 20 and 22

### 01/04/2018
* Added script and results for tumor data analysis

### 12/15/2017
* Added Venn diagrams within each timepoint and within each comparison    

From Dr. Kong:    
I discussed the RNA-seq data for UVB and UA and SFN study yesterday with Anne and Ran, so, here is what I like to see in the ms:    
1. Control no UVB   
2. Control with UVB    
3. UVB + UA    
4. UVB + SFN    
So, split the study into 2 ms - UVB+UA and UVB+SFN    
For UVB+UA - compare 4 sets separately and NOT at the same time:  
1. 2 wk vs 15 wks - initiation/promotion stage    
2. 15 wks vs 25 wks - progression stage    
3. 2 wk vs 25 wks - all the way to end stage - very high risk skin BUT no tumor    
4. 2 wks vs tumor - all the way to tumor   
Do all at the same time - BUT not sure what we get out of this?    
Repeat for UVB+SFN.    
 
### 12/14/2017
* Added Venn diagrams    
    
From Anne:   
For the venn2, could you specify the comparison between two groups, like (1) UVB vs Control (No UVB), (2) UVA+UA vs UVB, and (3) UVB+SFN vs UVB? So we need these three venn graphics similar to Venn2 with still four time-points each. Also, it's will be great if you can send me an excel file of shared genes with changes (value 1, value 2, log2 fold change, q value) at different time points.    

### 12/12/2017
* Repository created