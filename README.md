# Asssessment-of-RNA-seq-performance
R or Python code for data analysis of the manuscript entitled "Toward Best Practice in Identifying Subtle Differential Expression with RNA-seq: A Real-World Multi-Center Benchmarking Study Using Quartet and MAQC Reference Materials"

Here we provide scripts for assessing RNA-seq performance, inclduing data quality, accuracy of absolute and relative gene expression, and accuracy of differentially expressed genes.
Three types of reference datasets are used, inclduing Quartet reference datasets, Quartet TaqMan datasets, and MAQC TaqMan datsets. We have provided these files named "quartet_reference_datasets.txt", "quartet_taqMan_datasets.txt", and "maqc_taqMan_datasets.txt".

1. SNR for expression data quality assessment
   
   A expression matrix for Quartet and mixed samples is required. The SNR18 (SNR based on all 18 samples) and SNR17 (SNR based on any 17 out of 18 samples) could be calculated.

2. Accuracy of absolute gene expression

   A expression matrix for all 24 samples is required.

   absolute_ERCC.py will conduct a correlation analysis between absolute expression data and concentrations of 92 ERCC genes;

   Absolute_TaqMan.py will conduct a correlation analysis between absolute expression data and TaqMan gene abundance for Quartet and MAQC samples.

3. Accuracy of relative gene expression

   A expression matrix for all 24 samples is required.

   Relative_reference.py will conduct a correlation analysis between relative expression data and three types of reference datastes;

   ERCC_ratio.py will analyze the consitency between observed relative expression and expected ratios for 92 ERCC genes. This is based on the expected ratios of four subgroups of ERCC genes, 4:1, 3:2, 2:3, 1:1, respectively;

   Mixing_ratio.py will analyze the consitency between observed relative expression and expected relative expression in T1/D6 and T2/D6 for all detected genes. This is based on the mixing ratios in T1 and T2 samples, which are mixed from M8 and D6 with the ratios of 3:1 and 1:3;

4. Accuracy of DEG identification

   Read counts for 24 samples are required.

   DEG_analysis.R includes four types of differential expression analysis tools, including edgeR, limma, DESeq2, DEGseq, and EBSeq. Users can set paramater to select the tools. For each tool, a series of threshold for filtering low-expression genes is applied, then the accuracy of DEG (MCC) is calculated based on three types of reference datsets under different thresholds. The optimal MCC is determined and ouputed.

