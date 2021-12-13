# Overview
This pipeline is designed to analyze RNASeq data for differential expression across embryonic stage and sex, to estimate a gene regulatory network, and to map differential expression data onto the estimated network. The code can be tweaked to analyze more general cases of differential expression. 

# Pipeline order 
The scripts are designed to be executed in the following order: 

1. PCA_and_DiffExp.R

   Creates PCA plot (for a more customized plot of the pca results, use plotly and add traces for each group). 
   
   Runs differential expression analysis between two specified groups. 
   
   Code adapted from [https://github.com/victorya-richardson/RNASeq_pipeline](https://github.com/victorya-richardson/RNASeq_pipeline)


2. GRN-finder.py

   Runs ARACHNe-AP to estimate a gene regulatory network for the gene expression input matrix, then uses differential expression data to identify male and female subnetworks. 
   
   Analyzes log2 fold change and betweenness centrality of all nodes in both subnetworks.
   
## Before Running GRN-finder.py
1. Install ARACHNe-AP and cytoscape. ARACHNe-AP installation instructions are at https://github.com/califano-lab/ARACNe-AP. Cytoscape can be downloaded from https://cytoscape.org/.
2. Update all filepaths in GRN-finder.py and in run-ARACHNe-AP.sh to match program and file locations in your system.
3. Make sure the cytoscape app is open in the background.

# Individual Files
**PCA_and_DiffExp.R**    
   Creates PCA plot and runs differential expression analysis between two specified groups. Code adapted from [https://github.com/victorya-richardson/RNASeq_pipeline](https://github.com/victorya-richardson/RNASeq_pipeline)    

**GRN-finder.py**    
Runs ARACHNe-AP to estimate a gene regulatory network for the gene expression input matrix, then uses differential expression data to identify male and female subnetworks. Finally, analyzes log2 fold change and betweenness centrality of all nodes in both subnetworks.

**utils.py**    
Contains helper classes for GRN-finder.py. The DEgenes class is designed to store and manipulate differential expression data previously obtained from the DEseq2 program as part of my research. The ArachneData class is unused.

**run-ARACHNe-AP.sh**    
Estimates the gene regulatory network from a gene expression matrix using ARACHNe-AP. ARACHNe-AP is provided with a list of transcription factors to consider as hub genes when building the network.



