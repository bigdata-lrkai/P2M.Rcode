# Code used to perform online analysis at P2M.db

###### P2M.db is equipped with 13 online analysis tools designed to characterize the differences between primary and metastatic tumors. 

* Among these tools are eight methodologies based on **scRNA-seq data**, including differential gene expression analysis, cell interaction analysis, cell trajectory analysis, and cell functional state analysis. 
* Additionally, there are five tools focused on **bulk data**, such as protein interaction analysis, immune cell infiltration analysis, and survival analysis.
* Enhanced with visualization plugins like Echarts, Highcharts, and Plotly, **our database provides users with an immersive visual interface, facilitating data exploration and interpretation.**

###### We provide all the original code used for online analysis. Here is a detailed introduction to the online analysis tools supported by P2M.db.

---

<img src="md/Figure 3.jpg" alt="Figure 3" style="zoom:80%;" align="left"/>


**The detailed functions of the online tool for scRNA-seq data are described as follows**

1. ***Data Characteristics*.** This function shows the difference in the number of counts, the number of genes, and the percent of mitochondrial genes detected between the primary and metastatic samples, which can reflect the quality of the current data to a certain extent.
2. ***Cell Cluster and Cell type*.** This function demonstrates the difference between cell types in primary and metastatic samples, as well as the distribution and number of various cell types in two-dimensional space.
3. ***Differentially Expressed Genes*.** This function calculates significantly differentially expressed genes between primary and metastatic samples, including analysis of primary vs. metastasis as a whole and for each cell type individually.
4. ***Functional Characteristics*.** Functional enrichment analysis of differentially expressed genes in different cell types between primary and metastasis.
5.  ***Cell Interaction Characteristics*.** Compare the difference in cell interaction by describing the cell interaction relationship in primary and metastasis.
6. ***Cell Trajectory Analysis*.** Comparison of cell differentiation characteristics between primary and metastatic sites by cell trajectory analysis.
7. ***Gene Characteristics*.** This function shows the proportion and means of expression of genes with a high coefficient of variation in primary and metastasis.
8. ***Cell Functional State*.** Compare the differences in 18 cellular functional states, including G1/S, G2/M, M/G1, Angiogenesis, Differentiation, DNA damage, DNA repair, EMT, Hypoxia, Invasion, Metastasis, Proliferation, Quiescence, and Stemness, between primary and metastatic sites.

---

<img src="md/Figure 4.jpg" alt="Figure 4" style="zoom:80%;"  align="left"/>


**The detailed functions of the online tool for bulk data are described as follows**

1. ***Differentially Expressed Genes*.** This function calculates significantly differentially expressed genes between primary and metastatic samples, including analysis of primary vs. metastasis as a whole and for each cell type individually.
2. ***Functional Characteristics*.** Functional enrichment analysis of differentially expressed genes in different cell types between primary and metastasis.
3. ***Protein-protein Interaction*.** Comparing the differences in protein interactions between primary and metastatic samples through protein interaction analysis
4. ***Immune Infiltration Analysis*.** Comparing the difference in immune cell proportion between primary and metastatic samples through immune infiltration analysis
5.  ***Survival Analysis*.** Comparing the effects of different genes on the survival status of primary and metastatic samples

