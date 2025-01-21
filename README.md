Background: CpG methylation is an epigenetic marker that varies across tissue types. However, the methylation status of a single CpG site is unreliable as a biomarker due to errors introduced by bisulfite sequencing, sampling techniques, and biological variability.
Definition: Phased Methylation Pattern (PMP) is a unique set of coordinates that includes the DNA strand (‘f’ for forward (+) or ‘r’ for reverse (-)), the relative positions of three CpG sites on the same strand (e.g., x:y:z), and their methylation status (e.g., ‘000’ for all unmethylated or ‘111’ for all methylated). It represents a combined epigenetic signature across these CpG sites.
In this project we try to find out if : Phased methylation patterns (PMPs) can act as reliable biomarkers to differentiate tissue types, providing higher specificity compared to individual CpG sites.

The python codes written will give the following information from the coordinates and relative positions of 3 CpG sites and their methylation status in a .csv file 
1. Median and coefficient of variation (CV) for single CpG coverage in each tissue and coverage statistics plotting.
2. Identify PMPs with high specificity, minimizing false positives using a statistical approach to assign confidence
3. Identify PMPs with high specificity, minimizing false positives using ML approach
4. Identify PMPs with high specificity, minimizing false positives using statistical and ML approach
5. Calculate the mean variant read fraction (VRF) for each PMP in both tissues
6. Sequencing depth and specificity confidence computing
7. Validation of the hypothesis of whether the sequencing depth affect confidence by comparing the specificity of the top 10 PMPs against individual CpG sites
