This is an organization of data for the manuscript "Temporal dynamics and tissue-specific variations of the blueberry phyllosphere mycobiome", including post-processed data and the code used to examine said data. Raw sequencing data generated for this study can be found in the NCBI SRA archive as Bioproject PRJNA1161654. 
The code should be properly annotated and polished to help repeat analysis. 

As is, the folders are structured to be compatible with how the code was written, but as long as paths are right the code should work with alternative folder structures.

Thank you for your interest in this study :)

some major packages used are not available through CRAN
from CRAN, install devtools
then, run install.packages("BiocManager") to ease installation from Bioconductor
BiocManager::install("phyloseq")
BiocManager::install("metagenomeSeq")
BiocManager::install("SpiecEasi")
BiocManager::install("RCy3") #this one allows to connect to Cytoscape with cytoscapePing()
For this project, the Functional Guilds from Funguilds were used.
The one utilized in this study has been provided (coredata/FunctionalGuilds.csv)
General file organization:
/timecourse/ #base folder
/timecourse/code #contains the code to be used
/timecourse/coredata #contains core data used for assembling phyloseq data and for further analysis
/timecourse/Figures/ #storage for figures generated
/timecourse/Generated Tables/ #generated tables from parts of analysis largely go here
/timecourse/Generated Tables/decon #generated tables from microDecon. 
/timecourse/Generated Tables/SpiecEasi #generated tables for use with SpiecEasi
/timecourse/Generated Tables/SpiecEasi/RDSobj/ #generated RDS objects from SpiecEasi analysis for easier reloading. 
