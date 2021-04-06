# ArrayExpress-data-into-Bioconductor-object

Functionality overview

(1) Downloads microarray genomics data from ArrayExpress repository.  
(2) Preprocesses data into phenotype and assayData.  
(3) Creates eSet objects based on the type of experiment.  
(4) Built on latest libraries.

**Requirements**  
require("XML")  
require(RCurl)  
require(affy)  
require(Biobase)  
require(limma)  
require(oligo)  
require(CCl4)  
require(stringr)  


**Conversion**  
<img src="https://user-images.githubusercontent.com/37441690/113650737-59342e80-96ae-11eb-8df5-2ea001edd123.png" width="800" height="450">


**ExpressionMatrix to Phenotype data and featuredata**

![image](https://user-images.githubusercontent.com/37441690/113650871-95678f00-96ae-11eb-91d9-733f7ca5df19.png)
