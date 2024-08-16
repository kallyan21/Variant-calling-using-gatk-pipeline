This is the pipeline for variant calling for whole exome sequence using GATK pipeline. 
Prior to working with this pipeline make sure you have checked your raw sequenece, if quality is not good then do some trim and qc of your sequence. Download all the resources from https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?pli=1&prefix=&forceOnObjectsSortingFiltering=true.

Usages:
Use wes_new.sh for mapping and intitial gvcf calling for individual sample.
And final_variant_calling_new.sh for final variant calling.
