Implementation of the regression method from Cai, Ling, et al. "Genomic regression analysis of coordinated expression." Nature communications 8.1 (2017): 2187. to correct the effect of cnvs on rna seq expression values. Used with BRCA dataset from TCGA

Steps: 
1. Get the Segmented files. Download manifest from TCGA. We're currently working with Masked Copy Number Segment files. 

gdc-client download -m $MANIFESTDIR/manifest.txt -d $DOWNLOADDIR --log-file download-segmented-files-healthy.log --retry-amount 3

2. Put everything into one single segmented file.
find $DOWNLOADDIR -name '*.seg.v2.txt' -exec mv '{}' $DOWNLOADDIR \;
tail -q -n +2 *.v2.txt > $SEGMENTEDFILE

3. Run GISTIC2 using the run_gistic.sh script. (Copy the script to the GISTIC2 home) 
markersfile was obtained from GDC: https://gdc.cancer.gov/about-data/data-harmonization-and-generation/gdc-reference-files/
snp6.na35.remap.hg38.subset.txt.gz

4. Segmented files have aliquot IDs. Use the aliquote_submitterid_match.py to query GDC and get the submitter ID. Use cnv-id-submitterid-match.R to match aliquotes and submitter IDs.  

5. Expression matrix has our own IDs. use exp-id-submitterid-match.R to match submitter IDs. File exp-caseids-cancer.csv obtained through the GDC API.

6. Match expression and cnvs cases using exp-cnv-match.R. 

7. Perform correction using exp-cnv-correlation.R
