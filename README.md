# HCV_FitnessCost

## Fitness Cost Analysis of HCV1A


* Data directory contains sample information. SampleSheet_Mac251.csv has the most comprehensive sample information.
* Output directory contains all output files, except the bam files (due to the memory limit)
* Rscripts directory contains all scripts for the analysis. 


## Analysis Steps 

### Pre-step 1. Create bash files to process FASTQ files with CreateBashFiles.R 

	• Template files are in Data/template (raw FASTQ files are not available for privacy reasons).
	• Processed bam files are available at https://figshare.com/articles/dataset/HCV_Bam_Files/13239491 (download bam2_1a.zip)
	• Trimming/filtering was done in BBTools and mapping was done in BWA.
    • Need sam for the next step (not publicly provided) 
    
### Pre-step 2. Final quality check with 1.Reads.check.R   
    • Requires SAM files mapped to their own consensus sequences 
    • Eliminate potential chimeric or mismapped sequences 
    • Need to convert the output sam files to bam/ban.bai for next steps


### Step 1. Create frequency tables for each sample/population (Scripts 2-4)
	• Requires BAM/BAI files mapped to a consensus sequence of each file to start 
	• CSV files in Overview2 (unfiltered datasets) 
	
### Step 2. Create mutation frequency summaries and figures (Scripts 5.1-5.3)
	• The filtered (reads >100) overview files are Overview3.
	
### Step 3. Beta regression on mutation frequencies (Scripts 6.1-6.3)
	• Conduct beta regression to udnerstand factors affecting mutation frequencies at each site.
	
### Step 4. Selection coefficient (fitness cost) analysis (Scripts 7.1-7.7)
	• Summarize, analyze, run stats, run beta regression and create figures on selection coefficients

### Step 5. Additional analysis (Scripts 8.1-8.8)
	• 8.1. Estiamte mutation rates from stop (non-sense) mutations
	• 8.2. Compare in-vitro data from Geller et al (2016)
	• 8.3. Look at the known drug resistance sites
	• 8.4. Assess high-frequency and low-frequency regions identified in-vitro by Geller et al.
	• 8.5. Assess highly conserved sites across populations
	• 8.6. Assess correlation between within-host and between-host  
	• 8.7. Estimate gamma parameters for DFE (Distribution of Fitness Effects)
	• 8.8. Create a plot ordered by mut freq rankings

	

     
