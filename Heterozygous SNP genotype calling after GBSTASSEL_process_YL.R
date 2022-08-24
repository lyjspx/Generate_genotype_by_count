## Headers ===================================
## Purpose of script: Take allele counts and re-output 
##                    alleles according to self-defined criteria
## Author: Dr.Yuan Liu
## Date Created: 2020 
## Date Edited: 2022 Aug 
## Email: lydesert@live.com
## Input file: vcf file
## Output file: hmp file
## 
## Validated SessionInfo
## forcats_0.5.1   dplyr_1.0.7     purrr_0.3.4     readr_2.0.2
## tidyr_1.1.4     tibble_3.1.5    ggplot2_3.3.5   tidyverse_1.3.1 stringr_1.4.0   vcfR_1.13.0

#define variable ============================
input_file_path <- "/home/rstudio/rdata/triticale_Li/TRtoCSRYE_sorted(1).vcf"
output_simple_mode <- FALSE #true: heterzygous to 'H'; false: heterzygous to 'R,Y,S,etc.'


library('vcfR')
library('stringr')
library('tidyverse')
vcf.file <- read.vcfR(input_file_path) 

num.snp <- dim(vcf.file@gt)[1]
num.line <- dim(vcf.file@gt)[2]-1

generate_geno <- function(geno, output_simple_mode){
  IUPAC.code <- c('A/G'='R','G/A'='R',
                  'C/T'='Y','T/C'='Y',
                  'G/C'='S','C/G'='S',
                  'A/T'='W','T/A'='W',
                  'G/T'='K','T/G'='K',
                  'A/C'='R','C/A'='M')
  geno.striped <- str_trim(geno)
  names(geno.striped) <- names(geno)
  geno <- geno.striped
  mindepth <- 4 #minimum depth to call genotype
  mindepth_homo <- 4 #minimum depth to call homozygous for each genotype
  mindepth_het <- 4 #minimum depth to call heterozygous for each genotype
  minratio_het <- 0.1 #cutoff ratio of alleles to call heterozygous ie 0.1, 
  #if ratio higher than this will be heterozygous
  minratio_homo <- 0.01 #cutoff ratio of alleles to call homozygous ie 0.1,
  #if ratio lower than this will be homozygous
  minMinor <- 1 #minimum read depth of minor allele to call heterozygous
  
  #mindepth <= mindepth_homo <= mindepth_het
  if(geno['gt_DP'] < mindepth){
    return('N') #None or no enough reads
  }
  
  if(str_count(geno['gt_GT'],'/') > 1){
    return('N') #Multi-alleles
  }
  
  allele.count <- as.numeric(str_split(geno['gt_AD'],',')[[1]])
  allele <- str_split(geno['gt_GT_alleles'],'/')[[1]]
  
  if(max(allele.count) == geno['gt_DP']){ #simple homo
    return(allele[which.max(allele.count)])
  }
  
  if(sum(allele.count) >= mindepth_het){
    if(min(allele.count) < minMinor){
      return('N')
    }
    
    if(min(allele.count)/sum(allele.count) >= minratio_het){
      if(output_simple_mode){
        return('H')
      }else{
        return(IUPAC.code[geno['gt_GT_alleles']])
      }
    }else if(min(allele.count)/sum(allele.count) < minratio_homo){
      return(allele[which.min(allele.count)])
    }
  }
  return('N')
}

#Process snps in a batch of 1000 snps ======================
all.geno <- c()
for(i in 1:ceiling(num.snp/1000)){
  print(paste(c('processing',(i-1)*1000,'-',i*1000,'snps'),collapse = ' '))
  start.snp <- (i-1)*1000 + 1
  end.snp <- min(i*1000, num.snp)
  slice  <- extract_gt_tidy(
    vcf.file[start.snp:end.snp],
    format_fields = NULL,
    format_types = TRUE,
    dot_is_NA = TRUE,
    alleles = TRUE,
    allele.sep = "/",
    gt_column_prepend = "gt_",
    verbose = FALSE
  )
  slice <- slice %>% arrange(Key)
  temp <- apply(slice,1 ,generate_geno,output_simple_mode=output_simple_mode)
  all.geno <- c(all.geno, temp)
}

geno.matrix <- matrix(all.geno,ncol = num.line,byrow = T)
dim(geno.matrix)

colnames(geno.matrix) <- colnames(vcf.file@gt)[2:111]
rownames(geno.matrix) <- vcf.file@fix[,3]

geno.matrix[!geno.matrix %in% c('N','A','T','C','G','H','R','Y','S','W','K','M')] <- 'N'
geno.matrix[1:5,1:5]

## HAPMAP format =======================================
fix.col <- cbind.data.frame(vcf.file@fix[,3],
                            paste0(vcf.file@fix[,4],'/',vcf.file@fix[,5]),
                            vcf.file@fix[,1],
                            vcf.file@fix[,2])
colnames(fix.col) <- c('rs#', 'alleles', 'chrom','pos')
fix.col$strand <- '+'
fix.col$`assembly#` <- 'NA'
fix.col$center <- 'NA'
fix.col$protSID <- 'NA'
fix.col$assayLSID <- 'NA'
fix.col$panel <- 'NA'
fix.col$QCcode <- 'NA'

geno.matrix.hmp <- cbind(fix.col,geno.matrix)
geno.matrix.hmp[1:5,]



#Saving results ===============================
saveRDS(geno.matrix, file = 'processed_geno_grape_4_4_4_geno.matrix.Aug.rds')
saveRDS(all.geno, file = 'processed_geno_grape_4_4_4_all.geno.Aug.rds')
write.table(geno.matrix.hmp,file = 'processed_geno_grape_4_4_4_Aug9.hmp.txt',
            quote = F,sep = '\t',row.names = F)

