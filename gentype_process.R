library('vcfR')
library('stringr')
#library('tidyverse')
vcf.file <- read.vcfR("C:/Users/liu.yuan/Desktop/grapebypop/grape.keyfile_May18_GBS63.txtCov0.1.MAF0.01.vcf")

vcf.file@meta[1:5]
vcf.file@fix[1:5,]
vcf.file@gt[1:5,1:5]

dim(vcf.file@gt)


# all.DP <- apply(slice1,c(1,2),function(x) str_split(x,":")[[1]][3])
# 
# all.DP <- c()
# for(i in 1:98){
#   start.snp <- (i-1)*10000 + 1
#   end.snp <- min(i*10000, 972383)
#   slice <- vcf.file@gt[start.snp:end.snp,2:385]
#   slice.DP <- apply(slice,c(1,2),function(x) str_split(x,":")[[1]][3])
#   all.DP <- c(all.DP, slice.DP)
# }
# 
# dp.table <- table(as.numeric(all.DP))
# barplot(dp.table[1:20])

# vcf.long_format <- extract_gt_tidy(
#   vcf.file,
#   format_fields = NULL,
#   format_types = TRUE,
#   dot_is_NA = TRUE,
#   alleles = TRUE,
#   allele.sep = "/",
#   gt_column_prepend = "gt_",
#   verbose = TRUE
# )

#vcf.long_format$processd <- 'N'
#vcf.long_format[1:5,]




generate_geno <- function(geno){
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
      return('H')
    }else if(min(allele.count)/sum(allele.count) < minratio_homo){
      return(allele[which.min(allele.count)])
    }
  }
  
  return('N')
  
}

#((dim(vcf.file)[1] %/% 10000)+1)
all.geno <- matrix(ncol=dim(vcf.file)[3]-1)
for(i in 1:((dim(vcf.file)[1] %/% 10000)+1)){
  start.snp <- (i-1)*10000 + 1
  end.snp <- min(i*10000, dim(vcf.file)[1])
  print(paste(i,Sys.time()))
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
  
  temp <- apply(slice,1 ,generate_geno)
  
  all.geno <- rbind(all.geno, matrix(temp,ncol=dim(vcf.file)[3]-1,byrow = F))
  
  }


geno.matrix <- all.geno[-1,]
dim(geno.matrix)


sum(SNP.missing < 200)

vcf.file@gt[1:5,1:5]


colnames(geno.matrix) <- colnames(vcf.file@gt)[2:(dim(vcf.file)[3])]
rownames(geno.matrix) <- vcf.file@fix[,3]
geno.matrix[!geno.matrix  %in% c('N','A','T','C','G','H')] <- 'N'
geno.matrix[1:5,1:5]

#write.table(geno.matrix,file = 'processed_geno_grape_4_4_4.May31.txt',quote = F)
save(geno.matrix, file = 'processed_geno_grape_GBS63_4_4_4.May31.rdata')

#geno.matrix <- read.table('processed_geno_grape_4_4_4.May31.txt',header = F,skip = 1)

SNP.missing <- apply(geno.matrix,1, function(x) sum(x == 'N',na.rm = T))

par(ps=8)
hist(SNP.missing,breaks = 50,labels = T)

dim(geno.matrix[SNP.missing<362*0.8,])
dim(geno.matrix[SNP.missing<362*0.7,])

write.csv(geno.matrix[SNP.missing<362*0.8,],
          file = 'processed_geno_grape_GBS63_4_4_4_NA80.csv',
          quote = F)
write.csv(geno.matrix[SNP.missing<362*0.7,],
            file = 'processed_geno_GBS63_4_4_4_NA70.csv',
            quote = F)
write.csv(geno.matrix[SNP.missing<362*0.6,],
          file = 'processed_geno_GBS63_4_4_4_NA60.csv',
          quote = F)
write.csv(geno.matrix[SNP.missing<362*0.5,],
          file = 'processed_geno_GBS63_4_4_4_NA50.csv',
          quote = F)





# test - deprecated ------------------------
## HAPMAP format
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

write.table(geno.matrix.hmp,file = 'processed_geno_grape_4_4_4.hmp.txt',
            quote = F,sep = '\t',row.names = F)



gt <- extract.gt(vcf.file[1:10000,], return.alleles = TRUE)

gt[1:5,1:5]

gt.count <- extract.gt(vcf.file, return.alleles = F)

gt.count[1:15,1:15]

gt.dp <- extract.gt(vcf.file,element = "DP",as.numeric = T)
gt.dp[1:5,1:10]

gt.ad <- extract.gt(vcf.file,element = "AD",as.numeric = F)
gt.ad[1:15,1:5]

gt.ad <- extract.gt(vcf.file,element = "AD",as.numeric = F)
gt.ad[1:15,1:5]

