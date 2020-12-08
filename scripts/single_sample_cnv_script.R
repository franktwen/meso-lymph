# ExomeDepth RScript
# 12/8/2020
# Kelsey McNeely
#module load gcc/4.9.4 java-jdk/1.7.0_80 zlib/1.2.8 bzip2/1.0.6 xz/5.2.2 curl/7.51.0 pcre/8.39 R/3.1.2

#Run ExomeDepth on a single sample using an existing reference set
#If the individual being tested is related to any members of the reference set 
#they must be removed (Line 40)

#must provide file paths for bam, bai, exons.hg38, and reference set

# Get files to analyze
Bam <- c()

# Enter bam file path for sample to be analyzed
# One for each sample
Bam[[1]] <- 'input bam file path 1' 

BamIndex <- c()

# Enter bai (bam index) file paths for samples to be analyzed
# MUST BE IN THE SAME ORDER AS BAM FILES
BamIndex[[1]] <- 'input bai file path matching bam file 1'


# Load libraries and necessary reference files
library(ExomeDepth)
exons.hg38 <- readRDS("file path for exons.hg38 ")
my.ref <- readRDS('file path to exisiting reference set')


# Get bam counts for samples 
ExomeCount  <- getBamCounts(bed.frame = exons.hg38, bam.files = Bam,index.files = BamIndex, include.chr=FALSE,)  
ExomeCount.dafr <- as(ExomeCount [,colnames(ExomeCount)],'data.frame')

# Select test sample
# Number below needs to be modified depending on total number of control + test samples
my.test <- ExomeCount.dafr [,25]

#If individual is related to members of reference set remove related samples: 
#individual= id of individual to remove
#my.reference.set$individual <- NULL

#my.ref [,6:25] will change based on how many individuals were removed from the reference set
my.reference.set <- as.matrix(my.ref [,6:25])

# Select optimized reference subset from reference set
my.choice <- select.reference.set (test.counts = my.test, reference.counts = my.reference.set, bin.length = (ExomeCount.dafr$end - ExomeCount.dafr$start)/1000, n.bins.reduced=10000)
my.matrix <- as.matrix(my.reference.set[,my.choice$reference.choice], drop=FALSE)
my.reference.selected <- apply(X= my.matrix, MAR=1, FUN=sum)

# Call CNVs
all.exons <- new ('ExomeDepth', test=my.test, reference= my.reference.selected, formula = 'cbind(test,reference) ~ 1')
all.exons <- CallCNVs (x= all.exons,transition.probability= 10^-4, chromosome= ExomeCount.dafr$space, start =ExomeCount.dafr$start, end= ExomeCount.dafr$end, name = ExomeCount.dafr$names)

#print(head(all.exons@CNV.calls))

output.file <- 'output file path'

write.csv(file = output.file, x = all.exons@CNV.calls, row.names = FALSE)










