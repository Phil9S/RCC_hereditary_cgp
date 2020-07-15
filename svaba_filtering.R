library(vcfR)

# Filter parameters
depth <- 5
qual <- 10
AD <- 0.2
miss <- 0.2
maf <- 0.05

# Load svaba vcf outputs
vcf.indel <- read.vcfR(file = "cgp_sv_calls.svaba.indel.vcf")
vcf.sv <- read.vcfR(file = "cgp_sv_calls.svaba.sv.vcf")

# svaba indel calls
vcf.indel
## QUAL filtering
vcf.indel <- vcf.indel[which(vcf.indel@fix[,6] > qual),]

## AD filtering
ad.indel <- extract.gt(vcf.indel,element = "AD",as.numeric = T) / extract.gt(vcf.indel, element='DP', as.numeric=TRUE)
ad.indel[ad.indel < AD] <- NA
is.na(vcf.indel@gt[,-1][is.na(ad.indel)]) <- TRUE

## Depth filtering
dp.indel <- extract.gt(vcf.indel, element='DP', as.numeric=TRUE)
dp.indel[dp.indel < depth] <- NA
is.na(vcf.indel@gt[,-1][is.na(dp.indel)]) <- TRUE

## Missingness filtering
indel.maf <- as.data.frame(maf(vcf.indel))
indel.maf$use <- ifelse((indel.maf$`NA` / indel.maf$nAllele)/2 < miss,"TRUE","FALSE")
indel.maf$use <- indel.maf$use[is.na(indel.maf$use)] <- FALSE

## MAF filtering
indel.maf$use <- ifelse(indel.maf$Frequency < maf & indel.maf$Frequency > 0,"TRUE","FALSE")
indel.maf$use <- indel.maf$use[is.na(indel.maf$use)] <- FALSE

vcf.indel <- vcf.indel[which(indel.maf$use == TRUE),]
dim(vcf.indel)

# svaba sv calls
vcf.sv
## QUAL filtering
vcf.sv <- vcf.sv[which(vcf.sv@fix[,6] > qual),]

## AD filtering
ad.sv <- extract.gt(vcf.sv,element = "AD",as.numeric = T) / extract.gt(vcf.sv, element='DP', as.numeric=TRUE)
ad.sv[ad.sv < AD] <- NA
ad.sv[is.nan(ad.sv)] <- NA
is.na(vcf.sv@gt[,-1][is.na(ad.sv)]) <- TRUE

## Depth filtering
dp.sv <- extract.gt(vcf.sv, element='DP', as.numeric=TRUE)
dp.sv[dp.sv < depth] <- NA
is.na(vcf.sv@gt[,-1][is.na(dp.sv)]) <- TRUE

## Missingness filtering
sv.maf <- as.data.frame(maf(vcf.sv))
sv.maf$use <- ifelse((sv.maf$`NA` / sv.maf$nAllele)/2 < miss,"TRUE","FALSE")
sv.maf$use <- sv.maf$use[is.na(sv.maf$use)] <- FALSE

## MAF filtering
sv.maf$use <- ifelse(sv.maf$Frequency < maf & sv.maf$Frequency > 0,"TRUE","FALSE")
sv.maf$use <- sv.maf$use[is.na(sv.maf$use)] <- FALSE

vcf.sv <- vcf.sv[which(indel.maf$use == TRUE),]
dim(vcf.sv)
