#' predict_tumor_purity.R
##########################################
#' This script predicts tumor purity levels from WSH score matrices and a predefined selection of sites
#' provided through the argument *sel.sites*. You can choose one of the predictors from the directory
#' and should specify the corresponding WSH score data matrix
sel.sites <- "predictors/sel_sites_qfdrp.csv"
wsh.location <- "qfdrp.csv.gz"
wsh.annotation "qfdrp_annotation.RData"

library(GenomicRanges)
qfdrp <- read.csv(wsh.location)
qfdrp <- qfdrp[,-1]
qfdrp <- qfdrp[,unlist(lapply(colnames(qfdrp),function(x)grepl("EwS_T",x)))]
ffpe <- c("EwS_T122","EwS_T123","EwS_T124",
"EwS_T125","EwS_T126","EwS_T127","EwS_T129","EwS_T130","EwS_T131","EwS_T132","EwS_T133")
qfdrp <- qfdrp[,!(colnames(qfdrp)%in%ffpe)]
load(wsh.annotation)
sel.sites <- read.csv(sel.sites)
sel.anno <- makeGRangesFromDataFrame(sel.sites[-1,])
op <- findOverlaps(sel.anno,annotation,type="equal")

reg.coefs <- as.matrix(sel.sites$Regression.Coefficient[-1])
intercept <- as.numeric(sel.sites$Regression.Coefficient[1])
dat <- qfdrp[subjectHits(op),]
pred.purity <- t(as.matrix(dat))%*%reg.coefs+intercept

