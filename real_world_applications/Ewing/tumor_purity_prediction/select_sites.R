#' select_sites.R
############################################################
#' This scripts selects sites to be associated with tumor purity by employing multiple random
#' initializations of 10 fold cross validations. For each of run, the number of times a particular
#' CpG had a non-zero coefficient in an elastic net regression run is counted. We then select those
#' that, on median, occur in more than 5 folds of the cross-validation run.
#' NOTE: The sample annotation sheet cannot be provided through GitHub. Please contact the Michael
#' Scherer (mscherer[at]mpi-inf.mpg.de) for the sheet or for questions.
wsh.matrix <- "qfdrp.csv.gz"
sample.sheet <- "sample_annotation.csv"
wsh.annotation <- "qfdrp_annotation.RData"

library(RnBeads)
library(glmnet)
library(pheatmap)
library(RColorBrewer)
qfdrp <- read.csv(wsh.matrix)
qfdrp <- qfdrp[,-1]
qfdrp <- qfdrp[,unlist(lapply(colnames(qfdrp),function(x)grepl("EwS_T",x)))]
ffpe <- c("EwS_T122","EwS_T123","EwS_T124",
"EwS_T125","EwS_T126","EwS_T127","EwS_T129","EwS_T130","EwS_T131","EwS_T132","EwS_T133")
qfdrp <- qfdrp[,!(colnames(qfdrp)%in%ffpe)]

ages <- read.csv(sample.sheet)
ages <- ages[unlist(lapply(as.character(ages$sample_id),function(x)grepl("EwS_T",x))),]
ages <- ages[!(ages$sample_id%in%ffpe),]
ages <- ages[ages$sample_id!="EwS_T123",]
ages <- ages[match(colnames(qfdrp),ages$sample_id),]
row.names(ages) <- ages$sample_id
tp <- ages$tumor_purity
tp[tp %in% "not quantified"] <- NA
ages$tumor_purity <- as.numeric(as.character(tp))

na.qfdrp <- apply(qfdrp,1,function(x)any(is.na(x)))
qfdrp <- qfdrp[!na.qfdrp,]
load(wsh.annotation)
annotation <- annotation[!na.qfdrp,]

res <- ages$tumor_purity
dat <- qfdrp[,!is.na(res)]
res <- res[!is.na(res)]
parallel.setup(10)
all.res <- foreach (j = 1:10,.combine="cbind") %dopar% {
	count.cpgs <- rep(0,nrow(dat))
	nSamples <- length(res)
	k <- 10
	rand.init <- sample(1:nSamples,nSamples)
	res <- res[rand.init]
	dat <- dat[,rand.init]
	size <- nSamples%/%k
	count <- 1
	ret <- c()    
	for(i in seq(1, k * size, by = size)){
		    logger.info(as.character(count))
		    choose <- rep(FALSE, nSamples)
		    choose[i:(i + size - 1)] <- TRUE
		    testSet <- dat[, choose]
		    test.response <- res[choose]
		    notChosen <- !choose
		    trainSet <- dat[, notChosen]
		    train.response <- res[notChosen]
			cv.mod <- cv.glmnet(t(trainSet),train.response,alpha=1)
			lin.mod <- glmnet(t(trainSet),train.response,alpha=1,lambda=cv.mod$lambda.min)
			coefs <- as.numeric(coef(lin.mod))[-1]
			count.cpgs[which(coefs != 0)] <- count.cpgs[which(coefs != 0)]+1
			predictor <- function(test.set){
				return(predict(lin.mod,test.set))
			}
		    predicted.response <- predictor(t(testSet))
		    cor <- cor(predicted.response, test.response)
		    mean <- mean(abs(predicted.response - test.response))
		    mean <- round(mean, 2)
		    column <- c(cor, mean)
		    count <- count + 1
		    ret <- cbind(ret,column)
	}
	lastCol <- c(mean(ret[1, 1:10]), mean(ret[2, 1:10]))
	ret <- cbind(ret, lastCol)
	colnames(ret) <- c(paste("Fold", seq(1:10)), "Mean")
	row.names(ret) <- c("Correlation", "Mean")
	count.cpgs
}
count.cpgs <- count(rowMedians(all.res))

sel.cpgs <- annotation[rowMedians(all.res)>5]
write.csv(sel.cpgs,"sel_sites.csv")
