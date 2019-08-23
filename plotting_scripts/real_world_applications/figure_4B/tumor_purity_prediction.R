############# tumor_purity_prediction.R ######################################################
#' This script executes tumor purity prediction from any WSH score you input. It return a final
#' model consiting of consistently present sites in a LASSO regression and also stores evaluation
#' plots and further diagnostics.
#' @param score.path: path to the results in matrix form (e.g. qfdrp.csv.gz)
#' @param score.annotation: the path to the .RData file containing the genomic annotations of the rows
#'           of score.path
#' @param save.sites: path to a CSV file where we want to store the selected sites used for prediction
score.path <- "path_to_score.csv.gz"
score.annotation <- "path_to_score_annotation.RData"
save.sites <- "selected_sites.csv"

sample.annotation <- "../sample_annotation.csv"


### FUNCTIONS ####
# Execute cross validation for a given function
general.cv <- function(fitFunction,response,data, k = 10, alpha =1){
    nSamples <- length(response)
	rand.init <- sample(1:nSamples,nSamples)
	response <- response[rand.init]
	data <- data[,rand.init]
    size <- nSamples%/%k
    count <- 1
	ret <- c()    
	for(i in seq(1, k * size, by = size)){
            logger.info(as.character(count))
            choose <- rep(FALSE, nSamples)
            choose[i:(i + size - 1)] <- TRUE
            testSet <- data[, choose]
            test.response <- response[choose]
            notChosen <- !choose
            trainSet <- data[, notChosen]
            train.response <- response[notChosen]
            predictor <- fitFunction(t(trainSet),train.response,alpha)
            predicted.response <- predictor(t(testSet))
            cor <- cor(predicted.response, test.response)
            mean <- mean(abs(predicted.response - test.response))
            mean <- round(mean, 2)
            column <- c(cor, mean)
            count <- count + 1
            ret <- cbind(ret,column)
    }
    lastCol <- c(mean(ret[1, 1:10],na.rm=T), mean(ret[2, 1:10],na.rm=T))
    ret <- cbind(ret, lastCol)
    colnames(ret) <- c(paste("Fold", seq(1:10)), "Mean")
    row.names(ret) <- c("Correlation", "Mean")
    return(ret)
}

# returns an executable function
fit.fun <- function(x,y,alpha){
	cv.mod <- cv.glmnet(x,y,alpha=alpha)
	lin.mod <- glmnet(x,y,alpha=alpha,lambda=cv.mod$lambda.min)
	coefs <- as.numeric(coef(lin.mod))
	ret.fun <- function(test.set){
		return(predict(lin.mod,test.set))
	}
	return(ret.fun)
}


### SCRIPT #######
library(RnBeads)
library(glmnet)
library(pheatmap)
library(RColorBrewer)
qfdrp <- read.csv(score.path)
qfdrp <- qfdrp[,-1]
qfdrp <- qfdrp[,unlist(lapply(colnames(qfdrp),function(x)grepl("EwS_T",x)))]
ffpe <- c("EwS_T122","EwS_T123","EwS_T124",
"EwS_T125","EwS_T126","EwS_T127","EwS_T129","EwS_T130","EwS_T131","EwS_T132","EwS_T133")
qfdrp <- qfdrp[,!(colnames(qfdrp)%in%ffpe)]

sample.anno <- read.csv(sample.annotation)
sample.anno <- sample.anno[unlist(lapply(as.character(sample.anno$sample_id),function(x)grepl("EwS_T",x))),]
sample.anno <- sample.anno[!(sample.anno$sample_id%in%ffpe),]
sample.anno <- sample.anno[sample.anno$sample_id!="EwS_T123",]
sample.anno <- sample.anno[match(colnames(qfdrp),sample.anno$sample_id),]
row.names(sample.anno) <- sample.anno$sample_id
tp <- as.character(sample.anno$tumor_purity)
tp[tp %in% "not quantified"] <- NA
sample.anno$tumor_purity <- as.numeric(as.character(tp))

na.qfdrp <- apply(qfdrp,1,function(x)any(is.na(x)))
qfdrp <- qfdrp[!na.qfdrp,]
load(score.annotation)
annotation <- annotation[!na.qfdrp,]

res <- sample.anno$tumor_purity
dat <- t(qfdrp[,!is.na(res)])
res <- res[!is.na(res)]
parallel.setup(10)

# This evaluates the cross-validation errors with resepct to different alpha values
parallel.setup(10)
cv.errors <- foreach(alph = seq(0,1,by=0.1),.combine="c",.export=c("general.cv")) %dopar% {
	ret <- general.cv(fit.fun,res,t(dat),alpha=alph)
	ret[2,11]
}

# This evaluates overall performance for a given alpha value (here 1)
cv.errors <- foreach(i = 1:10,.combine="cbind") %dopar% {
	ret <- general.cv(fit.fun,res,t(dat),alpha=1)
	ret[,11]
}

# We now select those sites that are present, on average, in more than five of the cross-validation 
# folds when executing 10 different random initializations of the CV

res <- sample.anno$tumor_purity
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
sel.sites <- annotation[which(rowMedians(all.res)>5)]
write.csv(sel.sites,save.sites)

## Now we perform simple linear models on the selected sites to evaluate final performance of the model
op <- findOverlaps(sel.sites,annotation,type="equal")

dat <- qfdrp[subjectHits(op),]
res <- sample.anno$tumor_purity
dat <- dat[,!is.na(res)]
res <- res[!is.na(res)]
dat <- as.data.frame(t(dat))
lin.mod <- glm(res~.,data=dat)
vals.train <- predict(lin.mod,dat)
to.plot <- data.frame(Pred=vals.train,Actual=res)
colnames(to.plot)[1] <- "Predicted"
plot <- ggplot(to.plot,aes(x=Predicted,y=Actual))+geom_point()+geom_abline(slope=1,intercept=0)+theme_bw()

# CV
dat <- qfdrp[subjectHits(op),]
res <- sample.anno$tumor_purity
dat <- dat[,!is.na(res)]
res <- res[!is.na(res)]
nSamples <- length(res)
k <- 10
parallel.setup(10)
res.all <- foreach(j = 1:10,.combine="cbind") %dopar% {
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
			testSet <- as.data.frame(t(dat[, choose]))
			test.response <- res[choose]
			notChosen <- !choose
			trainSet <- as.data.frame(t(dat[, notChosen]))
			train.response <- res[notChosen]
			lin.mod <- lm(train.response~.,data=trainSet)
			predictor <- function(test.set){
				return(predict(lin.mod,test.set))
			}
			predicted.response <- predictor(testSet)
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
	ret[,11]
}

