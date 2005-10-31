library(MASS)
library(caMassClass)
library(grDevices)

loadFiles <- function(directory) {
  data <- NULL
  files <- dir(directory)
  for (file in files){
    tmpData <- read.csv(file.path(directory, file))
    data$int <- rbind(data$int, tmpData$Intensity)
  }
  return(data)
}

load.mzXML.directory <- function(dirname) {
  data <- NULL
  files <- dir(dirname)
  for (file in files) {
    tmpData <- load.mzXML.file(file.path(dirname, file))
    data$peaks <- rbind(data$peaks, tmpData$peaks)
    data$mass <- rbind(data$mass, tmpData$mass)
  }
  return(data)
}

load.mzXML.file <- function(filename) {
  data <- NULL
  mz <- read.mzXML(filename)
  for (scan in mz$scan) {
    data$int <- rbind(data$int, scan$peaks)
#    data$mass <- rbind(data$mass, scan$mass)
  }
  return(data)
}

normalize <- function(mat) {
  for (ixRow in seq(nrow(mat))) {
    rowMin <- min(mat[ixRow,])
    rowRange <- max(mat[ixRow,]) - rowMin
    mat[ixRow,] <- (mat[ixRow,] - rowMin) / rowRange
  }
  return(mat)
}

colMeans <- function(mat) {
  means <- rep(0, 1, ncol(mat))
  for (ixCol in seq(ncol(mat))) {
    means[ixCol] <- mean(mat[,ixCol])
  }
  return(means)
}

pcaByHand <- function(mat, cDropEigenvecs) {
  # set each instance's mean to zero
  means <- colMeans(mat)
  for (ix in seq(nrow(mat))) {
    mat[ix,] <- mat[ix,] - means
  }

  eig <- eigen(mat %*% t(mat))
  # signs may differ between eig$vectors and vec in prepareData of Matlab Q5
  eig$vectors <- t(mat) %*% eig$vectors

  for (ixVec in seq(ncol(eig$vectors))) {
    eig$vectors[,ixVec] <- eig$vectors[,ixVec] / max(svd(eig$vectors[,ixVec])$d)
  }

  # drop the eigenvectors we're not keeping
  eig$vectors <- eig$vectors[,1:(ncol(eig$vectors) - cDropEigenvecs)]
  
  mat <- mat %*% eig$vectors
  # in the Matlab version:
  #    mat[1,1] == Q(1, end) && mat[1,ncol(mat)] == Q(1,1)
  return(list(mat=mat, eigVecs=eig$vectors, eigVals=eig$values))
}

pcaProject <- function(mat, eigVecs){
  # set each instance's mean to zero
  means <- colMeans(mat)
  for (ix in seq(nrow(mat))) {
    mat[ix,] <- mat[ix,] - means
  }

  # project into the pca space
  mat <- mat %*% eigVecs

  return(mat)
}

q5Discriminant <- function(trainData, trainClasses, normalizeFunction, pcaFunction, ldaFunction, cDropEigenVecs) {
  trainData <- normalizeFunction(trainData)
  pcaRes <- pcaFunction(trainData, cDropEigenVecs)

  discriminant <- ldaFunction(pcaRes$mat, trainClasses)
  return(list(discriminant=discriminant, eigVecs=pcaRes$eigVecs, eigVals=pcaRes$eigVals))
}

classifyData <- function(testData, normalizeFunction, projectionFunction, discriminant, eigVecs) {
  testData <- normalizeFunction(testData)

  testData <- projectionFunction(testData, eigVecs)

  results <- predict(discriminant, testData)
  return(results)
}

q5 <- function(trainData, trainClasses, testData, testClasses, normalizeFunction, pcaFunction, ldaFunction, cDropEigenVecs, projectionFunction){
  # for trainClasses and testClasses, 1 indicates disease.
  discriminant <- q5Discriminant(trainData, trainClasses, normalizeFunction, pcaFunction, ldaFunction, cDropEigenVecs)
  testDataResults <- classifyData(testData, normalizeFunction, projectionFunction, discriminant$discriminant, discriminant$eigVecs)

  truePos <- length(intersect(which(testDataResults$class == 1), 
  which(testClasses == 1)))

  falsePos <- length(intersect(which(testDataResults$class == 1), which(testClasses != 1)))

  trueNeg <- length(intersect(which(testDataResults$class != 1), 
  which(testClasses != 1)))

  falseNeg <- length(intersect(which(testDataResults$class != 1), which(testClasses == 1)))
  
  positivePredictiveValue <- truePos / (truePos + falsePos)
  percentageCorrectlyClassified <- (truePos + trueNeg) / (truePos + falsePos + trueNeg + falseNeg)
  sensitivity <- truePos / (truePos + falseNeg)
  specificity <- trueNeg / (falsePos + trueNeg)

  return(list(discriminant=discriminant, positivePredictiveValue=positivePredictiveValue, percentageCorrectlyClassified=percentageCorrectlyClassified, sensitivity=sensitivity, specificity=specificity))
}

cDropEigenVecs = 10
#cdir <- 'C:\\Documents and Settings\\Paul K. Courtney\\Desktop\\Q5 Demo\\control'
#controlData <- loadFiles(cdir)
controlData <- load.mzXML.file('C:\\Documents and Settings\\Paul K. Courtney\\Desktop\\Q5 Demo\\control_mzXML\\control.mzXML')
caseData <- load.mzXML.file('C:\\Documents and Settings\\Paul K. Courtney\\Desktop\\Q5 Demo\\cancer_mzXML\\cancer.mzXML')
#ncdir <- 'C:\\Documents and Settings\\Paul K. Courtney\\Desktop\\Q5 Demo\\cancer'
#caseData <- loadFiles(ncdir)

# ixTraCon means training indices of control matrix.
ixTraCon <- 1:floor(nrow(controlData$int) / 2)
# ixTraCase means training indices of case matrix
ixTraCase <- 1:floor(nrow(caseData$int) / 2)

# ixTesCon means testing indices of control matrix.
ixTesCon <- (max(ixTraCon)+1):nrow(controlData$int)
# ixTesCase means testing indices of case matrix
ixTesCase <- (max(ixTraCase)+1):nrow(caseData$int)

trainData = rbind(controlData$int[ixTraCon,], caseData$int[ixTraCase,])
trainClasses = c(rep(0, length(ixTraCon)), rep(1, length(ixTraCase)))

testData = rbind(controlData$int[ixTesCon,], caseData$int[ixTesCase,])
testClasses = c(rep(0, length(ixTesCon)), rep(1, length(ixTesCase)))

q5Res <- q5(trainData, trainClasses, testData, testClasses, normalize, pcaByHand, lda, cDropEigenVecs, pcaProject)

png('eigenVals.png')
plot(q5Res$discriminant$eigVals, ylab="Eigenvalues")
dev.off()
