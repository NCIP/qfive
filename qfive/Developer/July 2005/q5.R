library(stats)
library(MASS)

loadFiles <- function(directory){
  data <- NULL
  files <- dir(directory)
  for (file in files){
    tmpData <- read.csv(file.path(directory, file))
    data$int <- rbind(data$int, tmpData$Intensity)
    data$mz <- rbind(data$mz, tmpData$M.Z)
  }
  return(data)
}

normalize <- function(mat) {
  for (ixRow in seq(nrow(mat))) {
    rowMin = min(mat[ixRow,])
    rowRange = max(mat[ixRow,]) - rowMin
    mat[ixRow,] = (mat[ixRow,] - rowMin) / rowRange
  }
  return(mat)
}

colMeans <- function(mat) {
  means = rep(0, 1, ncol(mat))
  for (ixCol in seq(ncol(mat))) {
    means[ixCol] = mean(mat[,ixCol])
  }
  return(means)
}

pcaByHand <- function(mat, cDropEigenvecs) {
  # set each instance's mean to zero
  means = colMeans(mat)
  for (ix in seq(nrow(mat))) {
    mat[ix,] = mat[ix,] - means
  }

  eig = eigen(mat %*% t(mat))
  # signs may differ between eig$vectors and vec in prepareData of Matlab Q5
  eig$vectors = t(mat) %*% eig$vectors

  for (ixVec in seq(ncol(eig$vectors))) {
    eig$vectors[,ixVec] <- eig$vectors[,ixVec] / max(svd(eig$vectors[,ixVec])$d)
  }

  # drop the eigenvectors we're not keeping
  eig$vectors = eig$vectors[,1:(ncol(eig$vectors) - cDropEigenvecs)]
  
  mat = mat %*% eig$vectors
  # mat[1,1] == Q(1, end) && mat[1,ncol(mat)] == Q(1,1)
  return(list(mat=mat, eigVecs=eig$vectors, colMeans=means))
}

q5 <- function(trainData, trainClasses, cDropEigenvecs) {
  trainData = normalize(trainData)
  pcaRes = pcaByHand(trainData, cDropEigenvecs)

  discriminant <- lda(pcaRes$mat, trainClasses)
  return(list(discriminant=discriminant, colMeans=pcaRes$colMeans, eigVecs=pcaRes$eigVecs))
}

classifyData <- function(testData, colMeans, eigVecs, discriminant) {
  testData <- normalize(testData)
  for (ix in seq(nrow(testData))) {
    testData[ix,] = testData[ix,] - colMeans
  }

  # project into the pca space
  testData <- testData %*% eigVecs

  results = predict(discriminant, testData)
  return(results)
}

cDropEigenvecs = 10
cdir <- 'data/080702Data/Control'
controlData <- loadFiles(cdir)
ncdir <- 'data/080702Data/OvarianCancer'
caseData <- loadFiles(ncdir)

# ixTraCon means training indices of control matrix.
ixTraCon <- 1:floor(nrow(controlData$int) / 2)
# ixTraCase means training indices of case matrix
ixTraCase <- 1:floor(nrow(caseData$int) / 2)

# ixTesCon means testing indices of control matrix.
ixTesCon <- (max(ixTraCon)+1):nrow(controlData$int)
# ixTesCase means testing indices of case matrix
ixTesCase <- (max(ixTraCase)+1):nrow(caseData$int)

traData = rbind(controlData$int[ixTraCon,], caseData$int[ixTraCase,])
traClasses = c(rep(1, length(ixTraCon)), rep(2, length(ixTraCase)))

q5Res <- q5(traData, traClasses, cDropEigenvecs)

tesData = rbind(controlData$int[ixTesCon,], caseData$int[ixTesCase,])

results <- classifyData(tesData, q5Res$colMeans, q5Res$eigVecs, q5Res$discriminant)
# got all of the controls right
all(which(results$class==1)==1:length(ixTesCon))
# got all of the cases right
all(which(results$class==2)==(length(ixTesCon))+1:(length(results$class) - length(ixTesCon)))
