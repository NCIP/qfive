##Testing Code for Q5
##Author: Ted Laderas
##Note: This script does not include the code used for testing SVM
##which was used for performance comparison.



TestTrainSplits <- function(group1, group2, numsplits=1, percent=0.50) {

        group1.samples <- nrow(group1)
        group2.samples <- nrow(group2)

	  #calculate training set sizes from each group
        group1.trainsize <- floor(percent * group1.samples)
        group2.trainsize <- floor(percent * group2.samples)

	  #calculate test set sizes from each group
        group1.testsize <- group1.samples - group1.trainsize
        group2.testsize <- group2.samples - group2.trainsize

        #splits <- list(length=numsplits)

	  #initialize splits object
        splits <- vector(mode="list", length=numsplits)

        for(i in 1:numsplits) {

		    #vector of all possible scans for group 1
                group1.vec <- c(1:group1.samples)
		    #vector of all possible scans for group 2
                group2.vec <- c(1:group2.samples)

		    #randomly sample scans for training group
                group1.trainindex <- sample(group1.vec, group1.trainsize)
                group2.trainindex <- sample(group2.vec, group2.trainsize)

		    #test group consists of those indices not in training group
                group1.testindex <- group1.vec[-group1.trainindex]
                group2.testindex <- group2.vec[-group2.trainindex]

		    #traingroups is training assignments
                traingroups <- c(rep(1, group1.trainsize), rep(2, group2.trainsize))
		    #concatenate training indices for group 1 and group 2
                trainindex <- c(group1.trainindex, group2.trainindex)

		    #test groups is testing assignments
                testgroups <- c(rep(1, group1.testsize), rep(2, group2.testsize))
		    #concatenate training indice for group 1 and group 2
                testindex <- c(group1.testindex, group2.testindex)

		    #save train assignments and indices as a data.frame
                train <- 
                        as.data.frame(cbind(trainindex, traingroups))

		    #save test assignments and indices as a matrix
                test <- 
                        as.data.frame(cbind(testindex, testgroups))

		    #save train and test into a slot of splits
                splits[[i]] <- list(train=train, test = test)
        }

	  #return test/train splits
        splits

}



run.Q5 <- function(group1, group2, splits)  {

    	  #initialize results objects
        statmatrix <- NULL
        discrimlist <- NULL

        resultlist <- vector(mode = "list", length = length(splits))
        for(i in 1:length(splits)) {

                #it <- paste("Train/Test Splits: ", i,  sep="")
                #print(it)
                trainassign <- splits[[i]]$train

		    #select those indices for the train set from group 1
                traingroup1 <- as.vector(subset(trainassign, traingroups == 1, 
                        select = (trainindex)))

		    #select those indices for the train set from group 2
                traingroup2 <- as.vector(subset(trainassign, traingroups == 2,
                        select = (trainindex)))

                testassign <- splits[[i]]$test

		    #select those indices for the test set from group 1
                testgroup1 <- as.vector(subset(testassign, testgroups ==1,
                        select = (testindex)))
		
		    #select those indices for the test set from group 2
                testgroup2 <- as.vector(subset(testassign, testgroups ==2,
                        select = (testindex)))

		    #grab assignments 
                trainclass <- as.vector(splits[[i]]$train$traingroups)
                testclass <- as.vector(splits[[i]]$test$testgroups)

		    #build the training and test matrices by subsetting group1 and 2
                trainset <- rbind(group1[traingroup1$trainindex,], 
                        group2[traingroup2$trainindex,])
                testset <- rbind(group1[testgroup1$testindex,], 
                        group2[testgroup2$testindex,])

                print("running Q5 test")
                print(i)

		    #run Q5 with subsetted matrices and assignments
                results <- q5(trainData = trainset, trainClasses = trainclass, 
                        testData = testset, testClasses = testclass)

		    #gather statistics into a list
                stats <- list(results$positivePredictiveValue,
                	results$percentageCorrectlyClassified,
                	results$sensitivity,
                	results$specificity)

                stats <- as.vector(stats, mode = "numeric")
                statmatrix <- rbind(statmatrix, stats)
                discrimlist <- c(discrimlist, results$discriminant)
        }
        colnames(statmatrix) <- c("PPV", "%CC", "Sens", "Spec")

	  #save statistics and discriminant
        resultlist <- list(statmatrix = statmatrix, discrimlist = discrimlist)
	  resultlist
}


multipleruns <- function(group1, group2, numsplits = 2, percent = 0.5) {

	  #generate test/train splits
        splits <- TestTrainSplits(group1, group2, numsplits = numsplits, percent)
	  #run Q5 on splits
        results <- run.Q5(group1, group2, splits)

	  #save results and splits in a list
        resultslist <- list(results = results, splits =splits)
        resultslist

}

#runRatios automates running 3 sets of train/test splits
#at 3 different ratios (0.5, 0.75, 0.9)
#main function used to conduct the analyses

runRatios <- function(group1=grp3, group2=grp4, numsplits = 100) {

        splitsize <- c(.5, .75, .9)

	  #initialize results object
        results <- vector(mode = "list", length = length(splitsize))

        for(i in 1:length(splitsize)) {

		    #run Q5 at a single train/test split size
                r1 <- multipleruns(group1, group2, numsplits, percent =
                        splitsize[i])
		    #save results into a slot
                results[[i]] <- r1
        }

        results
}


##following is the actual script used for testing
#open q5 library
library(Q5)
#declare split ratios
splitsize <- c(0.5, 0.75, 0.9)

#run comparisons on non winsorized data
#run grp3 vs grp 4 comparison
grp3.vs.grp4 <- runRatios(grp3, grp4, numsplits = 100)

#Note: Test-Train splits can be accessed for a single ratio
#as: grp3.vs.grp4[[1]]$splits <- gives all 100 test/train splits for train/test ratio of 0.5
#grp3.vs.grp4[[2]]$splits <- gives all 100 test/train splits for train/test ratio of 0.75
#grp3.vs.grp4[[3]]$splits <- gives all 100 test/train splits for train/test ratio of 0.9

#print statistics for grp 3 vs grp 4 comparison
for(i in 1:length(grp3.vs.grp4)) {
	#for each set of splits (0.5, 0.75, 0.9) report statistics
	#mean is reported first
	print(paste("Mean:", splitsize[i]))
	print("PPV\t%CC\tSens\tSpec")
	apply(grp3.vs.grp4[[i]]$results$statmatrix, 2, mean)
	#sd is reported next
	print(paste("SD:", splitsize[i]))
	print("PPV\t%CC\tSens\tSpec")
	apply(grp3.vs.grp4[[i]]$results$statmatrix, 2, sd)
}

#run group1 versus group 4 comparison on nonwinsorized data
grp1.vs.grp4 <- runRatios(grp1, grp4, numsplits = 100)

for(i in 1:length(grp1.vs.grp4)) {
	#for each set of splits (0.5, 0.75, 0.9) report statistics
	#mean is reported first
	print(paste("Mean:", splitsize[i]))
	print("PPV\t%CC\tSens\tSpec")
	apply(grp1.vs.grp4[[i]]$results$statmatrix, 2, mean)
	#sd is reported next
	print(paste("SD:", splitsize[i]))
	print("PPV\t%CC\tSens\tSpec")
	apply(grp1.vs.grp4[[i]]$results$statmatrix, 2, sd)
}

grp1.vs.grp3 <- runRatios(grp1, grp3, numsplits = 100)

for(i in 1:length(grp1.vs.grp3)) {
	#for each set of splits (0.5, 0.75, 0.9) report statistics
	#mean is reported first
	print(paste("Mean:", splitsize[i]))
	print("PPV\t%CC\tSens\tSpec")
	apply(grp1.vs.grp3[[i]]$results$statmatrix, 2, mean)
	#sd is reported next
	print("PPV\t%CC\tSens\tSpec")
	print(paste("SD:", splitsize[i]))
	apply(grp1.vs.grp3[[i]]$results$statmatrix, 2, sd)
}

#run comparisons on winsorized data
grp3.vs.grp4.winsor <- runRatios(grp3.winsor, grp4.winsor, numsplits = 100)

#run statistics
for(i in 1:length(grp3.vs.grp4.winsor)) {
	#for each set of splits (0.5, 0.75, 0.9) report statistics
	#mean is reported first
	print(paste("Mean:", splitsize[i]))
	print("PPV\t%CC\tSens\tSpec")
	apply(grp3.vs.grp4.winsor[[i]]$results$statmatrix, 2, mean)
	#sd is reported next
	print("PPV\t%CC\tSens\tSpec")
	print(paste("SD:", splitsize[i]))
	apply(grp3.vs.grp4.winsor[[i]]$results$statmatrix, 2, sd)
}

#run comparisons on winsorized data
grp1.vs.grp4.winsor <- runRatios(grp1.winsor, grp4.winsor, numsplits = 100)

#run statistics
for(i in 1:length(grp1.vs.grp4.winsor)) {
	#for each set of splits (0.5, 0.75, 0.9) report statistics
	#mean is reported first
	print(paste("Mean:", splitsize[i]))
	print("PPV\t%CC\tSens\tSpec")
	apply(grp1.vs.grp4.winsor[[i]]$results$statmatrix, 2, mean)
	#sd is reported next
	print(paste("SD:", splitsize[i]))
	print("PPV\t%CC\tSens\tSpec")
	apply(grp1.vs.grp4.winsor[[i]]$results$statmatrix, 2, sd)
}

#run comparisons on winsorized data
grp1.vs.grp3.winsor <- runRatios(grp1.winsor, grp3.winsor, numsplits = 100)

#run statistics
for(i in 1:length(grp1.vs.grp3.winsor)) {
	#for each set of splits (0.5, 0.75, 0.9) report statistics
	#mean is reported first
	print(paste("Mean:", splitsize[i]))
	print("PPV\t%CC\tSens\tSpec")
	apply(grp1.vs.grp3.winsor[[i]]$results$statmatrix, 2, mean)
	#sd is reported next
	print(paste("SD:", splitsize[i]))
	print("PPV\t%CC\tSens\tSpec")
	apply(grp1.vs.grp3.winsor[[i]]$results$statmatrix, 2, sd)
}
