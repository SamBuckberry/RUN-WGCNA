


library(WGCNA)
library(flashClust)
library(genefilter)
#disableWGCNAThreads()
options(stringsAsFactors = FALSE)

#' @param data1 data.frame of normalised gene expression values with samples as columns and rows as genes.
#' @param outFile character. A file path for the soft power plots as .pdf
#' @param propGenes A numeric value indicating the proportion of most variable genes to be retained in the analysis between 0 and 1.
#' @return This function does not return any value. A pdf file is generated with two plots.

determineSoftPowerWGCNA <- function(data1, outFile, propGenes=1){
        
        options(stringsAsFactors = FALSE)
        
        # Remove bad genes (missing values, ...)
        nGenesInput <- dim(data1)[1]
        data1 <- data1[goodSamplesGenes(datExpr=t(data1))$goodGenes, ]
        nGenesGood <- dim(data1)[1]
        nGenesRemoved <- nGenesInput - nGenesGood
        message(paste(nGenesRemoved, " genes filtered from dataset"))
        propGenes <- round(propGenes * dim(data1)[1])
        
        # Filter genes based in variance
        keepGenesExpr1 <- rank(-rowVars(data1)) <= propGenes
        data1 <- data1[keepGenesExpr1, ]
        genesRetained <- dim(data1)[1]
        message(paste(genesRetained, " genes retained in dataset"))
        
        # plot powers to work out soft-power 
        # Choose a set of soft-thresholding powers
        powers <- c(c(1:10), seq(from = 12, to=20, by=2))
        
        # Call the network topology analysis function
        sft <- pickSoftThreshold(t(data1), powerVector=powers, verbose=5)
        
        # Plot the results:
        png(outFile, width=180, height=125, units="mm", res=150)
        par(mfrow = c(1,2))
        cex1 <- 0.9
        
        # Scale-free topology fit index as a function of the soft-thresholding power
        plot(sft$fitIndices[,1],
             -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
             xlab="Soft Threshold (power)",
             ylab="Scale Free Topology Model Fit,signed R^2",
             type="n", main = paste("Scale independence"))
        
        text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
             labels=powers, cex=cex1, col="red")
        
        # this line corresponds to using an R^2 cut-off of h
        abline(h=0.80, col="red")
        
        # Mean connectivity as a function of the soft-thresholding power
        plot(sft$fitIndices[,1], sft$fitIndices[,5],
             xlab="Soft Threshold (power)",
             ylab="Mean Connectivity", type="n",
             main = paste("Mean connectivity"))
        
        text(sft$fitIndices[,1], sft$fitIndices[,5], 
             labels=powers, cex=cex1, col="red")
        
        dev.off()
}

#' @param data1 data.frame of normalised gene expression values with samples as columns and rows as genes.
#' @param propGenes A numeric value indicating the proportion of most variable genes to be retained in the analysis between 0 and 1.
#' @param softPower Integer. The soft thresholding power used to construct the WGCNA network. This can be determined using the determineSoftPowerWGCNA function.
#' @param signedNetwork Logical. Keep track of the sign of the correlation in network construction?
#' @return A list of 3. data1 is a data.frame of the input data after filtering. geneTreeA1 is the gene tree constructed by WGCNA. dissTOMA1 is the distance matirx.
#' @description A wrapper function for constructing a WGCNA network 

runWGCNA <- function(data1, propGenes, softPower=10, signedNetwork=TRUE){
        
        options(stringsAsFactors = FALSE)
        
        type <- ifelse(test=signedNetwork==TRUE, yes="signed", no="unsigned")
        
        # Remove bad genes (missing values, ...)
        nGenesInput <- dim(data1)[1]
        data1 <- data1[goodSamplesGenes(datExpr=t(data1))$goodGenes, ]
        nGenesGood <- dim(data1)[1]
        nGenesRemoved <- nGenesInput - nGenesGood
        message(paste(nGenesRemoved, " genes filtered from dataset"))
        propGenes <- round(propGenes * dim(data1)[1])
        
        # Filter genes based in variance
        keepGenesExpr1 <- rank(-rowVars(data1)) <= propGenes
        data1 <- data1[keepGenesExpr1, ]
        genesRetained <- dim(data1)[1]
        message(paste(genesRetained, " genes retained in dataset"))
        
        # Run WGCNA on the datasets
        datExprA1g <- data1
        adjacencyA1 <- adjacency(t(datExprA1g),power=softPower,type=type)
        diag(adjacencyA1) <- 0
        dissTOMA1 <- 1-TOMsimilarity(adjacencyA1, TOMType=type)
        
        geneTreeA1 <- flashClust(as.dist(dissTOMA1), method="average")
        
        # Return the relevant objects
        return(list(data1 = data1, geneTreeA1 = geneTreeA1,
                    dissTOMA1 = dissTOMA1))
}


plotModulesCut <- function(referenceDataset, outFile, minClusterSize=30){
        
        # Seperate list objects
        mColorh=NULL
        for (ds in 0:3){
                tree <- cutreeHybrid(dendro = referenceDataset[["geneTreeA1"]], 
                                     pamStage=FALSE, minClusterSize = minClusterSize,
                                     cutHeight = 0.99, deepSplit = ds,
                                     distM = referenceDataset[["dissTOMA1"]])
                
                mColorh <- cbind(mColorh, labels2colors(tree$labels))
        }
        
        pdf(outFile, height=10, width=25)
        plotDendroAndColors(referenceDataset[["geneTreeA1"]], mColorh,
                            paste("dpSplt =",0:3), main = "",
                            dendroLabels=FALSE)
        dev.off()
}

calculateModuleEigengenes <- function(referenceDataset, split, minClusterSize=30){
        # Seperate list objects 
        mColorh=NULL
        for (ds in 0:3){
                tree <- cutreeHybrid(dendro = referenceDataset[["geneTreeA1"]],
                                     pamStage=FALSE, minClusterSize = minClusterSize,
                                     cutHeight = 0.99, deepSplit = ds,
                                     distM = referenceDataset[["dissTOMA1"]])
                
                mColorh <- cbind(mColorh, labels2colors(tree$labels))
        }
        
        modulesA1 <- mColorh[ ,split] # (Chosen based on plot)
        
        PCs <- moduleEigengenes(t(referenceDataset$data1), colors=modulesA1)
        ME <- PCs$eigengenes
        rownames(ME) <- colnames(referenceDataset$data1)
        ME
}

geneModuleColors <- function(referenceDataset, split, minClusterSize=30){
        # Seperate list objects 
        mColorh=NULL
        for (ds in 0:3){
                tree <- cutreeHybrid(dendro = referenceDataset[["geneTreeA1"]],
                                     pamStage=FALSE, minClusterSize = minClusterSize,
                                     cutHeight = 0.99, deepSplit = ds,
                                     distM = referenceDataset[["dissTOMA1"]])
                
                mColorh <- cbind(mColorh, labels2colors(tree$labels))
        }
        
        modules <- mColorh[ ,split] # (Chosen based on plot)
        modules
}


moduleHubGenes <- function(referenceDataset, MEs, nGenes, split=1){
        
        message("ranking genes")
        kMEs <- signedKME(datExpr=t(referenceDataset$data1), datME=MEs)
        
        # rank the genes for each module on kMEs
        rankGenes <- function(x){
                kMErank <- rank(-kMEs[ ,x])
                genes <- rownames(kMEs)
                genes <- genes[order(kMErank)]
                genes[1:nGenes]
        }
        
        topGenes <- lapply(1:ncol(kMEs), rankGenes)
        
        # Get the top results in a data.frame
        topGenes <- do.call(cbind, topGenes)
        colnames(topGenes) <- substr(colnames(kMEs), start=4, stop=30)
        return(topGenes)
}

moduleMembership <- function(referenceDataset, MEs, split=1, kME_threshold=0.7){
        message("calculating kME's")
        
        kMEs <- signedKME(datExpr=t(referenceDataset$data1), datME=MEs)
        
        # Get gene names with kME > kME_threshold
        modGenes <- function(x){
                rownames(kMEs[abs(kMEs[ ,x]) > kME_threshold, ])
        }
        geneLists <- lapply(X=1:ncol(kMEs), FUN=modGenes)
        names(geneLists) <- substr(colnames(kMEs), start=4, stop=20)
        return(geneLists)
}

presStatsWGCNA <- function(referenceDataset, data2, colors, split=2, type="signed",
                           nPermutations=100, minClusterSize=30,
                           greyName="grey", qVal=FALSE){
        
        # Remove bad genes (missing values, ...)
        data2 <- data2[goodSamplesGenes(datExpr=t(data2))$goodGenes, ]
        
        #         # Seperate list objects 
        #         mColorh=NULL
        #         for (ds in 0:3){
        #                 tree <- cutreeHybrid(dendro = referenceDataset[["geneTreeA1"]],
        #                                      pamStage=FALSE, minClusterSize = minClusterSize,
        #                                      cutHeight = 0.99, deepSplit = ds,
        #                                      distM = referenceDataset[["dissTOMA1"]])
        #                 
        #                 mColorh <- cbind(mColorh, labels2colors(tree$labels))
        #         }
        #         
        #         modulesA1 <- mColorh[ ,split] # (Chosen based on plot)
        
        geneModules <- cbind(rownames(referenceDataset$data1), colors)
        
        ### Quantify module preservation
        
        data1 <- referenceDataset[["data1"]]
        
        multiExpr <- list(A1=list(data=t(data1)),
                          A2=list(data=t(data2)))
        
        multiColor <- list(A1 = colors)
        
        mp <- modulePreservation(multiData=multiExpr, multiColor=multiColor,
                                 referenceNetworks=1, verbose=3,
                                 calculateQvalue=qVal, 
                                 networkType=type,
                                 nPermutations=nPermutations,
                                 maxGoldModuleSize=100,
                                 greyName=greyName
                                 #maxModuleSize=400
        )
        
        stats <- mp$preservation$Z$ref.A1$inColumnsAlsoPresentIn.A2
        stats <- stats[order(-stats[,2]), ]
        print(stats[ ,1:3])
        return(list(mp = mp, geneModules = geneModules))
}

moduleTraitCorrelations <- function(MEs, trait){
        
        # Order traits based on clustering tree
        d <- dist(x=t(trait))
        traitTree <- hclust(d, method="a")
        plot(traitTree, xlab="", ylab="", main="", sub="")
        traitOrder <- traitTree$labels[traitTree$order]
        trait <- trait[ ,match(traitOrder, colnames(trait))]
        
        # Order modules based on clustering tree
        d <- dist(x=t(MEs))
        modTree <- hclust(d, method="a")
        plot(modTree, xlab="", ylab="", main="", sub="")
        modOrder <- modTree$labels[modTree$order]
        MEs <- MEs[ ,match(modOrder, colnames(MEs))]
        
        # Calculate module-trait correlations
        moduleTraitCor <- cor(x=trait, y=MEs, use="pairwise.complete")
        moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(MEs))
        
        return(list(moduleTraitCor=moduleTraitCor,
                    moduleTraitPvalue=moduleTraitPvalue))
}

plotModuleTraitCorrelations <- function(eigenTraitCor, pThresh, outFile){
        
        moduleTraitCor <- eigenTraitCor[["moduleTraitCor"]]
        moduleTraitPvalue <- eigenTraitCor[["moduleTraitPvalue"]]
        
        
        ## *** Experimental: organse traits by clustering correlation values
        d <- dist(x=moduleTraitCor)
        traitTree <- hclust(d, method="a")
        
        d <- dist(x=t(moduleTraitCor))
        moduleTree <- hclust(d, method="a")
        
        moduleTraitCor <- moduleTraitCor[traitTree$order, moduleTree$order]
        moduleTraitPvalue <- moduleTraitPvalue[traitTree$order, moduleTree$order]
        
        # Put the p-values in a text matrix for plotting
        textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                            signif(moduleTraitPvalue, 1), ")", sep = "")
        
        # Change text matrix to only contain significant p-values
        textMatrix <- signif(moduleTraitPvalue, 1)
        textMatrix[textMatrix >= pThresh] <- ""
        dim(textMatrix) <- dim(moduleTraitCor)
        
        # Reorder the the plot by dendrogram of correlation values
        
        # Display the correlation values within a heatmap plot
        myColRamp <- colorRampPalette(colors=c("#1F8A70", "#FFFFFF", "#FD7400"))
        
        # Graphical parameters
        par(mar = c(6, 8.5, 3, 3), cex.lab=0.5)
        
        modSymbols <- convertColorToLabel(substr(x=colnames(moduleTraitCor), 3, 30), prefix="M")
        
        pdf(outFile)
        labeledHeatmap(Matrix=moduleTraitCor,
                       xColorLabels=TRUE,
                       xLabels=colnames(moduleTraitCor),
                       yLabels=rownames(moduleTraitCor),
                       xSymbols=modSymbols,
                       textMatrix=textMatrix,
                       colors=myColRamp(50),
                       main=paste("Gene Module - Clinical Variable Correlations"),
                       cex.text=0.5,
                       cex.lab.y=0.5, 
                       cex.lab.x=0.7)
        
        dev.off()
        
}

# Convert module colors to module names (so thing look more professional)

convertColorToLabel <- function(colors, n=100, prefix="mod"){
        
        # Get the WGCNA colors
        cols <- standardColors(n)
        labels <- 1:n
        labels <- paste(prefix, labels, sep="")
        colToLabel <- data.frame(colors = cols, labels = labels)
        
        # Add grey as label 0
        grey <- c("grey", paste(prefix, "0", sep=""))
        colToLabel <- rbind(grey, colToLabel)
        
        # Function to return corresponding label for a given color
        convertFunc <- function(x){
                colToLabel[colToLabel$colors == x, ]$labels
        }
        
        # Apply the function to a list of colors
        unlist(lapply(colors, convertFunc))
}

