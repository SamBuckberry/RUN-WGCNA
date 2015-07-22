# set directory and source the WGCNA functions
setwd(dir="~/R/CTRgroup/bec_zincArray/coExpression/")
source("WGCNA_functions.R")
dat <- read.table("../steveMethod/zincdietarrays/zincDietExprEstimates.csv", sep=",", header=TRUE)
dat[1:6, 1:6]
# Remove genes with < 8 probes
dat <- dat[dat$nProbes > 7, ]

# Clean up for analyses
rownames(dat) <- dat$gene
dat <- dat[ ,-c(1:2, 31)]
dat[1:6, 1:6]

# Determine the soft thresholding power
determineSoftPowerWGCNA(data1=dat, outFile="softPowerPlots.pdf")

# Construct the network
#c1 <- runWGCNA(data1=dat, propGenes=1, softPower=6,
#signedNetwork=TRUE)
# Save the netowrk due to the long computation time
#save(c1, file="wgcna_network.Rda")
load("wgcna_network.Rda")

# Plot the dendrogram with different splits
plotModulesCut(referenceDataset=c1, outFile="clusterModulesPlot.pdf")

# Calculate the eigengenes
e1 <- calculateModuleEigengenes(referenceDataset=c1, split=2)
rownames(e1) <- colnames(dat)
colnames(e1) <- gsub(pattern="ME", replacement="", x=colnames(e1))

write.table(e1, file="moduleEigengenes.csv", quote=FALSE, sep=",",
            col.names=TRUE, row.names=TRUE)

# Get the module genes
mm <- moduleMembership(referenceDataset=c1, MEs=e1, split=2, kME_threshold=0.7)
## Make a matrix from the list, with shorter vectors filled out with ""s
n <- max(sapply(mm, length))
ll <- lapply(mm, function(X) {
        c(as.character(X), rep("", times = n - length(X)))
})
out <- do.call(cbind, ll)       

# Export the file
write.table(out, file="moduleGenes.csv", sep=",", quote=FALSE, col.names=TRUE)



# Get the group information and plot eigengenes
e1$group <- 




e1$group <- factor(substring(rownames(e1), first=1, last=1))
melted <- melt(data=e1, varnames="group")
pdf(file="eigengene_dietGroup_boxplots.pdf", width=12)
qplot(data=melted, y=value, x=group, facets=.~variable, geom=c("boxplot", "point"),
      ylab="Eignegene expression", colour=group)
dev.off()

# Test for differences in eigengene expression 
x <- e1$group
e1$group <- NULL

pValues <- c(rep(x="NA", times=ncol(e1)))

for(i in 1:ncol(e1)){
        y <- e1[ ,i]
        fit <- lm(formula=y ~ x)
        a <- anova(fit)
        pValues[i] <- a[1, 5]
}

## Limma differential expression
f <- factor(x)
design <- model.matrix(~0+f)
colnames(design) <- c("C","M","P")
fit <- lmFit(eset, design)
contrast.matrix <- makeContrasts(M-C, P-M, P-C, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
results <- decideTests(fit2)
pdf("limma_venn.pdf")
vennDiagram(results)
dev.off()


