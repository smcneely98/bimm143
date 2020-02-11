#' ---
#' title: "Class 05: Data exploration and visualization in R"
#' author: "Seth McNeely"
#' date: "2020-01-23"
#' ---

#Section 1
plot(1:10, col="blue", typ="o")

#Section 2
#2A
weight <- read.table("bimm143_05_rstats/weight_chart.txt",
                    header=TRUE)
plot(weight$Age, weight$Weight, col="blue", 
     typ="o", pch=15, cex=1.5, lwd=2, 
     ylim=c(2,10), xlab="Age (months)", 
     ylab="Weight (kg)", main="Weight by Age")

#2B
mouse <- read.table("bimm143_05_rstats/feature_counts.txt",
                          header=TRUE, sep="\t")
par(mar=c(3.1, 11.1, 4.1, 2))
barplot(mouse$Count, horiz=TRUE, 
        names.arg=mouse$Feature, 
        main="Number of Features in the Mouse GCRm38 genome",
        las=1, xlim=c(0,80000))

#2C
x <- c(rnorm(10000),rnorm(10000)+4)
hist(x, breaks=80)

# Section 3
#3A
mf <- read.delim("bimm143_05_rstats/male_female_counts.txt")
par(mar=c(6, 6, 3, 3))
barplot(mf$Count, names.arg=mf$Sample, 
        col=rainbow(nrow(mf)), 
        las=2, ylab="Counts")

barplot(mf$Count, names.arg=mf$Sample, 
        col=c("blue1", "pink"), 
        las=2, ylab="Counts")

#3B
genes <- read.delim("bimm143_05_rstats/up_down_expression.txt", 
                    header=TRUE)
table(genes$State)
palette(c("blue", "grey", "red"))
plot(genes$Condition1, genes$Condition2, 
     col=genes$State, xlab="Condition 1",
     ylab="Condition 2")

#3C
meth <- read.delim("bimm143_05_rstats/expression_methylation.txt")
inds <- meth$expression > 0
dcols <- densCols(meth$gene.meth[inds], meth$expression[inds])
plot(meth$gene.meth[inds], meth$expression[inds], 
     col=dcols, pch=20)

dcols.custom <- densCols(meth$gene.meth[inds], meth$expression[inds],
                         colramp=colorRampPalette(c("blue2",
                                                    "green2",
                                                    "red2",
                                                    "yellow")))
plot(meth$gene.meth[inds], meth$expression[inds], 
     col=dcols.custom, pch=20)
