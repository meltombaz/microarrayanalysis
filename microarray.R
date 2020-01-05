if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.10")
BiocManager::install("affy")
library(affy)
affy.data = ReadAffy()
eset.mas5 = mas5(affy.data)
exprSet.nologs = exprs(eset.mas5)
# List the column (chip) names
colnames(exprSet.nologs)
# Rename the column names if we want
colnames(exprSet.nologs) = c("high10-1", 
                             "high10-2", "high48-1","high48-2",
                             "low10-1",  "low10-2" , "low48-1",
                             "low48-2")
exprSet = log(exprSet.nologs, 2)
write.table(exprSet, file="estrogen_mas5_matrix.txt", quote=F, sep="\t")
# Run the Affy A/P call algorithm on the CEL files we processed above
data.mas5calls = mas5calls(affy.data)
# Get the actual A/P calls
data.mas5calls.calls = exprs(data.mas5calls)
# Print the calls as a matrix
write.table(data.mas5calls.calls, file="estrogen_mas5calls.txt", quote=F, sep="\t")
#Calculating log2 ratios
high10.mean = apply(exprSet[, c("high10-1", "high10-2")], 1, mean)
high48.mean = apply(exprSet[, c("high48-1", "high48-2")], 1, mean)
low10.mean = apply(exprSet[, c("low10-1", "low10-2")], 1, mean)
low48.mean = apply(exprSet[, c("low48-1", "low48-2")], 1, mean)
high48.to.high10 = high48.mean - high10.mean
low48.to.low10 = low48.mean - low10.mean
all.data = cbind(exprSet, high48.mean, high10.mean, low48.mean, low10.mean,
                 high48.to.high10, low48.to.low10)
# Check what data we have here
colnames(all.data)
write.table(all.data, file="Microarray_Analysis_data_1_SOLUTION.txt", quote=F, sep="\t")
exprSet = read.delim("estrogen_mas5_matrix.txt")
# Check how the chips are named
colnames(exprSet)
dataset.1 = exprSet[1, c("high48.1", "high48.2")]
dataset.2 = exprSet[1, c("high10.1", "high10.2")]
# or we can indicate columns by their numbers
dataset.1 = exprSet[1, c(1,2)]
dataset.2 = exprSet[1, c(3,4)]
t.test.gene.1 = t.test(dataset.1, dataset.2, "two.sided")
# Let's see what these data are
dataset.1
dataset.2
# Print just the p-value from the t-test
#install.packages("rt.test")
#library(rt.test)
t.test.gene.1$p.value
high.p.value.all.genes = apply(exprSet, 1, function(x) { t.test(x[1:2], x[3:4]) $p.value } )
low.p.value.all.genes = apply(exprSet, 1, function(x) { t.test(x[5:6], x[7:8]) $p.value } )
# Check the first few brain ones to make sure the first one agrees with our single-gene command
high.p.value.all.genes[1:5]
data.mas5calls.calls = read.delim("estrogen_mas5calls.txt")
# Example for one gene
AP.gene.1 = paste(data.mas5calls.calls[1,], collapse="")
# For all genes
AP = apply(data.mas5calls.calls, 1, paste, collapse="")
genes.present = names(AP[AP != "AAAAAAAA"])
# How many probetset/genes are present?
length(genes.present)
exprSet.present = exprSet[genes.present,]
high.raw.pvals.present = high.p.value.all.genes[genes.present]
low.raw.pvals.present = low.p.value.all.genes[genes.present]
high.fdr.pvals.present = p.adjust(high.raw.pvals.present, method="fdr")
low.fdr.pvals.present = p.adjust(low.raw.pvals.present, method="fdr")
high.fdr.pvals.present.sorted = 
  high.fdr.pvals.present[order(high.fdr.pvals.present)]
low.fdr.pvals.present.sorted = 
  low.fdr.pvals.present[order(low.fdr.pvals.present)]
# Look at the 10 lowest p-values
high.fdr.pvals.present.sorted[1:10]
low.fdr.pvals.present.sorted[1:10]
expression.plus.pvals = cbind(exprSet.present, high.raw.pvals.present, 
                              high.fdr.pvals.present, low.raw.pvals.present, low.fdr.pvals.present)
write.table(expression.plus.pvals, "estrogen_mas5_DE_analysis.txt", sep="\t", quote=F)
high.DE.probesets = names(high.raw.pvals.present[high.raw.pvals.present < 0.01])
low.DE.probesets = names(low.raw.pvals.present[low.raw.pvals.present < 0.01])
# Read the big file we created at the end of last class
all.data = read.delim("Microarray_Analysis_data_1_SOLUTION.txt")
# Get the log2 ratios for the probesets (rows we want)
high.DE.log2.ratios = all.data[high.DE.probesets, 
                                c("high48.to.high10", "low48.to.low10")]
low.DE.log2.ratios = all.data[low.DE.probesets, 
                                c("high48.to.high10", "low48.to.low10")]
write.table(high.DE.log2.ratios, "high.DE.log2.ratios.txt", sep="\t", quote=F)
write.table(low.DE.log2.ratios, "low.DE.log2.ratios.txt", sep="\t", quote=F)

colnames(all.data)
x.data = all.data[, "high48.mean"]
y.data = all.data[, "low48.mean"]
plot(x.data, y.data)
plot(x.data, y.data, main = "log2 expression in high 48h vs low 48 h",
     xlab="high48", ylab="low48", col="blue", cex=0.5)

x1.data = all.data[, "high10.mean"]
y1.data = all.data[, "low10.mean"]
plot(x.data, y.data, main = "log2 expression in high 10h vs low 10 h",
     xlab="high48", ylab="low48", col="blue", cex=0.5)


high = all.data[, "high48.mean"]
low = all.data[, "low48.mean"]
A = (high + low) / 2
M = high - low
plot(A, M)
plot(A, M, main="MA plot of low vs high", pch=19, cex=0.2, col="red")

expression.plus.pvals = read.delim("estrogen_mas5_DE_analysis.txt")
log2.ratios = expression.plus.pvals[, "high48.1"] -  expression.plus.pvals[, "low48.1"]
p.values = expression.plus.pvals[, "high.raw.pvals.present"]
plot(log2.ratios, -log(p.values, 10) )
