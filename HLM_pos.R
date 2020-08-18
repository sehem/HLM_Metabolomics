##
## Preparation of parameters
##

## Peak picking parameters are scanned from centWaveOptResults folder in line 47

## Change polarity due to the polarity of the ion source
polarity <- "positive"

## Enter the sequence order of your samples for batch correction using pooled QC samples! If you do
## not use batch correction set as NA
order <- c(4,5,6,7,8,16,17,18,19,20,10,11,12,13,14,1,2,3,9,15,21,22)
## This script requires tidyverse, xcms, CAMERA, ggrepel, Rtnse, and gplots!
## If not installed you can use the following lines to install them:
##install.packages(c("tidyverse", "xcms", "CAMERA", "ggrepel", "Rtnse", "gplots")
##source("https://bioconductor.org/biocLite.R")
##biocLite(c("xcms", "CAMERA"))

## Loading necessary packages
library("tidyverse")
library("xcms")
library("CAMERA")
library("ggrepel")
library("Rtsne")
library("gplots")
library("MASS")

## Load metaboLib.R
source("~/metaboLib_devel.R")


##
## Preprocessing
##


## Create directory for output files
dir.create("evaluateResults", showWarnings = FALSE)

## Load file paths and parameters
files <- list.files(path = getwd(), recursive = TRUE, full.names = TRUE, pattern = ".mzXML")
parameter <- scan("centWaveOptResults/parameter.csv", sep = ",")

## Preprocess files
set <- xcmsSet(files, method = "centWave", peakwidth = c(parameter[1], parameter[2]),
               ppm = parameter[3], snthresh = parameter[4], mzdiff = parameter[5],
               prefilter = c(parameter[6], parameter[7]))
set <- group.density(set, bw = parameter[8])
set <- retcor(set, method = "obiwarp", plottype = "none")
set <- group.density(set, bw = parameter[8])
set <- retcor(set, method = "obiwarp", plottype = "none")
set <- group.density(set, bw = parameter[8])
set <- fillPeaks(set)

## Annotation of isotopes and adducts and export peaklilst to csv
set.CAMERA <- CAMERA::annotate(set, polarity  = polarity)
write.csv(getPeaklist(set.CAMERA), file = "evaluateResults/peaklistCAMERA.csv")


##
## Start of evaluation
##

##load("evaluateResults/evaluate.RData")
pdf("evaluateResults/sampleStatistics.pdf", width = 6, height = 6)
Classes <- set$class
colours_classes <- c("Blank" = "blue", "High" = "red", "Low" = "green",
                     "QC" = "black")

## Extract feature abundances and replace 0 by lowest measured abundance (acc. to Wehrens 2016)
peaklist <- peakTable(set)
rownames(peaklist) <- groupnames(set)
class.levels <- levels(set$class)
matrix.unfiltered <- peaklist[,(8 + length(class.levels)):length(peaklist)]

## Log10 transform for further corrections
matrix.transformed <- log10.matrix(matrix.unfiltered)

## Perform batch correction on data
matrix.normalised <- norm.single.batch(matrix = matrix.transformed, class = set$class, order = order,
                                       type = "pooled", standard = (peaklist$QC == 7),
                                       poly = 5, plot = 50)


##
## Do Statistics
##


## Calculate significance by Volcano plot
significance.table <- t(matrix.normalised)[1:length(set$class[set$class != "QC"]),]
significance.result <- evalANOVA(features = significance.table,
                                 names = rownames(matrix.normalised),
                                 classes = set$class[set$class != "QC"],
                                 p.value = 0.01/nrow(peaklist),
                                 plot = TRUE,
                                 labels = TRUE)

## Filter by Volcano Plot
matrix <- column_to_rownames(filter(rownames_to_column(as.data.frame(matrix.normalised)),
                                    significance.result$Significant == "TRUE"))
matrix <- t(matrix)

## Perform multivariate statistics
if (length(matrix) > 1) {

    ## Plot heatmap with dendrogram
    heatmap.2(x = t(matrix), scale = "row", distfun = dist, margins = c(8,8), trace = "none")

    ## Perform t-SNE analysis and plot results
    tsne <- Rtsne(matrix, dim = 2, perplexity = 5, epoch = 100,
                  verbose = FALSE, pca_center = TRUE, pca_scale = TRUE)
    print(
        ggplot(data = as.data.frame(tsne$Y), aes(x = tsne$Y[,1], y = tsne$Y[,2], col = Classes)) +
        geom_point(size = 3) +
        scale_color_manual(values = colours_classes) +
        ggtitle("t-SNE") +
        xlab("X") +
        ylab("Y") +
        theme_bw()
    )
    
    ## Perform PCA and plot results
    pca <- evalPCA(matrix, plot = TRUE, annotate.scores = TRUE, annotate.loadings = 15,
                   classes = set$class, scale = FALSE, center = TRUE)
    
    ## Perform PC-LDA
    n.pc <- length(which(pca$sdev^2 > mean(pca$sdev^2)))
    if(n.pc < 3) {

        n.pc <- 3

    }
    pc.ldam <- lda(as.matrix(pca$x[,1:n.pc]), set$class)
    pc.lda.loadings <- pca$rotation[,1:n.pc] %*% pc.ldam$scaling

    ## Calculate prediction accuracy
    cv.ldam <- lda(as.matrix(pca$x)[,1:n.pc], set$class, CV = TRUE)
    pc.lda.accuracy <- sum(diag(prop.table(table(set$class, cv.ldam$class))))
    pc.lda.group.accuracy <- diag(prop.table(table(set$class, cv.ldam$class), 1))

    ## Predict model for test data
    p <- predict(pc.ldam, as.data.frame(pca$x[,1:n.pc]))

    ## Plot PC-DFA Scores without labels
    print(
        ggplot(data = as.data.frame(p$x), aes(x = LD1, y = LD2)) +
        geom_point(size = 3, aes(col = Classes)) +
        scale_color_manual(values = colours_classes) +
        geom_hline(yintercept = 0) +
        geom_vline(xintercept = 0) +
        ggtitle(paste0("PCs used: ", n.pc, ", ", "Accuracy: ",
                       round(pc.lda.accuracy*100, 0), "%")) +
        theme_bw()
    )
    
    ## Plot PC-DFA Loadings
    print(
        ggplot(data = as.data.frame(pc.lda.loadings), aes(x = LD1, y = LD2)) +
        geom_point(size = 3, aes(col = "red"), show.legend = FALSE) +
        geom_hline(yintercept = 0) +
        geom_vline(xintercept = 0) +
        geom_text_repel(point.padding = 0.2, data = as.data.frame(pc.lda.loadings),
                        aes(label = rownames(pc.lda.loadings))) +
        ggtitle(paste0("PCs used: ", n.pc, ", ", "Accuracy: ",
                       round(pc.lda.accuracy*100, 0), "%")) +
        theme_bw()
    )

}
dev.off()

if (length(matrix) > 1) {

    ## Plot Feature EICs and save to hard drive
    groupid <- significance.result$feature[significance.result$Significant == "TRUE"]
    eic <- getEIC(set, groupidx = groupid, rt = "corrected")
    pdf("evaluateResults/Chromatograms.pdf", width = 12, height = 8)
    plot(eic, set, groupidx = groupnames(eic))
    dev.off()
  
}

save.image(file = "evaluateResults/evaluate.RData")

## Summarise significant features for identification
significant.features <- dplyr::select(significance.result, feature, pvalue)
significant.features$mz <- rbind(peaklist[which(peaklist$QC == 7),], peaklist[which(peaklist$QC != 7),])$mz
significant.features$rt <- rbind(peaklist[which(peaklist$QC == 7),], peaklist[which(peaklist$QC != 7),])$rt
significant.features$adduct <- rbind(getPeaklist(set.CAMERA)[which(getPeaklist(set.CAMERA)$QC == 7),],
                                     getPeaklist(set.CAMERA)[which(getPeaklist(set.CAMERA)$QC != 7),])$adduct
significant.features$isotopes <- rbind(getPeaklist(set.CAMERA)[which(getPeaklist(set.CAMERA)$QC == 7),],
                                     getPeaklist(set.CAMERA)[which(getPeaklist(set.CAMERA)$QC != 7),])$isotopes
significant.features$pcgroup <- rbind(getPeaklist(set.CAMERA)[which(getPeaklist(set.CAMERA)$QC == 7),],
                                      getPeaklist(set.CAMERA)[which(getPeaklist(set.CAMERA)$QC != 7),])$pcgroup
significant.features <- significant.features[which(significance.result$Significant == TRUE),]
significant.features$abs.loading <- sqrt(abs(pc.lda.loadings[,1])^2 + abs(pc.lda.loadings[,2])^2)
significant.features <- arrange(significant.features, desc(abs.loading))
write.csv(significant.features, file = "evaluateResults/significant_features.csv")
