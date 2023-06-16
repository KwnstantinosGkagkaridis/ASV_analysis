library(dada2); packageVersion("dada2")

mythreads <- 4

path <- "./MiSeq_SOP" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Plot the per-base qualities
pdf(file = "seq_quality_scores.pdf", height = 4, width = 8, onefile = T)
for(i in 1:length(fnFs)){
  print(i)
  plot <- plotQualityProfile(c(fnFs[i],fnRs[i]))
  print(plot)
}
dev.off()

#Inspect read quality profiles
#plotQualityProfile(fnFs[1:2])
#plotQualityProfile(fnRs[1:2])

# Place filtered files in filtered/ subdirectory
filtFs <- file.path("filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, 
                     filtFs, 
                     fnRs, 
                     filtRs, 
                     truncLen = c(240,160),
                     maxN = 0, maxEE = c(2,2), 
                     truncQ = 2, 
                     rm.phix = TRUE,
                     compress = TRUE, 
                     multithread = mythreads) 
head(out)

#Learn the Error Rates using a machine learning algorithm
errF <- learnErrors(filtFs, multithread = mythreads)
errR <- learnErrors(filtRs, multithread = mythreads)
#plotErrors(errF, nominalQ=TRUE)

# dereplication of each one of the red pairs to unique sequences (collapsing of the identical sequences for each pair per sample)
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

## sample composition inference after read correction
dadaFs <- dada(derepFs, err=errF, multithread = mythreads)
dadaRs <- dada(derepRs, err=errR, multithread = mythreads)
#dadaFs[[1]]

# merge read pairs retaining the per sequence sample information
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)

# Inspect the merger data.frame from the first sample
#head(mergers[[1]])

#### construct the sequence table, remove the chimeras, and create a summary ----
# construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

#remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)


# Inspect distribution of sequence lengths
#table(nchar(getSequences(seqtab)))

# record the portion of good sequences out of the total prior the chimea removal
sum(seqtab.nochim)/sum(seqtab)

# track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

write.table(track, file = "readQC.txt", col.names = NA, sep = "\t", quote = FALSE)

#Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "./tax/silva_nr_v132_train_set.fa.gz", multithread = mythreads)

#taxa <- addSpecies(taxa, "./tax/silva_species_assignment_v132.fa.gz")














##### Perform the analysis with the transformed abundance matrix and agglomerate to genus ----
library(devtools)
#devtools::install_github("microEcology/pime")
library(pime)
library(phyloseq)
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")


samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender <- substr(subject,1,1)
subject <- substr(subject,2,999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out


ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps
ps_orig <- ps


#### write/read the final master phyloseq object ----
## write the phyloseq object on the disk with the following command
saveRDS(ps_orig, "my_ps_orig.rds")
## test and read from that file with the following
ps_orig <- readRDS("my_ps_orig.rds")

## If you inspect the taxonomy you will see plenty of NA values at the lower taxonomic ranks.
## These can be replaced with their higher rank annotations in order to eventually add this information in the plots 
## prepare a phyloseq object with taxonomy of genus level at best
## to do so you can implement the following for loop
ps_fin_ASV_best_tax <- ps_orig # assign the late phyloseq object to a new one
# run the loop
for(i in 2:ncol(tax_table(ps_fin_ASV_best_tax))){
  for(j in 1:nrow(tax_table(ps_fin_ASV_best_tax))){
    if(is.na(tax_table(ps_fin_ASV_best_tax)[j,i])){
      tax_table(ps_fin_ASV_best_tax)[j,i] <- tax_table(ps_fin_ASV_best_tax)[j,i-1]
    } 
  }
}

saveRDS(ps_fin_ASV_best_tax, file = "ps_fin_ASV_best_tax.RDS") # save the new phyloseq object
ps_fin_ASV_best_tax <- readRDS("ps_fin_ASV_best_tax.RDS") # read it and assign in a new variable


## depending on your question you might want to work at a genus level (genus level annotations are the most reliable ones for the 16S rRNA marker gene and its amplicons)
# prepare the genus level phyloseq object
ps_fin_genus <- tax_glom(physeq = ps_fin_ASV_best_tax,taxrank = "Genus") # agglomerate sequences of ASVs belonging in the same genus level
saveRDS(ps_fin_genus, file = "ps_fin_genus.RDS")
ps_fin_genus <- readRDS("ps_fin_genus.RDS")



#### write the master tables in a format that can be parsed with GUI software like Excel ----
ps_fin_ASV_ra <- transform_sample_counts(ps_fin_ASV_best_tax, function(x) 100 * x / sum(x) ) # convert ASV counts to percentages
#my_master_df_tax <- merge(data.frame(tax_table(ps_fin_ASV_ra)), data.frame(t(otu_table(ps_fin_ASV_ra))), by = "row.names") # combine taxonomy information with the relative abundance table
#mysmpldta <- data.frame(sample_data(ps_fin_ASV_ra))
#mysmpldta <- mysmpldta[order(mysmpldta$Day),] # reorder the metadata according to the day
#my_master_df_ra_tax <- my_master_df_tax[,c(colnames(my_master_df_tax)[1:8],row.names(mysmpldta))] # do the same for the data
#write.table(my_master_df_ra_tax, file = "000_master_table_rel_abund.txt", sep = "\t", quote = F, row.names = F) # write the tables
#write.table(mysmpldta, file = "000_master_saple_data.txt", sep = "\t", quote = F, col.names = NA) # write the tables

# write also the sequences
#write.table(data.frame(refseq(ps_fin_ASV_ra)), file = "000_master_seq.fasta", col.names = NA, quote = F, sep = "\t")




### prep the alpha diversity indices for setting up the plot xlim and/or ylim values ----
adiv <- plot_richness(ps_fin_ASV_best_tax, measures=c("Shannon","InvSimpson","Fisher","Observed","ACE"))
## calculate Good's coverage estimate as well
library(entropart)
good <- MetaCommunity(t(ps_fin_ASV_best_tax@otu_table))$SampleCoverage.communities # obtain the Good's coverage estimator
good_tbl <- data.frame(samples = names(good), variable = rep("coverage",length(good)), value=good)

alpha_long <- rbind(adiv$data[,c("samples", "variable", "value")], good_tbl)

# convert table formats from long to wide with a series of commands from the reshape 2 package and merge the various tables 
library(reshape2)
alpha_wide <- dcast(alpha_long, samples ~ variable, value.var="value") # convert the format
row.names(alpha_wide) <- alpha_wide$samples # assign the row.names
alpha_wide <- alpha_wide[,c("coverage","Shannon","InvSimpson","Fisher","Observed","ACE")] # select the alpha diversity indices only
colnames(alpha_wide) <- c("coverage","Shannon","Inv. Simpson","Fisher's α", "Observed", "ACE") # change the column names to add the proper Fisher's α

alpha_wide_fact <- merge(alpha_wide, data.frame(sample_data(ps_fin_ASV_best_tax)), by = "row.names") # merge the sample data with the diversity indices
row.names(alpha_wide_fact) <- alpha_wide_fact$Row.names # the row.names were moved to an integral column... add the sample names as row names
alpha_wide_fact <- alpha_wide_fact[-which(colnames(alpha_wide_fact)%in%"Row.names")] # remove the integral column to pretify this

alpha_wide_fact_final_collective <- alpha_wide_fact # assign the table to a final variable for downstream use


#### combine the sequencing QC, the design2, and the alpha diversity results ----
myreadqc <-read.table("readQC.txt", header = T, row.names = 1, sep = "\t")

myrqcdes <-merge(myreadqc, alpha_wide_fact_final_collective, by = "row.names")
row.names(myrqcdes) <- myrqcdes$Row.names

write.table(myrqcdes, file = "readQC_design_div.txt", quote = F, sep = "\t", row.names = F)



myrqcdes$When <- factor(myrqcdes$When, levels = c("Early","Late"))
boxplot(myrqcdes$Shannon~myrqcdes$When)





#### perform anova or equivalent for the alpha diversity indices with a for loop ----

mytestvars <- colnames(myrqcdes)[9:13] # make a collection of the names of variables to be tested

mytestfact <- "When" # assign the name of the independent variable to a variable
myheight <- 2.2 # this will be the height of the plot... it is useful to set it in the case of multiple plots in order to use this as a master switch

library("agricolae") # load the package with the anova-like commands of interest

cairo_pdf("alpha_div_plots.pdf", height = myheight, width = 10)
mystatsout <- list()
par(mfrow = c(1,5))
par(xpd = T)

for(mytestvar in mytestvars){
  
  # create the aov matrix
  myaovmatrix <- myrqcdes
  
  
  
  # run a shapiro wilk test to select parametric or non parametric analysis
  shap_out <- shapiro.test(myaovmatrix[,mytestvar])
  mystatsout[[paste(mytestvar, sep = " // ")]][["shap"]] <- shap_out
  
  # run the parametric or non-parametric analysis according to the shapiro.test results
  if(shap_out$p.value < 0.05){
    # non-parametric
    mykrusk <- kruskal(myaovmatrix[,mytestvar], myaovmatrix[,mytestfact], group = T, p.adj = "BH")
    
    mystatsout[[paste(mytestvar, sep = " // ")]][["krusk"]] <- mykrusk
    # prepare the barplot
    mytestvarord <- levels(factor(myaovmatrix[,mytestfact]))
    par(mar = c(4,5,4,2))
    barerrplt <- bar.err(mykrusk$means[mytestvarord[length(mytestvarord):1],], variation="SD", xlim=c(0,max(alpha_wide_fact_final_collective[,mytestvar], na.rm = T)),horiz=TRUE, bar=T, col="grey60", space=1, main= paste(mytestvar), las=1)
    
    # delete the significance group letters in the case that the anova was not significant
    par(xpd = T)
    if(mykrusk$statistics$p.chisq <= 0.05){
      text(x = mykrusk$means[mytestvarord[length(mytestvarord):1],1] + mykrusk$means[mytestvarord[length(mytestvarord):1],3], barerrplt$x,labels = mykrusk$groups[mytestvarord[length(mytestvarord):1],2], pos = 4, font = 3)
      par(xpd = F)
    }
  } else{
    # perform the parametric
    
    # select the alphadiv matrix and design rows
    myform <- as.formula(paste("`",mytestvar,"` ~ ",mytestfact, sep = ""))
    mymod <- aov(myform, data = myaovmatrix)
    mysumaov <- summary(mymod)
    
    mystatsout[[paste(mytestvar, sep = " // ")]][["ANOVA"]] <- mysumaov
    
    # order the matrices etc
    mytestvarord <- levels(factor(myaovmatrix[,mytestfact]))
    
    # run the Tukey test
    myHSDtest <- HSD.test(mymod, mytestfact, group=T)
    
    mystatsout[[paste(mytestvar, sep = " // ")]][["HSD test"]] <- myHSDtest
    
    # prepare the barplot
    par(mar = c(4,5,4,2))
    barerrplt <- bar.err(myHSDtest$means[mytestvarord[length(mytestvarord):1],], variation="SD", xlim=c(0,max(alpha_wide_fact_final_collective[,mytestvar], na.rm = T)),horiz=TRUE, bar=T, col="grey60", space=1, main= paste(mytestvar), las=1)
    
    # delete the significance group letters in the case that the anova was not significant
    par(xpd = T)
    if(mysumaov[[1]]$`Pr(>F)`[1] <= 0.05){
      text(x = myHSDtest$means[mytestvarord[length(mytestvarord):1],1] + myHSDtest$means[mytestvarord[length(mytestvarord):1],2], barerrplt$x,labels = myHSDtest$groups[mytestvarord[length(mytestvarord):1],2], pos = 4)
      par(xpd = F)
    }
  }
  
}

dev.off()
capture.output(mystatsout,file = "alpha_div_plots.txt")
















# b-diversity analysis

#### prepare a heatmap ----

# prepare a phyloseq object with the names being a combination of the phylum, genus and ASV levels
ps_fin_ASV_ra_heatmap <- ps_fin_ASV_ra # assign the relative abundance table to a fresh one
myphylum <- data.frame(tax_table(ps_fin_ASV_ra_heatmap))$Phylum # assign the ASV phyla to a variable... it is necessary to define the data.frame object class because the table is otherwise perceived as a vector
mygenus <- data.frame(tax_table(ps_fin_ASV_ra_heatmap))$Genus # do the same for the genus
myASV <- taxa_names(ps_fin_ASV_ra_heatmap) # assign also the ASV names to a variable

# combine all three into a final variable
mycmbndnm <- paste(myphylum,mygenus,myASV,sep = ":")
# save these in a table which is indexed by the ASV names
mycmbndnmtbl <- data.frame(mycmbndnm, row.names = myASV)

# replace the names of the ASVs with the generated lineage information
taxa_names(ps_fin_ASV_ra_heatmap) <- mycmbndnmtbl[taxa_names(ps_fin_ASV_ra_heatmap),1] # the "taxa_names(ps_fin_ASV_ra_heatmap)" in combination with the row.names of the table is used to assure that the right order of names is used. 




# select the x most abundant taxa in order to facilitate the plotting function
myXtaxNum <- 30
# obtain the mean values of the raxon relative abundances
mymeanras <- colMeans(data.frame(otu_table(ps_fin_ASV_ra_heatmap), check.names = F))
# get the taxon descending orders according to the taxa
mydescendingnames <- names(mymeanras[order(mymeanras, decreasing = T)])

# extract the top taxa
ps_fin_ASV_ra_heatmap_to30 <- prune_taxa(mydescendingnames[1:myXtaxNum],ps_fin_ASV_ra_heatmap)
# saanity check... remove empty samples in case any were generated after the removal of taxa # in this case it will not make a difference
ps_fin_ASV_ra_heatmap_to30 <- prune_samples(sample_sums(ps_fin_ASV_ra_heatmap_to30) > 0,ps_fin_ASV_ra_heatmap_to30)



# extract the microbial relative abundance table for plotting
mymicrobedata <- t(data.frame(otu_table(ps_fin_ASV_ra_heatmap_to30), check.names = F)) # extract the OTU data and transpose



## prepare a plot with two sample clusters
mycolors <- colorRampPalette(c("grey90","lightsteelblue4","darkorange2","brown","brown"))(n = 119) # set your color arrangement

library(pheatmap)

cairo_pdf(file = "pheatmap.pdf", width = 9, height = 7)
pheatmap(mymicrobedata, color = mycolors, cluster_rows = T, cutree_cols = 2, cutree_rows = 5, border_color = "grey30",cellwidth = 10, cellheight = 10)
dev.off()






#### perform a non-metric multidimentional scaling analysis ----

### NMDS and PERMANOVA hypothesis testing ----

library(plyr)

mynmds <- ordinate(ps_fin_ASV_ra_heatmap, method="NMDS", distance="bray") # run NMDS
library(vegan) # load the vegan library; necessary for extracting the site and species scores below
mynmdssit <- data.frame(scores(mynmds, display = "sites")) # retrieve the site scores
mynmdsspe <- data.frame(scores(mynmds, display = "species")) # retrieve the species scores

mymicrmatre <- ps_fin_ASV_ra_heatmap@otu_table # pass the OTU/ASV table in a new variable
domOTUnums <- 15 # select the number of dominant OTUs to plot 
mymicrmatredomOTUs <- names(sort(colMeans(mymicrmatre), decreasing = T))[1:domOTUnums] # extract the names of the selected dominant OTUs/ASVs
minradomOTUs <- round(colMeans(mymicrmatre)[mymicrmatredomOTUs[domOTUnums]],1) # calculate the minimum mean relative abundance of the selected OTUs/ASVs
maxradomOTUs <- round(colMeans(mymicrmatre)[1],1) # calculate the maximum mean relative abundance of the selected OTUs/ASVs
mynmdssit$fact <- factor(data.frame(sample_data(ps_fin_ASV_ra_heatmap))$When) # define the independent variables
meannms1 <- plyr::mapvalues(factor(mynmdssit$fact), from=levels(factor(mynmdssit$fact)), to=1:length(levels(factor(mynmdssit$fact)))) # convert the variables into numeric values in order to use for chossing colors later on
mynmdssitmeanssdpre <- aggregate(. ~ fact, mynmdssit[,1:3], function(x) c(mean = mean(x), sd = sd(x))) # calculate the means and standard deviations of the sample scores according to the post weaning period
mynmdssitmeanssd <- data.frame(mynmdssitmeanssdpre$fact, mynmdssitmeanssdpre$NMDS1, mynmdssitmeanssdpre$NMDS2) # combine these in a data.frame


## Perform the PERMANOVA 
mytbl <- data.frame(otu_table(ps_fin_ASV_ra_heatmap), stringsAsFactors = F) # extract the data table to be used with the permanova test
mydesign <- data.frame(sample_data(ps_fin_ASV_ra_heatmap)) # extract also the design information
library(vegan) # load the vegan library commands
mypermanova <- adonis2(mytbl ~ When, mydesign, permutations = 999, method = "bray",
                       sqrt.dist = FALSE, add = FALSE, by = "terms",
                       parallel = 4) # perform the PERMANOVA analysis 
write.table(data.frame(mypermanova, check.names = F), file = "PERMANOVA.txt", quote = F, col.names = NA) # save the permanova table


cairo_pdf("NMDS.pdf", height = 6, width = 7) # initiate a graphics device
myplotcols <- RColorBrewer::brewer.pal(n = 3, name = 'RdBu')[c(1,3)] # choose some colors out of the brewer palette
plot(mynmdssit[,1:2], bg = myplotcols[meannms1], frame = F, cex = 0, pch = 21, xlim = c(min(mynmdssit[,1]),max(1.3*mynmdssit[,1])), xlab = "NMDS1", ylab = "NMDS2") # create the scatter plot
vegan::ordiellipse(as.matrix(mynmdssit[,1:2]), factor(mynmdssit$fact), kind = "ehull", lty = 4, lwd=1) # use ellipses to visualize samples of the same 
points(mynmdssit[,1:2], bg = myplotcols[meannms1], cex = 1.5, pch = 21) # plot the sample points
title(main = bquote(atop("PERMANOVA: comm ~ When,  R"^2~.(round(100*mypermanova$R2[1],1))*"%, P"~.(round(mypermanova$`Pr(>F)`[1],3))))) # add the permanova information on top of the plot
par(adj = 1) # prepare for the addition (by moving the text to the far right by tweaking the adj global parameter) of the minimum and maximum mean percentage values of the top taxa to be plotted in this biplot
title(sub = paste(" (mean RA of presented ASVs: ",minradomOTUs,"-",maxradomOTUs,"%)", sep = "")) # add the relative abundance info
par(adj = .5) # return the adj parameter to the default value
arrows(0,0,mynmdsspe[mymicrmatredomOTUs,1] , mynmdsspe[mymicrmatredomOTUs,2], angle = 25, length = 0.15, col = rgb(153,153,153, max = 255, alpha = 100)) # add the arrows showing the microbial gradients
plotrix::thigmophobe.labels(1.2*mynmdsspe[mymicrmatredomOTUs,1], 1.2*mynmdsspe[mymicrmatredomOTUs,2], labels = mymicrmatredomOTUs, cex = .6, font = 2, col = rgb(153,153,153, max = 255, alpha = 175)) # add the microbial taxon labels
graphics::legend("topright",bty = "n", legend = levels(mynmdssit$fact), pch = 21, pt.bg = c(myplotcols[1:length(levels(mynmdssit$fact))],rep(rgb(0, 0, 0, max = 255, alpha = 0),4)), pt.cex = 1.5) # add the legend
# shut the graphics device down
dev.off()










#### core microbiome analysis ----
# more information about this analysis can be found at https://microbiome.github.io/tutorials/Core.html 
library("devtools")
install_github("microbiome/microbiome")


# create the post weaning period early and late phyloseq objects
early_ps_fin_ASV_ra_heatmap <- subset_samples(ps_fin_ASV_ra_heatmap, When == "Early") # extract the Early period samples
early_ps_fin_ASV_ra_heatmap <- prune_taxa(taxa_sums(early_ps_fin_ASV_ra_heatmap) > 0, early_ps_fin_ASV_ra_heatmap) # remove taxa that have 0 relative abundances in the remaining samples

late_ps_fin_ASV_ra_heatmap <- subset_samples(ps_fin_ASV_ra_heatmap, When == "Late") # extract the Late period samples
late_ps_fin_ASV_ra_heatmap <- prune_taxa(taxa_sums(late_ps_fin_ASV_ra_heatmap) > 0, late_ps_fin_ASV_ra_heatmap) # remove taxa that have 0 relative abundances in the remaining samples



library(microbiome)

early_prev_in_1perc <- prevalence(early_ps_fin_ASV_ra_heatmap, detection = 1, sort = TRUE) # obtain the prevalence for min abundance of 1 percent
early_core.taxa <- core_members(early_ps_fin_ASV_ra_heatmap, detection = 1, prevalence = 0.95) # obtain the taxa with 1 percent RA in at least 80% of the samples
early_core.abundance <- sample_sums(core(early_ps_fin_ASV_ra_heatmap, detection = 1, prevalence = 0.95)) # obtain the per-sample relative abundance of the core taxa given the set criteria
summary(early_core.abundance) # obtain the distribution of these relative abundances
early_ps_fin_ASV_ra_heatmap_core  <- prune_taxa(early_core.taxa,late_ps_fin_ASV_ra_heatmap) # obtain the core microbiome phyloseq object


# late samples
late_prev_in_1perc <- prevalence(late_ps_fin_ASV_ra_heatmap, detection = 1, sort = TRUE) # obtain the prevalence for min abundance of 1 percent
late_core.taxa <- core_members(late_ps_fin_ASV_ra_heatmap, detection = 1, prevalence = 50/100) # obtain the taxa with 1 percent RA in at least 80% of the samples
late_core.abundance <- sample_sums(core(late_ps_fin_ASV_ra_heatmap, detection = 1, prevalence = 0.95)) # obtain the per-sample relative abundance of the core taxa given the set criteria
summary(late_core.abundance) # obtain the distribution of these relative abundances
late_ps_fin_ASV_ra_heatmap_core  <- prune_taxa(late_core.taxa,late_ps_fin_ASV_ra_heatmap) # obtain the core microbiome phyloseq object

# plot the core microbiome sizes for microorganisms of relative abundance with prevalence of at least 95%
pdf("core_sizes_with_ra_1_in_at_least_95_of_exp_period_samples.pdf", width = 4, height = 4)
par(mar = c(4,6,4,4))
boxplot(c(early_core.abundance,late_core.abundance) ~ c(rep("Early",length(early_core.abundance)),rep("Late",length(late_core.abundance))), xlab = "", ylab = "core microbiome\nrelative abundance (%)", frame = FALSE)
dev.off()

### prepare also the heatmaps of the core microbiomes
library(RColorBrewer)
prevalences <- seq(.05, 1, .05) # set the prevalence range to be plotted
detections <- seq(0, 40, length = 21) # set also the detections
# both the above values also depend on the dataset

# prepare the core for the early samples
pdf("core_heatmap_early.pdf", width = 5.5, height = 2.5)
p <- plot_core(early_ps_fin_ASV_ra_heatmap_core, plot.type = "heatmap", 
               prevalences = prevalences,
               detections = detections,
               colours = rev(brewer.pal(5, "Spectral")),
               min.prevalence = 0, horizontal = FALSE) 
# change the x axis label
str(p) # will provide the information of the plot parameters
p$labels$x <- "Detection Threshold (%)"
print(p)
dev.off()
# prepare the core for the late samples samples
pdf("core_heatmap_late.pdf", width = 5.5, height = 2.5)
p <- plot_core(late_ps_fin_ASV_ra_heatmap_core, plot.type = "heatmap", 
               prevalences = prevalences,
               detections = detections,
               colours = rev(brewer.pal(5, "Spectral")),
               min.prevalence = 0, horizontal = FALSE) 
p$labels$x <- "Detection Threshold (%)"
print(p)
dev.off()


# prepare the core for all samples
pdf("core_heatmap_all.pdf", width = 5.5, height = 30)
p <- plot_core(ps_fin_ASV_ra_heatmap, plot.type = "heatmap", 
               prevalences = prevalences,
               detections = detections,
               colours = rev(brewer.pal(5, "Spectral")),
               min.prevalence = .2, horizontal = FALSE) 
p$labels$x <- "Detection Threshold (%)"
print(p)
dev.off()









#### perform the differential abundance analysis with the DESeq2 normalization approach ----
# this script is based on the approaches of the phyloseq extensions repository found at https://joey711.github.io/phyloseq-extensions/DESeq2.html
library("DESeq2") # load the necessary package... it is a bioconductor package
packageVersion("DESeq2") # check the version

# use the counts phyloseq object (instead of relative abundances or percentages), since the DESeq normalization approach requires a counts table 
diagdds <- phyloseq_to_deseq2(ps_fin_ASV_best_tax, ~ When) # prepare the object necessary for the DESeq analysis (provide the counts and the independent variable name) out of the phyloseq object 
diagdds <- DESeq(diagdds, test="Wald", fitType="parametric") # perform the Wald significance pairwise test


DEres <- results(diagdds) # extract the results table
alpha = 0.01
sigtab = DEres[which(DEres$padj < alpha), ] # select the significant ASVs
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps_fin_ASV_best_tax)[rownames(sigtab), ], "matrix")) # combine with the taxonomy information
head(sigtab) # take a look at the results



## prepare the plot
library("ggplot2") # use ggplot2 for the plot... ggplot2 plots are in general quite pretty and sophisticated, parted by layers of commands, but this means that manipulating them is difficult for beginners... one trick is to use the str() command in order to identify the parameter we want to manipulate (in case we cannot change it with another approach) and then change it accordingly
theme_set(theme_bw()) # use the black and white plot theme set
scale_fill_discrete <- function(palname = "Set1", ...) { # function that will help manipulating the plot colours
  scale_fill_brewer(palette = palname, ...)
}

# use the phylum and genus in the plot and order according to the log2 fold change values between the early and late groups
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
p <- ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

pdf("differential_abundance.pdf", width = 8, height = 7)
print(p)
dev.off()








#### perform the correlation network analysis ----
## Several nice tips for preparing the network plot at: https://www.biostars.org/p/285296/
library(Hmisc)
ps_fin_ASV_ra_network <- ps_fin_ASV_ra_heatmap
taxa_names(ps_fin_ASV_ra_network) <- gsub("^[A-Z,a-z]+:","",taxa_names(ps_fin_ASV_ra_heatmap))
# obtain the Spearman correlation values
correl <- rcorr(as.matrix(otu_table(ps_fin_ASV_ra_network)), type = "spearman")

# save the r and p values in objects for downstream manipulation according to the set cutoffs
r_vals<-correl$r
p_vals<-correl$P
# set the diagonal of the p_values object to one (to avoid self loops.... although this I remove also with the simplify command)
diag(p_vals)<-1

# adjust the p-values for multiple hypothesis testing
p_vals_arr <- array(p_vals) # convert to an array
p_vals_arr_adj <- p.adjust(p_vals_arr, method = "BH") # adjust the p-values according to the selected method
p_vals_arr_adj_mat <- matrix(p_vals_arr_adj, nrow = nrow(p_vals), ncol = ncol(p_vals)) # convert back to a matrix
colnames(p_vals_arr_adj_mat) <- row.names(p_vals_arr_adj_mat) <- row.names(p_vals) # set the row and column names

# set the r_vals according to decided adjusted p-value and r/rho cutoffs
r_vals_padj <- r_vals # assign the new r_values
r_vals_padj_mat <- ifelse(p_vals_arr_adj_mat > 0.05,0,r_vals_padj) # adjust by setting to zero non-significant values

# if all went well run the graph algorithm using the igraph package commands
library(igraph) # can install it from CRAN
## Several nice tips for preparing the network plot at: https://www.biostars.org/p/285296/
g<-graph.adjacency(r_vals_padj_mat
                   , weighted=TRUE # for the edge calculation using continuous values
                   , diag=FALSE # we have also set the r_values to zero further up... so not necessary, but better safe than sorry
                   , mode="undirected" # undirected is the mode of choice
)

# modify plotting parameters
E(g)[which(E(g)$weight<0)]$color <- "darkred" # Colour negative correlation edges as red
E(g)[which(E(g)$weight>0)]$color <- "darkgreen" # Colour positive correlation edges as green
E(g)$weight <- abs(E(g)$weight) # Convert edge weights to absolute values
V(g)$shape <- "sphere" # Change shape of graph vertices
V(g)$color <- "skyblue" # Change colour of graph vertices
V(g)$vertex.frame.color <- "white" # Change colour of vertex frames
vSizes <- 2 * colMeans(data.frame(otu_table(ps_fin_ASV_ra_network)), na.rm = T) # Calculate the size of the vertices to be proportional to the relative abundance of each taxon represented by each vertex
edgeweights <- E(g)$weight * 5 # Amplify or decrease the width of the edges according to the plotting needs (you might need to change the multiplier)

# Convert the graph adjacency object into a minimum spanning tree based on Prim's algorithm
mst <- mst(g, algorithm="prim")

# Identify sub-communities
#mst.communities <- edge.betweenness.community(mst, weights=NULL, directed=FALSE) # calculate the betweenness centrality for edges in order to use for the cluster calculations 
#mst.clustering <- make_clusters(mst, membership=mst.communities$membership) # calculate the memberships
#myvertcols <- as.factor(tax_table(ps_fin_ASV_ra_network)[,2]) # color vertices according to taxon phylum
#V(mst)$color <- myvertcols

# prepare the plot
# start a graphics device there
pdf("correlation_network_plot.pdf", height = 40, width = 40)
plot(
  mst
  , layout=layout.fruchterman.reingold
  , edge.curved=TRUE
  , vertex.size=vSizes,
  , vertex.label.dist = ifelse(vSizes/5 > 0.5, 0.5, vSizes/5)
  , vertex.label.color="black"
  , asp=FALSE
  , vertex.label.cex=1.5
  , edge.width=edgeweights
  , edge.arrow.mode=0
)


# close the graphics device
dev.off()







