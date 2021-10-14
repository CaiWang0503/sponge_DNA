#https://benjjneb.github.io/decontam/vignettes/decontam_intro.html
BiocManager::install("decontam")
library(decontam)
library(phyloseq); packageVersion("phyloseq")
library(ggplot2)
ps <- readRDS(system.file("extdata", "MUClite.rds", package="decontam"))
head(sample_data(ps))
df <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()

contamdf.freq <- isContaminant(ps, method="frequency", conc="quant_reading")#the default value of threshold = 0.1
head(contamdf.freq)
table(contamdf.freq$contaminant)
head(which(contamdf.freq$contaminant))
plot_frequency(ps, taxa_names(ps)[c(1,3)], conc="quant_reading") + 
  xlab("DNA Concentration (PicoGreen fluorescent intensity)")

set.seed(100)
plot_frequency(ps, taxa_names(ps)[sample(which(contamdf.freq$contaminant),3)], conc="quant_reading") +
  xlab("DNA Concentration (PicoGreen fluorescent intensity)")

ps.noncontam <- prune_taxa(!contamdf.freq$contaminant, ps)


####################
sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Control Sample"
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)
head(which(contamdf.prev$contaminant))

contamdf.prev05 <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05$contaminant)
ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Control Sample", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "True Sample", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

#########################
#Try all the methods to see the difference in detected contaminated sequences

contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg")
