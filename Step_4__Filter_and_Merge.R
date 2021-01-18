### Information ###

# It is important to filter out MOTUs that have a low read count, and are likely contamination.
# Samples with low read count should also be removed
# Non-prey will also be removed (such as host or parasitic nematodes for example)
# This script filters the Gillet and Zeale datasets separately, and then merges them into a single phyloseq dataset for downstream analysis


### Set up ###

library(phyloseq)
library(ggplot2)
library(dplyr)
library(vegan)
library(stringr)
library(knitr)
library(hilldiv)
library(kableExtra)


### 1) Load and filter MOTUs and Samples in Gillet dataset ###

load("../gillet_dataset/sumaclust98/phyloseq_object_clust_iden98_taxa_assigned_no_singletons.RData")
taxa_names(gillet.phylo) <- paste("gMOTU", seq(nrow(tax_table(gillet.phylo))), sep="_")

# remove samples with less than 1000 reads
gillet.phylo <- prune_samples(sample_sums(gillet.phylo) > 1000, gillet.phylo)

gillet.phylo


# Remove 'sample.' part of sample names that obitools annoyingly adds

chck <- sample_names(gillet.phylo)
chck <- str_remove(chck, "sample.")
sample_names(gillet.phylo) <- chck

### Examine the blanks?

n <- which(otu_table(gillet.phylo)[,"G_NEG"] > 0)

m <- which(otu_table(gillet.phylo)[,"G_Neg"] > 0)

l <- unique(c(n,m))
blnk.df <- as.data.frame(as.matrix(tax_table(gillet.phylo))[l,4:7])
blnk.df$total.reads <- taxa_sums(gillet.phylo)[l]
blnk.df$Neg1.reads <- otu_table(gillet.phylo)[l, "G_NEG"]
for (i in 1:nrow(blnk.df)) {blnk.df$Neg1.prop[i] <- (blnk.df$Neg1.reads[i] / blnk.df$total.reads[i]) * 100 }
blnk.df$Neg2.reads <- otu_table(gillet.phylo)[l, "G_Neg"]
for (i in 1:nrow(blnk.df)) {blnk.df$Neg2.prop[i] <- (blnk.df$Neg2.reads[i] / blnk.df$total.reads[i]) * 100 }
blnk.df$perc <- apply(blnk.df[,c("Neg1.prop","Neg2.prop")], 1, sum)
rownames(blnk.df) <- 1:nrow(blnk.df)

kable(blnk.df, caption = "MOTUs identified in the blanks")


### Remove taxa of which the blanks hold over 2% of the total reads for that MOTU

tab.nam <- blnk.df[,"perc"] > 2 # 2 is for 2%
tab.df <- blnk.df[tab.nam,]

removeTaxa <- rownames(tab.df) # Lists the MOTUs to remove
phy.obj <- subset_taxa(gillet.phylo, !(taxa_names(gillet.phylo) %in% removeTaxa))
phy.obj

### Visualise vertebrate amplification in samples

# Load colour palette
pal.o = c("#f0a3ff",
          "#0075dc",
          "#993f00",
          "#4c005c",
          "#191919",
          "#005c31",
          "#2bce48",
          "#ffcc99",
          "#808080",
          "#94ffb5",
          "#8f7c00",
          "#9dcc00",
          "#c20088",
          "blue",
          "#ffa405",
          "#ffa8bb",
          "#426600",
          "#ff0010",
          "#5ef1f2",
          "#00998f",
          "#e0ff66",
          "indianred",
          "#003380",
          "green",
          "khaki4",
          "darkred",
          "coral4",
          "violetred2",
          "#0075dc",
          "#993f00")


samples.phylo <- gillet.phylo

samples.phylo <- tax_glom(gillet.phylo, taxrank = "order")
n <- grep("Chiroptera", tax_table(samples.phylo)[,"order"])
m <- grep("Eulipotyphla", tax_table(samples.phylo)[,"order"])
tax_table(samples.phylo)[-c(n, m), 1:4] <- "Other"
unique(tax_table(samples.phylo)[,"order"])

mm.oth <- tax_glom(samples.phylo, taxrank = "order")
mm.oth
mm.oth.ra = transform_sample_counts(mm.oth, function(x) 100 * x/sum(x))

ra.samples.bar <- phyloseq::plot_bar(mm.oth.ra) # extracts information needed for barplots
ra.samples.bar.data <- ra.samples.bar$data
p1.2 <- ggplot(ra.samples.bar.data, aes(x= Sample, y=Abundance, fill = order))

p1.2 + geom_bar(stat="identity", color="black") +
  scale_fill_manual(values = pal.o) +
  facet_wrap(~ mammal, scale = "free_x") +
  theme_classic() +
  theme(legend.position = "right") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


### Check range of vertebrate amplification in samples

samples.phylo <- subset_samples(phy.obj, mammal != "BLANK")

n <- grep("Chordata", tax_table(samples.phylo)[,"phylum"])
tax_table(samples.phylo)[-n, 1:4] <- "Other"
unique(tax_table(samples.phylo)[,"phylum"])

mm.oth <- tax_glom(samples.phylo, taxrank = "phylum")
mm.oth
mm.oth.ra = transform_sample_counts(mm.oth, function(x) 100 * x/sum(x))

tax_table(mm.oth.ra)[,1:2]
otu_table(mm.oth.ra)

# Range of Vertebrate amplification across all samples
range(otu_table(mm.oth.ra)[1,])

# Range of Vertebrate amplification across GWTS
gwts.prop <- subset_samples(mm.oth.ra, mammal == "GWTS")
range(otu_table(gwts.prop)[1,])

# Range of Vertebrate amplification across Pygmies
pyg.prop <- subset_samples(mm.oth.ra, mammal == "Pygmy")
range(otu_table(pyg.prop)[1,])

# Range of Vertebrate amplification across Bats
bat.prop <- subset_samples(mm.oth.ra, mammal == "Bat")
range(otu_table(bat.prop)[1,])

#### Filtering ####

samples.phylo <- subset_samples(phy.obj, mammal != "BLANK") # remove blanks

# Remove non-prey
diet.prey <- subset_taxa(samples.phylo, !(class %in% c("Mammalia",
                                                       "none",
                                                       "Actinopteri",
                                                       "Bdelloidea",
                                                       "Udeonychophora", # velvet worms
                                                       "Merostomata", # horse shoe crabs
                                                       "Gammaproteobacteria", # bacteria
                                                       "Magnoliopsida", # plants
                                                       "Monogononta", # rotifers
                                                       "Dothideomycetes", # fungi
                                                       "Trebouxiophyceae", # green algae
                                                       "Chondrichthyes", # Cartilaginous fish
                                                       "Mucoromycetes", # fungi
                                                       "Phylum_Endomyxa", # micro things
                                                       "Eutardigrada", # tartigrades!!
                                                       "Elardia", # Amoebas
                                                       "Cephalopoda", # Cephalopods
                                                       "Amphibia", # Amphibians
                                                       "Aves", # Birds
                                                       "Chromadorea", # roundworms
                                                       "Hexanauplia",  # parasitic crustaceans
                                                       "Kingdom_Metazoa",
                                                       "Kingdom_",
                                                       "Phylum_Discosea", # amoebas
                                                       "Branchiopoda", # marine crustaceans
                                                       "Phylum_Nematoda")))

# remove samples with less than 1000 reads
sampl.filt <- prune_samples(sample_sums(diet.prey) > 1000, diet.prey)
otu.tab <- as.data.frame(otu_table(sampl.filt))
new.otu.tab <- copy_filt(otu.tab, 0.0001) # Remove MOTUs with less than 0.01% reads in each sample
new.otu.tab <- as.matrix(new.otu.tab)
otu_table(sampl.filt) <- otu_table(new.otu.tab, taxa_are_rows = TRUE)
sampl.filt

# Remove any remaining taxa with less than 5 reads in total from dataset
final.diet <- prune_taxa(taxa_sums(sampl.filt) > 4, sampl.filt)
final.diet

# Check range of read depth of samples
range(sample_sums(final.diet))

# Check average read depth
mean(sample_sums(final.diet))

# Check range of total reads per taxa
range(taxa_sums(final.diet))


# hill_div packages assessment of read depth per sample, according to a shannon diversity equivilent
depth_cov(new.otu.tab,
          qvalue = 1)


## Rarefaction analysis

Bat_G <- prune_samples(sample_data(final.diet)$mammal == "Bat", final.diet)
df1 <- as.data.frame(t(as.matrix(otu_table(Bat_G))))
gwts_G <- prune_samples(sample_data(final.diet)$mammal == "GWTS", final.diet)
df2 <- as.data.frame(t(as.matrix(otu_table(gwts_G))))
pyg_G <- prune_samples(sample_data(final.diet)$mammal == "Pygmy", final.diet)
df3 <- as.data.frame(t(as.matrix(otu_table(pyg_G))))

set.seed(57)
r1 <- rarecurve(df1[,])

r2 <- rarecurve(df2[,])

r3 <- rarecurve(df3[,])

out <- r1 # change to the rarefaction curve to plot (r1, r2 or r3 - see above)

par(mar = c(4.5, 4.5, 1, 1)) # bottom, left, top, right

plot(c(1, 15000), c(1, 120), xlab = "Reads",
     ylab = "MOTUs", type = "n", cex.axis = 2, cex.lab = 2, las = 1)
#abline(v = 1000)
#abline(v = 5000)
for (i in seq_along(out)) {
  N <- attr(out[[i]], "Subsample")
  lines(N, out[[i]], col = "black")
}

par(mar=c(5.1, 4.1, 4.1, 2.1)) # back to default plot parameters


# final phyloseq object for gillet dataset

final.diet.g <- final.diet

### Identify how many MOTUs are identified to Order, Family, Genus and Species level

df1 <- as.data.frame(as.matrix(tax_table(final.diet.g)))
df <- df1

df$pident <- as.character(df$pident)
df$pident <- as.numeric(df$pident)

###################
## Species
n <- which(df[,"pident"] > 97.9999)
gn <- grep("Genus_", df[n,"species"])
fm <- grep("Family_", df[n,"species"])
or <- grep("Order_", df[n,"species"])
cl <- grep("Class_", df[n,"species"])
ph <- grep("Phylum_", df[n,"species"])
kg <- grep("Kindgom_", df[n,"species"])
nn <- grep("none", df[n,"species"])
s.p <- grep("_sp._", df[n,"species"])
n.r <- grep("_nr._", df[n,"species"])
c.f <- grep("_cf._", df[n,"species"])
x <- length(c(gn,fm,or,cl,ph,kg,nn,s.p,n.r,c.f))
sp.y <- length(df[n,"species"]) - x
chck <- df[n[-c(gn,fm,or,cl,ph,kg,nn,s.p,n.r,c.f)],"species"]
un.sp.y <- length(unique(chck))

###################
## Genus
df <- df[-n,]
n <- which(df[,"pident"] > 94.9999)
gn <- grep("Genus_", df[n,"genus"])
fm <- grep("Family_", df[n,"genus"])
or <- grep("Order_", df[n,"genus"])
cl <- grep("Class_", df[n,"genus"])
ph <- grep("Phylum_", df[n,"genus"])
kg <- grep("Kindgom_", df[n,"genus"])
nn <- grep("none", df[n,"genus"])
s.p <- grep("_sp._", df[n,"genus"])
n.r <- grep("_nr._", df[n,"genus"])
c.f <- grep("_cf._", df[n,"genus"])
g.g <- grep("_gen._", df[n,"genus"])
x <- length(c(gn,fm,or,cl,ph,kg,nn,s.p,n.r,c.f,g.g))
gn.y <- length(df[n,"genus"]) - x

gn <- grep("Genus_", df1[,"genus"])
fm <- grep("Family_", df1[,"genus"])
or <- grep("Order_", df1[,"genus"])
cl <- grep("Class_", df1[,"genus"])
ph <- grep("Phylum_", df1[,"genus"])
kg <- grep("Kindgom_", df1[,"genus"])
nn <- grep("none", df1[,"genus"])
s.p <- grep("_sp._", df1[,"genus"])
n.r <- grep("_nr._", df1[,"genus"])
c.f <- grep("_cf._", df1[,"genus"])
g.g <- grep("_gen._", df1[,"genus"])
chck <- df1[-c(gn,fm,or,cl,ph,kg,nn,s.p,n.r,c.f),"genus"]
un.gn.y <- length(unique(chck))

###################
## Family
df <- df[-n,]
n <- which(df[,"pident"] > 92.999)
gn <- grep("Genus_", df[n,"family"])
fm <- grep("Family_", df[n,"family"])
or <- grep("Order_", df[n,"family"])
cl <- grep("Class_", df[n,"family"])
ph <- grep("Phylum_", df[n,"family"])
kg <- grep("Kindgom_", df[n,"family"])
nn <- grep("none", df[n,"family"])
s.p <- grep("_sp._", df[n,"family"])
n.r <- grep("_nr._", df[n,"family"])
c.f <- grep("_cf._", df[n,"family"])
g.g <- grep("_gen._", df[n,"family"])
x <- length(c(gn,fm,or,cl,ph,kg,nn,s.p,n.r,c.f,g.g))
fm.y <- length(df[n,"family"]) - x

gn <- grep("Genus_", df1[,"family"])
fm <- grep("Family_", df1[,"family"])
or <- grep("Order_", df1[,"family"])
cl <- grep("Class_", df1[,"family"])
ph <- grep("Phylum_", df1[,"family"])
kg <- grep("Kindgom_", df1[,"family"])
nn <- grep("none", df1[,"family"])
s.p <- grep("_sp._", df1[,"family"])
n.r <- grep("_nr._", df1[,"family"])
c.f <- grep("_cf._", df1[,"family"])
g.g <- grep("_gen._", df1[,"family"])
chck <- df1[-c(gn,fm,or,cl,ph,kg,nn,s.p,n.r,c.f),"family"]
un.fm.y <- length(unique(chck))

###################
## Order
df <- df[-n,]
n <- which(df[,"pident"] > 89.9999)
gn <- grep("Genus_", df[n,"order"])
fm <- grep("Family_", df[n,"order"])
or <- grep("Order_", df[n,"order"])
cl <- grep("Class_", df[n,"order"])
ph <- grep("Phylum_", df[n,"order"])
kg <- grep("Kindgom_", df[n,"order"])
nn <- grep("none", df[n,"order"])
s.p <- grep("_sp._", df[n,"order"])
n.r <- grep("_nr._", df[n,"order"])
c.f <- grep("_cf._", df[n,"order"])
g.g <- grep("_gen._", df[n,"order"])
x <- length(c(gn,fm,or,cl,ph,kg,nn,s.p,n.r,c.f,g.g))
or.y <- length(df[n,"order"]) - x

gn <- grep("Genus_", df1[,"order"])
fm <- grep("Family_", df1[,"order"])
or <- grep("Order_", df1[,"order"])
cl <- grep("Class_", df1[,"order"])
ph <- grep("Phylum_", df1[,"order"])
kg <- grep("Kindgom_", df1[,"order"])
nn <- grep("none", df1[,"order"])
s.p <- grep("_sp._", df1[,"order"])
n.r <- grep("_nr._", df1[,"order"])
c.f <- grep("_cf._", df1[,"order"])
g.g <- grep("_gen._", df1[,"order"])
chck <- df1[-c(gn,fm,or,cl,ph,kg,nn,s.p,n.r,c.f),"order"]
un.or.y <- length(unique(chck))

df <- df1

# Save results to table to compare to Zeale primers later

tabx <- data.frame(primer = NA,
                   order = NA, unique.orders = NA,
                   family = NA, unique.families = NA,
                   genus = NA, unique.genus = NA,
                   species = NA, unique.species = NA,
                   total.taxa = NA,
                   total.reads = NA)

tabx[1,] <- c("Gillet", or.y, un.or.y, fm.y, un.fm.y,
              gn.y, un.gn.y, sp.y, un.sp.y,
              ntaxa(final.diet.g), sum(sample_sums(final.diet.g)))


### 2) Load and filter MOTUs and Samples in Zeale dataset ###

load("../zeale_dataset/sumaclust98/phyloseq_object_zeale_clust_iden98_taxa_assigned_no_singletons.RData")
taxa_names(zeale.phylo) <- paste("zMOTU", seq(nrow(tax_table(zeale.phylo))), sep="_")

zeale.phylo

# Remove 'sample.' part of sample names that obitools adds

chck <- sample_names(zeale.phylo)
chck <- str_remove(chck, "sample.")
sample_names(zeale.phylo) <- chck

# Examine the blanks?

n <- which(otu_table(zeale.phylo)[,"Z_NEG"] > 0)

m <- which(otu_table(zeale.phylo)[,"Z_Neg"] > 0)

l <- unique(c(n,m))
blnk.df <- as.data.frame(as.matrix(tax_table(zeale.phylo))[l,4:7])
blnk.df$total.reads <- taxa_sums(zeale.phylo)[l]
blnk.df$Neg1.reads <- otu_table(zeale.phylo)[l, "Z_NEG"]
for (i in 1:nrow(blnk.df)) {blnk.df$Neg1.prop[i] <- (blnk.df$Neg1.reads[i] / blnk.df$total.reads[i]) * 100 }
blnk.df$Neg2.reads <- otu_table(zeale.phylo)[l, "Z_Neg"]
for (i in 1:nrow(blnk.df)) {blnk.df$Neg2.prop[i] <- (blnk.df$Neg2.reads[i] / blnk.df$total.reads[i]) * 100 }
blnk.df$perc <- apply(blnk.df[,c("Neg1.prop","Neg2.prop")], 1, sum)
rownames(blnk.df) <- 1:nrow(blnk.df)

kable(blnk.df, caption = "MOTUs identified in the blanks")


# Remove taxa of which the blanks hold over 2% of the total reads for that MOTU

tab.nam <- blnk.df[,"perc"] > 2 # 2 is for 2%
tab.df <- blnk.df[tab.nam,]

removeTaxa <- rownames(tab.df) # Lists the MOTUs to remove
phy.obj <- subset_taxa(zeale.phylo, !(taxa_names(zeale.phylo) %in% removeTaxa))
phy.obj

## Visualise any vertebrate amplification

samples.phylo <- zeale.phylo

samples.phylo <- tax_glom(zeale.phylo, taxrank = "order")
n <- grep("Chiroptera", tax_table(samples.phylo)[,"order"])
m <- grep("Eulipotyphla", tax_table(samples.phylo)[,"order"])
tax_table(samples.phylo)[-c(n, m), 1:4] <- "Other"
unique(tax_table(samples.phylo)[,"order"])

mm.oth <- tax_glom(samples.phylo, taxrank = "order")
mm.oth
mm.oth.ra = transform_sample_counts(mm.oth, function(x) 100 * x/sum(x))

ra.samples.bar <- phyloseq::plot_bar(mm.oth.ra) # extracts information needed for barplots
ra.samples.bar.data <- ra.samples.bar$data
p1.2 <- ggplot(ra.samples.bar.data, aes(x= Sample, y=Abundance, fill = order))

p1.2 + geom_bar(stat="identity", color="black") +
  scale_fill_manual(values = pal.o) +
  facet_wrap(~ mammal, scale = "free_x") +
  theme_classic() +
  theme(legend.position = "right") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# Check range of vertebrate amplification in samples

samples.phylo <- subset_samples(phy.obj, mammal != "BLANK")

n <- grep("Chordata", tax_table(samples.phylo)[,"phylum"])
tax_table(samples.phylo)[-n, 1:4] <- "Other"
unique(tax_table(samples.phylo)[,"phylum"])

mm.oth <- tax_glom(samples.phylo, taxrank = "phylum")
mm.oth
mm.oth.ra = transform_sample_counts(mm.oth, function(x) 100 * x/sum(x))

tax_table(mm.oth.ra)[,1:2]
otu_table(mm.oth.ra)

# Range of Vertebrate amplification across all samples
range(otu_table(mm.oth.ra)[1,])

# Range of Vertebrate amplification across GWTS
gwts.prop <- subset_samples(mm.oth.ra, mammal == "GWTS")
range(otu_table(gwts.prop)[1,])

# Range of Vertebrate amplification across Pygmies
pyg.prop <- subset_samples(mm.oth.ra, mammal == "Pygmy")
range(otu_table(pyg.prop)[1,])

# Range of Vertebrate amplification across Bats
bat.prop <- subset_samples(mm.oth.ra, mammal == "Bat")
range(otu_table(bat.prop)[1,])

#### Filtering ####

samples.phylo <- subset_samples(phy.obj, mammal != "BLANK") # Remove negative controls

# Remove non-prey
diet.prey <- subset_taxa(samples.phylo, !(class %in% c("Mammalia",
                                                       "none",
                                                       "Actinopteri",
                                                       "Bdelloidea",
                                                       "Udeonychophora", # velvet worms
                                                       "Merostomata", # horse shoe crabs
                                                       "Gammaproteobacteria", # bacteria
                                                       "Magnoliopsida", # plants
                                                       "Monogononta", # rotifers
                                                       "Dothideomycetes", # fungi
                                                       "Trebouxiophyceae", # green algae
                                                       "Chondrichthyes", # Cartilaginous fish
                                                       "Mucoromycetes", # fungi
                                                       "Phylum_Endomyxa", # micro things
                                                       "Eutardigrada", # tartigrades!!
                                                       "Elardia", # Amoebas
                                                       "Cephalopoda", # Cephalopods
                                                       "Amphibia", # Amphibians
                                                       "Aves", # Birds
                                                       "Chromadorea", # roundworms
                                                       "Hexanauplia",  # parasitic crustaceans
                                                       "Kingdom_Metazoa",
                                                       "Kingdom_",
                                                       "Phylum_Discosea", # amoebas
                                                       "Branchiopoda", # marine crustaceans
                                                       "Phylum_Nematoda")))

sampl.filt <- prune_samples(sample_sums(diet.prey) > 1000, diet.prey)
otu.tab <- as.data.frame(otu_table(sampl.filt))
new.otu.tab <- copy_filt(otu.tab, 0.0001) # Remove MOTUs with less than 0.01% reads in each sample
new.otu.tab <- as.matrix(new.otu.tab)
otu_table(sampl.filt) <- otu_table(new.otu.tab, taxa_are_rows = TRUE)
sampl.filt

# Remove any remaining taxa with less than 5 reads in total from dataset
final.diet <- prune_taxa(taxa_sums(sampl.filt) > 4, sampl.filt)
final.diet

# Check range of read depth of samples
range(sample_sums(final.diet))

# Check average read depth
mean(sample_sums(final.diet))

# Check range of total reads per taxa
range(taxa_sums(final.diet))


# hill_div packages assessment of read depth per sample, according to a shannon diversity equivilent
depth_cov(new.otu.tab,
          qvalue = 1)


## Rarefaction analysis

Bat_Z <- prune_samples(sample_data(final.diet)$mammal == "Bat", final.diet)
df1 <- as.data.frame(t(as.matrix(otu_table(Bat_Z))))
gwts_Z <- prune_samples(sample_data(final.diet)$mammal == "GWTS", final.diet)
df2 <- as.data.frame(t(as.matrix(otu_table(gwts_Z))))
pyg_Z <- prune_samples(sample_data(final.diet)$mammal == "Pygmy", final.diet)
df3 <- as.data.frame(t(as.matrix(otu_table(pyg_Z))))

set.seed(57)
r1 <- rarecurve(df1[,])

r2 <- rarecurve(df2[,])

r3 <- rarecurve(df3[,])

out <- r1 # change to the rarefaction curve to plot (r1, r2 or r3 - see above)

par(mar = c(4.5, 4.5, 1, 1)) # bottom, left, top, right

plot(c(1, 15000), c(1, 120), xlab = "Reads",
     ylab = "MOTUs", type = "n", cex.axis = 2, cex.lab = 2, las = 1)
#abline(v = 1000)
#abline(v = 5000)
for (i in seq_along(out)) {
  N <- attr(out[[i]], "Subsample")
  lines(N, out[[i]], col = "black")
}

par(mar=c(5.1, 4.1, 4.1, 2.1)) # back to default plot parameters

# Save Zeale dataset phyloseq object
final.diet.z <- final.diet


### Identify how many MOTUs are identified to Order, Family, Genus and Species level

df1 <- as.data.frame(as.matrix(tax_table(final.diet.z)))
df <- df1

df$pident <- as.character(df$pident)
df$pident <- as.numeric(df$pident)

###################
## Species
n <- which(df[,"pident"] > 97.9999)
gn <- grep("Genus_", df[n,"species"])
fm <- grep("Family_", df[n,"species"])
or <- grep("Order_", df[n,"species"])
cl <- grep("Class_", df[n,"species"])
ph <- grep("Phylum_", df[n,"species"])
kg <- grep("Kindgom_", df[n,"species"])
nn <- grep("none", df[n,"species"])
s.p <- grep("_sp._", df[n,"species"])
n.r <- grep("_nr._", df[n,"species"])
c.f <- grep("_cf._", df[n,"species"])
x <- length(c(gn,fm,or,cl,ph,kg,nn,s.p,n.r,c.f))
sp.y <- length(df[n,"species"]) - x
chck <- df[n[-c(gn,fm,or,cl,ph,kg,nn,s.p,n.r,c.f)],"species"]
un.sp.y <- length(unique(chck))

###################
## Genus
df <- df[-n,]
n <- which(df[,"pident"] > 94.9999)
gn <- grep("Genus_", df[n,"genus"])
fm <- grep("Family_", df[n,"genus"])
or <- grep("Order_", df[n,"genus"])
cl <- grep("Class_", df[n,"genus"])
ph <- grep("Phylum_", df[n,"genus"])
kg <- grep("Kindgom_", df[n,"genus"])
nn <- grep("none", df[n,"genus"])
s.p <- grep("_sp._", df[n,"genus"])
n.r <- grep("_nr._", df[n,"genus"])
c.f <- grep("_cf._", df[n,"genus"])
g.g <- grep("_gen._", df[n,"genus"])
x <- length(c(gn,fm,or,cl,ph,kg,nn,s.p,n.r,c.f,g.g))
gn.y <- length(df[n,"genus"]) - x

gn <- grep("Genus_", df1[,"genus"])
fm <- grep("Family_", df1[,"genus"])
or <- grep("Order_", df1[,"genus"])
cl <- grep("Class_", df1[,"genus"])
ph <- grep("Phylum_", df1[,"genus"])
kg <- grep("Kindgom_", df1[,"genus"])
nn <- grep("none", df1[,"genus"])
s.p <- grep("_sp._", df1[,"genus"])
n.r <- grep("_nr._", df1[,"genus"])
c.f <- grep("_cf._", df1[,"genus"])
g.g <- grep("_gen._", df1[,"genus"])
chck <- df1[-c(gn,fm,or,cl,ph,kg,nn,s.p,n.r,c.f),"genus"]
un.gn.y <- length(unique(chck))

###################
## Family
df <- df[-n,]
n <- which(df[,"pident"] > 92.999)
gn <- grep("Genus_", df[n,"family"])
fm <- grep("Family_", df[n,"family"])
or <- grep("Order_", df[n,"family"])
cl <- grep("Class_", df[n,"family"])
ph <- grep("Phylum_", df[n,"family"])
kg <- grep("Kindgom_", df[n,"family"])
nn <- grep("none", df[n,"family"])
s.p <- grep("_sp._", df[n,"family"])
n.r <- grep("_nr._", df[n,"family"])
c.f <- grep("_cf._", df[n,"family"])
g.g <- grep("_gen._", df[n,"family"])
x <- length(c(gn,fm,or,cl,ph,kg,nn,s.p,n.r,c.f,g.g))
fm.y <- length(df[n,"family"]) - x

gn <- grep("Genus_", df1[,"family"])
fm <- grep("Family_", df1[,"family"])
or <- grep("Order_", df1[,"family"])
cl <- grep("Class_", df1[,"family"])
ph <- grep("Phylum_", df1[,"family"])
kg <- grep("Kindgom_", df1[,"family"])
nn <- grep("none", df1[,"family"])
s.p <- grep("_sp._", df1[,"family"])
n.r <- grep("_nr._", df1[,"family"])
c.f <- grep("_cf._", df1[,"family"])
g.g <- grep("_gen._", df1[,"family"])
chck <- df1[-c(gn,fm,or,cl,ph,kg,nn,s.p,n.r,c.f),"family"]
un.fm.y <- length(unique(chck))

###################
## Order
df <- df[-n,]
n <- which(df[,"pident"] > 89.9999)
gn <- grep("Genus_", df[n,"order"])
fm <- grep("Family_", df[n,"order"])
or <- grep("Order_", df[n,"order"])
cl <- grep("Class_", df[n,"order"])
ph <- grep("Phylum_", df[n,"order"])
kg <- grep("Kindgom_", df[n,"order"])
nn <- grep("none", df[n,"order"])
s.p <- grep("_sp._", df[n,"order"])
n.r <- grep("_nr._", df[n,"order"])
c.f <- grep("_cf._", df[n,"order"])
g.g <- grep("_gen._", df[n,"order"])
x <- length(c(gn,fm,or,cl,ph,kg,nn,s.p,n.r,c.f,g.g))
or.y <- length(df[n,"order"]) - x

gn <- grep("Genus_", df1[,"order"])
fm <- grep("Family_", df1[,"order"])
or <- grep("Order_", df1[,"order"])
cl <- grep("Class_", df1[,"order"])
ph <- grep("Phylum_", df1[,"order"])
kg <- grep("Kindgom_", df1[,"order"])
nn <- grep("none", df1[,"order"])
s.p <- grep("_sp._", df1[,"order"])
n.r <- grep("_nr._", df1[,"order"])
c.f <- grep("_cf._", df1[,"order"])
g.g <- grep("_gen._", df1[,"order"])
chck <- df1[-c(gn,fm,or,cl,ph,kg,nn,s.p,n.r,c.f),"order"]
un.or.y <- length(unique(chck))

df <- df1

#un.gn.y <- length(unique(df1[,"genus"]))
#un.fm.y <- length(unique(df1[,"family"]))
#un.or.y <- length(unique(df1[,"order"]))

tabx[2,] <- c("Zeale", or.y, un.or.y, fm.y, un.fm.y,
              gn.y, un.gn.y, sp.y, un.sp.y,
              ntaxa(final.diet.z), sum(sample_sums(final.diet.z)))

# View table
landscape(knitr::kable(tabx))


## Visualise the differences between taxa identification between primers ##

library(reshape2)

tabx

df <- melt(data = tabx, id.vars = "primer",
           measure.vars = c("unique.orders",
                            "unique.families",
                            "unique.genus",
                            "unique.species"))
df$value <- as.numeric(df$value)
df$primer <- factor(df$primer, levels = c("Zeale", "Gillet"))

p <- ggplot(df, aes(x = variable, y=value)) +
  geom_bar(aes(fill = primer),
           stat = "identity",
           color = "black",
           position = position_dodge()) +
  theme_classic() +
  scale_fill_manual(values = c("darkblue", "pink")) +
  scale_x_discrete(labels=c("unique.orders" = "Identified\nOrders",
                            "unique.families" = "Identified\nFamilies",
                            "unique.genus" = "Identified\nGenera",
                            "unique.species" = "Identified\nSpecies")) +
  theme(axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5,
                                   face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold")) +
  labs(y = "Number Identified") +
  theme(legend.position = c(0.15,0.8)) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 20)) +
  #theme(legend.position = "bottom") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15, face = "bold"))

p



### 3) Merge The Gillet and Zeale Primer Datasets ###

final.diet.g
final.diet.z

sample_data(final.diet.g)$primer <- rep("Gillet", nsamples(final.diet.g))
sample_data(final.diet.z)$primer <- rep("Zeale", nsamples(final.diet.z))
mrg <- merge_phyloseq(final.diet.g, final.diet.z)
new.df <- read.csv("./samplesheet.csv")
rownames(new.df) <- new.df$id
new.df <- new.df[,-1]
sample_data(mrg) <- new.df

sample_data(mrg)$mammal_primer <- paste(sample_data(mrg)$mammal,
                                               sample_data(mrg)$primer, sep="_")

save(mrg, file = "./merged_primer_dataset.RData")
