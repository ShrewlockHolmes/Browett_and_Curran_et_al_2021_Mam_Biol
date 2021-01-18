### Information ###

# This script performs analyses and generates figures seen in manuscript

# This script requires the "./merged_primer_dataset.RData" file created in Step_4__Filter_and_Merge.R

# Analyses are:
# Alpha and Beta diversity
# NMDS
# hierarchical clustering
# Venn Diagrams
# Composition of diet
# Random Forest Classifier

#### Set Up ####

# setwd("~/Dropbox/Analyses")

library(stringr)
library(ggplot2)
library(phyloseq)
library(FSA)
library(microbiome)
library(vegan)
library(grid)
library(knitr)
library(spaa)
library(dplyr)
library(reshape2)
library(gridExtra)
library(dendextend)

#### Load Phyloseq object ####

load("./merged_primer_dataset.RData")

#### BOLD taxa ####

# Abundant MOTUs without taxonomic assignment from NCBI were manually blasted against the BOLD reference database
# These are the MOTUs with improved taxonomic resolution

final.diet <- mrg

tax_table(final.diet)["gMOTU_14", 4:9] <- c("Diptera",
                                            "Limoniidae",
                                            "Molophilus",
                                            "Molophilus_griseus",
                                            as.numeric(100.00), "BOLD")


tax_table(final.diet)["gMOTU_21", 4:9] <- c("Archaeognatha",
                                            "Machilidae",
                                            "Dilta",
                                            "Dilta_hibernica",
                                            as.numeric(99.24), "BOLD")

tax_table(final.diet)["gMOTU_26", 4:9] <- c("Isopoda",
                                            "Asellidae",
                                            "Family_Asellidae",
                                            "Family_Asellidae",
                                            as.numeric(94.87), "BOLD")

tax_table(final.diet)["gMOTU_51", 4:9] <- c("Diptera",
                                            "Limoniidae",
                                            "Austrolimnophila",
                                            "Austrolimnophila_ochracea",
                                            as.numeric(100.00), "BOLD")

tax_table(final.diet)["gMOTU_24", 4:9] <- c("Hymenoptera",
                                            "Formicidae",
                                            "Stenamma",
                                            "Stenamma_debile",
                                            as.numeric(100.00), "BOLD")

tax_table(final.diet)["gMOTU_61", 4:9] <- c("Hymenoptera",
                                            "Formicidae",
                                            "Myrmecina",
                                            "Myrmecina_graminicola",
                                            as.numeric(100.00), "BOLD")

tax_table(final.diet)["gMOTU_76", 4:9] <- c("Diptera",
                                            "Sciaridae",
                                            "Prosciara",
                                            "Prosciara_ungulata",
                                            as.numeric(100.00), "BOLD")

tax_table(final.diet)["gMOTU_77", 4:9] <- c("Lithobiomorpha",
                                            "Lithobiidae",
                                            "Lithobius",
                                            "Genus_Lithobius",
                                            as.numeric(95.45), "BOLD")

tax_table(final.diet)["gMOTU_96", 4:9] <- c("Isopoda",
                                            "Oniscidae",
                                            "Oniscus",
                                            "Oniscus_asellus",
                                            as.numeric(100.00), "BOLD")

tax_table(final.diet)["gMOTU_140", 4:9] <- c("Diptera",
                                             "Tipulidae",
                                             "Tipula",
                                             "Genus_Tipula",
                                             as.numeric(95.24), "BOLD")

tax_table(final.diet)["gMOTU_144", 4:9] <- c("Geophilomorpha",
                                             "Geophilidae",
                                             "Geophilus",
                                             "Geophilus_easoni",
                                             as.numeric(100.00), "BOLD")

tax_table(final.diet)["zMOTU_111", 3:9] <- c("Insecta", "Diptera",
                                             "Psychodidae",
                                             "Evandromyia",
                                             "Genus_Evandromyia",
                                             as.numeric(97.53), "BOLD")

tax_table(final.diet)["zMOTU_162", 4:9] <- c("Blattodea",
                                             "Ectobiidae",
                                             "Ectobius",
                                             "Ectobius_pallidus",
                                             as.numeric(100.00), "BOLD")

tax_table(final.diet)["zMOTU_173", 4:9] <- c("Lepidoptera",
                                             "Argyresthiidae",
                                             "Family_Argyresthiidae",
                                             "Family_Argyresthiidae",
                                             as.numeric(93.8), "BOLD")

tax_table(final.diet)["zMOTU_178", 4:9] <- c("Diptera",
                                             "Chloropidae",
                                             "Dasyopa",
                                             "Genus_Dasyopa",
                                             as.numeric(95.24), "BOLD")

tax_table(final.diet)["zMOTU_181", 4:9] <- c("Diptera",
                                             "Order_Diptera",
                                             "Order_Diptera",
                                             "Order_Diptera",
                                             as.numeric(92.5), "BOLD")

tax_table(final.diet)["zMOTU_196", 4:9] <- c("Blattodea",
                                             "Ectobiidae",
                                             "Ectobius",
                                             "Genus_Ectobius",
                                             as.numeric(96.83), "BOLD")


#### combining species ####

# This section merges MOTUs from the same species, using the percentage identity value

n <- which(as.numeric(tax_table(final.diet)[,"pident"]) > 98)
m <- taxa_names(final.diet)[n]
p <- taxa_names(final.diet)[-n]

chck <- prune_taxa(m, final.diet)
chck2 <- tax_glom(chck, taxrank="species")
chck2
chck3 <- prune_taxa(p, final.diet)
chck3
chck4 <- merge_phyloseq(chck2, chck3)
chck4

final.diet <- chck4


#### Include Dataset Using Both Primers ####

# Make samples using both primers combined

comb.phy <- final.diet

mrg.pr <- merge_samples(comb.phy, "sample")
sample_data(mrg.pr)$primer <- "Both"
k <- which(sample_data(mrg.pr)[,"mammal"] == 1)
sample_data(mrg.pr)[k,"mammal"] <- "Bat"
k <- which(sample_data(mrg.pr)[,"mammal"] == 2)
sample_data(mrg.pr)[k,"mammal"] <- "GWTS"
k <- which(sample_data(mrg.pr)[,"mammal"] == 3)
sample_data(mrg.pr)[k,"mammal"] <- "Pygmy"

sample_data(mrg.pr)$mammal_primer <- paste(sample_data(mrg.pr)$mammal, sample_data(mrg.pr)$primer, sep = "_")

# Combine phyloseq objects

final <- merge_phyloseq(final.diet, mrg.pr)

final

# Add extra column of data - bat vs shrew (combine GWTS and Pygmy)

sample_data(final)$habitat <- as.factor(c(rep("Air", 22), rep("Terrestrial", 22),
                                          rep("Air", 23), rep("Terrestrial", 15),
                                          rep("Air", 24), rep("Terrestrial", 22)))

#### Save New Phyloseq ####

save(final, file = "./final_phyloseq_object.RData")

#### Load Phyloseq for analyses #####

load("./final_phyloseq_object.RData")

#### load mean merge samples function ####

merge_samples_mean <- function(physeq, group){
  group_sums <- as.matrix(table(sample_data(physeq)[ ,group]))[,1]
  merged <- merge_samples(physeq, group)
  x <- as.matrix(otu_table(merged))
  if(taxa_are_rows(merged)){ x<-t(x) }
  out <- t(x/group_sums)
  out <- otu_table(out, taxa_are_rows = TRUE)
  otu_table(merged) <- out
  return(merged)
}


#### Load Colour Palette ####

pal.o = c("#f0a3ff", # Araneae  # Spare palette from https://graphicdesign.stackexchange.com/questions/3682/where-can-i-find-a-large-palette-set-of-contrasting-colors-for-coloring-many-d
          "#0075dc", # Arch * # Palette is from study to create most contrasting colours
          "#993f00", # Coleoptera
          "#4c005c", # Diptera *
          "#191919", # Ento *
          "#005c31", # Geo *
          "#2bce48", # Glom
          "#ffcc99", # Hap *
          "#808080", # Hemip
          "#94ffb5", # Hymen
          "#8f7c00", # Isopoda
          "#9dcc00", # Julida
          "#c20088", # Lepidop
          "blue", # Lithobiomorpha *
          "#ffa405", # Mesostig
          "#ffa8bb", # Neuropter *
          "#426600", # Oppiliones *
          "#ff0010", # Orbitida *
          "#5ef1f2",# Orthoptera *
          "#00998f", # Poly *
          "#e0ff66", # psocoptera
          "indianred",     # Sarc
          "#003380",    # stylo
          "green", # trom *
          "khaki4",
          "darkred",
          "coral4",
          "violetred2",
          "#0075dc",
          "#993f00")

#### Alpha Diversity Individuals ####

# Agglomerate MOTUs to 'species' level
sep.gn <- tax_glom(final, taxrank = "species")

# Alternative to not agglomerate MOTUs
#sep.gn <- final

min_lib <- min(sample_sums(sep.gn)) # establish smallest sampling depth
nsamp = nsamples(sep.gn)
trials = 100 # how many trials to run
richness <- matrix(nrow = nsamp, ncol = trials)
row.names(richness) <- sample_names(sep.gn)
evenness <- matrix(nrow = nsamp, ncol = trials)
row.names(evenness) <- sample_names(sep.gn)
shannon <- matrix(nrow = nsamp, ncol = trials)
row.names(shannon) <- sample_names(sep.gn)
pielou <- matrix(nrow = nsamp, ncol = trials)
row.names(pielou) <- sample_names(sep.gn)

set.seed(3)

# loop to rerun rarifying and diversity estimates
for (i in 1:trials) {
  #Subsample
  r <- rarefy_even_depth(sep.gn, sample.size = min_lib, verbose = FALSE, replace = TRUE)

  otu.tb <- t(as.matrix(otu_table(r)))

  #Calculate richness
  rich <- as.numeric(specnumber(otu.tb))
  richness[ ,i] <- rich

  #Calculate evenness
  even <- as.numeric(diversity(otu.tb, index = "invsimpson"))
  evenness[ ,i] <- even

  #Calculate shannon
  shan <- as.numeric(diversity(otu.tb, index = "shannon"))
  shannon[ ,i] <- shan

  #Calculate Pielou's Evenness J
  J <- shan / log(rich)
  pielou[ ,i] <- J
}

# Combine diversity estimates with Sample Data
data <- as.data.frame(as.matrix(sample_data(sep.gn)))
Sample <- row.names(richness)
rich.mean <- apply(richness, 1, mean)
rich.sd <- apply(richness, 1, sd)
#measure <- rep("Richness", nsamp)
even.mean <- apply(evenness, 1, mean)
even.sd <- apply(evenness, 1, sd)
shan.mean <- apply(shannon, 1, mean)
shan.sd <- apply(shannon, 1, sd)
piel.mean <- apply(pielou, 1, mean)
piel.sd <- apply(pielou, 1, sd)
#rich_stats <- data.frame(Sample, mean, sd, measure)
rich_stats <- cbind(data, rich.mean, rich.sd, even.mean, even.sd, shan.mean, shan.sd, piel.mean, piel.sd)

rich_stats$primer <- factor(rich_stats$primer,
                            levels = c("Zeale", "Both", "Gillet"))

rich_stats$mammal_primer <- factor(rich_stats$mammal_primer,
                                   levels = c("GWTS_Gillet", "GWTS_Both", "GWTS_Zeale",
                                              "Pygmy_Gillet", "Pygmy_Both", "Pygmy_Zeale",
                                              "Bat_Gillet", "Bat_Both", "Bat_Zeale"))

rich_stats$mammal <- factor(rich_stats$mammal,
                            levels = c("Bat", "Pygmy", "GWTS"),
                            labels = c("R. hipposideros", "S. minutus", "C. russula"))

#c("#1B9E77", "#D95F02", "#7570B3", "#E7298A")

G.1 <- ggplot(rich_stats, aes(mammal_primer, rich.mean, group = mammal_primer)) +
  #geom_boxplot(fill= c("#440154ff"), alpha = 0.25) +
  geom_boxplot(aes(fill = primer), alpha = 0.25) +
  geom_jitter(aes(colour = primer), position = position_jitter(width = 0.1), size=3) +
  stat_summary(fun = mean, geom="point", shape=20, size = 10) +
  coord_flip() +
  scale_colour_manual(values= c("#1B9E77", "#D95F02", "#7570B3")) + # Manually sets colours according to object "palette"
  scale_fill_manual(values= c("#1B9E77", "#D95F02", "#7570B3")) + # Manually sets colours according to object "palette"
  #scale_x_discrete(labels= c("Inside", "Middle", "Edge")) + # Manually labels groups using set described in "key" object
  # scale_y_continuous(limits = c(0, 0.55), expand = c(0, 0)) + # Sets Graph Value Limits
  theme_classic() + # Chooses Theme
  theme(axis.text = element_text(size = 10, face= "bold"), # Size of Axis measurements/Values Text
        axis.title=element_text(size=14), # Sets size of Axis Titles. face= "bold" will change the font to bold
        panel.grid.major = element_line(colour = "grey95"), # Sets colour of major grid lines
        panel.grid.minor = element_line(colour = "grey95") # Sets colour of minor grid lines
  ) +
  labs(y= "Species Richness", x = "") + # Labels Axis'
  #facet_wrap(~mammal, scale = "free_x") +
  facet_wrap(~ mammal, scale = "free_y", ncol = 1) +
  theme(strip.text.x = element_blank()) +
  theme(legend.position="none") +
  ylim(c(0,40))

G.1

G.2 <- ggplot(rich_stats, aes(mammal_primer, shan.mean, group = mammal_primer)) +
  #geom_boxplot(fill= c("#440154ff"), alpha = 0.25) +
  geom_boxplot(aes(fill = primer), alpha = 0.25) +
  geom_jitter(aes(colour = primer), position = position_jitter(width = 0.1), size=3) +
  stat_summary(fun = mean, geom="point", shape=20, size = 10) +
  coord_flip() +
  scale_colour_manual(values= c("#1B9E77", "#D95F02", "#7570B3")) + # Manually sets colours according to object "palette"
  scale_fill_manual(values= c("#1B9E77", "#D95F02", "#7570B3")) + # Manually sets colours according to object "palette"
  #scale_x_discrete(labels= c("Inside", "Middle", "Edge")) + # Manually labels groups using set described in "key" object
  # scale_y_continuous(limits = c(0, 0.55), expand = c(0, 0)) + # Sets Graph Value Limits
  theme_classic() + # Chooses Theme
  theme(axis.text = element_text(size = 10, face= "bold"), # Size of Axis measurements/Values Text
        axis.title=element_text(size=14), # Sets size of Axis Titles. face= "bold" will change the font to bold
        panel.grid.major = element_line(colour = "grey95"), # Sets colour of major grid lines
        panel.grid.minor = element_line(colour = "grey95") # Sets colour of minor grid lines
  ) +
  labs(y= "Shannon Diversity", x = "") + # Labels Axis'
  #facet_wrap(~mammal, scale = "free_x") +
  facet_wrap(~ mammal, scale = "free_y", ncol = 1) +
  theme(strip.text.x = element_blank()) +
  theme(legend.position="none")

G.2

grid.arrange(G.1, G.2 + theme(axis.text.y = element_blank()), nrow=1)

legd <- cowplot::get_legend(G.2 + theme(legend.position = "right",
                                        legend.title = element_blank(),
                                        legend.text = element_text(size = 20),
                                        legend.key.size = unit(1.5, "cm")))
grid.newpage()
grid.draw(legd)

##### Alpha Diversity By Group ####

#sep.gn <- tax_glom(final, taxrank = "species")
#mrg.pr = transform_sample_counts(sep.gn, function(x) 100 * x/sum(x))
#mrg.pr <- merge_samples_mean(mrg.pr, "mammal_primer")
mrg.pr <- merge_samples(sep.gn, "mammal_primer")
sample_data(mrg.pr)$mammal_primer <- rownames(sample_data(mrg.pr))
sample_data(mrg.pr)$mammal <- c(rep("Bat", 3), rep("GWTS", 3), rep("Pygmy", 3))
sample_data(mrg.pr)$primer <- rep(c("Both", "Gillet", "Zeale"), 3)

min_lib <- min(sample_sums(mrg.pr))
mrg.pr <- rarefy_even_depth(mrg.pr, sample.size = min_lib, verbose = FALSE, replace = TRUE)

all.otb <- t(as.data.frame(as.matrix(otu_table(mrg.pr))))
#all.otb <- all.otb[which(apply(all.otb,1,sum) > 0),]

ovrlp <- matrix(nrow = 5, ncol = ncol(all.otb))
ovrlp[1,] <- specnumber(t(all.otb))
ovrlp[2,] <- diversity(t(all.otb), index = "shannon")
for(i in 1:ncol(all.otb)){
  ovrlp[3,i] <- ovrlp[2,i] / log(ovrlp[1,i])
}
ovrlp[4,] <- as.numeric(niche.width(all.otb, method = "levins"))
for (i in 1:ncol(all.otb)){
  B <- ovrlp[4,i]
  n <- as.numeric(length(which(all.otb[,i] > 0)))
  Ba <- (B - 1) / (n - 1)
  ovrlp[5,i] <- Ba
}
rownames(ovrlp) <- c("Richness", "Shannon", "Pielou's", "levins", "stnd.levins")
colnames(ovrlp) <- colnames(all.otb)
ovrlp


#### Compare Alpha Diversity ####

# turn off e-values
options(scipen = 999)

### Kruskal Wallice

kruskal.test(rich.mean ~ mammal_primer, data = rich_stats)
kruskal.test(shan.mean ~ mammal_primer, data = rich_stats)

### Pairwise Wilcoxin Test

pairwise.wilcox.test(rich_stats$shan.mean,
                     rich_stats$mammal_primer,
                     p.adjust.method = "BH")

### Dunn Test - Pairwise comprison

dn.df <- dunnTest(rich_stats$shan.mean ~ rich_stats$mammal_primer,
                  method = "hochberg")

dn.df$res[dn.df$res[,4] < 0.05,]

### ANOVA and Tukey pairwise comparison

# ANOVA using the shannon diversity
a1 <- aov(rich_stats$shan.mean ~ rich_stats$mammal_primer)
tk.df <- TukeyHSD(x=a1, 'rich_stats$mammal_primer', conf.level=0.95)
tk.df$`rich_stats$mammal_primer`[tk.df$`rich_stats$mammal_primer`[,4] < 0.05,]

par(mar = c(3, 14, 2, 1))
plot(TukeyHSD(x=a1, 'rich_stats$mammal_primer', conf.level=0.95), las = 1)

# ANOVA using the species richness
a1 <- aov(rich_stats$rich.mean ~ rich_stats$mammal_primer)
tk.df <- TukeyHSD(x=a1, 'rich_stats$mammal_primer', conf.level=0.95)
tk.df$`rich_stats$mammal_primer`[tk.df$`rich_stats$mammal_primer`[,4] < 0.05,]
plot(TukeyHSD(x=a1, 'rich_stats$mammal_primer', conf.level=0.95), las = 1)

# return margins to default
par(mar=c(5.1, 4.1, 4.1, 2.1))

#### Agglom at Species Level ####

sep.gn <- tax_glom(final, taxrank= "species")

#### Agglom at genus level ####

sep.gn <- tax_glom(final, taxrank= "genus")

#### Agglom at family level ####

sep.gn <- tax_glom(final, taxrank= "family")

#### Agglom at order level ####

sep.gn <- tax_glom(final, taxrank= "order")


#### PERMANOVAs Zeale, Gillet and Both ####

# Make sure data has been agglomerated to species level
sep.gn <- tax_glom(final, taxrank= "species")

### Adonis

sep.ra = transform_sample_counts(sep.gn, function(x) 100 * x/sum(x))
otu <- abundances(sep.ra)
meta <- meta(sep.ra)

set.seed(59009)
permanova <- adonis2(t(otu) ~ mammal * primer,
                     data = meta,
                     permutations=10000,
                     method = "bray")
#method = "jaccard")

permanova

# Which taxa contributing most to split of beta diversity between mammals

coef <- coefficients(permanova)["mammal1",]
top.coef <- coef[rev(order(abs(coef)))[1:20]]
n <- names(top.coef)
m <- tax_table(sep.ra)[n,"order"]
#k <- tax_table(mrg.ra)[n,"genus"]
#names(top.coef) <- paste(m, k, sep = " ") # playing with labels
names(top.coef) <- m
par(mar = c(3, 14, 2, 1))
barplot(sort(top.coef), horiz = T, las = 1,
        main = "Top 20 taxa\nseperating Mammals")


## Pairwise Adonis

# Function
# courtesy of https://www.researchgate.net/post/How_can_I_do_PerMANOVA_pairwise_contrasts_in_R
pairwise.adonis <- function(x,factors, sim.function = 'vegdist', sim.method = 'bray', p.adjust.m ='bonferroni')
{
  library(vegan)

  co = combn(unique(as.character(factors)),2)
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()


  for(elem in 1:ncol(co)){
    if(sim.function == 'daisy'){
      library(cluster); x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
    } else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}

    ad = adonis(x1 ~ factors[factors %in% c(co[1,elem],co[2,elem])] );
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  sig = c(rep('',length(p.adjusted)))
  sig[p.adjusted <= 0.05] <-'.'
  sig[p.adjusted <= 0.01] <-'*'
  sig[p.adjusted <= 0.001] <-'**'
  sig[p.adjusted <= 0.0001] <-'***'

  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted,sig)
  print("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1")
  return(pairw.res)

}

pairwise.adonis(t(otu), meta$mammal_primer)



#### betadisperser Zeale, Gillet and Both ####

# Make sure data has been agglomerated to species level
sep.gn <- tax_glom(final, taxrank= "species")

sep.ra = transform_sample_counts(sep.gn, function(x) 100 * x/sum(x))
otu <- abundances(sep.ra)
meta <- meta(sep.ra)

options(scipen = 999)
dist_matrix <- phyloseq::distance(sep.ra,
                                  method = "bray")
#method = "UniFrac") # requires phylogenetic information
mod <- betadisper(dist_matrix, meta$mammal_primer)
boxplot(mod, main = "Dispersion Between Shrews")
anova(mod)
permutest(mod, pairwise = TRUE)
TukeyHSD(mod)
par(mar = c(5.1, 12, 2, 1)) # bottom, left, top, right
plot(TukeyHSD(mod), las = 1)
par(mar=c(5.1, 4.1, 4.1, 2.1))

dis.stat <- data.frame(group = mod$group, dist = mod$distances)

G.1 <- ggplot(dis.stat, aes(group, dist)) +
  geom_boxplot(aes(fill= group), alpha = 0.35) +
  geom_jitter(aes(colour = group), position = position_jitter(width = 0.1), size=3) +
  stat_summary(fun = mean, geom="point", shape=20, size = 10) +
  #scale_colour_manual(values= c("#440154ff", '#21908dff', "indianred", "#fde725ff")) +
  #scale_x_discrete(labels= c("C. russla\nGillet",
  #                           "C. russla\nZeale",
  #                           "S. minutus\nGillet",
  #                           "C. russla\nZeale")) + # Manually labels groups using set described in "key" object
  # scale_y_continuous(limits = c(0, 0.55), expand = c(0, 0)) + # Sets Graph Value Limits
  theme_classic() + # Chooses Theme
  theme(axis.text = element_text(size = 10, face= "bold"), # Size of Axis measurements/Values Text
        axis.title=element_text(size=14), # Sets size of Axis Titles. face= "bold" will change the font to bold
        panel.grid.major = element_line(colour = "grey95"), # Sets colour of major grid lines
        panel.grid.minor = element_line(colour = "grey95") # Sets colour of minor grid lines
  ) +
  labs(y= "Level of Dispersion", x = "") + # Labels Axis'
  #facet_wrap(~Shrew, scale = "free_x") +
  theme(legend.position="none")

G.1

var(dis.stat[dis.stat[,1] == "Bat_Zeale",2])
var(dis.stat[dis.stat[,1] == "Bat_Gillet",2])
var(dis.stat[dis.stat[,1] == "Bat_Both",2])
var(dis.stat[dis.stat[,1] == "GWTS_Zeale",2])
var(dis.stat[dis.stat[,1] == "GWTS_Gillet",2])
var(dis.stat[dis.stat[,1] == "GWTS_Both",2])
var(dis.stat[dis.stat[,1] == "Pygmy_Zeale",2])
var(dis.stat[dis.stat[,1] == "Pygmy_Gillet",2])
var(dis.stat[dis.stat[,1] == "Pygmy_Both",2])
sd(dis.stat[dis.stat[,1] == "Bat_Zeale",2])
sd(dis.stat[dis.stat[,1] == "Bat_Gillet",2])
sd(dis.stat[dis.stat[,1] == "Bat_Both",2])
sd(dis.stat[dis.stat[,1] == "GWTS_Zeale",2])
sd(dis.stat[dis.stat[,1] == "GWTS_Gillet",2])
sd(dis.stat[dis.stat[,1] == "GWTS_Both",2])
sd(dis.stat[dis.stat[,1] == "Pygmy_Zeale",2])
sd(dis.stat[dis.stat[,1] == "Pygmy_Gillet",2])
sd(dis.stat[dis.stat[,1] == "Pygmy_Both",2])



#### NMDS Zeale, Gillet and Both ####

sep.ra = transform_sample_counts(sep.gn, function(x) 100 * x/sum(x))
set.seed(3847)
ord_diet <- phyloseq::ordinate(sep.ra,
                               method = "NMDS",
                               distance = "bray",
                               k =6)

scrs <- scores(ord_diet, display='sites')
smpl.df <- as.data.frame(as.matrix(sample_data(sep.ra)))
scrs <- cbind(as.data.frame(scrs), smpl.df)
cent <- aggregate(cbind(NMDS1, NMDS2) ~ mammal_primer, data = scrs, FUN = "mean")
segs <- merge(scrs, setNames(cent, c("mammal_primer","oNMDS1","oNMDS2")),
              by = "mammal_primer", sort = FALSE)


pal_3 <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "black", "#A65628", "#F781BF", "#999999")

p.1.3 <- ggplot(scrs, aes(x = NMDS1, y = NMDS2, colour = mammal_primer)) +
  geom_segment(data = segs, alpha = 0.5,
               mapping = aes(xend = oNMDS1, yend = oNMDS2)) +
  scale_color_manual(values = pal_3) +
  geom_point(data = cent, aes(colour = mammal_primer), size = 5) +
  geom_point(data = scrs, aes(colour = mammal_primer), size=3, alpha=0.9) +
  theme_bw() +
  theme(axis.text = element_text(size = 14, face= "bold", colour = "black"),
        axis.title=element_text(size=18),
        legend.title= element_text(size= 18),
        legend.text = element_text(size = 14),
        legend.key = element_rect(fill = "white", colour = NA),
        legend.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "grey95"),
        panel.grid.minor = element_line(colour = "grey95")
  ) +
  guides(colour = guide_legend(override.aes = list(size=5)),
         colour = guide_legend("mammal_primer"), fill = FALSE) +
  theme(legend.position = "none")

p.1.3

legd <- cowplot::get_legend(p.1.3 + theme(legend.position = "bottom",
                                          legend.title = element_blank(),
                                          legend.text = element_text(size = 20),
                                          legend.key.size = unit(1.5, "cm")))
grid.newpage()
grid.draw(legd)


#### Barplots ####

top25ord <- final

# isolate taxa that we don't want to group as 'other'
othr.ord <- subset_taxa(final, !(order %in% c("Araneae", "Archaeognatha", "Coleoptera", "Diptera", "Entomobryomorpha",
                                              "Glomerida", "Haplotaxida", "Hymenoptera", "Isopoda", "Julida",
                                              "Lepidoptera", "Mesostigmata", "Neuroptera", "Opiliones", "Polydesmida",
                                              "Stylommatophora", "Trichoptera", "Trombidiformes")))

n <- taxa_names(othr.ord)
tax_table(top25ord)[n, 1:4] <- "Other"

top25ord <- tax_glom(top25ord, taxrank = "order")

ra.samples = transform_sample_counts(top25ord, function(x) 100 * x/sum(x))
ra.samples.bar <- phyloseq::plot_bar(ra.samples) # extracts information needed for barplots
ra.samples.bar.data <- ra.samples.bar$data
ra.samples.bar.data$primer <- factor(ra.samples.bar.data$primer,
                                     levels = c("Zeale", "Gillet", "Both"))
p1.1 <- ggplot(ra.samples.bar.data, aes(x= Sample, y=Abundance, fill = order)) +
  geom_bar(stat="identity", color="black") +
  scale_fill_manual(values = pal.o) +
  facet_wrap(~ mammal + primer, scale = "free_x") +
  theme_classic() + # appearance composition
  theme(strip.text.x = element_text(size = 15, face = "bold")) + # facet text
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.key.size = unit(1, "cm")) +
  labs(y = "Relative Abundance") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"))

p1.1

# Merged Sample Barplots

sampl.filt <- transform_sample_counts(top25ord, function(x) 100 * x/sum(x))
mrg.pr <- merge_samples_mean(sampl.filt, "mammal_primer")

sample_data(mrg.pr)$mammal_primer <- rownames(sample_data(mrg.pr))
sample_data(mrg.pr)$mammal <- c(rep("Bat", 3), rep("GWTS", 3), rep("Pygmy", 3))
sample_data(mrg.pr)$primer <- rep(c("Both", "Gillet", "Zeale"), 3)

ra.samples.bar <- phyloseq::plot_bar(mrg.pr) # extracts information needed for barplots
ra.samples.bar.data <- ra.samples.bar$data
ra.samples.bar.data$primer <- factor(ra.samples.bar.data$primer,
                                     levels = c("Zeale", "Gillet", "Both"))
p1.2 <- ggplot(ra.samples.bar.data, aes(x= Sample, y=Abundance, fill = order)) +
  geom_bar(stat="identity", color="black") +
  scale_fill_manual(values = pal.o) +
  facet_wrap(~ primer, scale = "free_x") +
  theme_classic() + # appearance composition
  theme(strip.text.x = element_text(size = 15, face = "bold")) + # facet text
  scale_x_discrete(labels=c("GWTS_Gillet" = "C.russula\nN=7",
                            "Pygmy_Gillet" = "S.minutus\nN=15",
                            "Bat_Gillet" = "H. hipposideros\nN=22",
                            "GWTS_Zeale" = "C.russula\nN=4",
                            "Pygmy_Zeale" = "S.minutus\nN=11",
                            "Bat_Zeale" = "H. hipposideros\nN=23",
                            "Pygmy_Both" = "S.minutus\nN=15",
                            "GWTS_Both" = "C.russula\nN=7",
                            "Bat_Both" = "H. hipposideros\nN=23")) +
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.8, "cm")) +
  labs(y = "Relative Abundance (%)") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20,
                                    face = "bold")) +
  theme(axis.text.x = element_text(size = 8 ,
                                   face = "bold",
                                   colour = "black",
                                   angle = 0,
                                   hjust = 0.5),
        axis.text.y = element_text(size = 12,
                                   face = "bold"))

p1.2

# Horizontal barplot
ra.samples.bar.data$Sample <- factor(ra.samples.bar.data$Sample,
                                     levels = c("GWTS_Gillet", "GWTS_Both", "GWTS_Zeale",
                                                "Pygmy_Gillet", "Pygmy_Both", "Pygmy_Zeale",
                                                "Bat_Gillet", "Bat_Both", "Bat_Zeale"))

ra.samples.bar.data$mammal <- factor(ra.samples.bar.data$mammal,
                                     levels = c("Bat", "Pygmy", "GWTS"),
                                     labels = c("H. hipposideros", "S. minutus", "C. russula"))


p1.2 <- ggplot(ra.samples.bar.data, aes(x= Sample, y=Abundance, fill = order)) +
  geom_bar(stat="identity", color="black") +
  scale_fill_manual(values = pal.o) +
  #facet_wrap(~ mammal, scale = "free", ncol = 1, switch = "y") +
  facet_wrap(~ mammal, scale = "free_y", ncol = 1, strip.position = "left") +
  coord_flip() +
  theme_classic() + # appearance composition
  theme(strip.text.x = element_text(size = 15, face = "bold")) + # facet text
  scale_x_discrete(labels=c("GWTS_Gillet" = "Gillet",
                            "Pygmy_Gillet" = "Gillet",
                            "Bat_Gillet" = "Gillet",
                            "GWTS_Zeale" = "Zeale",
                            "Pygmy_Zeale" = "Zeale",
                            "Bat_Zeale" = "Zeale",
                            "Pygmy_Both" = "Both",
                            "GWTS_Both" = "Both",
                            "Bat_Both" = "Both")) +
  theme(strip.text.y = element_text(size = 12, color = "black", face = "bold.italic")) +
  theme(legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.8, "cm")) +
  labs(y = "Relative Abundance (%)") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20,
                                    face = "bold")) +
  theme(axis.text.x = element_text(size = 8 ,
                                   #face = "bold",
                                   colour = "black",
                                   angle = 0,
                                   hjust = 0.5),
        axis.text.y = element_text(size = 12,
                                   face = "bold"))

p1.2

legd <- cowplot::get_legend(p1.2 + theme(legend.position = "bottom") + guides(fill = guide_legend(nrow = 2)))
grid.newpage()
grid.draw(legd)


#### hclust ####

fin.cr <- tax_glom(final, taxrank = "species")
fin.cr = transform_sample_counts(fin.cr, function(x) 100 * x/sum(x))
mrg.pr <- merge_samples_mean(fin.cr, "mammal_primer")
sample_data(mrg.pr)$mammal_primer <- rownames(sample_data(mrg.pr))
sample_data(mrg.pr)$mammal <- c(rep("Bat", 3), rep("GWTS", 3), rep("Pygmy", 3))
sample_data(mrg.pr)$primer <- rep(c("Both", "Gillet", "Zeale"), 3)

chck <- otu_table(mrg.pr)
dist.mat <- distance(chck,
                     method = "bray")

hier.cl <- hclust(dist.mat,
                  method = "average")

hier.cl.den <- as.dendrogram(hier.cl)

par(mar = c(5.1, 2, 2, 6))

hier.cl.den %>% rotate(rev(c(3, 2, 1, 4, 5, 6, 7, 8, 9))) %>%  plot(horiz = TRUE)

par(mar=c(5.1, 4.1, 4.1, 2.1))

#### Venn Diagram ####

library(VennDiagram)

sep.gn <- tax_glom(final, taxrank= "species")

sep.gn.shrew <- sep.gn
sample_data(sep.gn.shrew)$mammal_primer <- str_replace(sample_data(sep.gn.shrew)$mammal_primer, "Pygmy_Zeale", "Shrew_Zeale")
sample_data(sep.gn.shrew)$mammal_primer <- str_replace(sample_data(sep.gn.shrew)$mammal_primer, "Pygmy_Gillet", "Shrew_Gillet")
sample_data(sep.gn.shrew)$mammal_primer <- str_replace(sample_data(sep.gn.shrew)$mammal_primer, "GWTS_Zeale", "Shrew_Zeale")
sample_data(sep.gn.shrew)$mammal_primer <- str_replace(sample_data(sep.gn.shrew)$mammal_primer, "GWTS_Gillet", "Shrew_Gillet")

mrg.pr <- merge_samples(sep.gn.shrew, "mammal_primer")

sample_data(mrg.pr)$mammal_primer <- rownames(sample_data(mrg.pr))
#mrg.pr <- tax_glom(mrg.pr, taxrank = "species")

df <- as.data.frame(otu_table(mrg.pr))
colnames(df) <- tax_table(mrg.pr)[,"species"]

df <- df[,!grepl("Class_|Order_|Family_|Genus_|none|_sp._|_cf._", colnames(df))]

Shrew_Zeale <- colnames(df["Shrew_Zeale", apply(df["Shrew_Zeale",], MARGIN=2, function(x) any(x >0))])
Shrew_Gillet <- colnames(df["Shrew_Gillet", apply(df["Shrew_Gillet",], MARGIN=2, function(x) any(x >0))])
Bat_Zeale <- colnames(df["Bat_Zeale", apply(df["Bat_Zeale",], MARGIN=2, function(x) any(x >0))])
Bat_Gillet <- colnames(df["Bat_Gillet", apply(df["Bat_Gillet",], MARGIN=2, function(x) any(x >0))])

# number of species detected by Zeale
length(unique(c(Shrew_Zeale, Bat_Zeale)))
# number of species detected by Zeale in bats
length(unique(Bat_Zeale))
# number of species detected by Zeale in shrews
length(unique(Shrew_Zeale))
# number of species detected by Gillet
length(unique(c(Shrew_Gillet, Bat_Gillet)))
# number of species detected by Gillet in bats
length(unique(Bat_Gillet))
# number of species detected by Gillet in shrews
length(unique(Shrew_Gillet))

fig.ire <- venn.diagram(
  x = list(Shrew_Zeale, Bat_Zeale, Shrew_Gillet, Bat_Gillet),
  category.names = c("Shrew\nZeale", "Bat\nZeale","Shrew Gillet", "Bat Gillet"),
  filename = NULL,
  cex = 3,
  cat.cex = 2,
  alpha = 0.5,
  col= c("#1B9E77",
         "#1B9E77",
         "#1B9E77",
         "#1B9E77"),
  fill = c("#1B9E77",
           "#1B9E77",
           "#1B9E77",
           "#1B9E77"))



grid.newpage()
grid.draw(fig.ire)

mrg.pr.gn <- tax_glom(mrg.pr, taxrank = "genus")

df <- as.data.frame(otu_table(mrg.pr.gn))
colnames(df) <- tax_table(mrg.pr.gn)[,"genus"]

df <- df[,!grepl("Class_|Order_|Family_|Genus_|none|_sp._|_cf._", colnames(df))]

Shrew_Zeale <- colnames(df["Shrew_Zeale", apply(df["Shrew_Zeale",], MARGIN=2, function(x) any(x >0))])
Shrew_Gillet <- colnames(df["Shrew_Gillet", apply(df["Shrew_Gillet",], MARGIN=2, function(x) any(x >0))])
Bat_Zeale <- colnames(df["Bat_Zeale", apply(df["Bat_Zeale",], MARGIN=2, function(x) any(x >0))])
Bat_Gillet <- colnames(df["Bat_Gillet", apply(df["Bat_Gillet",], MARGIN=2, function(x) any(x >0))])

# number of genera detected by Zeale
length(unique(c(Shrew_Zeale, Bat_Zeale)))
# number of genera detected by Zeale in bats
length(unique(Bat_Zeale))
# number of genera detected by Zeale in shrews
length(unique(Shrew_Zeale))
# number of genera detected by Gillet
length(unique(c(Shrew_Gillet, Bat_Gillet)))
# number of genera detected by Gillet in bats
length(unique(Bat_Gillet))
# number of genera detected by Gillet in shrews
length(unique(Shrew_Gillet))

fig.ire <- venn.diagram(
  x = list(Shrew_Zeale, Bat_Zeale, Shrew_Gillet, Bat_Gillet),
  category.names = c("Shrew\nZeale", "Bat\nZeale","Shrew Gillet", "Bat Gillet"),
  filename = NULL,
  cex = 3,
  cat.cex = 2,
  alpha = 0.5,
  col= c("#D95F02",
         "#D95F02",
         "#D95F02",
         "#D95F02"),
  fill = c("#D95F02",
           "#D95F02",
           "#D95F02",
           "#D95F02"))



grid.newpage()
grid.draw(fig.ire)

mrg.pr.fm <- tax_glom(mrg.pr, taxrank = "family")

df <- as.data.frame(otu_table(mrg.pr.fm))
colnames(df) <- tax_table(mrg.pr.fm)[,"family"]

df <- df[,!grepl("Class_|Order_|Family_|Genus_|none|_sp._|_cf._", colnames(df))]

Shrew_Zeale <- colnames(df["Shrew_Zeale", apply(df["Shrew_Zeale",], MARGIN=2, function(x) any(x >0))])
Shrew_Gillet <- colnames(df["Shrew_Gillet", apply(df["Shrew_Gillet",], MARGIN=2, function(x) any(x >0))])
Bat_Zeale <- colnames(df["Bat_Zeale", apply(df["Bat_Zeale",], MARGIN=2, function(x) any(x >0))])
Bat_Gillet <- colnames(df["Bat_Gillet", apply(df["Bat_Gillet",], MARGIN=2, function(x) any(x >0))])

# number of familes detected by Zeale
length(unique(c(Shrew_Zeale, Bat_Zeale)))
# number of familes detected by Zeale in bats
length(unique(Bat_Zeale))
# number of familes detected by Zeale in shrews
length(unique(Shrew_Zeale))
# number of families detected by Gillet
length(unique(c(Shrew_Gillet, Bat_Gillet)))
# number of families detected by Gillet in bats
length(unique(Bat_Gillet))
# number of families detected by Gillet in shrews
length(unique(Shrew_Gillet))

fig.ire <- venn.diagram(
  x = list(Shrew_Zeale, Bat_Zeale, Shrew_Gillet, Bat_Gillet),
  category.names = c("Shrew\nZeale", "Bat\nZeale","Shrew Gillet", "Bat Gillet"),
  filename = NULL,
  cex = 3,
  cat.cex = 2,
  alpha = 0.5,
  col= c("#7570B3",
         "#7570B3",
         "#7570B3",
         "#7570B3"),
  fill = c("#7570B3",
           "#7570B3",
           "#7570B3",
           "#7570B3"))

grid.newpage()
grid.draw(fig.ire)

mrg.pr.or <- tax_glom(mrg.pr, taxrank = "order")

df <- as.data.frame(otu_table(mrg.pr.or))
colnames(df) <- tax_table(mrg.pr.or)[,"order"]

df <- df[,!grepl("Class_|Order_|Family_|Genus_|none|_sp._|_cf._", colnames(df))]

Shrew_Zeale <- colnames(df["Shrew_Zeale", apply(df["Shrew_Zeale",], MARGIN=2, function(x) any(x >0))])
Shrew_Gillet <- colnames(df["Shrew_Gillet", apply(df["Shrew_Gillet",], MARGIN=2, function(x) any(x >0))])
Bat_Zeale <- colnames(df["Bat_Zeale", apply(df["Bat_Zeale",], MARGIN=2, function(x) any(x >0))])
Bat_Gillet <- colnames(df["Bat_Gillet", apply(df["Bat_Gillet",], MARGIN=2, function(x) any(x >0))])

# number of orders detected by Zeale
length(unique(c(Shrew_Zeale, Bat_Zeale)))
# number of orders detected by Zeale in bats
length(unique(Bat_Zeale))
# number of orders detected by Zeale in shrews
length(unique(Shrew_Zeale))
# number of orders detected by Gillet
length(unique(c(Shrew_Gillet, Bat_Gillet)))
# number of orders detected by Gillet in bats
length(unique(Bat_Gillet))
# number of orders detected by Gillet in shrews
length(unique(Shrew_Gillet))

fig.ire <- venn.diagram(
  x = list(Shrew_Zeale, Bat_Zeale, Shrew_Gillet, Bat_Gillet),
  category.names = c("Shrew\nZeale", "Bat\nZeale","Shrew Gillet", "Bat Gillet"),
  filename = NULL,
  cex = 3,
  cat.cex = 2,
  alpha = 0.5,
  col= c("#E7298A",
         "#E7298A",
         "#E7298A",
         "#E7298A"),
  fill = c("#E7298A",
           "#E7298A",
           "#E7298A",
           "#E7298A"))



grid.newpage()
grid.draw(fig.ire)

#### Random Forrest Classifier ####

# https://rpubs.com/michberr/randomforestmicrobe

# A more robust method to determine what invert taxa are responsible for PERMANOVA differences (here look at primers)

library(randomForest)

sep.gn <- final
sep.gn <- subset_samples(final, mammal == "Bat") # what mammal to include
#sep.gn <- subset_samples(final, habitat == "Terrestrial")
sep.gn <- subset_samples(sep.gn, primer == "Both") # which primer set to include
sep.gn <- tax_glom(sep.gn, taxrank= "species")
sep.ra = transform_sample_counts(sep.gn, function(x) 100 * x/sum(x))

# Make a dataframe of training data with OTUs as column and samples as rows
predictors <- t(otu_table(sep.ra))
dim(predictors)

# Make one column for our outcome/response variable
response <- as.factor(sample_data(sep.ra)$mammal)

# Combine them into 1 data frame
rf.data <- data.frame(response, predictors)

set.seed(2)
rn.for.classify <- randomForest(response ~ ., data = rf.data, ntree = 10000)
print(rn.for.classify)
names(rn.for.classify)

# Make a data frame with predictor names and their importance
imp <- importance(rn.for.classify)
imp <- data.frame(predictors = rownames(imp), imp)
#imp$predictors <- paste(tax_table(sep.ra)[,"family"], imp$predictors, sep = ": ")

# Order the predictor levels by importance
imp.sort <- arrange(imp, desc(MeanDecreaseGini))
#imp.sort$predictors <- str_replace(imp.sort$predictors, "SV.", "SV ")
imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)

# Select the top 20 predictors
imp.20 <- imp.sort[1:20, ]

# ggplot
ggplot(imp.20, aes(x = predictors, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "indianred") +
  coord_flip() +
  ggtitle("Most important OTUs for classifying... ")

# What are those OTUs?
#otunames <- imp.20$predictors
#r <- rownames(tax_table(sep.gn)) %in% otunames
#tax_table(sep.gn)[r,2:7]

tax_table(sep.gn)[as.character(imp.20$predictors),2:7]

sep.gn <- final
sep.gn <- subset_samples(sep.gn, primer != "Both")
sep.gn <- tax_glom(sep.gn, taxrank= "species")

mrg.pr = transform_sample_counts(sep.gn, function(x) 100 * x/sum(x))
mrg.pr <- merge_samples_mean(mrg.pr, "primer")
sample_data(mrg.pr)$primer <- rownames(sample_data(mrg.pr))
taxa_sums(mrg.pr)[as.character(imp.20$predictors)]
otu_table(mrg.pr)[as.character(imp.20$predictors),]

imp.20 <- cbind(imp.20, otu_table(mrg.pr)[as.character(imp.20$predictors),])

mrg.pr = transform_sample_counts(sep.gn, function(x) 100 * x/sum(x))
mrg.pr <- merge_samples_mean(mrg.pr, "mammal")
sample_data(mrg.pr)$mammal <- rownames(sample_data(mrg.pr))
#taxa_sums(mrg.pr)[as.character(imp.20$predictors)]
otu_table(mrg.pr)[as.character(imp.20$predictors),]

imp.20 <- cbind(imp.20, otu_table(mrg.pr)[as.character(imp.20$predictors),])
