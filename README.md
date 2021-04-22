These are the scripts used to run the analyses and generate figures for the following publication:
Browett S.S., Curran T.G., O’Meara D.B., Harrington A.P., Sales N.G., Antwis R.E. O'Neill D. & McDevitt A.D. Primer biases in the molecular assessment of diet in multiple insectivorous mammals. Mamm Biol (2021). https://doi.org/10.1007/s42991-021-00115-4

The study uses DNA metabarcoding to determine the diet of 3 insectivorous species: Rhinolophus hipposideros (lesser horseshoe bat), Crocidura russula (greater white-toothed shrew) and Sorex minutus (pygmy shrew). The paper evaluates the performance of two popular primers to show that empirical tesing of primers is required for determining the diet of novel species. The primers evaluated here are the Zeale primers (Zeale MRK, Butlin RK, Barker GLA, Lees DC, Jones G (2011) Taxon- specific PCR for DNA barcoding arthropod prey in bat faeces. Mol Ecol Resour 11(2):236–244.) and the modified Gillet primers (Esnaola A, Arrizabalaga-Escudero A, González-Esteban J, Elosegi A, Aihartza J (2018) Determining diet from faeces: selection of metabarcoding primers for the insectivore Pyrenean desman (Galemys pyrenaicus). PLoS ONE 13(12):e0208986.)

Abstract:
Our understanding of trophic interactions of small insectivorous mammals has been drastically improved with the advent of DNA metabarcoding. The technique has continued to be optimised over the years, with primer choice repeatedly being a vital factor for dietary inferences. However, the majority of dietary studies examining the effect of primer choice often rely on in silico analyses or comparing between species that occupy an identical niche type. Here, we apply DNA metabarcoding to empirically compare the prey detection capabilities of two widely used primer sets when assessing the diets of a flying (lesser horseshoe bat; Rhinolophus hipposideros) and two ground-dwelling insectivores (greater white-toothed shrew; Cro- cidura russula and pygmy shrew; Sorex minutus). Although R. hipposideros primarily rely on two prey orders (Lepidoptera and Diptera), the unique taxa detected by each primer shows that a combination of primers may be the best approach to fully describe bat trophic ecology. However, random forest classifier analysis suggests that one highly degenerate primer set detected the majority of both shrews’ diet despite higher levels of host amplification. The wide range of prey consumed by ground-dwelling insectivores can therefore be accurately documented from using a single broad-range primer set, which can decrease cost and labour. The results presented here show that dietary inferences will differ depending on the primer or primer combination used for insectivores occupying different niches (i.e., hunting in the air or ground) and demonstrate the importance of performing empirical pilot studies for novel study systems.


Raw sequence data are available on DRYAD digital repository (https://doi.org/10.5061/dryad.7d7wm37ts).

Raw sequences were generated on an Illumina MiSeq, V2 300 cycle kit. Samples are multiplexed with MID tags on the primers. See full details on library prep in manuscript. 

Each step from processing FASTQ files to full analyses are included here (Steps 1 - 5). The other files are:

Bat_ngsfilter.tsv: tab seperated value file used to demultiplex bat samples during obitools pipeline
Shrew_ngsfilter.tsv: tab seperated value file used to demultiplex shrew samples during obitools pipeline
final_phyloseq_object.RData: The filtered dataset formatted as a phyloseq object in R. This is the file used for analyses.
merged_primer_dataset.RData: This is the phyloseq object after low read MOTUs were removed and primers and samples were merged into a single phyloseq object.
phyloseq_object_clust_iden98_taxa_assigned_no_singletons.RData: The phyloseq object with taxa detected using the Gillet primers. MOTUs have been assigned to taxa, but MOTUs have not yet been filtered according to read counts.
phyloseq_object_zeale_clust_iden98_taxa_assigned_no_singletons.RData: The phyloseq object with taxa detected using the Zeale primers. MOTUs have been assigned to taxa, but MOTUs have not yet been filtered according to read counts.
samplesheet.csv: This file contains the details of each sample.
species_identified_in_pygmy.csv: The MOTUs identified to species level in pygmy shrews in this study
species_identified_in_gwts.csv: The MOTUs identified to species level in greater white-toothed shrews in this study
species_identified_in_bats.csv: The MOTUs identified to species level in bats in this study

If there are any questions regarding this data and/or scripts, please do not hesistate to contact Dr Sam Browett at browett.sam@gmail.com

