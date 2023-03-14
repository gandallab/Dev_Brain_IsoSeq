#!/usr/bin/env Rscript
#$ -cwd
#$ -N ISO-2a
#$ -V
#$ -o joblogs/2a_ISO_wgcna.$JOB_ID
#$ -j y
#$ -l h_rt=4:00:00,h_data=32G
#$ -pe shared 16
library(WGCNA)

if(FALSE) {
  suppressPackageStartupMessages({
    library(tidyverse)
    library(edgeR)
    library(WGCNA)
  })
  
  myDesign = tribble(
    ~sampleID, ~condition,~donor,
    "VZ_209", "VZ","209",
    "VZ_334", "VZ","334",
    "VZ_336", "VZ","336",
    "CP_209", "CP","209",
    "CP_334", "CP","334",
    "CP_336", "CP","336",
  ) %>%
    dplyr::mutate(
      dplyr::across(condition, as_factor)
    )
  
  cts = read_tsv("data/cp_vz_0.75_min_7_recovery_talon_abundance_filtered.tsv.gz")
  
  datExpr.counts = as.data.frame(cts[,12:35])
  rownames(datExpr.counts) = cts$annot_transcript_id
  datMeta = data.frame(sample=colnames(datExpr.counts))
  datMeta$Region = substr(datMeta$sample, 7,9)
  datMeta$Subject = substr(datMeta$sample, 1,3)
  datMeta$batch = substr(datMeta$sample, 5,5)
  
  
  talonSwitchList = readRDS("data/working/talonSwitchList_preFilter.rds")
  talon_filtered_isoforms = talonSwitchList$isoformFeatures$isoform_id
  localAnnoation <- unique(as.data.frame(talonSwitchList$exons@elementMetadata[,c("gene_id", "isoform_id")]))
  rm(talonSwitchList)
  
  celltypemarkers <- openxlsx::read.xlsx('https://www.cell.com/cms/10.1016/j.neuron.2019.06.011/attachment/508ae008-b926-487a-b871-844a12acc1f8/mmc5', sheet='Cluster enriched genes') %>% as_tibble()
  celltypemarkers_expressed = read_tsv("/Users/mgandal/GandalLab Dropbox/Michael Gandal/datasets/fetal_hub/data/polioudakis_neuron2020/TableS4-clusterExpressedTidy.tsv")
  
  celltypemarkers_tableS5 = openxlsx::read.xlsx('/Users/mgandal/GandalLab Dropbox/Michael Gandal/datasets/fetal_hub/data/polioudakis_neuron2020/TableS5_mmc6.xlsx',sheet = 2)
  
  
  # Use filtered isoforms from isoformSwitchAnalyer
  count_offset = 0.5
  
  datExpr.cpm = cpm(calcNormFactors(DGEList(counts=datExpr.counts[talon_filtered_isoforms,], group = datMeta$Region)))
  datExpr.cpm = log2(datExpr.cpm + count_offset)
  plot(density(datExpr.cpm[,1]))
  for(i in 2:24) lines(density(datExpr.cpm[,i]))
  
  datExpr.cpm <- as.data.frame(limma::removeBatchEffect(x = datExpr.cpm, batch = datMeta$batch,design=model.matrix(~datMeta$Region + datMeta$Subject)))
  
  isoformRepExpression <- 2^datExpr.cpm - count_offset
  isoformRepExpression[which(isoformRepExpression < 0, arr.ind = TRUE)] <- 0
  
  datExpr.localIF <- isoformToIsoformFraction(isoformRepExpression = isoformRepExpression, isoformGeneAnnotation = localAnnoation, quiet = F)
  datExpr.localIF <- round(datExpr.localIF, digits = 3)
  datExpr.localIF[is.na((datExpr.localIF))] = 0
  
  rm(isoformRepExpression)
  rm(datExpr.cpm)
  
  rm(datExpr.counts)
  
  softPower = 14
  maxNetSize=103000
  networkType='signed'
  fileBase = paste0("data/working/WGCNA/WGCNA_isoformFraction_top", maxNetSize/1000, "k_", networkType, "_sft", softPower)
  
  
  
  datIsoforms = data.frame(transcript_id=rownames(datExpr.localIF))
  datIsoforms <- left_join(datIsoforms, cts %>% 
                             dplyr::select(transcript_id=annot_transcript_id, gene_id =annot_gene_id, gene_name=annot_gene_name))
  datIsoforms$gene_id = substr(datIsoforms$gene_id,0,15)
  net=readRDS(paste0(fileBase, "_net.rds"))
  save.image(file='data/working/WGCNA/final/datExpr.localIF_batchCorrected_103k_image.RData')
}

setwd("/u/project/gandalm/gandalm/Github/Fetal_bulk_IsoSeq")
load('/u/project/gandalm/gandalm/Github/Fetal_bulk_IsoSeq/data/working/WGCNA/final/datExpr.localIF_batchCorrected_103k_image.RData')

fileBase= '/u/project/gandalm/gandalm/Github/Fetal_bulk_IsoSeq/data/working/WGCNA/WGCNA_isoformFraction_top103k_signed_sft14'

load(paste0(fileBase, "-block.1.RData"))
net$tree = hclust(1-as.dist(TOM), method="average")
net$colors.initial = net$colors
net$genes = rownames(datExpr.localIF)

# Recut the TOM 
ds = 4; minModSize = 100; dthresh = 0.05; pam = TRUE;
net$cut1 = cutreeHybrid(dendro = net$tree, pamStage=pam, minClusterSize= minModSize, cutHeight =  0.999999, 
                        deepSplit=ds, distM=as.matrix(1-TOM))
net$cut1$merged = mergeCloseModules(t(datExpr.localIF), colors = net$cut1$labels,cutHeight = dthresh)
names(net$cut1$labels) = rownames(datExpr.localIF)

# Recut the TOM 
ds = 4; minModSize = 100; dthresh = 0.05; pam = FALSE;
net$cut2 = cutreeHybrid(dendro = net$tree, pamStage=pam, minClusterSize= minModSize, cutHeight =  0.999999, 
                        deepSplit=ds, distM=as.matrix(1-TOM))
net$cut2$merged = mergeCloseModules(t(datExpr.localIF), colors = net$cut2$labels,cutHeight = dthresh)
names(net$cut2$labels) = rownames(datExpr.localIF)

saveRDS(net,file=paste0(fileBase, "_net.rds"))