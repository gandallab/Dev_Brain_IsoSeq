library(rtracklayer)
library(tidyverse)

geneAnno = rtracklayer::import("ref/gencode.v33lift37.annotation.gtf.gz") %>% as_tibble() %>% filter(type=='gene')
geneAnno$ensg = substr(geneAnno$gene_id,1,15)
geneAnno = geneAnno[match(unique(geneAnno$ensg), geneAnno$ensg),]

#Binary vector 0 for all protein coding genes, NA for non-coding; will replace with 1 for rare variant disease risk gene
protein_coding_bg = rep(NA, times=nrow(geneAnno))
protein_coding_bg[geneAnno$gene_type=="protein_coding"] = 0
names(protein_coding_bg) = geneAnno$ensg


#Start with compiled list of risk genes from SFARI, Kaplanis,
risk_genes =read.csv("ref/ASD+SCZ+DDD_2022.csv")
rareVar.binary =data.frame(ASD.SFARI.1 = protein_coding_bg)
rownames(rareVar.binary)= geneAnno$ensg
rareVar.binary$ASD.SFARI.1[geneAnno$gene_name %in% risk_genes$Gene[risk_genes$Set=="ASD (SFARI score 1)"]] = 1
rareVar.binary$ASD.SFARI.1or2 = protein_coding_bg
rareVar.binary$ASD.SFARI.1or2[geneAnno$gene_name %in% risk_genes$Gene[risk_genes$Set=="ASD (SFARI score 1or2)"]] = 1
rareVar.binary$ASD.SFARI.S = protein_coding_bg
rareVar.binary$ASD.SFARI.S[geneAnno$gene_name %in% risk_genes$Gene[risk_genes$Set=="ASD (SFARI syndromic)"]] = 1
rareVar.binary$DDD.kaplanis = protein_coding_bg
rareVar.binary$DDD.kaplanis[geneAnno$gene_name %in% risk_genes$Gene[risk_genes$Set=="DDD (Kaplanis et al. 2019)"]] = 1

epilepsy_helbig = openxlsx::read.xlsx('http://epilepsygenetics.net/wp-content/uploads/2023/01/Channelopathist_genes_internal_2023_v2.xlsx')
rareVar.binary$EPI.helbig = protein_coding_bg
rareVar.binary$EPI.helbig[geneAnno$gene_name %in% epilepsy_helbig$Gene] = 1


# Fu et al., Nat Genetics 2022 -- ASD, NDD TADA
fu=openxlsx::read.xlsx(('https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-022-01104-0/MediaObjects/41588_2022_1104_MOESM3_ESM.xlsx'),'Supplementary Table 11')
fu$p_TADA_ASD[fu$p_TADA_ASD==0] = min(fu$p_TADA_ASD[fu$p_TADA_ASD >0],na.rm=T)
fu$p_TADA_NDD[fu$p_TADA_NDD==0] = min(fu$p_TADA_NDD[fu$p_TADA_NDD >0],na.rm=T)

rareVar.logit = rareVar.binary
rareVar.logit$ASD.fuTADA= -log10(fu$p_TADA_ASD)[match(geneAnno$ensg, fu$gene_id)]
rareVar.binary$ASD.fuTADA= as.numeric(fu$FDR_TADA_ASD[match(geneAnno$ensg, fu$gene_id)] < .1)
rareVar.logit$NDD.fuTADA = -log10(fu$p_TADA_NDD)[match(geneAnno$ensg, fu$gene_id)]
rareVar.binary$NDD.fuTADA = as.numeric(fu$FDR_TADA_NDD[match(geneAnno$ensg, fu$gene_id)] < .1)

SCZ.schema = read_tsv('ref/risk_genes/SCHEMA_gene_results.tsv')
rareVar.logit$SCZ.schema = -log10(SCZ.schema$`P meta`[match(geneAnno$ensg, SCZ.schema$gene_id)])
rareVar.binary$SCZ.schema = as.numeric(SCZ.schema$`Q meta`[match(geneAnno$ensg, SCZ.schema$gene_id)] < .1)

BIP.bipex = read_tsv('ref/risk_genes/BipEx_gene_results.tsv') %>% filter(group=="Bipolar Disorder")
rareVar.logit$BIP.bipex = -log10(BIP.bipex$ptv_fisher_gnom_non_psych_pval[match(geneAnno$ensg,BIP.bipex$gene_id)])
rareVar.binary$BIP.bipex = as.numeric(BIP.bipex$ptv_fisher_gnom_non_psych_pval[match(geneAnno$ensg,BIP.bipex$gene_id)] < 0.01 )

EPI.epi25 = read_tsv('ref/risk_genes/Epi25_gene_results.tsv') %>% filter(group=="EPI")
rareVar.logit$EPI.epi25 = -log10(EPI.epi25$pval[match(geneAnno$ensg,EPI.epi25$gene_id)])
rareVar.binary$EPI.epi25 = as.numeric(EPI.epi25$pval[match(geneAnno$ensg,EPI.epi25$gene_id)] < .01)

apply(rareVar.binary,2, table)
rareVar.logit %>% as_tibble()
