## GENCODE GTF

Obtained from <https://www.gencodegenes.org/human/release_33lift37.html>.

## PolyASite Atlas

See:

https://polyasite.unibas.ch/atlas

https://academic.oup.com/nar/article/48/D1/D174/5588346

```
wget https://polyasite.unibas.ch/download/atlas/2.0/GRCh38.96/atlas.clusters.2.0.GRCh38.96.bed.gz
wget ftp://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/GRCh38_to_GRCh37.chain.gz
CrossMap.py bed GRCh38_to_GRCh37.chain.gz atlas.clusters.2.0.GRCh38.96.bed.gz atlas.clusters.2.0.GRCh37.96.bed
sed -i '/_PATCH/!s/^/chr/' atlas.clusters.2.0.GRCh37.96.bed
```

## IsoAnnotLite/tappAS reference GFF

```
wget http://app.tappas.org/resources/downloads/gffs/Homo_sapiens_GRCh38_Ensembl_86.zip
unzip Homo_sapiens_GRCh38_Ensembl_86.zip
rm Homo_sapiens_GRCh38_Ensembl_86.zip
#CrossMap.py gff GRCh38_to_GRCh37.chain.gz Homo_sapiens_GRCh38_Ensembl_86.gff3 Homo_sapiens_GRCh37_Ensembl_86.gff3

wget ftp://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/GRCh37_to_GRCh38.chain.gz
CrossMap.py gff ref/GRCh37_to_GRCh38.chain.gz cp_vz_0.75_min_7_recovery_talon.gff2 cp_vz_0.75_min_7_recovery_talon.GRCh38.gff2

conda activate .snakemake/conda/59fab4b877169aaff69a189e3f1daa5b
python3 SQANTI3/utilities/IsoAnnotLite_SQ3.py \
    -o cp_vz_0.75_min_7_recovery_talon.isoannot.gff3 \
    results/cp_vz_0.75_min_7_recovery_talon.SQANTI3/cp_vz_0.75_min_7_recovery_talon_corrected.gtf \
    results/cp_vz_0.75_min_7_recovery_talon.SQANTI3/cp_vz_0.75_min_7_recovery_talon_classification.txt \
    results/cp_vz_0.75_min_7_recovery_talon.SQANTI3/cp_vz_0.75_min_7_recovery_talon_junctions.txt \
    -gff3 ref/Homo_sapiens_GRCh38_Ensembl_86.gff3 \
    -stdout cp_vz_0.75_min_7_recovery_talon.isoannot.log

```

```
conda activate .snakemake/conda/59fab4b877169aaff69a189e3f1daa5b
talon_abundance \
    --db talon.db.bak \
    --build GRCh37 \
    --annot gencode.v33lift37 \
    --whitelist talon_whitelist_0.75_min_7_recovery.txt \
    --o cp_vz_0.75_min_7_recovery

(head -n1 cp_vz_0.75_min_7_recovery_talon_abundance_filtered.tsv | cut -f 12-; tail -n+2 cp_vz_0.75_min_7_recovery_talon_abundance_filtered.tsv | cut -f 4,12-) > cp_vz_0.75_min_7_recovery_talon.fl_counts.tsv
```
