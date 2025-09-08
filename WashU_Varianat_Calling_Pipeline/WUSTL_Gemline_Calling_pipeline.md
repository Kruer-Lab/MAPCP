# Germline Calling, Joint Genotyping, and Quality Control Workflow

This document provides a step-by-step overview of the pipeline for **whole-exome sequencing (WES)** analysis, including germline variant calling, joint genotyping through GATK GenomicsDBImport, and downstream quality control analyses. The commands are shown with their corresponding Docker images to ensure reproducibility across computational environments.

> Sources: `Pbrun_Germline_Functional_Equivalent_CMD.sh`, `Gatk_jointDBImport_CMD.sh`, `Quality_control_CMD.sh`.

---

## 1. Germline Variant Calling (WES)

The first stage processes paired-end FASTQ files to generate gVCFs per sample. This step includes alignment, base quality score recalibration, and variant calling.

### Step 1A. Alignment (Parabricks fq2bam)

**Docker image:** `nvcr.io/nvidia/clara/clara-parabricks:4.0.0-1`

```bash
/usr/local/parabricks/pbrun fq2bam --x4 \
  --ref Homo_sapiens_assembly38.fasta \
  --in-fq Sample1_R1.fastq.gz Sample1_R2.fastq.gz \
  --knownSites snvs.hq.vcf.gz \
  --knownSites indels1.hq.vcf.gz \
  --knownSites indels2.hq.vcf.gz \
  --interval-file Exome-IDT_V1V2_span50bp.bed \
  --out-bam Sample1_wes_germline.cram \
  --out-recal-file Sample1_wes_germline_report.txt \
  --out-duplicate-metrics Sample1_wes_dup_metrics.txt \
  --bwa-options='-Y -K 100000000' \
  --read-group-pl ILLUMINA \
  --read-group-sm Sample1 \
  --tmp-dir germ_Sample1_Temp \
  --num-gpus 1 \
  --memory-limit 280GB
```

### Step 1B. Base Quality Score Recalibration (BQSR)

**Docker image:** `broadinstitute/gatk:4.2.0.0`

```bash
/gatk/gatk ApplyBQSR \
  -R Homo_sapiens_assembly38.fasta \
  -I Sample1_wes_germline.cram \
  --bqsr-recal-file Sample1_wes_germline_report.txt \
  -O Sample1_wes_germline_applyBQSR.cram
```

### Step 1C. Variant Calling (HaplotypeCaller)

**Docker image:** `broadinstitute/gatk:4.2.0.0`

```bash
/gatk/gatk HaplotypeCaller \
  -I Sample1_wes_germline_applyBQSR.cram \
  -O Sample1_wes_germline.g.vcf.bgz \
  -R Homo_sapiens_assembly38.fasta \
  -L Exome-IDT_V1V2_span50bp.bed \
  -ERC BP_RESOLUTION
```

> For gVCF mode, replace `-ERC BP_RESOLUTION` with `-ERC GVCF`.

---

## 2. Joint Genotyping with GATK

The second stage aggregates gVCFs from multiple samples using **GATK GenomicsDBImport** and performs joint genotyping, variant quality score recalibration (VQSR), and filtering.

### Step 2A. GenomicsDBImport (Create New Workspace)

**Docker image:** `broadinstitute/gatk:4.2.0.0`

```bash
/gatk/gatk GenomicsDBImport \
  --java-options "-Xmx288GB -Xms288GB" \
  -R Homo_sapiens_assembly38.fasta \
  -L Exome-IDT_V1V2_span50bp.bed \
  --merge-input-intervals \
  --genomicsdb-workspace-path DBI_workspace \
  --batch-size 50 \
  --reader-threads 24 \
  --genomicsdb-shared-posixfs-optimizations true \
  --tmp-dir DBI_Temp \
  --sample-name-map Test_DBImap.txt
```

### Step 2B. Joint Genotyping

**Docker image:** `broadinstitute/gatk:4.2.0.0`

```bash
/gatk/gatk GenotypeGVCFs \
  --java-options "-Xmx288GB -Xms288GB" \
  -R Homo_sapiens_assembly38.fasta \
  -V gendb://DBI_workspace \
  -O Project_BP_genotype.vcf.gz \
  --tmp-dir GEN_Temp
```

### Step 2C. Variant Recalibration (VQSR)

**Docker image:** `broadinstitute/gatk:4.2.0.0`

```bash
/gatk/gatk --java-options "-Xmx288GB -Xms288GB" VariantRecalibrator \
  -V Project_BP_genotype.vcf.gz \
  -O Project_joint_vqsr.recal \
  --tranches-file Project_joint_vqsr.tranches \
  --resource:omni,known=false,training=true,truth=true,prior=12.0 1000G_omni2.5.hg38.vcf.bgz \
  --resource:1000g,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.high_confidence.hg38.vcf.bgz \
  --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 Homo_sapiens_assembly38.dbsnp138.vcf.gz \
  --resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.hg38.vcf.gz \
  --resource:mills,known=false,training=true,truth=true,prior=12.0 Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP -an MQ \
  -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 \
  -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 \
  -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
  --mode BOTH
```

### Step 2D. Apply VQSR

**Docker image:** `broadinstitute/gatk:4.2.0.0`

```bash
/gatk/gatk --java-options "-Xmx288GB -Xms288GB" ApplyVQSR \
  -R Homo_sapiens_assembly38.fasta \
  -V Project_BP_genotype.vcf.gz \
  --recal-file Project_joint_vqsr.recal \
  --tranches-file Project_joint_vqsr.tranches \
  -O Project_joint_vqsr.vcf.gz \
  --mode BOTH
```

### Step 2E. Extract PASS Variants

**Docker image:** `broadinstitute/gatk:4.2.0.0`

```bash
/gatk/gatk --java-options "-Xmx288GB -Xms288GB" SelectVariants \
  -R Homo_sapiens_assembly38.fasta \
  -V Project_joint_vqsr.vcf.gz \
  -O Project_joint_vqsr_gs.vcf.gz \
  -select 'vc.isNotFiltered()'
```

### Step 2F. Normalize and Index with BCFtools

**Docker image:** `dreammaerd/genomic-tools:v1`

```bash
/opt/conda/bin/bcftools norm -m-both -f Homo_sapiens_assembly38.fasta -O z Project_joint_vqsr_gs.vcf.gz \
  | /opt/conda/bin/bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' -O z -o Project_joint_vqsr_gs_sp_ln_reID.vcf.gz

/opt/conda/bin/bcftools index -t -f Project_joint_vqsr_gs_sp_ln_reID.vcf.gz
```

---

## 3. Quality Control (QC)

The final stage evaluates the quality of sequencing data and performs kinship, sex check, and population structure analysis.

### Step 3A. BamMetrics

**Docker image:** `dreammaerd/genomic-tools:v1`

```bash
/opt/bamMetrics -t 4 -r Homo_sapiens_assembly38.fasta -b Exome-IDT_V1V2_span50bp.bed -o Sample1_bamMetrics.tsv Sample1.cram
```

### Step 3B. Kinship Analysis

**Docker image:** `dreammaerd/genomic-tools:v1`

```bash
/bin/plink --geno 0.01 --genome --hwe 0.001 --maf 0.05 \
  --snps-only --allow-extra-chr --vcf-half-call m \
  --vcf Test.vcf.gz --out Project_kinship
```

### Step 3C. Sex Check

**Docker image:** `dreammaerd/genomic-tools:v1`

```bash
/bin/plink --double-id --allow-extra-chr --vcf-half-call m \
  --vcf Test.vcf.gz --out sexcheck_s1

/bin/plink -split-x hg38 --make-bed --allow-extra-chr \
  --bfile sexcheck_s1 --out sexcheck_s2

/bin/plink --impute-sex ycount --make-bed --allow-extra-chr \
  --bfile sexcheck_s2 --out sexcheck_s3
```

### Step 3D. PCA (Population Structure)

**Docker image:** `dreammaerd/genomic-tools:v1`

```bash
vcf2geno --inVcf Test_pca_input.vcf.bgz --out Test_pca_input

trace -p my_parameterfile -k 10 -K 20 \
  -s Test_pca_input.geno \
  -g HGDP_938.geno \
  -c HGDP_938.RefPC.coord \
  -o Project_tracePCA
```

**Visualization (R PCA Plotting)** **Docker image:** `r-base`

```bash
Rscript hdgp_PCA.R -i Project_tracePCA.ProPC.coord
```

---

# 4. In-house Post-processing & Annotation (Hail)

### Step 4A. VCF → Variant Table

**Script:** `VCF_to_VariantTable_with_ParentalGT.py`
**Docker image:** `dreammaerd/hail_vep_gnomad:0.2.120_104_0.6.4`

```bash
python VCF_to_VariantTable_with_ParentalGT.py Test.vcf Test.ped
```

**Inputs**

* `Test.vcf` — joint called VCF file.
* `Test.ped` — pedigree file (tab-delimited).&#x20;

**PED file example (`Test.ped`)**

famID  ProbandID  MotherID  FatherID
fam1   Proband1   Mother1   Father1

### Step 4B. Hail Annotation

**Script:** `Hail_annotation.py` (uses `my_hail_function_general.py`, `my_hail_function_variantTable_annotation.py`)
**Docker image:** `dreammaerd/hail_vep_gnomad:0.2.120_104_0.6.4`

```bash
python Hail_annotation.py \
  -i VariantTable_Sample1.txt \
  -m all \
  -p Phenotype.tsv \
  -c CADD_web.tsv.bgz
```

**Manifest file example (********`Phenotype.tsv`****\*\*\*\*)**

```
ID    Relation_to_Proband    Sex    Phenotype
P001  Proband                M      Affected
P002  Mother                 F      Unaffected
P003  Father                 M      Unaffected
```

**CADD web file example (********`CADD_web.tsv.bgz`****\*\*\*\*)**

```
Chrom  Pos     Ref     Alt     RawScore    PHRED
1       123456  G       A       1.234       12.3
```

---

## Notes

- File names shown here (e.g., `Sample1_R1.fastq.gz`, `Test.vcf.gz`) are illustrative examples. Replace them with actual input files.
- The pipeline is modular: each stage (germline calling, joint genotyping, QC) can be executed independently.
- Docker images are listed explicitly to facilitate reproducibility across HPC and cloud environments.
- Memory requirements (e.g., `288GB` for GATK) are set according to script defaults but may be adjusted depending on hardware availability.

---

This workflow provides a complete, containerized solution for WES germline analysis, extending from raw FASTQ files through joint genotyped VCFs, QC assessments, and Hail-based variant annotation.

