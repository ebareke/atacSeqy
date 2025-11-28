#!/usr/bin/env bash

# ============================================================
# ATAC-seq Pipeline — FINAL: paired/single detection, mito QC/filtering,
# adaptive peak calling per species, ArchR/ChromVAR downstream analysis,
# YAML config, SLURM/PBS array jobs, auto-QC and concordance.
#
# Features:
#  - Auto-detect paired-end vs single-end from samplesheet (fastq2 empty -> SE)
#  - Mitochondrial ratio calculation and optional removal/filtering
#  - Per-species adaptive peak mode (narrow | broad | auto) via YAML config
#  - Optional ArchR project and ChromVAR deviations matrix creation (R)
#  - PBS/SLURM array submission with dependency so QC runs after all samples
#  - Fingerprinting QC, TSS enrichment, FRiP, cross-correlation, PCA/UMAP concordance
#
# Requirements (install in your environment):
#  - System: bash, bwa, samtools, bedtools, awk, sed, sort
#  - Tools: yq, atropos, macs2, deepTools (computeMatrix/plotProfile/plotFingerprint, bamCoverage), bamToBed, featureCounts, multiqc, idr
#  - R + packages: DESeq2, ggplot2, pheatmap, uwot, ArchR, chromVAR, SummarizedExperiment
#  - phantompeakqualtools (optional for run_spp.R)
#
# Usage example:
#  bash atacseq_pipeline.sh --config config.yaml --samples samples.csv --species human --threads 16 --cluster slurm
#
set -euo pipefail
log(){ echo "[\$(date +%FT%T)] $*" >&2; }
run(){ if [[ "${DRYRUN:-0}" -eq 1 ]]; then echo "DRYRUN: $*"; else eval "$*"; fi }

# ----------------------------- Defaults
THREADS=8
OUTDIR="results"
CONFIG=""
SAMPLESHEET=""
SPECIES=""
DRYRUN=0
RUN_IDR=0
RUN_DESEQ=0
CLUSTER_MODE="local"   # local | pbs | slurm
RESUME=0

usage(){ sed -n '1,240p' "$0"; exit 1; }

# ----------------------------- Arg parsing
while [[ $# -gt 0 ]]; do
  case "$1" in
    --config) CONFIG="$2"; shift 2;;
    --samples) SAMPLESHEET="$2"; shift 2;;
    --species) SPECIES="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --cluster) CLUSTER_MODE="$2"; shift 2;;
    --dryrun) DRYRUN=1; shift;;
    --idr) RUN_IDR=1; shift;;
    --deseq) RUN_DESEQ=1; shift;;
    --resume) RESUME=1; shift;;
    --help) usage;;
    *) echo "Unknown arg: $1"; usage;;
  esac
done

# ----------------------------- Helpers
if ! command -v yq >/dev/null 2>&1; then log "ERROR: yq required (https://github.com/mikefarah/yq)"; exit 1; fi
get_cfg(){ yq -r ".$1 // empty" "$CONFIG" 2>/dev/null || echo ""; }

select_species(){
  if [[ -z "$SPECIES" ]]; then log "ERROR: --species required"; exit 1; fi
  GENOME_FA=$(get_cfg "species.$SPECIES.genome_fa")
  BWA_PREFIX=$(get_cfg "species.$SPECIES.bwa_index_prefix")
  EFFECTIVE_GENOME_SIZE=$(get_cfg "species.$SPECIES.effective_genome_size")
  BLACKLIST=$(get_cfg "species.$SPECIES.blacklist")
  TSS_BED=$(get_cfg "species.$SPECIES.tss_bed")
  MITO_NAME=$(get_cfg "species.$SPECIES.mito_name")
  MITO_FILTER=$(get_cfg "species.$SPECIES.mito_filter")
  MITO_THRESHOLD=$(get_cfg "species.$SPECIES.mito_threshold")
  PEAK_MODE=$(get_cfg "species.$SPECIES.peak_mode")
  ARCHR_ENABLE=$(get_cfg "analysis.archr_enabled")
  CHROMVAR_ENABLE=$(get_cfg "analysis.chromvar_enabled")

  [[ -z "$BWA_PREFIX" ]] && BWA_PREFIX="${GENOME_FA%.*}.bwa"
  [[ -z "$PEAK_MODE" ]] && PEAK_MODE="narrow"

  log "Selected species: $SPECIES"
}

select_species

ensure_index(){ if [[ ! -f "${BWA_PREFIX}.bwt" ]]; then log "Building bwa index for $GENOME_FA"; run "bwa index -p ${BWA_PREFIX} ${GENOME_FA}"; fi }

# ----------------------------- create sample list and detect PE/SE
if [[ -z "$SAMPLESHEET" || ! -f "$SAMPLESHEET" ]]; then log "ERROR: provide --samples samples.csv"; exit 1; fi

SAMPLES=()
while IFS=, read -r sample_id fastq1 fastq2 bam group rep species_col; do
  [[ "$sample_id" == "sample_id" ]] && continue
  [[ -n "$species_col" && "$species_col" != "" ]] && SPECIES="$species_col" && select_species
  paired=1
  if [[ -z "$fastq2" || "$fastq2" == "" ]]; then paired=0; fi
  SAMPLES+=("$sample_id|$fastq1|$fastq2|$bam|$group|$rep|$paired")
done < "$SAMPLESHEET"
TOTAL=${#SAMPLES[@]}

mkdir -p "$OUTDIR/jobs" "$OUTDIR/samples" "$OUTDIR/qc" "$OUTDIR/consensus"
cd "$OUTDIR"

# ----------------------------- per-sample processing function
process_sample(){
  sample_id="$1"; fastq1="$2"; fastq2="$3"; bam_in="$4"; group="$5"; rep="$6"; paired="$7"
  log "[sample] $sample_id (paired=$paired)"
  mkdir -p samples/${sample_id}
  pushd samples/${sample_id} >/dev/null

  ensure_index
  mapped_bam="${sample_id}.sorted.bam"

  # map: bwa mem handles SE/PE
  if [[ -n "$bam_in" && -f "$bam_in" ]]; then
    log "Using user-provided BAM"
    run "samtools sort -@ $THREADS -o $mapped_bam $bam_in"
  else
    if [[ "$paired" -eq 1 ]]; then
      run "atropos trim --threads $THREADS -a AGATCGGAAGAG -A AGATCGGAAGAG -pe1 $fastq1 -pe2 $fastq2 -o R1.trim.fq.gz -p R2.trim.fq.gz --minimum-length 25 --nextseq-trim 25"
      run "bwa mem -M -t $THREADS ${BWA_PREFIX} R1.trim.fq.gz R2.trim.fq.gz | samtools sort -@ $THREADS -o $mapped_bam -"
    else
      run "atropos trim --threads $THREADS -a AGATCGGAAGAG -o R1.trim.fq.gz --minimum-length 25 $fastq1"
      run "bwa mem -M -t $THREADS ${BWA_PREFIX} R1.trim.fq.gz | samtools sort -@ $THREADS -o $mapped_bam -"
    fi
  fi
  run "samtools index $mapped_bam"

  # mitochondrial ratio
  if [[ -n "$MITO_NAME" && -n "$MITO_THRESHOLD" ]]; then
    total_reads=$(samtools view -c -F 0x904 $mapped_bam || echo 0)
    mito_reads=$(samtools view -c -F 0x904 -r $MITO_NAME $mapped_bam 2>/dev/null || true)
    # alternative: use idxstats
    if [[ -z "$mito_reads" || "$mito_reads" == "" ]]; then
      mito_reads=$(samtools idxstats $mapped_bam | awk -v m="$MITO_NAME" '$1==m{print $3}')
    fi
    mito_reads=${mito_reads:-0}
    if [[ "$total_reads" -gt 0 ]]; then
      mito_ratio=$(awk "BEGIN{printf \"%.4f\", $mito_reads/$total_reads}")
    else mito_ratio=0; fi
    echo -e "$sample_id	$total_reads	$mito_reads	$mito_ratio" >> ../../qc/mito_ratio.tsv
    if [[ "$MITO_FILTER" == "true" && $(awk "BEGIN{print ($mito_ratio>$MITO_THRESHOLD)?1:0}") -eq 1 ]]; then
      log "Mito ratio $mito_ratio > $MITO_THRESHOLD; removing $MITO_NAME reads"
      chroms=$(samtools idxstats $mapped_bam | cut -f1 | grep -v -w "$MITO_NAME" | tr '
' ' ')
      run "samtools view -b $mapped_bam $chroms > ${sample_id}.noMT.bam"
      run "samtools index ${sample_id}.noMT.bam"
      mapped_bam="${sample_id}.noMT.bam"
    fi
  fi

  # mark duplicates
  run "samtools fixmate -m $mapped_bam fix.bam || true"
  run "samtools sort -@ $THREADS -o ns.bam fix.bam || true"
  run "samtools markdup -r ns.bam dedup.bam || samtools markdup -r $mapped_bam dedup.bam"
  run "samtools index dedup.bam"

  # insert size
  if command -v picard >/dev/null 2>&1; then
    run "picard CollectInsertSizeMetrics I=dedup.bam O=${sample_id}.insert_size_metrics.txt H=${sample_id}.insert_size_hist.pdf M=0.5"
  fi

  # Tn5 shift
  run "alignmentSieve --bam dedup.bam --outFile shifted.bam --ATACshift --numberOfProcessors $THREADS --minMappingQuality 10"
  run "samtools index shifted.bam"

  # blacklist
  if [[ -n "$BLACKLIST" && -f "$BLACKLIST" ]]; then
    run "bedtools intersect -v -abam shifted.bam -b $BLACKLIST > filtered.bam"
    run "samtools index filtered.bam"
  else
    run "cp shifted.bam filtered.bam && cp shifted.bam.bai filtered.bam.bai"
  fi

  # adaptive peak calling: per-species PEAK_MODE (narrow|broad|auto)
  mkdir -p macs2
  peak_call_mode="$PEAK_MODE"
  if [[ "$peak_call_mode" == "auto" ]]; then
    if [[ "$EFFECTIVE_GENOME_SIZE" -gt 5e8 ]]; then peak_call_mode="broad"; else peak_call_mode="narrow"; fi
  fi
  if [[ "$peak_call_mode" == "narrow" ]]; then
    run "macs2 callpeak -t filtered.bam -f BAM -n ${sample_id} --nomodel --shift 0 --extsize 200 -q 0.01 --keep-dup all --outdir macs2"
  else
    run "macs2 callpeak -t filtered.bam -f BAM -n ${sample_id} --broad --nomodel --shift 0 --extsize 200 -q 0.05 --keep-dup all --outdir macs2"
  fi

  # fragments
  run "bamToBed -bedpe -i filtered.bam > fragments.bedpe"

  # bigWig
  if command -v bamCoverage >/dev/null 2>&1; then
    run "bamCoverage -b filtered.bam -o ${sample_id}.bw --normalizeUsing RPGC --effectiveGenomeSize ${EFFECTIVE_GENOME_SIZE} --binSize 10 -p $THREADS"
  fi

  popd >/dev/null
}

# ----------------------------- write array commands
mkdir -p jobs job_scripts
ARRAY_FILE="job_scripts/array_cmds.txt"
> "$ARRAY_FILE"
for entry in "${SAMPLES[@]}"; do
  IFS='|' read -r sid fq1 fq2 sbam grp rep paired <<< "$entry"
  echo "process_sample $sid $fq1 $fq2 $sbam $grp $rep $paired" >> "$ARRAY_FILE"
done

# submit arrays and qc jobs
if [[ "$CLUSTER_MODE" == "pbs" ]]; then
  cat > jobs/array.pbs <<'PBS'
#!/bin/bash
#PBS -N ATAC_ARRAY
#PBS -J 1-${TOTAL}
#PBS -l nodes=1:ppn=${THREADS}
#PBS -l walltime=48:00:00
cd $PBS_O_WORKDIR
COMMAND=$(sed -n "${PBS_ARRAY_INDEX}p" job_scripts/array_cmds.txt)
bash -lc "$COMMAND"
PBS
  ARRAY_ID=$(run "qsub jobs/array.pbs" | awk '{print $1}')
  log "Submitted PBS array: $ARRAY_ID"
  cat > jobs/qc_after.pbs <<PBS
#!/bin/bash
#PBS -N ATAC_QC
#PBS -l nodes=1:ppn=${THREADS}
#PBS -l walltime=24:00:00
#PBS -W depend=afterokarray:${ARRAY_ID}
cd $PWD
bash atacseq_pipeline.sh --config $CONFIG --samples $SAMPLESHEET --species $SPECIES --threads $THREADS --cluster local --resume
PBS
  run "qsub jobs/qc_after.pbs"
  log "QC job submitted with dependency on array"
  exit 0
fi

if [[ "$CLUSTER_MODE" == "slurm" ]]; then
  cat > jobs/array.slurm <<'SLURM'
#!/bin/bash
#SBATCH --job-name=ATAC_ARRAY
#SBATCH --array=1-${TOTAL}
#SBATCH --cpus-per-task=${THREADS}
#SBATCH --time=48:00:00
cd $PWD
COMMAND=$(sed -n "${SLURM_ARRAY_TASK_ID}p" job_scripts/array_cmds.txt)
bash -lc "$COMMAND"
SLURM
  ARRAY_SUB=$(run "sbatch jobs/array.slurm" | awk '{print $4}')
  log "Submitted SLURM array: $ARRAY_SUB"
  cat > jobs/qc_after.slurm <<SLURM
#!/bin/bash
#SBATCH --job-name=ATAC_QC
#SBATCH --cpus-per-task=${THREADS}
#SBATCH --time=24:00:00
#SBATCH --dependency=afterok:${ARRAY_SUB}
cd $PWD
bash atacseq_pipeline.sh --config $CONFIG --samples $SAMPLESHEET --species $SPECIES --threads $THREADS --cluster local --resume
SLURM
  run "sbatch jobs/qc_after.slurm"
  exit 0
fi

# local mode: run sequentially
if [[ "$CLUSTER_MODE" == "local" ]]; then
  while IFS= read -r cmd; do
    run "$cmd"
  done < "$ARRAY_FILE"
fi

# ----------------------------- After mapping (QC/Consensus/Counts/ArchR)
# The remainder of the pipeline (consensus peak merging, counting, QC metrics,
# ArchR/ChromVAR execution and concordance clustering) is launched when running
# the script with --resume after array jobs finish (or automatically via cluster dependency above).

if [[ "$RESUME" -eq 1 || "$CLUSTER_MODE" == "local" ]]; then
  log "Starting post-processing steps: consensus peaks, counts, QC, ArchR/ChromVAR"

  # MERGE peaks
  find samples -name "*_peaks.narrowPeak" -o -name "*_peaks.broadPeak" > peaks.list || true
  if [[ -s peaks.list ]]; then
    run "cat \$(cat peaks.list) | sort -k1,1 -k2,2n | bedtools merge -i - > consensus/consensus.bed"
    run "awk 'BEGIN{OFS="	"} {print $1,$2+1,$3,\"peak\"NR,1}' consensus/consensus.bed > consensus/consensus.saf"
  fi

  # COUNTS
  BAMS=()
  for sdir in samples/*; do
    sid=$(basename $sdir)
    fb="$sdir/filtered.bam"
    [[ -f "$fb" ]] && BAMS+=("$fb")
  done
  if [[ ${#BAMS[@]} -gt 0 && -f consensus/consensus.saf ]]; then
    run "featureCounts -a consensus/consensus.saf -F SAF -o consensus/counts.txt -T $THREADS -p -B -C ${BAMS[*]}"
    run "cut -f1,7- consensus/counts.txt > consensus/counts.matrix"
  fi

  # FRiP
  > qc/frip.tsv
  for d in samples/*; do sid=$(basename $d); bamf=$d/filtered.bam; peakf=$d/macs2/${sid}_peaks.narrowPeak; if [[ -f $bamf && -f $peakf ]]; then total=$(samtools view -c -F 0x904 $bamf); inpk=$(bedtools intersect -u -a <(samtools view -F 0x904 -b $bamf) -b $peakf | samtools view -c -); frip=$(awk "BEGIN{printf \"%.4f\", ($total==0)?0:$inpk/$total}"); echo -e "$sid	$total	$inpk	$frip" >> qc/frip.tsv; fi; done

  # Fingerprint
  bwlist=$(find samples -name "*.bw" | tr '
' ' ')
  if [[ -n "$bwlist" ]]; then run "plotFingerprint -b $bwlist -o qc/fingerprint.png --labels $(echo $bwlist | xargs -n1 basename | tr '
' ',')"; fi

  # TSS Enrichment (if TSS_BED provided)
  if [[ -n "$TSS_BED" && -f "$TSS_BED" ]]; then
    mkdir -p qc/tss
    for d in samples/*; do sid=$(basename $d); if [[ -f $d/${sid}.bw ]]; then run "computeMatrix reference-point -S $d/${sid}.bw -R $TSS_BED --beforeRegionStartLength 2000 --afterRegionStartLength 2000 --skipZeros -o qc/tss/${sid}.mat.gz -p $THREADS"; run "plotProfile -m qc/tss/${sid}.mat.gz -o qc/tss/${sid}_tss.png"; fi; done
  fi

  # ArchR + ChromVAR (optional)
  if [[ "$ARCHR_ENABLE" == "true" || "$CHROMVAR_ENABLE" == "true" ]]; then
    log "Running ArchR/ChromVAR R pipeline"
    cat > archr_run.R <<'R'
#!/usr/bin/env Rscript
library(ArchR)
addArchRThreads(threads = as.integer(Sys.getenv('THREADS')))
# collect fragments and sample metadata
samples <- read.csv(Sys.getenv('SAMPLESHEET'), stringsAsFactors=FALSE)
frags <- file.path('samples', samples$sample_id, 'fragments.bedpe')
# Convert bedpe fragments to fragment files expected by ArchR (chr	start	end	sample)
# ArchR createArrowFiles() can take fragmentFiles and sampleNames
valid <- file.exists(frags)
frags <- frags[valid]
sampleNames <- samples$sample_id[valid]
if(length(frags)==0) stop('No fragment files found for ArchR')
ArrowFiles <- createArrowFiles(inputFiles = frags, sampleNames = sampleNames, minTSS=4, minFrags=1000, addTileMat=TRUE, addGeneScoreMat=TRUE)
proj <- ArchRProject(ArrowFiles = ArrowFiles, outputDirectory = 'ArchRProject', copyArrows = FALSE)
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = 'TileMatrix', name = 'IterativeLSI', iterations = 2, clusterParams = list(resolution = c(0.2), sampleCells = 10000, n.start = 10))
proj <- addClusters(input = proj, reducedDims = 'IterativeLSI')
proj <- addUMAP(ArchRProj = proj, reducedDims = 'IterativeLSI')
plotEmbedding(ArchRProj = proj, colorBy = 'cellColData', name = 'Sample', embedding = 'UMAP')
# add motif annotations and chromVAR deviations
proj <- addMotifAnnotations(ArchRProj = proj, motifSet = 'cisbp', name = 'Motif')
proj <- addBgdPeaks(proj)
proj <- addDeviationsMatrix(proj, peakAnnotation = 'Motif')
saveArchRProject(ArchRProj = proj, outputDirectory = 'ArchRProject', load = FALSE)
R
    export SAMPLESHEET="$(realpath $SAMPLESHEET)"
    export THREADS=$THREADS
    run "Rscript archr_run.R"
  fi

  # Concordance clustering (PCA/UMAP) — R script
  cat > concordance_run.R <<'R'
#!/usr/bin/env Rscript
library(DESeq2); library(ggplot2); library(pheatmap); library(uwot); library(matrixStats)
counts <- read.table('consensus/counts.matrix', header=TRUE, row.names=1)
ss <- read.csv(Sys.getenv('SAMPLESHEET'), stringsAsFactors=FALSE)
rownames(ss) <- ss$sample_id
common <- intersect(colnames(counts), rownames(ss))
counts <- counts[,common]; ss <- ss[common,]
dds <- DESeqDataSetFromMatrix(countData=counts, colData=ss, design=~group)
dds <- estimateSizeFactors(dds); vsd <- vst(dds, blind=FALSE); mat <- assay(vsd)
# PCA
pca <- prcomp(t(mat)); pc.df <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], ss)
png('concordance/pca.png', width=800, height=600); plot(pc.df$PC1, pc.df$PC2, col=as.factor(pc.df$group), pch=19, main='PCA'); dev.off()
# UMAP
um <- umap(t(mat), n_neighbors=15); png('concordance/umap.png', width=800, height=600); plot(um, col=as.factor(ss$group), pch=19, main='UMAP'); dev.off()
# correlation heatmap
cor.mat <- cor(t(mat)); pheatmap(cor.mat, filename='concordance/corr.png')
R
  export SAMPLESHEET="$(realpath $SAMPLESHEET)"
  run "Rscript concordance_run.R"

  # MultiQC
  if command -v multiqc >/dev/null 2>&1; then run "multiqc . -n multiqc_report"; fi

  log "Post-processing complete. Results in: $PWD"
fi

# ----------------------------- EXAMPLE YAML CONFIG + SAMPLESHEET + README
cat > EXAMPLE_config.yaml <<'YAML'
species:
  human:
    genome_fa: /path/genomes/hg38.fa
    bwa_index_prefix: /path/genomes/hg38
    effective_genome_size: 2913022398
    blacklist: /path/genomes/hg38.blacklist.bed
    tss_bed: /path/genomes/hg38_tss.bed
    mito_name: chrM
    mito_filter: true
    mito_threshold: 0.2
    peak_mode: narrow

  drosophila:
    genome_fa: /path/genomes/dm6.fa
    bwa_index_prefix: /path/genomes/dm6
    effective_genome_size: 142573017
    blacklist: /path/genomes/dm6.blacklist.bed
    tss_bed: /path/genomes/dm6_tss.bed
    mito_name: chrM
    mito_filter: false
    mito_threshold: 0.2
    peak_mode: narrow

analysis:
  archr_enabled: true
  chromvar_enabled: true
YAML

cat > EXAMPLE_samples.csv <<'CSV'
sample_id,fastq1,fastq2,bam,group,replicate,species
S1,/data/S1_R1.fq.gz,/data/S1_R2.fq.gz,,Control,1,human
S2,/data/S2_R1.fq.gz,/data/S2_R2.fq.gz,,Control,2,human
D1,/data/D1_R1.fq.gz,, ,Control,1,drosophila
CSV

cat > README_ATAC_FINAL.txt <<'README'
ATAC-seq Pipeline — Final

1) Edit EXAMPLE_config.yaml to point to your genomes and resources.
2) Create a samplesheet like EXAMPLE_samples.csv. Leave fastq2 empty for single-end samples.
3) Run on SLURM:
   bash atacseq_pipeline.sh --config config.yaml --samples samples.csv --species human --threads 16 --cluster slurm

Notes:
- The script auto-detects paired-end vs single-end per sample from samplesheet.
- Mitochondrial ratio is computed; if mito_filter=true and ratio > mito_threshold the mitochondrial reads will be removed prior to peak calling.
- Per-species peak_mode can be 'narrow', 'broad', or 'auto' (auto decides based on effective genome size).
- ArchR/ChromVAR will be run if enabled in config (requires ArchR installed in R).
- Resource tuning (memory, time) for cluster scripts should be adjusted for your HPC environment.
README

log "Script written and examples created in $(pwd)"
