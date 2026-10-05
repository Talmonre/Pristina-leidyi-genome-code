# Observations on the genome of selectively asexual organism *Pristina leidyi*

Code for the genome assembly, annotation and comparative analyses in:

> Achilleos K\*, Álvarez-Campos P\*, Bailie MA, Null R, Garcia-Castro H, Solana J, Özpolat BD, Kenny NJ. *Observations on the genome of selectively asexual organism* Pristina leidyi. *Genome Biology and Evolution*. \*Co-first authors.

The sequencing data, assembly and annotation are available from NCBI under the accessions given in the paper. This page contains the code only.

The analyses were run on the NeSI and University of Otago SLURM clusters. In the code below, cluster-specific settings (job headers, module and conda environment loading, absolute paths) have been removed. Input and output files are given generic names, set as variables at the top of each block. Unless stated otherwise, programs were run with default settings.

## Code authorship and contributions

This repository was developed by MAB and KA. MAB contributed to Sections 1.1-1.5 and 1.11-1.13. KA contributed to Sections 1.6-1.10, 1.14-1.16 and 2-3.

## Contents

- [Software](#software)
- [1. Genome assembly and assessment](#1-genome-assembly-and-assessment)
  - [1.1 Read quality control](#11-read-quality-control)
  - [1.2 k-mer profiling](#12-k-mer-profiling)
  - [1.3 Genome assembly with Hifiasm](#13-genome-assembly-with-hifiasm)
  - [1.4 Genome assemby polishing](#14-genome-assembly-polishing)
  - [1.5 Long-read scaffolding](#15-long-read-scaffolding)
  - [1.6 Assembly statistics](#16-assembly-statistics)
  - [1.7 DIAMOND BLASTx](#17-diamond-blastx)
  - [1.8 Read coverage](#18-read-coverage)
  - [1.9 Taxify](#19-taxify)
  - [1.10 BlobTools, BUSCO and snail plot (BlobToolkit)](#110-blobtools-busco-and-snail-plot-blobtoolkit)
  - [1.11 Hi-C scaffolding, Alignment and YaHS](#111-hi-c-scaffolding-alignment-and-yahs)
  - [1.12 Scaffold size filter](#112-scaffold-size-filter)
  - [1.13 Contact map](#113-contact-map)
  - [1.14 Retain longest scaffolds](#114-retain-longest-scaffolds)
  - [1.15 BLASTn](#115-blastn)
  - [1.16 Rename headers](#116-rename-headers)
- [2. Genome content and functional annotation](#2-genome-content-and-functional-annotation)
  - [2.1 Repeat annotation and masking](#21-repeat-annotation-and-masking)
    - [2.1.1 Repeat annotation using EDTA](#211-repeat-annotation-using-edta)
    - [2.1.2 Repeat annotation using RepeatMasker](#212-repeat-annotation-using-repeatmasker)
    - [2.1.3 Merge masking](#213-merge-masking)
    - [2.1.4 Comparison of masking tools](#214-comparison-of-masking-tools)
  - [2.2 Gene prediction](#22-gene-prediction)
    - [2.2.1 BRAKER3 using short RNA reads](#221-braker3-using-short-RNA-reads)
    - [2.2.2 BRAKER3 using protein and short RNA reads](#222-braker3-using-protein-and-short-RNA-reads)
    - [2.2.3. TSEBRA](#223-tsebra)
  - [2.3 Proteome taxonomic assessment](#23-proteome-taxonomic-assessment)
  - [2.4. Scaffold removal and final genome assembly](#24-scaffold-removal-and-final-genome-assembly)
  - [2.5. Reclassifying repeats ](#25-reclassifying-repeats)
  - [2.6. Telomeric and centromeric regions ](#26-telomeric-and-centromeric-regions)
  - [2.7. Functional annotation](#27-functional-annotation)
- [3. Gene family evolution in phylogenetic context](#3-gene-family-evolution-in-phylogenetic-context)
  - [3.1 Phylogenomic analysis](#31-phylogenomic-analysis)
     - [3.1.1 Orthofinder](#311-orthofinder)
     - [3.1.2 IQTREE and ASTRAL](#312-iqtree-and-astral)
     - [3.1.3 Phylogenetic tree](#313-phylogenetic-tree)
  - [3.2 Gene Family Evolution](#32-gene-familyevolution)
	- [3.2.1 Ultrametric tree](#321-ultrametric-tree)
    - [3.2.2 CAFE5 Analysis](#322-cafe5-analysis)
    - [3.2.3 Gene Family Expansions and Contractions](#323-gene-family-expansions-and-contractions)
    - [3.2.4 Annotation of Significant Familiess](#324-annotation-of-significant-familiess)

## Software

| Software | Version | Section |
|---|---|---|
| NanoQC / NanoComp | 0.9.4 / 1.20.0 | 1 |
| FastQC | – | 1 |
| FastK / Histex | 1.1.0 | 1 |
| PloidyPlot / Smudgeplot | 0.3.0 | 1 |
| GeneScopeFK | 1.0.0 | 1 |
| Hifiasm | 0.19.5 | 2 |
| Minimap2 | 2.28 (sections 2, 4); 2.22 (section 6) | 2, 4, 6 |
| Racon | 1.5.0 | 2 |
| LongStitch (Tigmint, ntLink, ARKS/LINKS) | 1.0.5 (LINKS 1.8.7) | 2 |
| BWA-MEM | 0.7.19 | 3 |
| SAMtools | 1.23 | 3 |
| YaHS | 1.2.2 | 3 |
| Python / NumPy / Matplotlib | 3.10.14 / 1.26.4 / 3.9.1 | 3 |
| QUAST | 5.3.0 | 4 |
| DIAMOND | 2.0.13 | 4, 7 |
| BlobTools | 1.0 | 4 |
| BlobToolKit | 2.0 | 4 |
| BUSCO (metazoa_odb10) | 5.8.2 | 4, 6 |
| EDTA | 2.0.0 | 5 |
| RepeatModeler / RepeatMasker | 2.0.6 / 4.1.8 | 5 |
| BEDTools | 2.30.0 | 5 |
| R VennDiagram | – | 5 |
| Entrez Direct | – | 6 |
| StringTie2 | – | 6 |
| BRAKER3 | 3.0.8 | 6 |
| TSEBRA | 1.1.2.5 | 6 |
| gffread | – | 6 |
| R dplyr / readr / stringr | – | 7 |
| eggNOG-mapper | 2.1.3 | 8 |
| InterProScan | 5.76-107.0 | 8 |
| OrthoFinder | 2.5.5 | 9 |
| MAFFT / ClipKIT / AMAS | 7.525 / 2.7.0 / 1.0 | 9 |
| IQ-TREE | 3.0.1 | 9 |
| GNU parallel | – | 9 |
| ASTER (ASTRAL-IV) | – | 9 |
| Newick Utilities | 1.6 | 9 |
| R ape | 5.8-1 | 10 |
| CAFE5 | – | 10 |
| CafePlotter | 0.2.0 | 10 |

– version not recorded.

## 1. Genome assembly and assessment

### 1.1 Read quality control

```bash
HIFI=hifi_ccs.fastq
OMNIC_R1=omnic_R1.fastq.gz
OMNIC_R2=omnic_R2.fastq.gz

nanoQC "$HIFI" -o NanoQC

NanoComp --fastq "$HIFI" "$OMNIC_R1" \
         --names HiFi_CCS OmniC \
         --outdir nanocomp_comparison

fastqc "$OMNIC_R1" "$OMNIC_R2"
```

### 1.2 k-mer profiling

```bash
HIFI=hifi_ccs.fastq

FastK -v -t32 -k21 -M90 -T4 "$HIFI" -NFastK_Table
PloidyPlot -e12 -k -v -T4 -okmerpairs FastK_Table
smudgeplot.py plot -t Pristina -o pristina kmerpairs_text.smu

Histex -G FastK_Table > pristina.hist
GeneScopeFK.R -i pristina.hist -o . -k 21
```

### 1.3 Genome assembly with Hifiasm

HiFi reads assembled in Hi-C integrated mode with the Omni-C reads.

```bash
HIFI=hifi_ccs.fastq
OMNIC_R1=omnic_R1.fastq.gz
OMNIC_R2=omnic_R2.fastq.gz

hifiasm -o pristina_ccs.asm -t 32 \
        --h1 "$OMNIC_R1" \
        --h2 "$OMNIC_R2" \
        "$HIFI"

# primary contigs to FASTA
awk '/^S/ {print ">"$2; print $3}' pristina_ccs.asm.hic.p_ctg.gfa > pristina_p_ctg.fasta
```

### 1.4 Genome assemby polishing

Three rounds of Racon, each using minimap2 alignments of the HiFi reads to the previous round.

```bash
HIFI=hifi_ccs.fastq
ASSEMBLY=pristina_p_ctg.fasta
THREADS=8

for i in 1 2 3; do
    minimap2 -ax map-pb -t "$THREADS" "$ASSEMBLY" "$HIFI" > alignment_iter${i}.sam
    racon -t "$THREADS" "$HIFI" alignment_iter${i}.sam "$ASSEMBLY" > racon_iter${i}.fasta
    ASSEMBLY=racon_iter${i}.fasta
done
```

### 1.5 Long-read scaffolding

Two rounds of LongStitch: Tigmint-long + ntLink, then Tigmint-long + ntLink + ARKS on the first-round output. LongStitch expects the draft as `<name>.fa` and the reads as `<name>.fq`.

```bash
ln -s racon_iter3.fasta racon_iter3.fa
ln -s hifi_ccs.fastq long_pacbio.ccs.reads.fq

# round 1: Tigmint-long + ntLink
longstitch tigmint-ntLink \
    draft=racon_iter3 \
    reads=long_pacbio.ccs.reads \
    G=600000 \
    longmap=hifi

# round 2: Tigmint-long + ntLink + ARKS
longstitch tigmint-ntLink-arks \
    draft=racon_iter3.cut250.tigmint.fa.k32.w100.z1000.ntLink.scaffolds \
    reads=long_pacbio.ccs.reads \
    G=600000 \
    longmap=hifi \
    c=4

cp racon_iter3.cut250.tigmint.fa.k32.w100.z1000.ntLink.scaffolds.cut250.tigmint.fa.k32.w100.z1000.ntLink.scaffolds_c4_m8-10000_cut250_k20_r0.05_e30000_z1000_l4_a0.3.scaffolds.fa \
   pristina.scaffolds.collapsed.fa
```

### 1.6 Assembly statistics

Assembly quality was assessed using QUAST

```bash
FASTA=pristina.scaffolds.collapsed.fa

quast.py "$FASTA" -o quast_results
```
Scaffold length statistics were also calculated from the final assembly, including total assembly length, scaffold count, average and largest scaffold length, and N50–N100 values.

```bash
FASTA=pristina.scaffolds.collapsed.fa
WRK_DIR=Working_directory

BASENAME=$(basename "$FASTA" .fasta)

# Create a lengths file: scaffold name + length, sorted descending
awk '/^>/{if(seqlen){print header"\t"seqlen}; header=$0; seqlen=0; next}{seqlen+=length($0)} END{print header"\t"seqlen}' "$FASTA" \
  | sort -k2,2nr > "$WRK_DIR/${BASENAME}.lengths.sorted.tsv

LENGTHS=$(awk '{print $2}' "$WRK_DIR/${BASENAME}.lengths.sorted.tsv")

# Total number of scaffolds
total_scaffolds=$(echo "$LENGTHS" | wc -l)

# Total genome size (sum)
total_length=$(echo "$LENGTHS" | paste -sd+ - | bc)

# Largest scaffold
largest=$(echo "$LENGTHS" | sort -nr | head -n1)

# Average length (with decimals)
average=$(awk -v sum=$total_length -v n=$total_scaffolds 'BEGIN {printf "%.2f", sum/n}')

# Median length
median_length=$(echo "$LENGTHS" | sort -n | awk -v n=$total_scaffolds '
  NR == int((n+1)/2) {med1=$1}
  NR == int((n+2)/2) {med2=$1}
  END {if (n%2==1) print med1; else print int((med1+med2)/2)}')

# Function to calculate Nx values
calc_Nx() {
  local x=$1
  awk -v x=$x -v total=$total_length '
  {len[NR]=$1}
  END {
    half = total * (x / 100)
    cumulative = 0
    for (i=1; i<=NR; i++) {
      cumulative += len[i]
      if (cumulative >= half) {
        print len[i], i
        exit
      }
    }
  }' < <(echo "$LENGTHS")
}

# Make sure lengths are sorted descending for Nx calculation
LENGTHS_SORTED=$(echo "$LENGTHS" | sort -nr)

# Calculate N50, N60, N70, N80, N90, N100
read N50 n50_n <<< $(calc_Nx 50 <<< "$LENGTHS_SORTED")
read N60 n60_n <<< $(calc_Nx 60 <<< "$LENGTHS_SORTED")
read N70 n70_n <<< $(calc_Nx 70 <<< "$LENGTHS_SORTED")
read N80 n80_n <<< $(calc_Nx 80 <<< "$LENGTHS_SORTED")
read N90 n90_n <<< $(calc_Nx 90 <<< "$LENGTHS_SORTED")
read N100 n100_n <<< $(calc_Nx 100 <<< "$LENGTHS_SORTED")

# N_count: count of scaffolds above 1000 bp (adjust as needed)
N_count=$(echo "$LENGTHS" | awk '$1 > 1000' | wc -l)

# Gaps: count of scaffolds of length exactly 1 bp (adjust if needed)
Gaps=$(echo "$LENGTHS" | awk '$1 == 1' | wc -l)

# Output to txt file
{
echo "Scaffold Length Summary Statistics"
echo "----------------------------------"
echo "Sum = $total_length, n = $total_scaffolds, ave = $average, largest = $largest"
echo ""
echo "N50 = $N50, n = $n50_n"
echo "N60 = $N60, n = $n60_n"
echo "N70 = $N70, n = $n70_n"
echo "N80 = $N80, n = $n80_n"
echo "N90 = $N90, n = $n90_n"
echo "N100 = $N100, n = $n100_n"
echo ""
echo "N_count = $N_count"
echo "Gaps = $Gaps"
} > "$WRK_DIR/${BASENAME}.stats.txt"
```

### 1.7 DIAMOND BLASTx

The assembly was split into six parts, with each part searched in a separate array job. Six array jobs were specified using #SBATCH --array=0-5 and submitted with sbatch --array=0-5 diamond.blastx.sh.

```bash
QUERY=pristina.scaffolds.collapsed.fa
DIAMOND_PATH=diamond
DB_PATH=nr.dmnd
OUTPUT_DIR=Output_directory
SPLIT_DIR=split_queries
SPLIT_OUT=split_output
PREFIX="pristina"

mkdir -p "$OUTPUT_DIR" "$SPLIT_DIR" "$SPLIT_OUT"

echo "----------------------------------------"
echo " STEP 1: Split query file (only once) "
echo "----------------------------------------"

if [ "$SLURM_ARRAY_TASK_ID" -eq 0 ]; then
    echo "Splitting query file..."
    num_files=6
    total_seqs=$(grep -c '^>' "$QUERY")
    seq_per_file=$(( total_seqs / num_files ))
    remainder=$(( total_seqs % num_files ))

    awk -v num_files="$num_files" -v seq_per_file="$seq_per_file" -v remainder="$remainder" -v outdir="$SPLIT_DIR" '
        BEGIN { file_idx=1; seq_count=0 }
        /^>/ {
            if (seq_count >= seq_per_file + (file_idx <= remainder ? 1 : 0)) {
                file_idx++
                seq_count=0
            }
            seq_count++
        }
        { print > (outdir "/query_part_" file_idx ".fasta") }
    ' "$QUERY"

    ls "$SPLIT_DIR"/query_part_*.fasta > "$SPLIT_DIR/parts_list.txt"
    echo "Splitting complete: $(wc -l < "$SPLIT_DIR/parts_list.txt") files created."
fi

# Wait a bit for file creation to propagate
sleep 60

echo "----------------------------------------"
echo " STEP 2: Run Diamond blastx "
echo "----------------------------------------"

PART_LIST_FILE="$SPLIT_DIR/parts_list.txt"
QUERY_FILE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$PART_LIST_FILE")
OUTFILE="$SPLIT_OUT/${PREFIX}_part${SLURM_ARRAY_TASK_ID}.out"

echo "------------------------------------------------"
echo "DIAMOND BLASTx job started for task $SLURM_ARRAY_TASK_ID on $(date '+%Y-%m-%d %H:%M:%S')"
echo "------------------------------------------------"

"$DIAMOND_PATH" blastx
--query "$QUERY_FILE" \
--db "$DB_PATH" \
--evalue 1e-5 \
--max-hsps 1 \
--max-target-seqs 1 \
--outfmt 6 qseqid bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
--out "$OUTFILE" 2 \
--more-sensitive \
--threads 32

echo "------------------------------------------------"
echo "Task $SLURM_ARRAY_TASK_ID completed: $OUTFILE"
echo "Finished at: $(date '+%Y-%m-%d %H:%M:%S')"
echo "------------------------------------------------"

```

### 1.8 Read coverage

```bash
QUERY=pristina.scaffolds.collapsed.fa
READS_FQ=hifi_ccs.fastq
OUTPUT_DIR=Output_directory

mkdir "$OUTPUT_DIR"
cd "$OUTPUT_DIR" || exit

# Extract base name of query file
BASENAME=$(basename "$QUERY" .fasta)

minimap2 -ax map-hifi "$QUERY" "$READS_FQ" -t 32 -o "$OUTPUT_DIR/${BASENAME}.sam" 2>&1 | tee ./alignments.log

# Convert SAM to BAM, sort and index
samtools view -S -b "$OUTPUT_DIR/${BASENAME}.sam" > "$OUTPUT_DIR/${BASENAME}.bam"
samtools sort "$OUTPUT_DIR/${BASENAME}.bam" -@ 32 -o "$OUTPUT_DIR/${BASENAME}.sorted.bam"
samtools index "$OUTPUT_DIR/${BASENAME}.sorted.bam"

echo "Alignment completed on $(date '+%Y-%m-%d %H:%M:%S')!"
echo "Results are stored in \"$OUTPUT_DIR\""
```

### 1.9 Taxify

Taxify was run on the DIAMOND results to assign taxonomy to the assembled sequences. The resulting taxonomic assignments were then used as input for BlobTools.

```bash
BLAST_OUTPUT=blastx.tsv
DB=prot.accession2taxid
OUTPUT_DIR=Output_directory

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR" || exit

blobtools taxify -f "$BLAST_OUTPUT" -c 2 -m "$DB" -s 1 -t 2

echo "Taxify job completed!"
```

### 1.10 BlobTools, BUSCO and snail plot (BlobToolkit)

DIAMOND BLASTx taxonomic assignments and read coverage were combined in BlobTools to assess the assembly. Assembly completeness was assessed using BUSCO, and key assembly metrics were summarised in a BlobToolKit snail plot.

```bash
QUERY=pristina.scaffolds.collapsed.fa
TAXIFY_OUTPUT=taxified.out
DB=nodesDB.txt
SORTED_BAM=genome.sorted.bam
OUTPUT_DIR=Output_directory
BUSCO_DIR=Busco_output_directory
BUSCO_DB=metazoa_odb10
SNAILPLOT_DIR=Snailplot_directory
BLOB_DIR=BlobDir_directory

mkdir "$OUTPUT_DIR" "$BUSCO_DIR" "$BLOB_DIR"
cd "$OUTPUT_DIR" || exit

echo "------------------------------------------------"
echo " 1.BLOBTOOLS PLOTTING (v1)"
echo "------------------------------------------------"

# Extract base name of query file
BASENAME=$(basename "$QUERY" .fasta)

echo "Creating blobplot"
blobtools create -i "$QUERY" -t "$TAXIFY_OUTPUT" --db "$DB" -b "$SORTED_BAM" --title ${BASENAME}.blobplot -o "$OUTPUT_DIR"/${BASENAME}.blobplot

echo "Generating tabular summary"
blobtools view -i "${OUTPUT_DIR}/${BASENAME}.blobplot.blobDB.json" -r phylum -o "${OUTPUT_DIR}/blobplot.Phylum_cov" --cov

echo "Plotting: Single plot per phylum"
blobtools plot -i "${OUTPUT_DIR}/${BASENAME}.blobplot.blobDB.json" -r phylum -o "${OUTPUT_DIR}/blobplot.Phylum_comb" --format pdf

echo "Plotting: Multiplot"
blobtools plot -i "${OUTPUT_DIR}/${BASENAME}.blobplot.blobDB.json" -r phylum -o "${OUTPUT_DIR}/blobplot.Phylum_multi" --format pdf --multiplot

echo "Blobplots completed!"

echo "------------------------------------------------"
echo "2. BUSCO ANALYSIS"
echo "------------------------------------------------"

busco -i "$QUERY" -l "$BUSCO_DB" -o "${BASENAME}_BUSCO" -m geno -c 32
mv "${BASENAME}_BUSCO" "$BUSCO_DIR"

echo "------------------------------------------------"
echo "3. SNAILPLOT USING BLOBTOOLKIT (v2)"
echo "------------------------------------------------"

cd "$SNAILPLOT_DIR"

blobtools create --fasta "$QUERY" "$BLOB_DIR"
blobtools add --busco "$BUSCO_DIR/${BASENAME}_BUSCO/run_metazoa_odb10/full_table.tsv" "$BLOB_DIR"
blobtools view --plot --view snail "$BLOB_DIR"

echo "Snail plot created in $SNAILPLOT_DIR"

echo "All analyses completed successfully on $(date '+%Y-%m-%d %H:%M:%S')"
```

### 1.11 Hi-C scaffolding, Alignment and YaHS

Omni-C read files were aligned independently, then merged. Mate information was added and duplicates removed before scaffolding with YaHS.

```bash
ASSEMBLY=pristina.scaffolds.collapsed.fa
HIC_R1=omnic_R1.fastq.gz
HIC_R2=omnic_R2.fastq.gz
OUTDIR=yahs_output
THREADS=32

mkdir -p "$OUTDIR"
bwa index "$ASSEMBLY"
samtools faidx "$ASSEMBLY"

bwa mem -t "$THREADS" -5SP "$ASSEMBLY" "$HIC_R1" | samtools view -@ "$THREADS" -bS - > "$OUTDIR/hic_R1.bam"
bwa mem -t "$THREADS" -5SP "$ASSEMBLY" "$HIC_R2" | samtools view -@ "$THREADS" -bS - > "$OUTDIR/hic_R2.bam"

samtools merge   -@ "$THREADS" "$OUTDIR/hic_combined.bam" "$OUTDIR/hic_R1.bam" "$OUTDIR/hic_R2.bam"
samtools sort    -@ "$THREADS" -n -o "$OUTDIR/hic_sorted.bam" "$OUTDIR/hic_combined.bam"
samtools fixmate -@ "$THREADS" -m "$OUTDIR/hic_sorted.bam" "$OUTDIR/hic_fixmate.bam"
samtools sort    -@ "$THREADS" -o "$OUTDIR/hic_positionsort.bam" "$OUTDIR/hic_fixmate.bam"
samtools markdup -@ "$THREADS" -r "$OUTDIR/hic_positionsort.bam" "$OUTDIR/hic_markdup.bam"
samtools index   "$OUTDIR/hic_markdup.bam"

yahs -o "$OUTDIR/yahs_out" "$ASSEMBLY" "$OUTDIR/hic_markdup.bam"
samtools faidx "$OUTDIR/yahs_out_scaffolds_final.fa"
```

### 1.12 Scaffold size filter

```bash
SCAFFOLDS=yahs_output/yahs_out_scaffolds_final.fa
MIN_SIZE=50000
OUT=pristina_scaffolds_filtered_50kb.fa

awk -v min="$MIN_SIZE" '$2 >= min {print $1}' "${SCAFFOLDS}.fai" > scaffolds_to_keep.txt
samtools faidx "$SCAFFOLDS" $(cat scaffolds_to_keep.txt) > "$OUT"
samtools faidx "$OUT"
```

### 1.13 Contact map

Contacts from the Hi-C alignments were placed on the scaffolds through the YaHS AGP file and binned at 100 kb for the 21 scaffolds longer than 10 Mb.

```python
#!/usr/bin/env python3
import subprocess

import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np

FAI = "pristina_scaffolds_filtered_50kb.fa.fai"
AGP = "yahs_output/yahs_out_scaffolds_final.agp"
BAM = "yahs_output/hic_markdup.bam"
N_SCAFFOLDS = 21
RESOLUTION = 100_000
MAX_READS = 20_000_000
OUT_PNG = "pristina_hic_contact_map_21_chromosomes.png"

# Scaffold order and genome-wide offsets
boundaries, scaffold_order, cumulative = {}, [], 0
with open(FAI) as f:
    for i, line in enumerate(f):
        if i >= N_SCAFFOLDS:
            break
        name, size = line.split("\t")[:2]
        boundaries[name] = (cumulative, cumulative + int(size))
        scaffold_order.append(name)
        cumulative += int(size)
scaffold_set = set(scaffold_order)

# Contig-to-scaffold coordinates from the AGP
contig_mapping = {}
with open(AGP) as f:
    for line in f:
        if line.startswith("#"):
            continue
        p = line.strip().split("\t")
        if len(p) < 9 or p[4] != "W":
            continue
        contig_mapping.setdefault(p[5], []).append(
            {"scaffold": p[0], "start": int(p[1]), "old_start": int(p[6]), "old_end": int(p[7])}
        )


def map_position(contig, pos):
    for m in contig_mapping.get(contig, []):
        if m["old_start"] <= pos <= m["old_end"]:
            return m["scaffold"], m["start"] + pos - m["old_start"]
    return None, None


# Contact matrix
n_bins = int(np.ceil(cumulative / RESOLUTION))
contact_matrix = np.zeros((n_bins, n_bins))

proc = subprocess.Popen(f"samtools view {BAM} | head -{MAX_READS}",
                        shell=True, stdout=subprocess.PIPE, text=True)
for line in proc.stdout:
    fields = line.strip().split("\t")
    if len(fields) < 8:
        continue
    contig1, pos1 = fields[2], int(fields[3])
    contig2 = contig1 if fields[6] == "=" else fields[6]
    try:
        pos2 = int(fields[7])
    except ValueError:
        continue
    if pos2 == 0:
        continue

    scaffold1, spos1 = map_position(contig1, pos1)
    scaffold2, spos2 = map_position(contig2, pos2)
    if scaffold1 not in scaffold_set or scaffold2 not in scaffold_set:
        continue

    bin1 = int((boundaries[scaffold1][0] + spos1 - 1) / RESOLUTION)
    bin2 = int((boundaries[scaffold2][0] + spos2 - 1) / RESOLUTION)
    if 0 <= bin1 < n_bins and 0 <= bin2 < n_bins:
        contact_matrix[bin1, bin2] += 1
        contact_matrix[bin2, bin1] += 1
proc.wait()

# Plot
contact_matrix[contact_matrix == 0] = np.nan
fig, ax = plt.subplots(figsize=(14, 12))
im = ax.imshow(contact_matrix, cmap="YlOrRd",
               norm=colors.LogNorm(vmin=1, vmax=np.nanmax(contact_matrix)),
               interpolation="nearest", origin="upper")

edges = [0] + [boundaries[s][1] / RESOLUTION for s in scaffold_order]
for pos in edges:
    ax.axhline(y=pos, color="blue", linewidth=1, alpha=0.4)
    ax.axvline(x=pos, color="blue", linewidth=1, alpha=0.4)

ticks = [(edges[i] + edges[i + 1]) / 2 for i in range(N_SCAFFOLDS)]
labels = [str(i + 1) for i in range(N_SCAFFOLDS)]
ax.set_xticks(ticks)
ax.set_xticklabels(labels, fontsize=10)
ax.set_yticks(ticks)
ax.set_yticklabels(labels, fontsize=10)
ax.set_xlabel("Chromosome", fontsize=13, fontweight="bold")
ax.set_ylabel("Chromosome", fontsize=13, fontweight="bold")
ax.set_title("Pristina leidyi Hi-C Contact Map\n21 Main Chromosomes (>10Mb)",
             fontsize=15, fontweight="bold", pad=15)

cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
cbar.set_label("Hi-C contacts (log₁₀)", rotation=270, labelpad=20, fontsize=11)

plt.tight_layout()
plt.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
```

### 1.14 Retain longest scaffolds

The 26 longest scaffolds were retained for downstream analysis. 

```bash
awk '/^>/{n++} n<=26' pristina_scaffolds_filtered_50kb.fa > genome.fasta
```

### 1.15 BLASTn

The assembly was reassessed for potential contamination against the NCBI nucleotide (nt) database. Steps 1.8 and 1.10 were repeated before downstream analysis.

```bash
QUERY=genome.fasta
DB=nt
OUTPUT_DIR=Output_directory

cd "$OUTPUT_DIR" || exit

echo "Starting blastn at $(date)"
blastn \
-query "$QUERY" \
-db "$DB" \
-evalue 1e-5 \
-max_target_seqs 1 \
-outfmt "6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen std" \
-out "$OUTPUT_DIR/blastn.tsv" \
-num_threads 32
```

### 1.16 Rename headers

The scaffold headers were renamed for simplicity before downstream analysis. A mapping file was generated to retain the correspondence between the new and original scaffold names.

```bash
input=genome.fasta
output=genome.sc26.fasta
mapfile=genome.sc26.fasta.idmap.tsv

awk 'BEGIN{FS="\n"; OFS="\n"} 
/^>/ {
  old=$0;
  gsub(/^>/,"",old);
  i++;
  new="scaffold_" i;
  print ">"new >> "'"$output"'";
  print new "\t" old >> "'"$mapfile"'";
  next
}
{
  print $0 >> "'"$output"'"
}' "$input"
```
## 2. Genome content and functional annotation

### 2.1 Repeat annotation and masking

#### 2.1.1. Repeat annotation using EDTA

```bash
QUERY=genome.sc26.fasta

EDTA.pl --genome "$QUERY" --sensitive 1 --anno 1 --threads 32 2>&1 | tee EDTA.log
```

#### 2.1.2. Repeat annotation using RepeatMasker

Repeat annotation using RepeatMasker and a *de novo* repeat library generated with RepeatModeler.

```bash
# Define paths
QUERY=genome.sc26.fasta
OUTPUT_DIR=Output_directory

mkdir "$OUTPUT_DIR"
cd "$OUTPUT_DIR" || exit 1

echo  "Building a Database from Genome Assembly"
BuildDatabase -name db "$QUERY"

echo "Now running RepeatModeler"
RepeatModeler -database db -threads 32 -LTRStruct 2>&1 | tee RepeatModeler.log

echo "---------------------------------------------------------"
echo "RepeatModeler completed on $(date '+%Y-%m-%d %H:%M:%S')"
echo "Results are stored in \"$OUTPUT_DIR\""
echo "---------------------------------------------------------"
```

```bash
REPEATMASKER_DB=db-families.fa
QUERY=genome.sc26.fasta
OUTPUT_DIR=Output_directory

# Extract base name of query file
BASENAME=$(basename "$QUERY" .fasta)
PREFIX="${BASENAME}_RepeatMasker_RepeatMod"

mkdir $OUTPUT_DIR" 
cd "$OUTPUT_DIR" 

echo "---------------------------------------------"
echo "Job started on $(date '+%Y-%m-%d %H:%M:%S')"
echo "Processing genome: $BASENAME"
echo "---------------------------------------------"

# Run RepeatMasker with databases
RepeatMasker -lib "$REPEATMASKER_DB" "$QUERY" -dir "$OUTPUT_DIR" 2>&1 | tee "$OUTPUT_DIR/RepeatMa_RepeatMod_run.log"

echo "---------------------------------------------"
echo "Renaming output files to use prefix: $PREFIX"
echo "---------------------------------------------"

# Loop through expected outputs
for ext in out tbl masked cat.gz; do
    # Find first matching file (RepeatMasker should only output one set)
    file=$(ls *.$ext 2>/dev/null | head -n 1 || true)

    # Only proceed if file exists
    if [[ -n "$file" && -f "$file" ]]; then
        new_name="${PREFIX}.${ext}"

        # Double-check no accidental overwrite
        if [[ -e "$new_name" ]]; then
            echo "Skipping rename for $file — target $new_name already exists."
        else
            echo "Renaming: $file → $new_name"
            mv -- "$file" "$new_name"
        fi
    fi
done

echo "---------------------------------------------"
echo "RepeatMasker using de novo RepeatModeler database completed on $(date '+%Y-%m-%d %H:%M:%S')"
echo "Results stored in: $OUTPUT_DIR"
echo "---------------------------------------------"
```

Repeat annotation using RepeatMasker and the Dfam 3.9 database (partitions 15 and 16).

```bash
REPEATMASKER_DB=Dfam.lib
QUERY=genome.sc26.fasta
OUTPUT_DIR=Output_directory

# Extract base name of query file
BASENAME=$(basename "$QUERY" .fasta)
PREFIX="${BASENAME}_RepeatMasker_Dfam"

mkdir $OUTPUT_DIR" 
cd "$OUTPUT_DIR" || exit 1

echo "-----------------------------------------------------------------"
echo "RepeatMasker using Dfam started on $(date '+%Y-%m-%d %H:%M:%S')"
echo "Processing genome: $BASENAME"
echo "-----------------------------------------------------------------"

# Run RepeatMasker with databases
RepeatMasker -lib "$REPEATMASKER_DB" "$QUERY" -dir "$OUTPUT_DIR" 2>&1 | tee "$OUTPUT_DIR/RepeatMa_Dfam_run.log"

echo "----------------------------------------------"
echo "Renaming output files to use prefix: $PREFIX"
echo "----------------------------------------------"

# Loop through expected outputs
for ext in out tbl masked cat.gz; do
    # Find first matching file (RepeatMasker should only output one set)
    file=$(ls *.$ext 2>/dev/null | head -n 1 || true)

    # Only proceed if file exists
    if [[ -n "$file" && -f "$file" ]]; then
        new_name="${PREFIX}.${ext}"

        # Double-check no accidental overwrite
        if [[ -e "$new_name" ]]; then
            echo "Skipping rename for $file — target $new_name already exists."
        else
            echo "Renaming: $file → $new_name"
            mv -- "$file" "$new_name"
        fi
    fi
done

echo "---------------------------------------------"
echo "RepeatMasker using Dfam database completed on $(date '+%Y-%m-%d %H:%M:%S')"
echo "Results stored in: $OUTPUT_DIR"
echo "---------------------------------------------"
```

#### 2.1.3. Merge masking

The three outputs from the masking steps were merged into a single non-overlapping masked genome.
The inputs required to generate a Venn diagram comparing the three tools were also produced within this script.
Note: BEDTools is required to run this script.

```bash
exec > >(tee -a merged_masking.log) 2>&1

REPEATMOD_OUT=RepeatMasker_RepeatMod.out
DFAM_OUT=RepeatMasker_Dfam.out
EDTA_OUT=EDTA.TEanno.out
QUERY=genome.sc26.fasta
OUTPUT_DIR=Output_directory

# Extract base name of query file
NAME=$(basename "$QUERY" .fasta.mod)

echo "Genome name detected: $NAME"

mkdir "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

echo "All outputs will be written to: $OUTPUT_DIR"

# ==== STEP 1. Convert to BED format ====
echo "[1/5] Converting outputs to BED..."
awk 'BEGIN {OFS="\t"} $1 !~ /^SW|^perc|^score|^sequence|^ID/ && NF > 1 {print $5, $6-1, $7, $10}' "$REPEATMOD_OUT" > "${NAME}_RepeatMod.bed"
awk 'BEGIN {OFS="\t"} $1 !~ /^SW|^perc|^score|^sequence|^ID/ && NF > 1 {print $5, $6-1, $7, $10}' "$DFAM_OUT" > "${NAME}_dfam.bed"
awk 'BEGIN{OFS="\t"} !/^#/ {print $5, $6-1, $7, $11}' "$EDTA_OUT" > "${NAME}_edta.bed"

# ==== STEP 2. Merge overlapping regions ====
echo "[2/5] Merging overlapping regions..."
bedtools sort -i "${NAME}_RepeatMod.bed" | bedtools merge -i - -c 4 -o distinct > "${NAME}_RepeatMod_merged.bed"
bedtools sort -i "${NAME}_dfam.bed" | bedtools merge -i - -c 4 -o distinct > "${NAME}_dfam_merged.bed"
bedtools sort -i "${NAME}_edta.bed" | bedtools merge -i - -c 4 -o distinct > "${NAME}_edta_merged.bed"

# ==== STEP 3. Compare overlaps ====
echo "[3/5] Running bedtools multiinter..."
bedtools multiinter -i "${NAME}_RepeatMod_merged.bed" "${NAME}_dfam_merged.bed" "${NAME}_edta_merged.bed" > "${NAME}_repeats_overlap.bed"

# ==== STEP 4. Summarize for Venn ====
echo "[4/5] Summarizing overlaps for Venn..."
mkdir -p Summaries_Venn_Diagram
cd Summaries_Venn_Diagram

awk '$4 == 1 && $5 == "1"' ../"${NAME}_repeats_overlap.bed" > "${NAME}_Repeatmodeler_only.bed"
awk '$4 == 1 && $5 == "2"' ../"${NAME}_repeats_overlap.bed" > "${NAME}_dfam_only.bed"
awk '$4 == 1 && $5 == "3"' ../"${NAME}_repeats_overlap.bed" > "${NAME}_edta_only.bed"
awk '$4 == 3' ../"${NAME}_repeats_overlap.bed" > "${NAME_RAW}_Repeatmodeler.dfam.edta.bed"
awk '$4 == 2 && $5 ~ /1/ && $5 ~ /2/' ../"${NAME}_repeats_overlap.bed" > "${NAME}_Repeatmodeler_dfam.bed"
awk '$4 == 2 && $5 ~ /1/ && $5 ~ /3/' ../"${NAME}_repeats_overlap.bed" > "${NAME}_Repeatmodeler_edta.bed"
awk '$4 == 2 && $5 ~ /2/ && $5 ~ /3/' ../"${NAME}_repeats_overlap.bed" > "${NAME}_dfam_edta.bed"

cd ..

# ==== STEP 5. Final merged masked genome ====
echo "[5/5] Building merged masked genome..."

cat "${NAME}_RepeatMod_merged.bed" "${NAME}_dfam_merged.bed" "${NAME}_edta_merged.bed" | bedtools sort -i - > all_sorted.bed

bedtools merge -i all_sorted.bed -c 4 -o distinct > "${NAME}.merged_masked_genome.bed"
bedtools maskfasta -fi "$QUERY" -bed "${NAME}.merged_masked_genome.bed" -fo "${NAME}.merged.hardmasked_genome.fasta"
bedtools maskfasta -fi "$QUERY" -bed "${NAME}.merged_masked_genome.bed" -fo "${NAME}.merged.softmasked_genome.fasta" -soft

echo "Completed on $(date '+%Y-%m-%d %H:%M:%S'). All results are inside: $OUTPUT_DIR"

mv ../merged_masking.log  "$OUTPUT_DIR"
```

#### 2.1.4. Comparison of masking tools

A Venn diagram was generated in R statistical software to compare the outputs from the three tools.

```r
install.packages("VennDiagram")
library(VennDiagram)
library(grid)
library(dplyr)

# Count regions
only_modeler <- length(readLines("Repeatmodeler_only.bed"))
only_dfam <- length(readLines("dfam_only.bed"))
only_edta <- length(readLines("edta_only.bed"))

modeler_dfam <- length(readLines("Repeatmodeler_dfam.bed"))
modeler_edta <- length(readLines("Repeatmodeler_edta.bed"))
dfam_edta <- length(readLines("dfam_edta.bed"))

all_three <- length(readLines("Repeatmodeler.dfam.edta.bed"))

# Original count labels as characters (used to match regions)
original_labels <- as.character(c(
  only_modeler,
  only_dfam,
  modeler_dfam,
  only_edta,
  modeler_edta,
  dfam_edta,
  all_three
))

total <- sum(c(only_modeler, only_dfam, modeler_dfam, only_edta, modeler_edta, dfam_edta, all_three))

# Custom formatted labels
custom_labels <- sprintf("%d\n(%.1f%%)",
                         c(only_modeler, only_dfam, modeler_dfam, only_edta, modeler_edta, dfam_edta, all_three),
                         100 * c(only_modeler, only_dfam, modeler_dfam, only_edta, modeler_edta, dfam_edta, all_three) / total)

# Draw the Venn diagram and hide the default text
grid.newpage()
venn_plot <- draw.triple.venn(
  area1 = only_modeler + modeler_dfam + modeler_edta + all_three,
  area2 = only_dfam + modeler_dfam + dfam_edta + all_three,
  area3 = only_edta + modeler_edta + dfam_edta + all_three,
  n12 = modeler_dfam + all_three,
  n13 = modeler_edta + all_three,
  n23 = dfam_edta + all_three,
  n123 = all_three,
  category = c("RepeatModeler", "Dfam", "EDTA"),
  fill = c("skyblue", "pink1", "palegreen"),
  #cex = 2,        # we will override these
  cat.cex = 1.5
)

# Replace labels by matching the original values
for (i in seq_along(venn_plot)) {
  if (grepl("text", venn_plot[[i]]$name) &&
      venn_plot[[i]]$label %in% original_labels) {
    index <- match(venn_plot[[i]]$label, original_labels)
    venn_plot[[i]]$label <- custom_labels[[index]]
    venn_plot[[i]]$gp <- gpar(fontsize = 24)
  }
}

# Redraw with white background
grid.newpage()
grid.rect(gp = gpar(fill = "white", col = NA))  # <- Add white rectangle

grid.draw(gList(venn_plot))
```

### 2.2. Gene prediction
#### 2.2.1. BRAKER3 using short RNA reads

Short-read RNA-sequencing data from Álvarez-Campos (2024) were aligned to the genome for use as evidence in BRAKER3.

```bash
echo "-----------------------------------------"
echo " SHORT RNA READS ALIGNMENT TO GENOME "
echo "-----------------------------------------"

QUERY=genome.sc26.softmasked.fasta
INPUT_DIR=Short_RNA_directory
OUTPUT_DIR=Output_directory

mkdir -p "$OUTPUT_DIR"
cd "$INPUT_DIR" || exit

# Array to store all sorted BAM paths for merging
SORTED_BAMS=()

# Loop through all _1.fastq files
for R1 in *_1.fastq; do
  R2="${R1/_1.fastq/_2.fastq}"

  if [[ -f "$R2" ]]; then
    SAMPLE_NAME="${R1%%_1.fastq}"
    SAM="${OUTPUT_DIR}/${SAMPLE_NAME}.sam"
    BAM="${OUTPUT_DIR}/${SAMPLE_NAME}.bam"
    SORTED_BAM="${OUTPUT_DIR}/${SAMPLE_NAME}.sorted.bam"
    LOG="${OUTPUT_DIR}/${SAMPLE_NAME}.output.log"

    echo "Aligning $R1 and $R2 to genome..."
    minimap2 -ax splice -t 32 "$QUERY" "$R1" "$R2" -o "$SAM" 2>&1 | tee "$LOG"

    echo "Converting SAM to BAM for $SAMPLE_NAME..."
    samtools view -@ 32 -bS "$SAM" > "$BAM"

    echo "Sorting BAM file for $SAMPLE_NAME..."
    samtools sort -@ 32 "$BAM" -o "$SORTED_BAM"

    echo "Indexing sorted BAM file for $SAMPLE_NAME..."
    samtools index "$SORTED_BAM"

    # Add to array for merging
    SORTED_BAMS+=("$SORTED_BAM")

    echo "Done with $SAMPLE_NAME. Final file: $SORTED_BAM"
  else
    echo "Warning: No matching _2.fastq file found for $R1"
  fi
done

# Merge all sorted BAMs
if [[ ${#SORTED_BAMS[@]} -gt 0 ]]; then
  MERGED_BAM="${OUTPUT_DIR}/merged.bam"
  MERGED_SORTED_BAM="${OUTPUT_DIR}/merged.sorted.bam"

  echo "Merging all sorted BAM files..."
  samtools merge -@ 32 "$MERGED_BAM" "${SORTED_BAMS[@]}"

  echo "Sorting merged BAM file..."
  samtools sort -@ 32 "$MERGED_BAM" -o "$MERGED_SORTED_BAM"

  echo "Indexing merged and sorted BAM file..."
  samtools index "$MERGED_SORTED_BAM"

  echo "Cleaning up intermediate files (.sam and .bam)..."
  find "$OUTPUT_DIR" -type f \( -name "*.sam" -o -name "*.bam" \) ! -name "*.sorted.bam" -delete

  echo "Cleanup complete. Kept only sorted BAMs and merged outputs."
  echo "Completed on $(date '+%Y-%m-%d %H:%M:%S'). Final merged and indexed BAM: $MERGED_SORTED_BAM"
else
  echo "No BAM files to merge."
fi

echo "Completed on $(date '+%Y-%m-%d %H:%M:%S'). All files are inside: $OUTPUT_DIR"
```

Short-read RNA-seq data were used as evidence in BRAKER3. BUSCO was used to assess the completeness of the predicted proteome.

```bash
AUGUSTUS_CONFIG_PATH=config
SINGULARITY=braker3.sif
QUERY=genome.sc26.softmasked.fasta
BAM_FILES=merged.sorted.bam
OUTPUT_DIR=Output_directory
BUSCO_DIR="$OUTPUT_DIR/BUSCO"			
BUSCO_DB=metazoa_odb10

mkdir -p "$OUTPUT_DIR" "$BUSCO_DIR"
cd "$OUTPUT_DIR" || exit
BUSCO_OUTPUT_NAME="$(basename $OUTPUT_DIR)_BUSCO"

echo "-----------------------------------------"
echo " 1. BRAKER3 USING RNA EVIDENCE"
echo "-----------------------------------------"

echo "Starting BRAKER3 at $(date)"
singularity exec --env AUGUSTUS_CONFIG_PATH="$AUGUSTUS_CONFIG_PATH" "$SINGULARITY" braker.pl --genome="$QUERY" --bam="$BAM_FILES" --softmasking --threads 32

echo "Gene prediction completed at $(date). See results in $OUTPUT_DIR/braker"

echo "-----------------------------------------"
echo " 2. CLEAN UP PROTEIN FILE FOR BUSCO"
echo "-----------------------------------------"

RAW_PROTEOME="$OUTPUT_DIR/braker/braker.aa"
CLEAN_PROTEOME="$OUTPUT_DIR/braker/braker_clean.aa"

# Remove asterisks from the protein sequences to run BUSCO
sed 's/\*//g' "$RAW_PROTEOME" > "$CLEAN_PROTEOME"

echo "Cleaned protein file created: $CLEAN_PROTEOME"

echo "-----------------------------------------"
echo " 3. BUSCO ANALYSIS"
echo "-----------------------------------------"

busco -i "$CLEAN_PROTEOME" -l "$BUSCO_DB" -o "$BUSCO_OUTPUT_NAME" -m prot -c 32

mv "${BASENAME}_BUSCO" "$BUSCO_DIR"

echo "-----------------------------------------"
echo " All analyses completed successfully on $(date '+%Y-%m-%d %H:%M:%S')"
echo "-----------------------------------------"
```

#### 2.2.2. BRAKER3 using protein and short RNA reads

Protein sequences from the class Clitellata (phylum: Annelida) were retrieved from the NCBI database for use as evidence in BRAKER3.

```bash
SEARCH_TERM='Clitellata[Organism]'
BATCH_SIZE=500
OUTFILE=clitellata_proteins.fasta

# -----------------------------------------
# Step 1: Get all protein IDs
# -----------------------------------------

echo "Searching NCBI for: $SEARCH_TERM"

esearch -db protein -query "$SEARCH_TERM" | \
    efetch -format uid > "$ID_FILE"

TOTAL_IDS=$(wc -l < "$ID_FILE")

echo "Found $TOTAL_IDS protein records."

# -----------------------------------------
# Step 2: Fetch protein sequences in batches
# -----------------------------------------

Create the output file
> "$OUTFILE"

# Split the protein IDs into batches
split -l "$BATCH_SIZE" "$ID_FILE" ids_batch_

# Fetch each batch and append it to the output FASTA
for f in ids_batch_*; do

    echo "Fetching batch: $f"

    BATCH_IDS=$(paste -sd, "$f")

    efetch \
        -db protein \
        -id "$BATCH_IDS" \
        -format fasta > "$TMP_FILE"

    cat "$TMP_FILE" >> "$OUTFILE"

    # Pause between NCBI requests
    sleep 1

done

# -----------------------------------------
# Step 3: Clean up temporary files
# -----------------------------------------

rm ids_batch_* "$TMP_FILE" "$ID_FILE"
echo "Download complete: $OUTFILE"
```

Short-read RNA-seq data from the previous step and protein sequences were used together as evidence in BRAKER3. BUSCO was used to assess the completeness of the predicted proteome.

```bash
AUGUSTUS_CONFIG_PATH=config
SINGULARITY=braker3.sif
QUERY=genome.sc26.softmasked.fasta
BAM_FILES=merged.sorted.bam
PROT=clitellata_proteins.fasta
OUTPUT_DIR=Output_directory
BUSCO_DIR="$OUTPUT_DIR/BUSCO"			
BUSCO_DB=metazoa_odb10

mkdir -p "$OUTPUT_DIR" "$BUSCO_DIR"
cd "$OUTPUT_DIR" || exit
BUSCO_OUTPUT_NAME="$(basename $OUTPUT_DIR)_BUSCO"

echo "-----------------------------------------"
echo " 1. BRAKER3 USING RNA & PROTEIN EVIDENCE"
echo "-----------------------------------------"

echo "Starting BRAKER3 at $(date)"
singularity exec --env AUGUSTUS_CONFIG_PATH="$AUGUSTUS_CONFIG_PATH" "$SINGULARITY" braker.pl --genome="$QUERY" --bam="$BAM_FILES" --prot_seq="$PROT" --softmasking --threads 32 

echo "Gene prediction completed at $(date). See results in $OUTPUT_DIR/braker"

echo "-----------------------------------------"
echo " 2. CLEAN UP PROTEIN FILE FOR BUSCO"
echo "-----------------------------------------"

RAW_PROTEOME="$OUTPUT_DIR/braker/braker.aa"
CLEAN_PROTEOME="$OUTPUT_DIR/braker/braker_clean.aa"

# Remove asterisks from the protein sequences to run BUSCO
sed 's/\*//g' "$RAW_PROTEOME" > "$CLEAN_PROTEOME"

echo "Cleaned protein file created: $CLEAN_PROTEOME"

echo "-----------------------------------------"
echo " 3. BUSCO ANALYSIS"
echo "-----------------------------------------"

busco -i "$CLEAN_PROTEOME" -l "$BUSCO_DB" -o "$BUSCO_OUTPUT_NAME" -m prot -c 32

mv "${BASENAME}_BUSCO" "$BUSCO_DIR"

echo "-----------------------------------------"
echo " All analyses completed successfully on $(date '+%Y-%m-%d %H:%M:%S')"
echo "-----------------------------------------"
```
#### 2.2.3. TSEBRA

The outputs from BRAKER3 were combined with long-read sequencing data in TSEBRA, with priority given to RNA-sequencing evidence.

Long reads were aligned to the genome for use as evidence in TSEBRA.

```bash
echo "-----------------------------------------"
echo " LONG RNA READS ALIGNMENT TO GENOME "
echo "-----------------------------------------"

QUERY=genome.sc26.softmasked.fasta
READS_FQ=long_RNA_reads.fasta
OUTPUT_DIR=Output_directory

# Extract base name of query file
BASENAME=$(basename "$QUERY" .fasta)

mkdir "$OUTPUT_DIR"
cd "$OUTPUT_DIR" || exit

minimap2 -ax map-pb "$QUERY" "$READS_FQ" -t 32 -o "$OUTPUT_DIR/${BASENAME}.sam" 2>&1 | tee "$OUTPUT_DIR/${BASENAME}.log"

# Convert SAM to BAM, sort and index
samtools view -S -b "$OUTPUT_DIR/${BASENAME}.sam" > "$OUTPUT_DIR/${BASENAME}.bam"
samtools sort "$OUTPUT_DIR/${BASENAME}.bam" -@ 32 -o "$OUTPUT_DIR/${BASENAME}.sorted.bam"
samtools index "$OUTPUT_DIR/${BASENAME}.sorted.bam"

echo "Alignment completed on $(date '+%Y-%m-%d %H:%M:%S')!"

echo "Creating gtf for TSEBRA"

stringtie "$OUTPUT_DIR/${BASENAME}.sorted.bam" -p 32 -o "$OUTPUT_DIR/${BASENAME}.gtf"

echo "Results are stored in \"$OUTPUT_DIR\""2
```

The tsebra.conf configuration file was used to run TSEBRA.

```text
# Weight for each hint source
P 0.9
E 1
C 0.7
M 0.4
# Required fraction of supported introns or supported start/stop codons for a transcript
intron_support 0.8
stasto_support 1
# Allowed difference for each feature
e_1 0.0
e_2 0.5
e_3 0.096
e_4 0.02
e_5 0.18
e_6 0.18
```

TSEBRA

```bash
echo "--------------------------------------------------"
echo " TSEBRA USING BRAKER3 OUTPUT AND LONG READ RNA"
echo "--------------------------------------------------"

GTF_1=shortRNA.braker.gtf	
GTF_2=shortRNA.prot.braker.gtf

HINTFILE_GFF_1=shortRNA.hintsfile.gff
HINTFILE_GFF_2=shortRNA.prot.hintsfile.gff"

RNA_GTF=longRNA.gtf
CONF_FILE=tsebra.con

OUTPUT_DIR=Output_directory

PREFIX="Pristina"
OUTPUT_FILE="$OUTPUT_DIR/TSEBRA_${PREFIX}.gtf"

mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR" || exit

# Check if output file exists and create a unique name
if [[ -e "$OUTPUT_FILE" ]]; then
    i=1
    while [[ -e "${OUTPUT_DIR}/TSEBRA_${PREFIX}_${i}.gtf" ]]; do
        ((i++))
    done
    OUTPUT_FILE="${OUTPUT_DIR}/TSEBRA_${PREFIX}_${i}.gtf"
    echo "Output file exists. Using new filename: $OUTPUT_FILE"
fi

# Add the long read RNA using the -k option
echo "Starting TSEBRA at $(date)"
echo "Using GTFs: $GTF_1 and $GTF_2"
echo "Using hints: $HINTFILE_GFF_1 and $HINTFILE_GFF_2"
echo "Using long read RNA GTF: $RNA_GTF"
echo "Using config: $CONF_FILE"

tsebra.py --gtf "$GTF_1","$GTF_2" --hintfiles "$HINTFILE_GFF_1","$HINTFILE_GFF_2" -k "$RNA_GTF" --cfg "$CONF_FILE" --out "$OUTPUT_DIR/TSEBRA_${PREFIX}.gtf" -v 2

echo "-------------------------------------------------------------------"
echo "TSEBRA gene prediction completed on $(date '+%Y-%m-%d %H:%M:%S')"
echo "Results stored in: $OUTPUT_DIR"
echo "-------------------------------------------------------------------"
```

Formatting the GTF produced by TSEBRA, converting GTF files to FASTA, and assessing the proteome using BUSCO.

```bash
TSEBRA_GTF=tsebra.output.gtf
TSEBRA_FORMATTED_GTF=tsebra.format.gtf
QUERY=genome.sc26.softmasked.fasta
OUTPUT_DIR=Output_directory
BUSCO_DIR="$OUTPUT_DIR/BUSCO"			
BUSCO_DB=metazoa_odb10
BUSCO_OUTPUT_NAME="$(basename $OUTPUT_DIR)_BUSCO"

PREFIX="Pristina"

mkdir -p "$OUTPUT_DIR" "$BUSCO_DIR"

echo "-------------------------------------------------------------------"
echo "1. FORMATTING GTF PRODUCED BY TSEBRA"
echo "-------------------------------------------------------------------"

echo "Job starting on on $(date '+%Y-%m-%d %H:%M:%S')"
echo "Using gtf: $TSEBRA_GTF"
echo "Using genome: $QUERY"
echo "New gtf: $TSEBRA_FORMATTED_GTF"

awk 'BEGIN{OFS="\t"} 
    $3 == "gene" { print $1, $2, $3, $4, $5, $6, $7, $8, "gene_id \""$9"\";" }
    $3 == "transcript" { 
        split($9, a, /\./); 
        gene_id = a[1]; 
        print $1, $2, $3, $4, $5, $6, $7, $8, "gene_id \""gene_id"\"; transcript_id \""$9"\";" 
    }
    $3 != "gene" && $3 != "transcript" { print }
' "$TSEBRA_GTF" > "$TSEBRA_FORMATTED_GTF"

echo "-------------------------------------------------------------------"
echo "2. CONVERTING GTF TO FASTA"
echo "-------------------------------------------------------------------"

gffread "$TSEBRA_FORMATTED_GTF" -g "$GENOME" -y "$OUTPUT_DIR/TSEBRA_${PREFIX}_proteins.fasta"
gffread "$TSEBRA_FORMATTED_GTF" -g "$GENOME" -x "$OUTPUT_DIR/TSEBRA_${PREFIX}_cds.fasta"

echo "Proteins and cds saved in $OUTPUT_DIR"

echo "-------------------------------------------------------------------"
echo "3. BUSCO ANALYSIS"
echo "-------------------------------------------------------------------"

PROTEOME=TSEBRA_${PREFIX}_proteins.fasta

busco -i "$PROTEOME" -l "$BUSCO_DB" -o "$BUSCO_OUTPUT_NAME" -m prot -c 32

mv "$BUSCO_OUTPUT_NAME" "$BUSCO_DIR"

echo "-------------------------------------------------------------------"
echo "All analyses completed successfully on $(date '+%Y-%m-%d %H:%M:%S')"
echo "-------------------------------------------------------------------"
```

### 2.3. Proteome taxonomic assessment

DIAMOND BLASTp against the NCBI nr database was performed to further assess the taxonomic composition of the proteome.

```bash
QUERY=proteins.fasta
DIAMOND_PATH=diamond
DB_PATH=nr.dmnd
OUTPUT_DIR=Output_directory

PREFIX="Pristina"

echo "Starting blastp at $(date)"

"$DIAMOND_PATH" blastp\
 --query "$QUERY"\
 --db "$DB_PATH"\
 --evalue 1e-5\
 --max-target-seqs 1\
 --outfmt 6 qseqid sseqid bitscore pident length qlen slen evalue\
 --more-sensitive \
 --out "$OUTPUT_DIR/${PREFIX}_blastp.out" \
 --threads 32

echo "-------------------------------------------------------------------------------"
echo " blastp completed."
echo " See results in $OUTPUT_DIR"
echo " All analyses completed successfully on $(date '+%Y-%m-%d %H:%M:%S')"
echo "-------------------------------------------------------------------------------"
```

The protein BLAST results were mapped back to scaffolds by linking annotated proteins to their corresponding scaffolds and calculating query and subject coverage. The three best-scoring hits for each scaffold were then selected to assess its taxonomic identity.

```r 
library(dplyr)
library(readr)
library(stringr)

## ---- Read GTF (AUGUSTUS-style) ----
gtf <- read_tsv(
  gtf_file,
  col_names = c(
    "scaffold", "source", "type",
    "start", "end", "score",
    "strand", "phase", "gene_id"
  ),
  comment = "#"
)

## ---- Keep only gene-level info ----
genes <- gtf %>%
  filter(type == "transcript") %>%
  mutate(
    gene_id = str_extract(gene_id, '(?<=transcript_id ")[^"]+')
  ) %>%
  select(scaffold, gene_id)
  
## ---- Read BLASTP output ----
# Make sure the col order is correct
blast <- read_tsv(
  blast_file,
  col_names = c(
    "qseqid",  
    "sseqid",
    "bitscore",
    "pident",
    "length",
    "qlen",
    "slen",
    "evalue"
  )
)

## ---- Merge ----
genes_with_blast <- genes %>%
  left_join(blast, by = "qseqid")
genes_with_blast <- na.omit(genes_with_blast)

## ---- Merge and sort by scaffold + BLAST hit ----
genes_with_blast_sorted <- genes_with_blast %>%
  arrange(scaffold, sseqid)
  
## ---- Write output ----
write_csv(genes_with_blast_sorted, "out_file.csv")

## ---- Compute query and subject coverage ----
genes_with_blast <- genes_with_blast %>%
  mutate(
    query_cov = ifelse(!is.na(length) & !is.na(qlen), length / qlen, NA_real_),
    subject_cov = ifelse(!is.na(length) & !is.na(slen), length / slen, NA_real_),
    combined_score = bitscore * query_cov * subject_cov
  )
  
write_csv(genes_with_blast, "compute_query.csv")

## ---- Select top three hits per scaffold ----
top_n <- 3  

top_hits_per_scaffold <- genes_with_blast %>%
  filter(!is.na(combined_score)) %>%
  group_by(scaffold) %>%
  arrange(desc(combined_score)) %>%
  slice_head(n = top_n) %>%
  ungroup()

write_csv(top_hits_per_scaffold, "bestblast_out_file.csv")
```

Taxify was also run as described in Section 1.9 to obtain taxonomic assignments for the protein BLAST hits. The taxonomic assignments and the top three BLAST hits for each scaffold were then inspected manually.
Scaffold 24 was identified as potential contamination and was removed.

### 2.4. Scaffold removal and final genome assembly

Scaffold 24 was removed from the genome assembly.

```bash
seqkit grep -v -p "scaffold_24" genome.sc26.fasta > genome.sc25.fasta
```

The best BLAST hit for each scaffold was retained based on the combined score calculated in Section 2.3. Section 1.10 was then repeated to obtain the final genome assembly metrics.
The masked genomes were merged again as described in Section 2.1.3, with an additional step added to the script to remove scaffold 24. This step was inserted between steps 1 and 3 of the previous script. Section 2.1.4 was also repeated to generate the final Venn diagram.

```bash
# ==== STEP 2. REMOVE scaffold 24 from ALL BEDs ====
echo "[2/6] Filtering scaffold 24..."

for bed in "${NAME_RAW}_RepeatMod.bed" \
           "${NAME_RAW}_dfam.bed" \
           "${NAME_RAW}_edta.bed"
do
    awk '$1!="scaffold_24"' "$bed" > tmp.bed
    mv tmp.bed "$bed"
done
```

To remove the transcripts corresponding to scaffold 24 from the proteome, the transcript IDs associated with scaffold 24 were first extracted and saved to a TXT file. A map linking scaffolds to transcripts was also created.

```bash
awk '$3=="transcript" {
    match($0, /transcript_id "[^"]+"/)
    transcript=substr($0, RSTART+15, RLENGTH-16)
    print $1 "\t" transcript
}' tsebra.format.gtf > scaffold24_transcripts_map.tsv
```

```bash
awk '$1=="scaffold_24" {print $2}' scaffold24_transcripts_map.tsv > scaffold24_transcripts.txt
```

Scaffold 24 was removed from the GTF produced by TSEBRA in Section 2.2.3.

```bash
awk '$1 != "scaffold_24"' tsebra.format.gtf > tsebra.format.sc25.gtf
```

Transcripts corresponding to scaffold 24 were removed from the proteome and CDS files.

```bash
seqkit grep -v -f scaffold24_transcripts.txt proteins.fasta > proteins.sc25.fasta
```

```bash
seqkit grep -v -f scaffold24_transcripts.txt cds.fasta > cds.sc25.fasta
```

### 2.5. Reclassifying repeats 

As each masking tool classifies repeats using different names, the repeats in the final merged masked genome were reclassified into major categories to characterise the repeat content of the genome.
Each interval in the merged mask was assigned a repeat class based on the overlapping annotations, with the following priority: EDTA > RepeatModeler > Dfam. The script produces the final repeat BED/GFF3 files and summaries by repeat class and repeat group.

First, scaffold 24 was removed from the input files. 

EDTA GFF3

```bash
awk '$1 != "scaffold_24"' EDTA.TEanno.gff3 > sc25.EDTA.TEanno.gff3
```

RepeatMasker Dfam output

```bash
awk '$5 != "scaffold_24"' RepeatMasker_Dfam.out > sc25_RepeatMasker_Dfam.out
```

RepeatMasker RepeatModeler output

```bash
awk '$5 != "scaffold_24"' RepeatMasker_RepeatMod.out > sc25_RepeatMasker_RepeatMod.out
```

Final repeat annotation and summary

```bash
MERGED_BED=merged_masked_genome.bed
GENOME_FA=genome.sc25.fasta
EDTA_GFF=sc25.EDTA.TEanno.gff3
RM_REPEATMOD=sc25_RepeatMasker_RepeatMod.out
RM_DFAM=sc25_RepeatMasker_Dfam.out
OUTPUT_DIR=Output_directory

mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

# -------------------------------------------------------------------
# 1. Convert EDTA GFF3 to BED
# -------------------------------------------------------------------

echo "Converting EDTA annotation to BED file..."

# Convert GFF3 1-based coordinates to BED 0-based coordinates
awk 'BEGIN{OFS="\t"} !/^#/ {name="NA"; class="NA"; split($9,a,";"); for(i in a){if(a[i] ~ /^Name=/){name=a[i]; sub("Name=","",name)} if(a[i] ~ /^classification=/){class=a[i]; sub("classification=","",class)}}
print $1,$4-1,$5,name,class}' "$EDTA_GFF" > EDTA.bed

# -------------------------------------------------------------------
# 2. Convert RepeatMasker outputs to BED
# -------------------------------------------------------------------

echo "Converting RepeatMasker annotations to BED file..."

# Extract sequence, start, end, repeat name and repeat class
awk 'BEGIN{OFS="\t"} $1 !~ /^SW|^score|^perc|^query|^position/ && NF>10 {print $5,$6-1,$7,$10,$11}' "$RM_REPEATMOD" > RepeatModeler.bed
awk 'BEGIN{OFS="\t"} $1 !~ /^SW|^score|^perc|^query|^position/ && NF>10 {print $5,$6-1,$7,$10,$11}' "$RM_DFAM" > Dfam.bed

# -------------------------------------------------------------------
# 3. Add ID to final merged mask
# -------------------------------------------------------------------

echo "Preparing merged mask..."

# Add an ID to each interval so overlapping annotations can be linked back to the corresponding final masked interval.
awk 'BEGIN{OFS="\t"} {print $0,"mask_"NR}' "$MERGED_BED" > merged_mask_ID.bed

# -------------------------------------------------------------------
# 4. Find overlaps
# -------------------------------------------------------------------

echo "Finding overlaps..."

bedtools intersect \
    -a merged_mask_ID.bed \
    -b EDTA.bed \
    -wa -wb \
    > merged_EDTA.tsv

bedtools intersect \
    -a merged_mask_ID.bed \
    -b RepeatModeler.bed \
    -wa -wb \
    > merged_RepeatModeler.tsv

bedtools intersect \
    -a merged_mask_ID.bed \
    -b Dfam.bed \
    -wa -wb \
    > merged_Dfam.tsv
	
# -------------------------------------------------------------------
# 5. Assign repeat classification with priority
# -------------------------------------------------------------------

echo "Assigning repeat classifications..."

# Classification priority: EDTA > RepeatModeler > Dfam.
# If an interval overlaps annotations from multiple sources, the first available classification in this order is assigned.

awk 'BEGIN{OFS="\t"}
FILENAME=="merged_EDTA.tsv" {id=$5; class[id]=$10; next}
FILENAME=="merged_RepeatModeler.tsv" {id=$5; if(!(id in class)) class[id]=$10; next}
FILENAME=="merged_Dfam.tsv" {id=$5; if(!(id in class)) class[id]=$10; next}
FILENAME=="merged_mask_ID.bed" {id=$5; if(id in class) print $1,$2,$3,class[id]; else print $1,$2,$3,"Unclassified"}
' merged_EDTA.tsv merged_RepeatModeler.tsv merged_Dfam.tsv merged_mask_ID.bed | sort -k1,1 -k2,2n > final_mask_annotated.bed

# -------------------------------------------------------------------
# 6. Create GFF3
# -------------------------------------------------------------------

echo "Creating GFF3..."

awk 'BEGIN{OFS="\t"; print "##gff-version 3"} \
{len=$3-$2; print $1,"FinalMask","repeat_region",$2+1,$3,".",".",".","ID=repeat_"NR";classification="$4";length="len}' \
final_mask_annotated.bed > final_mask_annotated.gff3

# -------------------------------------------------------------------
# 7. Calculate genome size
# -------------------------------------------------------------------

echo "Calculating genome size..."

GENOME_SIZE=$(awk '/^>/{if(length(seq)>0) total+=length(seq); seq=""; next}
{seq=seq $0} END{total+=length(seq); print total}' "$GENOME_FA")

echo "Genome size: $GENOME_SIZE bp"

# -------------------------------------------------------------------
# 8. Calculate repeat class summaries
# -------------------------------------------------------------------

echo "Creating detailed repeat summaries..."

# Calculate BP per class from the final annotated BED
awk 'BEGIN{OFS="\t"} {bp[$4]+=$3-$2} END{print "Repeat_class\tMasked_bp"; for(i in bp) print i,bp[i]}' final_mask_annotated.bed \
| sort -k2,2nr > final_repeat_class_summary_bp.txt

# Calculate percentage of genome for each class
awk -v genome="$GENOME_SIZE" 'BEGIN{OFS="\t"} NR==1{print $0,"Percent_of_genome"; next}
{printf "%s\t%.0f\t%.4f\n",$1,$2,($2/genome)*100}' final_repeat_class_summary_bp.txt \
> final_repeat_class_summary_genome_percent.txt

# -------------------------------------------------------------------
# 9. Calculate total masked statistics
# -------------------------------------------------------------------

echo "Calculating total masked statistics..."

# Calculate total masked BP from the annotated BED
TOTAL_MASKED_BP=$(awk '{sum+=$3-$2} END{print sum}' final_mask_annotated.bed)
TOTAL_MASKED_PERCENT=$(awk -v masked="$TOTAL_MASKED_BP" -v genome="$GENOME_SIZE" 'BEGIN{
printf "%.4f", (masked/genome)*100}')

# Verify sum of classes
SUM_OF_CLASSES=$(awk '{sum+=$2} END{print sum}' final_repeat_class_summary_bp.txt)

echo "Total masked BP from annotated BED: $TOTAL_MASKED_BP"
echo "Sum of individual class BPs: $SUM_OF_CLASSES"

if [ "$TOTAL_MASKED_BP" != "$SUM_OF_CLASSES" ]; then
    DIFF=$((TOTAL_MASKED_BP - SUM_OF_CLASSES))
    echo "WARNING: Difference of $DIFF bp between total and sum of classes!"
fi

# Write total repeat content
echo -e "Total_masked_bp\t$TOTAL_MASKED_BP" > total_repeat_content.txt
echo -e "Percent_of_genome\t$TOTAL_MASKED_PERCENT" >> total_repeat_content.txt
echo -e "Genome_size_bp\t$GENOME_SIZE" >> total_repeat_content.txt
echo -e "Sum_of_class_bp\t$SUM_OF_CLASSES" >> total_repeat_content.txt
echo -e "Class_sum_percent\t$(awk -v sum="$SUM_OF_CLASSES" -v genome="$GENOME_SIZE" 'BEGIN{printf "%.4f", (sum/genome)*100}')" >> total_repeat_content.txt

echo "Genome_size_bp	$GENOME_SIZE" > genome_size.txt

# -------------------------------------------------------------------
# 10. Collapsed repeat landscape
# -------------------------------------------------------------------

echo "Creating collapsed repeat summary..."

# Extract all unique classes for mapping
tail -n +2 final_repeat_class_summary_bp.txt | cut -f1 | sort -u > all_classes.txt

# Create a mapping of class -> group
awk 'BEGIN{OFS="\t"; print "Class\tGroup"}
{
    class = $1
    
    # Pattern matching
    if (class ~ /LTR/ || class ~ /^LTR/) {
        group = "LTR_retrotransposons"
    } else if (class ~ /LINE/ || class ~ /^LINE/) {
        group = "LINE_retrotransposons"
    } else if (class ~ /SINE/ || class ~ /^SINE/) {
        group = "SINE_retrotransposons"
    } else if (class ~ /PLE/ || class == "Retroposon") {
        group = "PLE_retrotransposons"
    } else if (class ~ /Helitron/ || class ~ /^RC\// || class ~ /^DNA\/Helitron/) {
        group = "Helitrons"
    } else if (class ~ /DNA/ || class ~ /MITE/) {
        group = "DNA_transposons"
    } else if (class == "Simple_repeat" || class ~ /Simple_repeat/) {
        group = "Simple_repeats"
    } else if (class == "Low_complexity" || class ~ /Low_complexity/) {
        group = "Low_complexity"
    } else if (class ~ /Satellite/) {
        group = "Satellite"
    } else if (class == "rRNA" || class == "tRNA" || class == "snRNA" || class == "srpRNA" || class ~ /RNA/) {
        group = "RNA_genes"
    } else if (class == "Unknown" || class == "unknown" || class == "Unclassified") {
        group = "Unclassified"
    } else if (class == "ARTEFACT") {
        group = "Artifact"
    } else {
        group = "Other"
    }
    
    print class, group
}
' all_classes.txt > class_mapping_debug.txt

# Now create the collapsed summary using the mapping
awk -v genome="$GENOME_SIZE" '
BEGIN{OFS="\t"; while ((getline < "class_mapping_debug.txt") > 0) {if ($1 != "Class") class_to_group[$1]=$2}
close("class_mapping_debug.txt")}
{
    class = $1
    bp = $2
    
    # Skip header
    if (class == "Repeat_class") next
    }
    
    # Get the group for this class
    if (class in class_to_group) {
        group = class_to_group[class]
    } else {
        # Fallback: try to determine group on the fly
        if (class ~ /DNA/ || class ~ /MITE/) {
            group = "DNA_transposons"
        } else if (class ~ /LTR/) {
            group = "LTR_retrotransposons"
        } else if (class ~ /LINE/) {
            group = "LINE_retrotransposons"
        } else if (class ~ /SINE/) {
            group = "SINE_retrotransposons"
        } else if (class ~ /Unknown/ || class ~ /unknown/) {
            group = "Unclassified"
        } else {
            group = "Other"
        }
    }
    
    # Sum BP by group
    group_bp[group] += bp
}
END {
    print "Repeat_group\tMasked_bp\tPercent_of_genome"
    total = 0
    for (g in group_bp) {
        total += group_bp[g]
        printf "%s\t%.0f\t%.4f\n", g, group_bp[g], (group_bp[g]/genome)*100
    }
    printf "Total_masked\t%.0f\t%.4f\n", total, (total/genome)*100
}
' final_repeat_class_summary_bp.txt | sort -k2,2nr > collapsed_repeat_summary.txt

# -------------------------------------------------------------------
# 11. Final statistics and verification
# -------------------------------------------------------------------

echo ""
echo "-------------------------------------"
echo "FINAL VERIFICATION"
echo "-------------------------------------"

# Extract values
TOTAL_FROM_TOTAL=$(grep "Total_masked_bp" total_repeat_content.txt | cut -f2)
TOTAL_FROM_COLLAPSED=$(grep "^Total_masked" collapsed_repeat_summary.txt | cut -f2)

echo "Total masked BP comparison:"
echo "  From total_repeat_content.txt: $TOTAL_FROM_TOTAL"
echo "  From collapsed_repeat_summary.txt: $TOTAL_FROM_COLLAPSED"
echo "  From final_mask_annotated.bed: $TOTAL_MASKED_BP"

if [ "$TOTAL_FROM_TOTAL" == "$TOTAL_FROM_COLLAPSED" ]; then
    echo "✓ All totals match!"
else
    DIFF=$((TOTAL_FROM_TOTAL - TOTAL_FROM_COLLAPSED))
    echo "✗ WARNING: Totals do not match! Difference: $DIFF bp"
    echo ""
    echo "Check class_mapping_debug.txt to see how each class was categorized"
    echo ""
    echo "Classes in 'Other' category (these need to be categorized):"
    if [ -f class_mapping_debug.txt ]; then
        grep -E "Other$" class_mapping_debug.txt || echo "  None found"
    fi
    echo ""
    echo "Classes in 'Artifact' category (PCR artifacts):"
    if [ -f class_mapping_debug.txt ]; then
        grep -E "Artifact$" class_mapping_debug.txt || echo "  None found"
    fi
fi

echo ""
echo "-------------------------------------"
echo "Repeat Content Summary"
echo "-------------------------------------"
echo "Genome size: $(printf "%'d" $GENOME_SIZE) bp"
echo "Total masked: $(printf "%'d" $TOTAL_MASKED_BP) bp ($TOTAL_MASKED_PERCENT%)"
echo "Number of repeat regions: $(wc -l < final_mask_annotated.bed)"
echo ""

# Show collapsed summary
echo "Collapsed repeat summary:"
cat collapsed_repeat_summary.txt

echo ""
echo "-------------------------------------"
echo "Output files:"
echo "-------------------------------------"
ls -lh *.txt *.bed *.gff3 2>/dev/null | awk '{print "  " $9}'
echo ""

echo "For more details, check:"
echo "  - all_classes.txt: All unique repeat classes"
echo "  - class_mapping_debug.txt: How each class was categorized"
echo "  - final_repeat_class_summary_bp.txt: Detailed class breakdown"
echo "  - final_repeat_class_summary_genome_percent.txt: Class percentages"
```

### 2.6. Telomeric and centromeric regions

The quarTeT pipeline was used to identify telomeric and centromeric regions.
TE annotations generated using EDTA were provided as input for CentroMiner, as recommended by the quarTeT workflow.

```bash
QUERY=genome.sc25.fasta
TE_ANNOTATION=EDTA.TEanno.gff3
SPECIES_TYPE="animal"
OUTPUT_DIR=Output_directory

PREFIX="Pristina"

mkdir "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

echo "[Step 1/2] Running TeloExplorer for telomere identification..."

quartet TeloExplorer \
    -i "$QUERY" \
    -c "$SPECIES_TYPE" \
    -p "${PREFIX}_telo"

echo "TeloExplorer complete."

echo "[Step 2/2] Running CentroMiner for centromere prediction..."

quartet CentroMiner \
        -i "$QUERY" \
        --TE "$TE_ANNOTATION" \
        -p "${PREFIX}_centro"
		
echo "CentroMiner complete."
```

### 2.7. Functional annotation

InterProScan

```bash
QUERY=proteins.sc25.fasta
OUTPUT_DIR=Output_directory
PREFIX="Pristina"

mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR" || exit

echo "Starting Interproscan at $(date)"
interproscan.sh \
--input "$QUERY" \
--seqtype p \
--goterms \
--pathways \
--iprlookup \
--formats tsv \
--outfile "$OUTPUT_DIR/${PREFIX}_interpro.tsv" 
--cpu 32
```

eggNOG-mapper

```bash
QUERY=proteins.sc25.fasta
EGGNOG_DB=eggnog_db
EGGNOG_DIAMOND_DB=eggnog_proteins.dmnd
OUTPUT_DIR=Output_directory
PREFIX="Pristina"

mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR" || exit

echo "Starting functional annotation using eggNOG and Diamond at $(date)"
emapper.py \
  -i "$QUERY" \
  --itype proteins \
  -m diamond \
  --cpu 32 \
  --data_dir "$EGGNOG_DB" \
  --dmnd_db "$EGGNOG_DIAMOND_DB" \
  --pfam_realign realign \
  --go_evidence all \
  --target_orthologs all \
  --dbmem \
  --report_orthologs \
  --decorate_gff yes \
  --output_dir "$OUTPUT_DIR" \
  --output "$PREFIX"_eggnog_diamond \
  --override
  ```

## 3. Gene family evolution in phylogenetic context

### 3.1. Phylogenomic analysis

#### 3.1.1. Orthofinder

Prior to Orthofinder, protein IDs were made uniquely identifiable by species for downstream comparative analyses. The below script adds the .faa filename (species name) as a prefix to each protein FASTA header.

```bash
input_dir=Protein_Data
set -euo pipefail

files=("$input_dir"/*.faa)

for file in "${files[@]}"; do
    species=$(basename "$file")
    species="${species%.*}"
    echo "Processing: $file -> ${species}_scaffold"
    tmp="${file}.tmp"
    awk -v sp="$species" '
        /^>/ {print ">" sp "_" substr($0,2); next}
        {print}
    ' "$file" > "$tmp"
    mv "$tmp" "$file"
done

echo "Done."
```

```bash
PROT_DIR=Protein_Data
OUTPUT_DIR=Output_directory
PREFIX="Pristina"

mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR" || exit

echo "Starting Orthofinder on $(date '+%Y-%m-%d %H:%M:%S')"
orthofinder -f "$PROT_DIR" -t 32 -a 32 -o "$OUTPUT_DIR/${PREFIX}orthofinder_out_$(date +%Y%m%d_%H%M%S)"

echo "----------------------------------------------------------------------"
echo "Orthofinder completed. See results in $OUTPUT_DIR"
echo "All analyses completed successfully on $(date '+%Y-%m-%d %H:%M:%S')"
echo "----------------------------------------------------------------------"
```

Prior to downstream analysis, a mapping file was created linking individual gene sequences to their respective species for species-tree inference with ASTRAL.

```bash
INPUT_DIR=Protein_Data
OUTPUT="mapping.txt"

awk '/^>/ {
    h=$1; sub(/^>/,"",h)
    if (match(h,/^(.*)_[^_]+$/,m)) print h "\t" m[1]
}' "$INPUT_DIR"/*.faa > "$OUTPUT"
```

#### 3.1.2. IQTREE and ASTRAL
A phylogenetic tree was constructed using IQ-TREE and ASTRAL

```bash
ORTHOFINDER_RESULTS=Orthofinder_results_directory
SINGLE_CP_ORTH_TXT=Orthogroups_SingleCopyOrthologues.txt
SINGLE_CP_ORTH_FASTA_DIR="$ORTHOFINDER_RESULTS/Orthogroup_Sequences"
THREADS=32
OUTGROUP_NAME="Pecten_maximus"
WORKDIR=Working_directory
PREFIX="Pristina"

mkdir -p "$WORKDIR_DIR"/{Single_Copy_Orthologues,Aligned_SCOs,Trimmed_SCOs,IQTREE/tmp_genes,IQTREE/gene_trees,IQTREE/model_summary}

# ---------------------------------------------
# STEP 1: COPY SINGLE COPY ORTHOLOGS
# ---------------------------------------------
echo "Step 1: Copy SCOs"

while read og; do
    cp "$SINGLE_CP_ORTH_FASTA_DIR/${og}.fa" "$WORKDIR/Single_Copy_Orthologues/"
done < "$SINGLE_CP_ORTH_TXT"

# ---------------------------------------------
# STEP 2: ALIGNMENT (MAFFT)
# ---------------------------------------------

echo "Step 2: MAFFT alignment"

for f in "$WORKDIR"/Single_Copy_Orthologues/*.fa; do
    base=$(basename "$f" .fa)
    mafft --auto --thread "$THREADS" "$f" > "$WORKDIR/Aligned_SCOs/${base}.aln.fa"
done

# ---------------------------------------------
# STEP 3: TRIMMING (CLIPKIT)
# ---------------------------------------------

echo "Step 3: ClipKIT trimming"

for aln in "$WORKDIR"/Aligned_SCOs/*.aln.fa; do
    base=$(basename "$aln" .aln.fa)
    clipkit "$aln" -m smart-gap -o "$WORKDIR/Trimmed_SCOs/${base}.trim.fa"
done

# ---------------------------------------------
# STEP 4: CONCATENATION (AMAS)
# ---------------------------------------------

echo "Step 4: AMAS concatenation"

AMAS.py concat -i "$WORKDIR"/Trimmed_SCOs/*.trim.fa \
               -f fasta -d aa \
               -p "$WORKDIR/${PREFIX}_partitions.txt" \
               -t "$WORKDIR/${PREFIX}_concatenated_alignment.fa"

# ---------------------------------------------
# STEP 5: IQ-TREE PREP
# ---------------------------------------------

echo "Step 5: Prepare IQ-TREE commands"

TRIM_DIR="$WORKDIR/Trimmed_SCOs"
TMP_DIR="$WORKDIR/IQTREE/tmp_genes"
TREE_DIR="$WORKDIR/IQTREE/gene_trees"
MODEL_DIR="$WORKDIR/IQTREE/model_summary"

CMD_FILE="$WORKDIR/IQTREE/iqtree_commands.txt"
MODEL_OUT="$MODEL_DIR/all_models.txt"

rm -f "$CMD_FILE" "$MODEL_OUT"

for f in "$TRIM_DIR"/*.trim.fa; do
    base=$(basename "$f" .trim.fa)

    TOTAL=$(grep -v "^>" "$f" | tr -d '\n' | wc -c)
    if (( TOTAL == 0 )); then
        echo "Skipping $base (empty alignment)"
        continue
    fi

    TMP_FASTA="$TMP_DIR/${base}.fasta"
    cp "$f" "$TMP_FASTA"

    echo "iqtree -s $TMP_FASTA -m MFP -nt 1 --safe -pre $TMP_DIR/${base} > /dev/null && \
cp $TMP_DIR/${base}.treefile $TREE_DIR/${base}.treefile && \
model=\$(grep 'Best-fit model' $TMP_DIR/${base}.iqtree | awk -F': ' '{print \$2}') && \
flock $MODEL_OUT echo -e '${base}\t'\$model >> $MODEL_OUT && \
rm -f $TMP_FASTA $TMP_DIR/${base}.*" >> "$CMD_FILE"
done

# ---------------------------------------------
# STEP 6: RUN IQ-TREE PARALLEL
# ---------------------------------------------

echo "Step 6: IQ-TREE parallel run"

parallel -j 32 < "$CMD_FILE"

# ---------------------------------------------
# STEP 7: CONCAT TREESET
# ---------------------------------------------

echo "Step 7: Concatenate gene trees"

cat "$TREE_DIR"/*.treefile > "$TREE_DIR/all_gene_trees.tre"

# ---------------------------------------------
# STEP 8: SPECIES TREE
# ---------------------------------------------

echo "Step 8: ASTRAL species tree"

astral4 -a "$WORKDIR/mapping.txt" -o "$WORKDIR/${PREFIX}.tre" "$TREE_DIR/all_gene_trees.tre"

# ---------------------------------------------
# STEP 9: REROOT SPECIES TREE
# ---------------------------------------------

echo "Step 9: Reroot species tree"

nw_reroot "$WORKDIR/${PREFIX}.tre" "$OUTGROUP_NAME" > "$WORKDIR/${PREFIX}_rooted.tre"

echo "---------------------------------------------"
echo "PIPELINE COMPLETE: $(date)"
echo "Output: $WORKDIR"
echo "---------------------------------------------"
```

#### 3.1.3. Phylogenetic tree

```r 

library(ape)
library(ggtree)
library(ggplot2)

tree <- read.tree("rooted.tre")

ggtree(tree) +
  geom_tiplab(size = 3) +
    geom_text2(aes(label = ifelse(is.na(label),
                                "",
                                paste0(round(as.numeric(label) * 100, 0), "%"))),
             size = 3,
             hjust = -0.02) +
  
  theme_tree2(base_size = 10) +
  theme(
    axis.text.x = element_text(size = 8),
    axis.title.x = element_blank(),
  ) +
  scale_y_reverse()
```

### 3.2. Gene Family Evolution

#### 3.2.1 Ultrametric tree

The phylogenetic tree was converted into an ultrametric tree. Divergence time estimates were obtained from the TimeTree database 

```r 
library(ape)
library(phytools)

my_tree <- read.tree("rooted.tre")

# ------------------------------------------------------------------
# FIND THE NODE NUMBERS FOR CALIBRATIONS
# ------------------------------------------------------------------

# Get the node numbers for the Most Recent Common Ancestor (MRCA) of the calibration species pairs

# Calibration 1: Root node (Mollusca vs Annelida)
# The common ancestor of Pecten maximus and Lumbricus terrestris
root_node <- getMRCA(my_tree, tip = c("Pecten_maximus", "Lumbricus_terrestris"))
print(paste("Root node number:", root_node))

# Calibration 2: Clitellata node (Earthworm vs Leech)
# The common ancestor of Lumbricus terrestris and Hirudinaria manillensis
clitellata_node <- getMRCA(my_tree, tip = c("Lumbricus_terrestris", "Hirudinaria_manillensis"))
print(paste("Clitellata node number:", clitellata_node))

# ------------------------------------------------------------------
# CREATE THE CALIBRATION DATA FRAME
# ------------------------------------------------------------------

# Create vectors for nodes and ages
node <- c(root_node, clitellata_node)

# From TimeTree:
# Calibration 1 (Root): Pecten maximus vs Lumbricus terrestris
#   CI: 447.8 - 650.0 MYA
# Calibration 2 (Clitellata): Lumbricus terrestris vs Hirudinaria manillensis
#   Range: 295.9 - 345.0 MYA
age.min <- c(447.8, 295.9)
age.max <- c(650.0, 345.0)

soft.bounds <- c(FALSE, FALSE)

calibration_df <- data.frame(
    node = node,
    age.min = age.min,
    age.max = age.max,
    soft.bounds = soft.bounds
)

print("Calibration data frame:")
print(calibration_df)

# ------------------------------------------------------------------
# VISUALIZE THE CALIBRATIONS ON YOUR TREE 
# ------------------------------------------------------------------

plot(my_tree, show.tip.label = TRUE, cex = 0.8)

# Add calibration nodes
nodelabels(node = root_node, pch = 19, col = "red", cex = 1.5)
nodelabels(node = clitellata_node, pch = 19, col = "blue", cex = 1.5)

# ------------------------------------------------------------------
# RUN THE DATING ANALYSIS WITH chronos
# ------------------------------------------------------------------

dated_tree <- chronos(
    my_tree,
    calibration = calibration_df,
    model = "correlated",          # Recommended for most datasets
    quiet = FALSE,                 
    control = chronos.control()    # Default settings
)

# ------------------------------------------------------------------
# CHECK THE RESULTS
# ------------------------------------------------------------------

# Check if the tree is now ultrametric
is_ultrametric <- is.ultrametric(dated_tree)
print(paste("Is the tree ultrametric?", is_ultrametric))

# View the new branch lengths
print(dated_tree$edge.length)

# Check root-to-tip distances (should all be equal)
dist_root_to_tip <- diag(vcv(dated_tree))
print("Root-to-tip distances (should all be roughly equal):")
print(dist_root_to_tip)

# ------------------------------------------------------------------
# SAVE ULTRAMETRIC TREE
# ------------------------------------------------------------------

# Save as Newick file for CAFE5
write.tree(dated_tree, file = "ultrametric.tre")
```

#### 3.2.2 CAFE5 Analysis

Prepare data for input into CAFE5. The script removes the Total column from the OrthoFinder Orthogroups.GeneCount.tsv file, filters out families with more than 100 copies in any species, adds the required Desc column, and generates summary logs.

```bash
awk 'BEGIN {OFS = "\t"; kept = 0; removed = 0; total = 0}

NR == 1 {
    # Remove the Total column from header (last column)
    header = ""
    for (i = 1; i < NF; i++) {
        if (i > 1) header = header OFS
        header = header $i
    }
    # Print header with Desc column
    print "Desc" OFS header > "Orthogroups_filtered.tsv"
    print "Desc" OFS header > "Orthogroups_removed.tsv"
    next
}
{
    total++
    
    # Remove the Total column from data (last column)
    data = ""
    large = 0
    for (i = 1; i < NF; i++) {
        if (i > 1) {
            data = data OFS
            # Check if this value is >=100 (skip first column = family ID)
            if ($i >= 100) {
                large = 1
            }
        }
        data = data $i
    }
    
    if (large == 0) {
        # Keep this family
        kept++
        print "A" kept OFS data > "Orthogroups_filtered.tsv"
    } else {
        # Remove this family
        removed++
        print "A" total OFS data > "Orthogroups_removed.tsv"
    }
}
END {
    print "Total families:", total
    print "Families kept:", kept
    print "Families removed:", removed
}
' Orthogroups.GeneCount.tsv
```

CAFE5

```bash
WORKDIR=Working_directory
ORTH_TSV=Orthogroups_filtered.tsv
TREE=ultrametric.tre
ERROR_MODEL="$WORKDIR/Base_error_model.txt"
FIGURE_DIR="$WORKDIR/cafe5_plots"

mkdir -p "$WORKDIR"/{results,logs}

# Estimate the CAFE5 error model

cafe5 -i "$ORTH_TSV" \
      -t "$TREE" \
      -o "$WORKDIR" \
      -p \
      -e \
      > "$WORKDIR/logs/error_estimation.log" 2>&1

echo "Error model created"

# Run CAFE5 with k=1 using the estimated error model

cafe5 -i "$ORTH_TSV" \
      -t "$TREE" \
      -o "$WORKDIR/results" \
      -p \
      -e"$ERROR_MODEL" \
      -k 1 \
      > "$WORKDIR/logs/k1.log" 2>&1
	  
# Run CAFEPLOTTER to visualise the CAFE5 Results

cafeplotter -i "$WORKDIR/results" \
            -o "$FIGURE_DIR" \
            --format png \
            --fig_width 12 \
            --fig_height 0.7 \
            --leaf_label_size 10 \
            --count_label_size 7 \
            --dpi 600 \
            --expansion_color red \
            --contraction_color blue \
            --ignore_branch_length
```

#### 3.2.3 Gene Family Expansions and Contractions

Extract significant gene family sequences. Work in the same working directory used for the CAFE5 analysis. Orthogroup sequences are located in the OrthoFinder output described in Section 3.1.1.

Significant expansions: Extracts significantly expanded gene families (p < 0.05) for the target species, as well as families expanded specifically in the target species.

The same script can be used to extract significant contractions by making the following changes:

In Step 2 – Change the target-species filter from `$col > 0` to `$col < 0`.
In Step 3 – Change the target-species filter from $col > 0 to $col < 0, and change the filter for other species from $i > 0 to $i < 0.
Change significant_expansions to significant_contractions in the output directory name.
Update the corresponding expansions/contractions wording in the script messages and summary.

All other steps, including significance filtering and sequence extraction, remain unchanged.

```bash
SPECIES="Pristina_leidyi"

WORKDIR=Working_directory
CHANGE_FILE="$WORKDIR/results/Base_change.tab"
REPORT_FILE="$WORKDIR/results/Base_report.cafe"
ORTHOGROUP_DIR=Orthogroup_Sequences

PVAL_THRESHOLD="0.05"
OUTPUT_DIR="$WORKDIR/${SPECIES}_significant_expansions"
LOG_FILE="$OUTPUT_DIR/significant_extraction.log"

mkdir -p "$OUTPUT_DIR"/{all_significant,specific}

exec > >(tee -a "$LOG_FILE") 2>&1

# ------------------------------------------------------------------
# STEP 1: FIND SPECIES COLUMN
# ------------------------------------------------------------------
echo "Step 1: Finding column for $SPECIES..."

COLUMN_NUM=$(head -1 "$CHANGE_FILE" | awk -v species="$SPECIES" '
{
    for(i=1; i<=NF; i++) {
        if(index($i, species) == 1) {
            print i
            exit
        }
    }
}')

echo "  Found $SPECIES in column $COLUMN_NUM"
echo ""

# ------------------------------------------------------------------
# STEP 2: EXTRACT SIGNIFICANT EXPANSIONS IN TARGET SPECIES
# ------------------------------------------------------------------
echo "Step 2: Extracting significant expansions in $SPECIES (p < $PVAL_THRESHOLD)..."

# First, get all significant families from report
grep "^OG" "$REPORT_FILE" | awk -v pval="$PVAL_THRESHOLD" '
{
    for(i=2; i<=NF; i++) {
        if($i ~ /^[0-9.]+$/ && $i < pval && $i >= 0) {
            print $1
            break
        }
    }
}' > "$OUTPUT_DIR/all_significant_ids.txt"

TOTAL_SIG=$(wc -l < "$OUTPUT_DIR/all_significant_ids.txt")
echo "  Found $TOTAL_SIG total significant families"

# Filter for expansions in target species (column > 0)
# For contractions, change "$col > 0" below to "$col < 0"

awk -v col="$COLUMN_NUM" '
    NR==FNR {sig[$1]=1; next} 
    $1 in sig && $col > 0 {
        print $1
    }
' "$OUTPUT_DIR/all_significant_ids.txt" "$CHANGE_FILE" > "$OUTPUT_DIR/all_significant_families.txt"

SIG_COUNT=$(wc -l < "$OUTPUT_DIR/all_significant_families.txt")
echo "  Found $SIG_COUNT significantly EXPANDED families in $SPECIES"

# Copy to all_significant directory
cp "$OUTPUT_DIR/all_significant_families.txt" "$OUTPUT_DIR/all_significant/families.txt"

if [ "$SIG_COUNT" -eq 0 ]; then
    echo "  ERROR: No significant expansions found in $SPECIES."
    exit 1
fi

# ------------------------------------------------------------------
# STEP 3: FIND SPECIES-SPECIFIC EXPANSIONS
# ------------------------------------------------------------------
echo ""
echo "Step 3: Finding $SPECIES-specific significant expansions..."

# For contractions:
#   change "$col > 0" to "$col < 0"
#   change "$i > 0" to "$i < 0"

awk -v col="$COLUMN_NUM" '
    NR==FNR {sig[$1]=1; next} 
    $1 in sig && $col > 0 {
        for(i=2; i<=NF; i++) {
            if(i != col && $i > 0) { next }
        }
        print $1
    }
' "$OUTPUT_DIR/all_significant_families.txt" "$CHANGE_FILE" > "$OUTPUT_DIR/specific_families.txt"

SPECIFIC_COUNT=$(wc -l < "$OUTPUT_DIR/specific_families.txt")
echo "  Found $SPECIFIC_COUNT $SPECIES-specific significant expansions"

# Copy to specific directory
cp "$OUTPUT_DIR/specific_families.txt" "$OUTPUT_DIR/specific/families.txt"

# ------------------------------------------------------------------
# STEP 4: EXTRACT SEQUENCES
# ------------------------------------------------------------------
echo ""
echo "Step 4: Extracting sequences..."

extract_seqs() {
    local list="$1" out="$2"
    [ ! -s "$list" ] && { echo "0"; return; }
    
    > "$out"
    while read family; do
        [ -f "$ORTHOGROUP_DIR/${family}.fa" ] || continue
        awk -v og="$family" -v species="$SPECIES" '
            $0 ~ "^>" species {print $0 "|OG:" og; flag=1; next}
            /^>/ && flag {flag=0}
            flag {print}
        ' "$ORTHOGROUP_DIR/${family}.fa" >> "$out"
    done < "$list"
    
    echo "$(grep -c "^>" "$out" 2>/dev/null || echo 0)"
}

# Extract all significant expansions
ALL_SEQ=$(extract_seqs "$OUTPUT_DIR/all_significant_families.txt" "$OUTPUT_DIR/all_significant/sequences.fa")

# Extract species-specific expansions
SPEC_SEQ=$(extract_seqs "$OUTPUT_DIR/specific_families.txt" "$OUTPUT_DIR/specific/sequences.fa")

# ------------------------------------------------------------------
# SUMMARY
# ------------------------------------------------------------------
echo ""
echo "=========================================="
echo "SUMMARY"
echo "=========================================="
echo "All significant expansions:   $SIG_COUNT (sequences: $ALL_SEQ)"
echo "$SPECIES-specific expansions: $SPECIFIC_COUNT (sequences: $SPEC_SEQ)"
echo ""
echo "Output: $OUTPUT_DIR"
echo "Log:    $LOG_FILE"
echo "=========================================="
```

#### 3.2.4 Annotation of Significant Families

The sequences.fa files generated in Section 3.2.2 are located in the significant expansion/contraction output directories, within either the all_significant or specific directory.

The script was run separately for each sequences.fa file requiring annotation.

Functional annotation outputs generated in Section 2.7 were used as input to extract annotations for the significant gene families.

```bash
SEQUENCE_FILE=sequences.fa
EGGNOG_FILE=eggnog.tsv
INTERPRO_FILE=interproscan.tsv
OUTPUT_DIR=Output_directory
LOG_FILE="$OUTPUT_DIR/CAFE5_annotation.log"

mkdir -p "$OUTPUT_DIR"

exec > >(tee -a "$LOG_FILE") 2>&1

# ------------------------------------------------------------------
# STEP 1: EXTRACT GENE IDS FROM FASTA HEADERS
# ------------------------------------------------------------------
echo "Extracting gene IDs from FASTA headers..."

# Convert FASTA headers to the gene ID format used in the eggNOG and InterPro files
grep "^>" "$SEQUENCE_FILE" | \
    sed 's/^>//' | \
    cut -d'|' -f1 | \
    sed 's/^.*_anno/anno/' > "$OUTPUT_DIR/gene_ids.txt"

GENE_COUNT=$(wc -l < "$OUTPUT_DIR/gene_ids.txt")
echo "Found $GENE_COUNT gene IDs"

# ------------------------------------------------------------------
# STEP 2: EXTRACT EGGNOG ANNOTATIONS
# ------------------------------------------------------------------
echo ""
echo "Extracting eggNOG annotations..."

grep -F -f "$OUTPUT_DIR/gene_ids.txt" "$EGGNOG_FILE" | grep -v "^#" > "$OUTPUT_DIR/eggnog_full.tsv"

EGGNOG_COUNT=$(wc -l < "$OUTPUT_DIR/eggnog_full.tsv" 2>/dev/null || echo 0)
echo "  Found $EGGNOG_COUNT eggNOG annotations"

# ------------------------------------------------------------------
# STEP 3: EXTRACT INTERPRO ANNOTATIONS
# ------------------------------------------------------------------
echo ""
echo "Extracting InterPro annotations..."

grep -F -f "$OUTPUT_DIR/gene_ids.txt" "$INTERPRO_FILE" 2>/dev/null > "$OUTPUT_DIR/interpro_full.tsv"

INTERPRO_COUNT=$(wc -l < "$OUTPUT_DIR/interpro_full.tsv" 2>/dev/null || echo 0)
echo "  Found $INTERPRO_COUNT InterPro annotation lines"

# ------------------------------------------------------------------
# STEP 4: EXTRACT UNIQUE GO TERMS
# ------------------------------------------------------------------
echo ""
echo "Extracting unique GO terms..."

# From eggNOG (column 10)
cut -f10 "$OUTPUT_DIR/eggnog_full.tsv" 2>/dev/null | tr ',' '\n' | grep -E "^GO:[0-9]+" > "$OUTPUT_DIR/go_eggnog.txt"

# From InterPro (column 12)
cut -f12 "$OUTPUT_DIR/interpro_full.tsv" 2>/dev/null | tr '|' '\n' | grep -E "^GO:[0-9]+" > "$OUTPUT_DIR/go_interpro.txt"

# Combine and deduplicate
cat "$OUTPUT_DIR/go_eggnog.txt" "$OUTPUT_DIR/go_interpro.txt" 2>/dev/null | sort -u > "$OUTPUT_DIR/go_terms_unique.txt"

GO_COUNT=$(wc -l < "$OUTPUT_DIR/go_terms_unique.txt" 2>/dev/null || echo 0)
echo "  Unique GO terms: $GO_COUNT"

# ------------------------------------------------------------------
# STEP 5: CLEAN UP
# ------------------------------------------------------------------
rm -f "$OUTPUT_DIR/gene_ids.txt" "$OUTPUT_DIR/go_eggnog.txt" "$OUTPUT_DIR/go_interpro.txt"

# ------------------------------------------------------------------
# SUMMARY
# ------------------------------------------------------------------
echo ""
echo "=========================================="
echo "DONE!"
echo "=========================================="
echo "Gene IDs:        $GENE_COUNT"
echo "eggNOG lines:   $EGGNOG_COUNT"
echo "InterPro lines: $INTERPRO_COUNT"
echo "Unique GO terms: $GO_COUNT"
echo ""
echo "Output files:"
echo "  1. eggnog_full.tsv    - all eggNOG columns"
echo "  2. interpro_full.tsv  - all InterPro columns"
echo "  3. go_terms_unique.txt - unique GO IDs"
echo ""
echo "Location: $OUTPUT_DIR"
echo "=========================================="
```

## License

MIT License

Copyright (c) 2026 Katerina Achilleos, Patricia Álvarez-Campos, Marc A. Bailie, Ryan Null, Helena Garcia-Castro, Jordi Solana, B. Duygu Özpolat and Nathan J. Kenny

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
