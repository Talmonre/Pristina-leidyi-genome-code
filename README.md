# Observations on the genome of selectively asexual organism *Pristina leidyi*

Code for the genome assembly, annotation and comparative analyses in:

> Achilleos K\*, Álvarez-Campos P\*, Bailie MA, Null R, Garcia-Castro H, Solana J, Özpolat BD, Kenny NJ. *Observations on the genome of selectively asexual organism* Pristina leidyi. *Genome Biology and Evolution*. \*Co-first authors.

The sequencing data, assembly and annotation are available from NCBI under the accessions given in the paper. This page contains the code only.

The analyses were run on the NeSI and University of Otago SLURM clusters. In the code below, cluster-specific settings (job headers, module and conda environment loading, absolute paths) have been removed. Input and output files are given generic names, set as variables at the top of each block. Unless stated otherwise, programs were run with default settings.

## Code authorship and contributions

This repository was developed by MAB and KA. MAB contributed to Sections 1.1-1.5 and 1.11-1.13. KA contributed to Sections 1.6-1.10, 1.14-1.16 and 2-5.

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
  - [1.10 BlobTools, BUSCO and snail plot (BlobToolkit)](#110-blobTools,-BUSCO-and-snail-plot-(BlobToolkit))
  - [1.11 Hi-C scaffolding, Alignment and YaHS](#111-hi-c-scaffolding,-alignment-and-YaHS)
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

- [5. Repeat annotation and masking](#5-repeat-annotation-and-masking)
  - [5.1 Repeat annotation using EDTA](#51-edta)
  - [5.2 RepeatModeler and RepeatMasker](#52-repeatmodeler-and-repeatmasker)
  - [5.3 Merged masking](#53-merged-masking)
  - [5.4 Venn diagram](#54-venn-diagram)
  - [5.5 Repeat classification and summary](#55-repeat-classification-and-summary)
- [6. Gene prediction](#6-gene-prediction)
  - [6.1 Protein evidence](#61-protein-evidence)
  - [6.2 Short-read RNA-seq alignment](#62-short-read-rna-seq-alignment)
  - [6.3 Iso-Seq transcripts](#63-iso-seq-transcripts)
  - [6.4 BRAKER3](#64-braker3)
  - [6.5 TSEBRA](#65-tsebra)
  - [6.6 Gene models, proteins and BUSCO](#66-gene-models-proteins-and-busco)
- [7. Proteome taxonomic screen](#7-proteome-taxonomic-screen)
- [8. Functional annotation](#8-functional-annotation)
- [9. Phylogenomics](#9-phylogenomics)
  - [9.1 Proteome headers and species mapping](#91-proteome-headers-and-species-mapping)
  - [9.2 OrthoFinder](#92-orthofinder)
  - [9.3 Gene trees and species tree](#93-gene-trees-and-species-tree)
- [10. Gene family evolution](#10-gene-family-evolution)
  - [10.1 Ultrametric tree](#101-ultrametric-tree)
  - [10.2 CAFE5 input](#102-cafe5-input)
  - [10.3 CAFE5](#103-cafe5)
  - [10.4 CafePlotter](#104-cafeplotter)
  - [10.5 Significant expansions and contractions](#105-significant-expansions-and-contractions)
  - [10.6 Annotation of significant families](#106-annotation-of-significant-families)
- [License](#license)

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
| BlobTools | – | 4 |
| BlobToolKit | – | 4 |
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
## 2. Gene prediction and assessment

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





Each interval of the merged mask was given a repeat class from the overlapping annotations (priority EDTA > RepeatModeler > Dfam). The script writes the final repeat BED/GFF3 and summaries per class and per repeat group.

```bash
MERGED_BED=genome.merged_masked_genome.bed
GENOME_FA=genome.fasta
EDTA_GFF=genome.EDTA.TEanno.gff3
RM_REPEATMOD=genome.RepeatMod.out
RM_DFAM=genome.Dfam.out

# ---- Repeat annotations to BED (chrom, start, end, name, class) ----
awk 'BEGIN {OFS = "\t"}
     !/^#/ {
         name = "NA"; class = "NA"
         n = split($9, a, ";")
         for (i = 1; i <= n; i++) {
             if (a[i] ~ /^Name=/)           { name  = a[i]; sub("Name=", "", name) }
             if (a[i] ~ /^classification=/) { class = a[i]; sub("classification=", "", class) }
         }
         print $1, $4 - 1, $5, name, class
     }' "$EDTA_GFF" > EDTA.bed

rm_to_bed() {
    awk 'BEGIN {OFS = "\t"}
         $1 !~ /^SW|^score|^perc|^query|^position/ && NF > 10 {print $5, $6 - 1, $7, $10, $11}' "$1"
}
rm_to_bed "$RM_REPEATMOD" > RepeatModeler.bed
rm_to_bed "$RM_DFAM"      > Dfam.bed

# ---- Classify each merged-mask interval (priority: EDTA > RepeatModeler > Dfam) ----
awk 'BEGIN {OFS = "\t"} {print $0, "mask_" NR}' "$MERGED_BED" > merged_mask_ID.bed
for source in EDTA RepeatModeler Dfam; do
    bedtools intersect -a merged_mask_ID.bed -b "${source}.bed" -wa -wb > "merged_${source}.tsv"
done

awk 'BEGIN {OFS = "\t"}
     FILENAME == "merged_EDTA.tsv"          { class[$5] = $10; next }
     FILENAME == "merged_RepeatModeler.tsv" { if (!($5 in class)) class[$5] = $10; next }
     FILENAME == "merged_Dfam.tsv"          { if (!($5 in class)) class[$5] = $10; next }
     { print $1, $2, $3, (($5 in class) ? class[$5] : "Unclassified") }' \
    merged_EDTA.tsv merged_RepeatModeler.tsv merged_Dfam.tsv merged_mask_ID.bed \
    | sort -k1,1 -k2,2n > final_mask_annotated.bed

awk 'BEGIN {OFS = "\t"; print "##gff-version 3"}
     { print $1, "FinalMask", "repeat_region", $2 + 1, $3, ".", ".", ".",
             "ID=repeat_" NR ";classification=" $4 ";length=" ($3 - $2) }' \
    final_mask_annotated.bed > final_mask_annotated.gff3

# ---- Summaries per repeat class and per repeat group ----
GENOME_SIZE=$(awk '!/^>/ {total += length($0)} END {print total}' "$GENOME_FA")

awk 'BEGIN {OFS = "\t"}
     { bp[$4] += $3 - $2 }
     END { print "Repeat_class\tMasked_bp"; for (c in bp) print c, bp[c] }' \
    final_mask_annotated.bed | sort -k2,2nr > final_repeat_class_summary_bp.txt

awk -v genome="$GENOME_SIZE" '
    function group(c) {
        if (c ~ /LTR/)                      return "LTR_retrotransposons"
        if (c ~ /LINE/)                     return "LINE_retrotransposons"
        if (c ~ /SINE/)                     return "SINE_retrotransposons"
        if (c ~ /PLE/ || c == "Retroposon") return "PLE_retrotransposons"
        if (c ~ /Helitron/ || c ~ /^RC\//)  return "Helitrons"
        if (c ~ /DNA/ || c ~ /MITE/)        return "DNA_transposons"
        if (c ~ /Simple_repeat/)            return "Simple_repeats"
        if (c ~ /Low_complexity/)           return "Low_complexity"
        if (c ~ /Satellite/)                return "Satellite"
        if (c ~ /RNA/)                      return "RNA_genes"
        if (c == "Unknown" || c == "unknown" || c == "Unclassified") return "Unclassified"
        if (c == "ARTEFACT")                return "Artifact"
        return "Other"
    }
    $1 != "Repeat_class" { group_bp[group($1)] += $2 }
    END {
        print "Repeat_group\tMasked_bp\tPercent_of_genome"
        for (g in group_bp) {
            total += group_bp[g]
            printf "%s\t%.0f\t%.4f\n", g, group_bp[g], 100 * group_bp[g] / genome
        }
        printf "Total_masked\t%.0f\t%.4f\n", total, 100 * total / genome
    }' final_repeat_class_summary_bp.txt | sort -k2,2nr > collapsed_repeat_summary.txt
```

## 6. Gene prediction

### 6.1 Protein evidence

```bash
SEARCH_TERM='Clitellata[Organism]'
BATCH_SIZE=500
OUTFILE=clitellata_proteins.fasta

esearch -db protein -query "$SEARCH_TERM" | efetch -format uid > clitellata_ids.txt
split -l "$BATCH_SIZE" clitellata_ids.txt ids_batch_

> "$OUTFILE"
for f in ids_batch_*; do
    efetch -db protein -id "$(paste -sd, "$f")" -format fasta >> "$OUTFILE"
    sleep 1
done
rm ids_batch_* clitellata_ids.txt
```

### 6.2 Short-read RNA-seq alignment

```bash
GENOME=softmasked_genome.fasta
THREADS=32

sorted_bams=()
for R1 in *_1.fastq; do
    R2="${R1/_1.fastq/_2.fastq}"
    [[ -f "$R2" ]] || continue
    sample="${R1%%_1.fastq}"

    minimap2 -ax splice -t "$THREADS" "$GENOME" "$R1" "$R2" -o "${sample}.sam"
    samtools view -@ "$THREADS" -bS "${sample}.sam" > "${sample}.bam"
    samtools sort -@ "$THREADS" "${sample}.bam" -o "${sample}.sorted.bam"
    samtools index "${sample}.sorted.bam"
    sorted_bams+=("${sample}.sorted.bam")
done

samtools merge -@ "$THREADS" merged.bam "${sorted_bams[@]}"
samtools sort  -@ "$THREADS" merged.bam -o merged.sorted.bam
samtools index merged.sorted.bam
```

### 6.3 Iso-Seq transcripts

```bash
GENOME=softmasked_genome.fasta
ISOSEQ=isoseq_hq.fasta
NAME=isoseq

minimap2 -ax map-pb "$GENOME" "$ISOSEQ" -t 32 -o "${NAME}.sam"
samtools view -S -b "${NAME}.sam" > "${NAME}.bam"
samtools sort "${NAME}.bam" -@ 32 -o "${NAME}.sorted.bam"
samtools index "${NAME}.sorted.bam"

stringtie "${NAME}.sorted.bam" -p 32 -o "${NAME}.gtf"
```

### 6.4 BRAKER3

BRAKER3 was run twice: with RNA-seq evidence only, and with RNA-seq plus Clitellata proteins.

```bash
GENOME=softmasked_genome.fasta
BAM=merged.sorted.bam
PROTEINS=clitellata_proteins.fasta
SIF=braker3.sif
export AUGUSTUS_CONFIG_PATH=augustus_config

# RNA-seq only
singularity exec --env AUGUSTUS_CONFIG_PATH="$AUGUSTUS_CONFIG_PATH" "$SIF" \
    braker.pl --genome="$GENOME" --bam="$BAM" --softmasking --threads 64 \
              --workingdir=braker3_rna

# RNA-seq + proteins
singularity exec --env AUGUSTUS_CONFIG_PATH="$AUGUSTUS_CONFIG_PATH" "$SIF" \
    braker.pl --genome="$GENOME" --bam="$BAM" --prot_seq="$PROTEINS" --softmasking --threads 32 \
              --workingdir=braker3_rna_prot

# protein BUSCO of each run
for run in braker3_rna braker3_rna_prot; do
    sed 's/\*//g' "$run/braker.aa" > "$run/braker_clean.aa"
    busco -i "$run/braker_clean.aa" -l metazoa_odb10 -o "${run}_BUSCO" -m prot -c 32
done
```

### 6.5 TSEBRA

The two BRAKER3 gene sets were combined, keeping the Iso-Seq transcript models (`-k`).

```bash
tsebra.py --gtf braker3_rna/braker.gtf,braker3_rna_prot/braker.gtf \
          --hintfiles braker3_rna/hintsfile.gff,braker3_rna_prot/hintsfile.gff \
          -k isoseq.gtf \
          --cfg tsebra.conf \
          --out TSEBRA.gtf \
          -v 2
```

`tsebra.conf`:

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

### 6.6 Gene models, proteins and BUSCO

```bash
GTF=TSEBRA.gtf
FIXED_GTF=TSEBRA_fixed.gtf
GENOME=softmasked_genome.fasta
PREFIX=Pristina

# add gene_id / transcript_id attributes to gene and transcript lines
awk 'BEGIN {OFS = "\t"}
     $3 == "gene"       { print $1, $2, $3, $4, $5, $6, $7, $8, "gene_id \"" $9 "\";" }
     $3 == "transcript" { split($9, a, /\./)
                          print $1, $2, $3, $4, $5, $6, $7, $8, "gene_id \"" a[1] "\"; transcript_id \"" $9 "\";" }
     $3 != "gene" && $3 != "transcript" { print }' "$GTF" > "$FIXED_GTF"

gffread "$FIXED_GTF" -g "$GENOME" -y "TSEBRA_${PREFIX}_proteins.fasta"
gffread "$FIXED_GTF" -g "$GENOME" -x "TSEBRA_${PREFIX}_cds.fasta"

busco -i "TSEBRA_${PREFIX}_proteins.fasta" -l metazoa_odb10 -o "TSEBRA_${PREFIX}_BUSCO" -m prot -c 32
```

## 7. Proteome taxonomic screen

Predicted proteins were searched against NCBI nr and the best hits were mapped back to their scaffolds.

```bash
QUERY=TSEBRA_Pristina_proteins.fasta
DB=nr.dmnd

diamond blastp \
    --query "$QUERY" \
    --db "$DB" \
    --evalue 1e-5 \
    --max-target-seqs 1 \
    --outfmt 6 qseqid sseqid bitscore pident length qlen slen evalue \
    --more-sensitive \
    --threads 32 \
    --out Pristina_blastp.out
```

```r
library(dplyr)
library(readr)
library(stringr)

gtf_file   <- "TSEBRA_fixed.gtf"
blast_file <- "Pristina_blastp.out"
top_n      <- 3

gtf <- read_tsv(gtf_file, comment = "#",
                col_names = c("scaffold", "source", "type", "start", "end",
                              "score", "strand", "phase", "attributes"))

genes <- gtf %>%
  filter(type == "transcript") %>%
  mutate(gene_id = str_extract(attributes, '(?<=transcript_id ")[^"]+')) %>%
  select(scaffold, gene_id)

blast <- read_tsv(blast_file,
                  col_names = c("gene_id", "sseqid", "bitscore", "pident",
                                "length", "qlen", "slen", "evalue"))

genes_with_blast <- genes %>%
  left_join(blast, by = "gene_id") %>%
  na.omit() %>%
  mutate(query_cov      = length / qlen,
         subject_cov    = length / slen,
         combined_score = bitscore * query_cov * subject_cov)

write_csv(arrange(genes_with_blast, scaffold, sseqid), "genes_scaffold_blast.csv")

top_hits_per_scaffold <- genes_with_blast %>%
  group_by(scaffold) %>%
  arrange(desc(combined_score)) %>%
  slice_head(n = top_n) %>%
  ungroup()

write_csv(top_hits_per_scaffold, "bestblast_per_scaffold.csv")
```

## 8. Functional annotation

```bash
PROTEINS=TSEBRA_Pristina_proteins.fasta
PREFIX=Pristina

emapper.py \
    -i "$PROTEINS" \
    --itype proteins \
    -m diamond \
    --cpu 32 \
    --data_dir eggnog_data \
    --dmnd_db eggnog_proteins.dmnd \
    --pfam_realign realign \
    --go_evidence all \
    --target_orthologs all \
    --dbmem \
    --report_orthologs \
    --decorate_gff yes \
    --output_dir eggnog \
    --output "${PREFIX}_eggnog_diamond"

interproscan.sh \
    --input "$PROTEINS" \
    --seqtype p \
    --goterms \
    --pathways \
    --iprlookup \
    --formats tsv \
    --outfile "${PREFIX}_interpro.tsv" \
    --cpu 32
```

## 9. Phylogenomics

### 9.1 Proteome headers and species mapping

Protein headers were prefixed with the species name (file name), and a gene-to-species map was built for ASTRAL.

```bash
PROT_DIR=Prot_Data

for file in "$PROT_DIR"/*.faa; do
    species=$(basename "$file" .faa)
    awk -v sp="$species" '/^>/ {print ">" sp "_" substr($0, 2); next} {print}' "$file" > "${file}.tmp"
    mv "${file}.tmp" "$file"
done

awk '/^>/ {
         h = substr($1, 2)
         if (match(h, /^(.*)_[^_]+$/, m)) print h "\t" m[1]
     }' "$PROT_DIR"/*.faa > mapping.txt
```

### 9.2 OrthoFinder

```bash
orthofinder -f Prot_Data -t 32 -a 32 -o orthofinder_out
```

### 9.3 Gene trees and species tree

```bash
RESULTS=orthofinder_out/Results
OUTGROUP=Pecten_maximus
PREFIX=Pristina
THREADS=32

mkdir -p SCOs Aligned_SCOs Trimmed_SCOs IQTREE/tmp IQTREE/gene_trees

# single-copy orthogroups
while read -r og; do
    cp "$RESULTS/Orthogroup_Sequences/${og}.fa" SCOs/
done < "$RESULTS/Orthogroups/Orthogroups_SingleCopyOrthologues.txt"

# align, trim, concatenate
for f in SCOs/*.fa; do
    base=$(basename "$f" .fa)
    mafft --auto --thread "$THREADS" "$f" > "Aligned_SCOs/${base}.aln.fa"
    clipkit "Aligned_SCOs/${base}.aln.fa" -m smart-gap -o "Trimmed_SCOs/${base}.trim.fa"
done

AMAS.py concat -i Trimmed_SCOs/*.trim.fa -f fasta -d aa \
               -p "${PREFIX}_partitions.txt" \
               -t "${PREFIX}_concatenated_alignment.fa"

# one IQ-TREE run (ModelFinder) per non-empty trimmed alignment, in parallel
for f in Trimmed_SCOs/*.trim.fa; do
    base=$(basename "$f" .trim.fa)
    [ "$(grep -v '^>' "$f" | tr -d '\n' | wc -c)" -gt 0 ] || continue
    echo "iqtree -s $f -m MFP -nt 1 --safe -pre IQTREE/tmp/${base} > /dev/null && cp IQTREE/tmp/${base}.treefile IQTREE/gene_trees/"
done > iqtree_commands.txt
parallel -j "$THREADS" < iqtree_commands.txt

# best-fit model per gene
for f in IQTREE/tmp/*.iqtree; do
    printf '%s\t%s\n' "$(basename "$f" .iqtree)" "$(grep 'Best-fit model' "$f" | awk -F': ' '{print $2}')"
done > all_models.txt

# species tree
cat IQTREE/gene_trees/*.treefile > all_gene_trees.tre
astral4 -a mapping.txt -o "${PREFIX}.tre" all_gene_trees.tre
nw_reroot "${PREFIX}.tre" "$OUTGROUP" > "${PREFIX}_rooted.tre"
```

## 10. Gene family evolution

### 10.1 Ultrametric tree

```r
library(ape)
library(phytools)

my_tree <- read.tree("Pristina_rooted.tre")

# calibrations (TimeTree): P. maximus-L. terrestris and L. terrestris-H. manillensis
root_node       <- getMRCA(my_tree, tip = c("Pecten_maximus", "Lumbricus_terrestris"))
clitellata_node <- getMRCA(my_tree, tip = c("Lumbricus_terrestris", "Hirudinaria_manillensis"))

calibration_df <- data.frame(
  node        = c(root_node, clitellata_node),
  age.min     = c(447.8, 295.9),
  age.max     = c(650.0, 345.0),
  soft.bounds = c(FALSE, FALSE)
)

dated_tree <- chronos(my_tree,
                      calibration = calibration_df,
                      model       = "correlated",
                      quiet       = FALSE,
                      control     = chronos.control())

if (!is.ultrametric(dated_tree)) dated_tree <- force.ultrametric(dated_tree, method = "nnls")
if (!is.binary.tree(dated_tree)) dated_tree <- multi2di(dated_tree, random = FALSE)

write.tree(dated_tree, file = "Pristina_ultrametric_tree_cafe5.tre")
write.nexus(dated_tree, file = "Pristina_ultrametric_tree_cafe5.nex")

node_ages <- branching.times(dated_tree)
write.csv(data.frame(node = names(node_ages), age = node_ages),
          file = "divergence_times.csv", row.names = FALSE)
```

### 10.2 CAFE5 input

`Orthogroups.GeneCount.tsv` from OrthoFinder was reformatted for CAFE5 (Desc column added, Total column removed). Families with ≥ 100 copies in any species were excluded.

```bash
awk 'BEGIN { OFS = "\t" }
     NR == 1 {
         header = $1
         for (i = 2; i < NF; i++) header = header OFS $i
         print "Desc", header > "Orthogroups_filtered.tsv"
         print "Desc", header > "Orthogroups_removed.tsv"
         next
     }
     {
         total++
         data = $1; large = 0
         for (i = 2; i < NF; i++) { data = data OFS $i; if ($i >= 100) large = 1 }
         if (large) { removed++; print "A" total, data > "Orthogroups_removed.tsv" }
         else       { kept++;    print "A" kept,  data > "Orthogroups_filtered.tsv" }
     }
     END { printf "Families kept: %d\nFamilies removed: %d\nTotal families: %d\n", kept, removed, total > "filtering_summary.txt" }' \
    Orthogroups.GeneCount.tsv
```

### 10.3 CAFE5

An error model was estimated first, then used in the base-model run.

```bash
INPUT=Orthogroups_filtered.tsv
TREE=Pristina_ultrametric_tree_cafe5.tre

cafe5 -i "$INPUT" -t "$TREE" -p -e -o error_estimation

cafe5 -i "$INPUT" -t "$TREE" -p \
      -eerror_estimation/Base_error_model.txt \
      -k 1 \
      -o results/k1
```

### 10.4 CafePlotter

```bash
cafeplotter -i results/k1 \
            -o cafe5_plots_k1 \
            --format png \
            --fig_width 12 \
            --fig_height 0.7 \
            --dpi 600 \
            --leaf_label_size 10 \
            --count_label_size 7 \
            --expansion_color red \
            --contraction_color blue \
            --ignore_branch_length
```

### 10.5 Significant expansions and contractions

Families with a significant change (p < 0.05) were split into expansions (positive change on the *P. leidyi* branch) and contractions (negative change). Sequences were then extracted for each set, both for all significant families and for those changed only in *P. leidyi*. Set `DIRECTION=contraction` for contractions.

```bash
SPECIES=Pristina_leidyi
DIRECTION=expansion             # expansion | contraction
CHANGE_FILE=results/k1/Base_change.tab
REPORT_FILE=results/k1/Base_report.cafe
OG_DIR=orthofinder_out/Results/Orthogroup_Sequences
PVAL=0.05
OUT="${SPECIES}_significant_${DIRECTION}s"

mkdir -p "$OUT"/{all_significant,specific}
[ "$DIRECTION" = expansion ] && SIGN=1 || SIGN=-1

# column of the species in the change table
col=$(head -1 "$CHANGE_FILE" | awk -v sp="$SPECIES" '{for (i = 1; i <= NF; i++) if (index($i, sp) == 1) {print i; exit}}')

# families with p < PVAL
grep "^OG" "$REPORT_FILE" \
    | awk -v p="$PVAL" '{for (i = 2; i <= NF; i++) if ($i ~ /^[0-9.]+$/ && $i < p && $i >= 0) {print $1; break}}' \
    > "$OUT/all_significant_ids.txt"

# significant families changed in the given direction in the species
awk -v col="$col" -v s="$SIGN" 'NR == FNR {sig[$1]; next} ($1 in sig) && $col * s > 0 {print $1}' \
    "$OUT/all_significant_ids.txt" "$CHANGE_FILE" > "$OUT/all_significant/families.txt"

# ... and changed in that direction only in the species
awk -v col="$col" -v s="$SIGN" '
    NR == FNR {sig[$1]; next}
    ($1 in sig) && $col * s > 0 {
        for (i = 2; i <= NF; i++) if (i != col && $i * s > 0) next
        print $1
    }' "$OUT/all_significant/families.txt" "$CHANGE_FILE" > "$OUT/specific/families.txt"

# species sequences of each family
extract_seqs() {
    [ -s "$1" ] || return
    > "$2"
    while read -r og; do
        [ -f "$OG_DIR/${og}.fa" ] || continue
        awk -v og="$og" -v sp="$SPECIES" '
            $0 ~ "^>" sp {print $0 "|OG:" og; keep = 1; next}
            /^>/ {keep = 0}
            keep' "$OG_DIR/${og}.fa" >> "$2"
    done < "$1"
}
extract_seqs "$OUT/all_significant/families.txt" "$OUT/all_significant/sequences.fa"
extract_seqs "$OUT/specific/families.txt"        "$OUT/specific/sequences.fa"
```

### 10.6 Annotation of significant families

```bash
SEQUENCES=sequences.fa
EGGNOG=Pristina_eggnog_diamond.emapper.annotations
INTERPRO=Pristina_interpro.tsv
OUT=annotation

mkdir -p "$OUT"

grep "^>" "$SEQUENCES" | sed 's/^>//' | cut -d'|' -f1 | sed 's/^.*_anno/anno/' > "$OUT/gene_ids.txt"

grep -F -f "$OUT/gene_ids.txt" "$EGGNOG"   | grep -v "^#" > "$OUT/eggnog_full.tsv"
grep -F -f "$OUT/gene_ids.txt" "$INTERPRO"                > "$OUT/interpro_full.tsv"

# unique GO terms from both annotations
{
    cut -f10 "$OUT/eggnog_full.tsv"   | tr ',' '\n'
    cut -f12 "$OUT/interpro_full.tsv" | tr '|' '\n'
} | grep -E "^GO:[0-9]+" | sort -u > "$OUT/go_terms_unique.txt"

rm "$OUT/gene_ids.txt"
```

## License

MIT License

Copyright (c) 2026 Katerina Achilleos, Patricia Álvarez-Campos, Marc A. Bailie, Ryan Null, Helena Garcia-Castro, Jordi Solana, B. Duygu Özpolat and Nathan J. Kenny

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
