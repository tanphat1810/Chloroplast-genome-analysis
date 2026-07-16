# Chloroplast Genome Analysis

## Overview

This repository documents a complete workflow for chloroplast genome analysis, including:

1. Raw-read quality assessment and preprocessing.
2. De novo chloroplast genome assembly.
3. Genome annotation and graphical visualization.
4. Structural and compositional feature analysis.
5. Phylogenetic analysis.
6. Relative synonymous codon usage analysis.
7. Selection pressure analysis.
8. Nucleotide diversity analysis.
9. Identification of simple and long sequence repeats.
10. Comparative chloroplast genome visualization using mVISTA.

The workflow is designed for paired-end Illumina sequencing data but may be adapted to other sequencing platforms where appropriate.

![Chloroplast genome analysis workflow](Images/Workflow.svg)

## Workflow Summary

| Step | Analysis | Main tools |
|---|---|---|
| 1 | Read quality control and preprocessing | FastQC, MultiQC, fastp |
| 2 | Chloroplast genome assembly | GetOrganelle, Bandage |
| 3 | Genome annotation and genome-map visualization | GeSeq, CPGAVAS2, OGDRAW |
| 4 | Genome feature analysis | Custom Python scripts |
| 5 | Phylogenetic analysis | MAFFT, MEGA |
| 6 | Relative synonymous codon usage | CodonW |
| 7 | Selection pressure analysis | MAFFT, PAL2NAL, KaKs_Calculator |
| 8 | Nucleotide diversity analysis | DnaSP |
| 9 | SSR and LSR identification | MISA, Vmatch |
| 10 | Comparative genome visualization | mVISTA |

> **Note**
>
> All paths enclosed in angle brackets, such as `<path/to/input.fasta>`, are placeholders and must be replaced with actual file or directory paths.

# 1. Data Preprocessing

Raw paired-end FASTQ files are initially evaluated using **FastQC v0.12.1**. Individual FastQC reports are subsequently summarized using **MultiQC v1.33** to facilitate comparison among sequencing files and samples.

The raw reads are then processed using **fastp v1.0.1** to:

- Remove low-quality bases.
- Trim bases from the beginning of reads where necessary.
- Detect and remove adapter sequences.
- Remove reads containing excessive ambiguous nucleotides.
- Generate cleaned paired-end FASTQ files.
- Produce HTML and JSON quality-control reports.

The cleaned reads are used as input for chloroplast genome assembly.

## 1.1. Quality Control Using FastQC

### Single-end Data

```bash
fastqc \
  <path/to/input.fastq.gz> \
  -o <path/to/fastqc_output_directory> \
  -t 16
```

### Paired-end Data

```bash
fastqc \
  <path/to/input_R1.fastq.gz> \
  <path/to/input_R2.fastq.gz> \
  -o <path/to/fastqc_output_directory> \
  -t 16
```

### Command Explanation

| Option | Description |
|---|---|
| `<path/to/input.fastq.gz>` | Input FASTQ file |
| `-o` | Output directory for FastQC reports |
| `-t 16` | Number of processor threads used |

FastQC generates one HTML report and one compressed result archive for each input FASTQ file.

Important quality metrics include:

- Per-base sequence quality.
- Per-sequence quality scores.
- Per-base nucleotide composition.
- GC-content distribution.
- Sequence duplication levels.
- Adapter content.
- Overrepresented sequences.
- Read-length distribution.

## 1.2. Generate a Summary Report Using MultiQC

```bash
multiqc \
  <path/to/fastqc_output_directory> \
  -o <path/to/multiqc_output_directory>
```

### Command Explanation

| Option | Description |
|---|---|
| `<path/to/fastqc_output_directory>` | Directory containing FastQC results |
| `-o` | Output directory for the MultiQC report |

MultiQC scans the specified input directory and combines results from multiple samples into a single interactive HTML report.

## 1.3. Read Filtering and Trimming Using fastp

```bash
fastp \
  -i <path/to/input_R1.fastq.gz> \
  -I <path/to/input_R2.fastq.gz> \
  -o <path/to/output_R1.fastq.gz> \
  -O <path/to/output_R2.fastq.gz> \
  --trim_front1 10 \
  --trim_front2 10 \
  --detect_adapter_for_pe \
  --n_base_limit 5 \
  --html <path/to/fastp_report.html> \
  --json <path/to/fastp_report.json> \
  --thread 8
```

### Command Explanation

| Option | Description |
|---|---|
| `-i` | Input forward-read file |
| `-I` | Input reverse-read file |
| `-o` | Output cleaned forward-read file |
| `-O` | Output cleaned reverse-read file |
| `--trim_front1 10` | Removes 10 bases from the 5′ end of forward reads |
| `--trim_front2 10` | Removes 10 bases from the 5′ end of reverse reads |
| `--detect_adapter_for_pe` | Automatically detects adapters in paired-end reads |
| `--n_base_limit 5` | Removes reads containing more than five ambiguous `N` bases |
| `--html` | Path to the HTML quality-control report |
| `--json` | Path to the JSON quality-control report |
| `--thread 8` | Uses eight processing threads |

The trimming parameters should be adjusted according to the FastQC results and the characteristics of the sequencing library.

## 1.4. Post-trimming Quality Assessment

FastQC and MultiQC should be run again on the cleaned FASTQ files to confirm that:

- Low-quality bases have been removed.
- Adapter contamination has been reduced.
- Read quality is suitable for genome assembly.
- Excessive read loss has not occurred during filtering.

# 2. Chloroplast Genome Assembly

The cleaned paired-end reads are used for de novo chloroplast genome assembly using **GetOrganelle v1.7.7.1**.

GetOrganelle performs organelle-read recruitment, iterative extension, de Bruijn graph construction, and organelle-genome identification. The `embplant_pt` database is used for assembling plastid genomes from embryophytes.

## 2.1. Activate the Conda Environment

```bash
conda activate getorg
```

The environment name may differ depending on the local installation.

## 2.2. De Novo Assembly Using GetOrganelle

```bash
get_organelle_from_reads.py \
  -1 <path/to/cleaned_R1.fastq.gz> \
  -2 <path/to/cleaned_R2.fastq.gz> \
  -o <path/to/assembly_output_directory> \
  -t 8 \
  -F embplant_pt \
  -R 15 \
  -k 21,45,65,85,105,127
```

### Command Explanation

| Option | Description |
|---|---|
| `-1` | Cleaned forward-read FASTQ file |
| `-2` | Cleaned reverse-read FASTQ file |
| `-o` | Assembly output directory |
| `-t 8` | Number of processor threads |
| `-F embplant_pt` | Uses the embryophyte plastid-genome database |
| `-R 15` | Maximum number of read-extension rounds |
| `-k` | K-mer sizes used for de Bruijn graph assembly |

## 2.3. Main GetOrganelle Outputs

GetOrganelle generates several intermediate and final output files. The most important files for this workflow are:

| File type | Purpose |
|---|---|
| FASTA | Contains one or more candidate assembled chloroplast genome sequences |
| GFA | Contains the assembly graph |
| Log files | Record assembly parameters and execution information |

A FASTA file containing `complete` in its filename generally represents a candidate complete chloroplast genome assembly. However, this output should still be inspected and validated before being accepted as the final genome.

## 2.4. Assembly Graph Inspection Using Bandage

The assembly graph in GFA format is visualized using **Bandage v0.9.0**.

The graph should be inspected to determine whether:

- The chloroplast genome forms a complete circular structure.
- The typical quadripartite structure is present.
- The large single-copy region is identifiable.
- The small single-copy region is identifiable.
- Two inverted-repeat regions are present.
- Alternative paths or unresolved branches remain.
- Contaminating contigs may be present.

The typical chloroplast genome structure consists of:

```text
LSC — IRb — SSC — IRa
```

where:

- `LSC` represents the large single-copy region.
- `SSC` represents the small single-copy region.
- `IRa` and `IRb` represent the two inverted-repeat regions.

## 2.5. Assembly Validation

Before downstream analysis, the assembled chloroplast genome should be checked for:

- Expected genome size.
- Circularity.
- Complete inverted-repeat regions.
- Consistency between the FASTA sequence and assembly graph.
- Absence of obvious contamination.
- Absence of unresolved nucleotide stretches.
- Appropriate sequencing-read coverage.
- Correct orientation of the single-copy regions.

Because chloroplast genomes are circular, different assemblies may begin at different genomic coordinates without representing biologically different structures. Genome orientation and starting coordinates should therefore be standardized before comparative analysis.

# 3. Chloroplast Genome Annotation and Visualization

The complete chloroplast genome sequence in FASTA format is annotated using **GeSeq v2.03**.

Alternatively, **CPGAVAS2** may be used for chloroplast genome annotation.

GeSeq is designed for rapid annotation of organellar genomes, particularly plastid genomes. CPGAVAS2 is an integrated platform for chloroplast genome annotation, visualization, and comparative analysis.

## 3.1. Genome Annotation Using GeSeq

GeSeq web server:

```text
https://chlorobox.mpimp-golm.mpg.de/geseq.html
```

The main input is the assembled chloroplast genome in FASTA format.

Where possible, closely related and well-annotated chloroplast genomes should be selected as reference sequences.

Typical GeSeq outputs include:

- GenBank-format annotation files.
- GFF annotation files.
- Protein and nucleotide FASTA files.
- Annotation tables.
- Graphical genome maps.
- Files for manual review and correction.

## 3.2. Alternative Annotation Using CPGAVAS2

CPGAVAS2 web server:

```text
http://47.96.249.172:16019/analyzer/home
```

CPGAVAS2 can be used as an alternative annotation platform or to compare annotation results obtained from GeSeq.

## 3.3. Annotation Validation

Automated annotation results should be manually inspected before downstream analysis.

The following features should be checked:

- Gene names.
- Gene coordinates.
- Start and stop codons.
- Intron and exon boundaries.
- Trans-spliced genes.
- Gene orientation.
- Duplicated genes in the inverted-repeat regions.
- Partial or pseudogene annotations.
- tRNA predictions.
- rRNA annotations.
- Overlapping genes.
- Unusually short or long coding sequences.

Annotation errors can influence gene counts, codon-usage calculations, selection-pressure analysis, and comparative genomics.

## 3.4. Genome Map Visualization Using OGDRAW

The annotated GenBank file is used to generate a graphical chloroplast genome map using **OGDRAW v1.3.1**.

OGDRAW may be accessed directly or through the Chlorobox platform associated with GeSeq.

Recommended output formats include:

- SVG for further editing.
- PDF for publication-quality figures.
- PNG for presentations and web display.

A typical chloroplast genome map displays:

- Gene names and functional groups.
- Gene orientation.
- LSC, SSC, IRa, and IRb regions.
- GC-content variation.
- Genome coordinates.

# 4. Chloroplast Genome Feature Analysis

Following genome annotation, the GenBank files are processed using custom Python scripts to summarize the structural, functional, and compositional characteristics of each chloroplast genome.

## 4.1. Genome Structure

The following structural measurements are extracted:

- Total chloroplast genome length.
- Length of the LSC region.
- Length of the SSC region.
- Length of the IRa region.
- Length of the IRb region.

For a typical quadripartite chloroplast genome:

```math
L_{\mathrm{genome}}
=
L_{\mathrm{LSC}}
+
L_{\mathrm{SSC}}
+
L_{\mathrm{IRa}}
+
L_{\mathrm{IRb}}
```

Because the two inverted-repeat regions are normally similar in length:

```math
L_{\mathrm{genome}}
\approx
L_{\mathrm{LSC}}
+
L_{\mathrm{SSC}}
+
2L_{\mathrm{IR}}
```

## 4.2. Coding Sequence Characteristics

The following coding-sequence characteristics are calculated:

- Total CDS length.
- Proportion of CDS relative to the complete genome.
- Number of protein-coding genes.
- Number of unique protein-coding genes.
- Number of duplicated protein-coding genes located in the IR regions.

The CDS proportion is calculated as:

```math
\mathrm{CDS\ proportion}\ (\%)
=
\frac{\mathrm{Total\ CDS\ length}}
{\mathrm{Genome\ length}}
\times 100
```

It is important to specify whether duplicated gene copies in the inverted repeats are counted once as unique genes or counted separately as physical gene copies.

## 4.3. RNA Gene Characteristics

The following RNA-gene features are summarized:

- Number of tRNA genes.
- Number of unique tRNA genes.
- Number of rRNA genes.
- Number of unique rRNA genes.
- Number of duplicated RNA genes in the inverted-repeat regions.

## 4.4. Nucleotide Composition

GC content is calculated for:

- The complete chloroplast genome.
- The LSC region.
- The SSC region.
- The IR regions.

GC content is calculated as:

```math
\mathrm{GC}\ (\%)
=
\frac{G+C}{A+T+G+C}
\times 100
```

The inverted-repeat regions commonly have higher GC content than the single-copy regions because they contain duplicated ribosomal RNA genes with relatively high GC content.

## 4.5. Codon-position Nucleotide Composition

For all protein-coding sequences, nucleotide composition is calculated at:

- The first codon position.
- The second codon position.
- The third codon position.

The following values may be reported:

- AT1.
- AT2.
- AT3.
- GC1.
- GC2.
- GC3.

Codon-position composition is useful for evaluating nucleotide-composition bias and synonymous codon usage.

## 4.6. Intron-containing Genes

Detailed information is extracted for genes containing introns, including:

- Gene name.
- Strand orientation.
- Number of exons.
- Number of introns.
- Length of each exon.
- Length of each intron.
- Start and end coordinates.
- Total gene length.
- Genomic region in which the gene is located.
- Cis-spliced or trans-spliced status, where applicable.

The results are summarized in tabular format for downstream comparative analysis.

## 4.7. Recommended Output Tables

The custom scripts may generate the following tables:

| Output table | Main contents |
|---|---|
| Genome summary | Genome length, LSC, SSC, IRa, IRb, and GC content |
| Gene summary | Protein-coding genes, tRNAs, rRNAs, and duplicated genes |
| CDS summary | CDS length and CDS proportion |
| Codon-position summary | AT and GC content at codon positions 1, 2, and 3 |
| Intron summary | Intron-containing genes and exon/intron statistics |

# 5. Phylogenetic Analysis

Complete chloroplast genome sequences generated in this study are analyzed together with publicly available chloroplast genome sequences from related taxa, including representatives of the same genus, species, subspecies, or closely related groups.

## 5.1. Sequence Collection

Reference chloroplast genomes should be obtained from reliable public databases such as GenBank.

For each reference sequence, record:

- Species name.
- Subspecies or variety, where applicable.
- GenBank accession number.
- Sequence length.
- Annotation status.
- Publication or source information.

At least one appropriate outgroup should be included to root the phylogenetic tree.

## 5.2. Multiple-sequence Alignment Using MAFFT

The complete chloroplast genome sequences are aligned using **MAFFT v7.490**.

Example command:

```bash
mafft \
  --auto \
  --thread 8 \
  <path/to/input_genomes.fasta> \
  > <path/to/aligned_genomes.fasta>
```

### Command Explanation

| Option | Description |
|---|---|
| `--auto` | Automatically selects an appropriate alignment strategy |
| `--thread 8` | Uses eight processing threads |
| `<path/to/input_genomes.fasta>` | Input multi-FASTA file |
| `>` | Redirects the alignment to an output file |

The alignment should be inspected for:

- Large unaligned regions.
- Incorrect genome orientations.
- Inverted SSC regions.
- Excessive gaps.
- Misaligned inverted-repeat boundaries.
- Poor-quality or incomplete reference sequences.

## 5.3. Phylogenetic Tree Reconstruction Using MEGA

The aligned chloroplast genome sequences are imported into **MEGA 12**.

Phylogenetic trees may be reconstructed using an appropriate method such as:

- Maximum Likelihood.
- Neighbor Joining.
- Maximum Parsimony.

For publication-quality analysis, Maximum Likelihood is generally preferred when an appropriate nucleotide-substitution model is selected.

Branch support is assessed using:

```text
1,000 bootstrap replicates
```

The following information should be reported:

- Tree-building method.
- Nucleotide-substitution model.
- Treatment of gaps and missing data.
- Number of bootstrap replicates.
- Outgroup used.
- Alignment length.
- Number of taxa.
- Software version.

Bootstrap values indicate the proportion of resampled datasets supporting a particular branch. They should not be interpreted as direct probabilities that a clade is correct.

# 6. Relative Synonymous Codon Usage Analysis

## 6.1. Definition

**Relative Synonymous Codon Usage**, abbreviated as **RSCU**, measures how frequently a synonymous codon is used relative to the expected frequency under equal use of all synonymous codons encoding the same amino acid.

RSCU is useful for examining:

- Codon-usage bias.
- Nucleotide-composition bias.
- Preferred and underrepresented codons.
- Evolutionary patterns in protein-coding genes.
- Differences in codon usage among genomes or taxa.

## 6.2. General Formula

```math
RSCU
=
\frac{\mathrm{Observed\ frequency\ of\ a\ codon}}
{\mathrm{Expected\ frequency\ under\ equal\ synonymous\ codon\ usage}}
```

If amino acid \(i\) is encoded by \(n_i\) synonymous codons, the RSCU of codon \(j\) is calculated as:

```math
RSCU_{ij}
=
\frac{X_{ij}}
{\frac{1}{n_i}\sum_{j=1}^{n_i}X_{ij}}
```

where:

- \(X_{ij}\) is the observed number of occurrences of codon \(j\) encoding amino acid \(i\).
- \(n_i\) is the number of synonymous codons encoding amino acid \(i\).

## 6.3. Interpretation

- **RSCU = 1**: the codon is used at the frequency expected under equal synonymous codon usage.
- **RSCU > 1**: the codon is used more frequently than expected.
- **RSCU < 1**: the codon is used less frequently than expected.
- **RSCU > 1.6**: the codon is often considered strongly preferred.
- **RSCU < 0.6**: the codon is often considered strongly underrepresented.

The thresholds of 1.6 and 0.6 are commonly used descriptive criteria rather than universal statistical cutoffs.

Methionine and tryptophan are each encoded by only one codon. Their RSCU values are therefore expected to equal 1.

Stop codons should normally be excluded from analyses of synonymous codon usage unless they are being investigated separately.

## 6.4. CDS Extraction

Protein-coding sequences are extracted from the annotated GenBank file.

Before analysis, each CDS should be checked to ensure that:

- The sequence is in the correct 5′–3′ orientation.
- The length is divisible by three.
- The reading frame is correct.
- No internal stop codons are present.
- Partial or pseudogene sequences are excluded where appropriate.
- Duplicate copies in the inverted repeats are handled consistently.
- The terminal stop codon is treated consistently among all sequences.

## 6.5. RSCU Calculation Using CodonW

RSCU values are calculated using **CodonW v1.4.4**.

```bash
codonw \
  <path/to/cds.fasta> \
  <path/to/RSCU.out> \
  <path/to/RSCU.blk> \
  -cutot \
  -nomenu \
  -silent
```

The same command may be written on one line:

```bash
codonw <path/to/cds.fasta> <path/to/RSCU.out> <path/to/RSCU.blk> -cutot -nomenu -silent
```

### Command Explanation

| Argument or option | Description |
|---|---|
| `<path/to/cds.fasta>` | FASTA file containing protein-coding sequences |
| `<path/to/RSCU.out>` | Standard CodonW output file |
| `<path/to/RSCU.blk>` | Bulk-output file containing total codon-usage results |
| `-cutot` | Calculates total codon usage across the input dataset |
| `-nomenu` | Disables the interactive menu |
| `-silent` | Runs without interactive prompts |

The `-cutot` option pools codon counts from all input CDSs and calculates codon usage for the complete dataset.

Therefore, the resulting RSCU values represent overall codon usage across all included coding sequences rather than gene-specific RSCU values.

## 6.6. Recommended RSCU Outputs

The final table should include:

| Column | Description |
|---|---|
| Amino acid | Amino acid encoded by the codon |
| Codon | DNA or RNA codon |
| Codon count | Number of codon occurrences |
| RSCU | Relative synonymous codon usage |
| Classification | Preferred, neutral, or underrepresented |

# 7. Selection Pressure Analysis

## 7.1. Definition

Selection pressure analysis evaluates the evolutionary constraints acting on protein-coding genes.

It is commonly performed using the ratio of nonsynonymous substitutions to synonymous substitutions, expressed as:

- Ka/Ks.
- dN/dS.
- \(\omega\).

## 7.2. Formula

```math
\omega
=
\frac{K_a}{K_s}
=
\frac{dN}{dS}
```

where:

- **Ka or dN** is the rate of nonsynonymous substitutions that alter the encoded amino acid.
- **Ks or dS** is the rate of synonymous substitutions that do not alter the encoded amino acid.
- **ω** is the Ka/Ks or dN/dS ratio used to characterize the dominant selection regime.

## 7.3. Interpretation

| Ka/Ks value | Selection regime | Interpretation |
|---|---|---|
| **Ka/Ks < 1** | Purifying or negative selection | Amino-acid-changing substitutions are preferentially removed, indicating functional conservation |
| **Ka/Ks ≈ 1** | Neutral evolution | Nonsynonymous and synonymous substitutions accumulate at approximately similar relative rates |
| **Ka/Ks > 1** | Positive selection | Amino-acid-changing substitutions may have been favored during evolution |

## 7.4. Important Interpretation Considerations

A Ka/Ks value greater than 1 should not automatically be considered definitive evidence of positive selection.

Results may be unstable when:

- Ks is zero or close to zero.
- Only one or a few substitutions are present.
- The gene is short.
- Sequences are highly similar.
- Sequences are incorrectly aligned.
- Nonhomologous genes are compared.
- Pseudogenes or partial CDSs are included.
- Only a single pairwise comparison is available.
- Statistical significance has not been evaluated.

A high Ka/Ks value based on one nonsynonymous substitution and no synonymous substitutions should be interpreted cautiously.

## 7.5. Input Preparation

Homologous protein-coding genes are extracted from each chloroplast genome.

For each gene:

1. Nucleotide CDSs are translated into amino-acid sequences.
2. Protein sequences are aligned using MAFFT.
3. The protein alignment is converted into a codon-aware nucleotide alignment using PAL2NAL.
4. The codon alignment is analyzed using KaKs_Calculator.

Codon-aware alignment is essential because direct nucleotide alignment may introduce frameshifted or biologically invalid codon comparisons.

## 7.6. Ka/Ks Calculation

Ka/Ks values are calculated using a custom pipeline incorporating **KaKs_Calculator v3.0**, MAFFT, and PAL2NAL.

```bash
python kaksPipeline.py \
  -i <path/to/input_directory> \
  -o <path/to/output_directory> \
  --threads 8 \
  --mafft-threads 1 \
  --pal2nal pal2nal.pl \
  --mafft-options=--auto \
  --keep-intermediate \
  --verbose \
  --log-file <path/to/pipeline.log>
```

### Command Explanation

| Option | Description |
|---|---|
| `-i` | Directory containing homologous CDS input files |
| `-o` | Output directory |
| `--threads 8` | Maximum number of parallel processing threads |
| `--mafft-threads 1` | Number of threads assigned to each MAFFT process |
| `--pal2nal pal2nal.pl` | Path to the PAL2NAL script |
| `--mafft-options=--auto` | Uses the MAFFT automatic strategy |
| `--keep-intermediate` | Retains intermediate alignment files |
| `--verbose` | Prints detailed execution information |
| `--log-file` | Writes pipeline messages to a log file |

Because `kaksPipeline.py` is a custom script, its exact input structure and internal behavior depend on the implementation included in this repository.

## 7.7. Z-score Calculation

Standardized Z-scores may be calculated to compare gene-level Ka/Ks values across the dataset.

The general formula is:

```math
Z
=
\frac{x-\mu}{\sigma}
```

where:

- \(x\) is the observed value for a gene.
- \(\mu\) is the mean of the analyzed values.
- \(\sigma\) is the standard deviation.

Example command:

```bash
python zscore.py \
  -i <path/to/input.tsv> \
  -o <path/to/output_directory> \
  --input-format long \
  --mode gene \
  --cap-value 5 \
  --inf-policy cap \
  --ks-zero-policy na \
  --min-non-na 2 \
  --z-cap 3
```

### Command Explanation

| Option | Description |
|---|---|
| `-i` | Input Ka/Ks table |
| `-o` | Output directory |
| `--input-format long` | Specifies long-format input data |
| `--mode gene` | Calculates values at the gene level |
| `--cap-value 5` | Caps extreme Ka/Ks values at 5 where applicable |
| `--inf-policy cap` | Replaces infinite values using the selected cap |
| `--ks-zero-policy na` | Converts cases with Ks = 0 to missing values |
| `--min-non-na 2` | Requires at least two valid observations |
| `--z-cap 3` | Caps extreme Z-scores at ±3 |

The exact behavior of these options depends on the implementation of the custom `zscore.py` script.

## 7.8. Recommended Output Fields

| Field | Description |
|---|---|
| Gene | Gene name |
| Sequence pair | Two sequences compared |
| Ka | Nonsynonymous substitution rate |
| Ks | Synonymous substitution rate |
| Ka/Ks | Selection-pressure ratio |
| Z-score | Standardized gene-level value |
| Interpretation | Purifying, neutral, or candidate positive selection |
| Warning | Low Ks, missing value, saturation, or insufficient variation |

# 8. Nucleotide Diversity Analysis

## 8.1. Definition

Nucleotide diversity, denoted by **π**, measures the average number of nucleotide differences per site between all possible pairs of sequences in a dataset.

DnaSP calculates nucleotide diversity and other population-genetic statistics from a multiple-sequence alignment.

Nucleotide diversity may be calculated for:

- The complete chloroplast genome.
- Individual protein-coding genes.
- Intergenic spacer regions.
- LSC, SSC, and IR regions.
- Sliding windows across the complete genome.

In this workflow, nucleotide diversity is calculated for:

- Each chloroplast gene.
- The complete chloroplast genome.

## 8.2. Simplified Formula

For \(n\) aligned sequences with an analyzed alignment length of \(L\):

```math
\pi =
\frac{\sum_{i<j} d_{ij}}
{\binom{n}{2} \times L}
```

where:

- \(d_{ij}\) is the number of nucleotide differences between sequences \(i\) and \(j\).
- \(\binom{n}{2}\) is the total number of possible sequence pairs.
- \(L\) is the number of nucleotide sites analyzed.
- \(\pi\) is the average number of pairwise nucleotide differences per site.

The exact calculation may depend on how gaps, missing data, and ambiguous nucleotides are treated in the selected DnaSP settings.

## 8.3. Interpretation

- **π = 0**: all analyzed sequences are identical at the included positions.
- **Low π**: the sequence region is relatively conserved.
- **High π**: the sequence region contains comparatively greater nucleotide variation.

There is no universal threshold for classifying π as high or low.

Values should be interpreted by comparison among:

- Genes.
- Genomic regions.
- Taxonomic groups.
- Populations.
- Species.
- Datasets generated using consistent methods.

## 8.4. Relevant DnaSP Statistics

| Statistic | Meaning |
|---|---|
| **Pi or π** | Average number of nucleotide differences per site among sequence pairs |
| **S** | Number of segregating or polymorphic sites |
| **Eta or η** | Estimated total number of mutations |
| **Theta-W or θw** | Watterson’s nucleotide-diversity estimator based on segregating sites |
| **Hd** | Haplotype diversity |
| **K** | Average number of nucleotide differences between sequence pairs |

A region may contain many polymorphic sites but still have a moderate π value if most variants occur in only one or a few sequences.

## 8.5. Gene-level Nucleotide Diversity Analysis

For gene-level analysis:

1. Homologous gene sequences are extracted from all genomes.
2. Each gene is aligned separately.
3. Alignment quality is manually checked.
4. The aligned gene is imported into DnaSP.
5. Nucleotide diversity is calculated.
6. Gene-level π values are exported and compared.

Genes with relatively high π values may be useful candidates for:

- DNA barcoding.
- Molecular marker development.
- Population-genetic analysis.
- Phylogenetic reconstruction.
- Taxonomic discrimination.

## 8.6. Whole-genome Nucleotide Diversity Analysis

Complete chloroplast genome sequences are aligned and imported into DnaSP.

Before calculation, the following should be checked:

- Genome orientation.
- Starting coordinates.
- SSC orientation.
- IR duplication.
- Alignment gaps.
- Ambiguous sites.
- Large insertion/deletion regions.
- Missing or incomplete sequence regions.

## 8.7. Sliding-window Analysis

DnaSP may calculate nucleotide diversity using a sliding-window approach.

Example parameters:

```text
Window length: 600 bp
Step size: 200 bp
```

The windows would be evaluated as:

```text
1–600 bp
201–800 bp
401–1000 bp
601–1200 bp
...
```

A smaller window provides higher spatial resolution but may produce noisier estimates. A larger window produces smoother estimates but may conceal short highly variable regions.

The window and step sizes must be reported in the Methods section and figure legends.

## 8.8. Recommended Output Table

| Column | Description |
|---|---|
| Gene or region | Name of the analyzed sequence region |
| Alignment length | Number of aligned nucleotide sites |
| Number of sequences | Number of sequences included |
| S | Number of segregating sites |
| Eta | Estimated number of mutations |
| Pi | Nucleotide diversity |
| Hd | Haplotype diversity, where applicable |
| Notes | Missing data, gaps, or quality-control comments |

# 9. Identification of Simple Sequence Repeats and Long Sequence Repeats

## 9.1. Overview

Repetitive sequences contribute to chloroplast genome variation and may be useful for:

- Molecular marker development.
- Population-genetic studies.
- Species identification.
- Genome-comparison studies.
- Investigation of structural variation.
- Examination of possible recombination-associated regions.

Two major repeat categories are examined:

- Simple sequence repeats using MISA.
- Long sequence repeats using Vmatch.

# 9.2. Simple Sequence Repeats

## 9.2.1. Definition

**Simple sequence repeats**, abbreviated as **SSRs**, are also known as microsatellites.

They consist of short nucleotide motifs repeated consecutively in tandem.

Examples include:

| SSR class | Motif length | Example |
|---|---:|---|
| Mononucleotide | 1 nucleotide | `AAAAAAAAAA` |
| Dinucleotide | 2 nucleotides | `ATATATATAT` |
| Trinucleotide | 3 nucleotides | `AATAATAATAAT` |
| Tetranucleotide | 4 nucleotides | `AAATAAATAAAT` |
| Pentanucleotide | 5 nucleotides | `AAAATAAAATAAAAT` |
| Hexanucleotide | 6 nucleotides | `AAAAATAAAAATAAAAAT` |

## 9.2.2. SSR Identification Using MISA

SSRs are identified using **MISA v2.1**, the MIcroSAtellite identification tool.

MISA detects:

- Perfect microsatellites.
- Compound microsatellites.
- Motif type.
- Repeat number.
- Start and end coordinates.
- SSR length.

## 9.2.3. Example MISA Configuration

A commonly used repeat-threshold configuration is:

```text
definition(unit_size,min_repeats): 1-10 2-5 3-4 4-3 5-3 6-3
interruptions(max_difference_for_2_SSRs): 100
```

This configuration means:

| Motif class | Minimum number of repeats |
|---|---:|
| Mononucleotide | 10 |
| Dinucleotide | 5 |
| Trinucleotide | 4 |
| Tetranucleotide | 3 |
| Pentanucleotide | 3 |
| Hexanucleotide | 3 |

These thresholds are analysis parameters rather than universal biological definitions. The selected thresholds must be reported when comparing SSR results among studies.

The interruption parameter determines the maximum distance between two adjacent SSRs for them to be classified as a compound SSR.

## 9.2.4. Run MISA

Depending on the installation, MISA may be run using:

```bash
perl misa.pl <path/to/genome.fasta>
```

or:

```bash
misa.pl <path/to/genome.fasta>
```

Typical output files include:

```text
<genome.fasta>.misa
<genome.fasta>.statistics
```

The `.misa` file generally contains individual SSR records, whereas the statistics file summarizes SSR counts and motif classes.

## 9.2.5. Recommended SSR Output Fields

| Field | Description |
|---|---|
| SSR ID | Unique identifier |
| SSR type | Perfect or compound SSR |
| Motif | Repeated nucleotide unit |
| Repeat number | Number of motif repetitions |
| Size | Total SSR length |
| Start | Start coordinate |
| End | End coordinate |
| Region | LSC, SSC, IRa, or IRb |
| Annotation | CDS, intron, intergenic spacer, tRNA, or rRNA |
| Gene or region name | Associated annotated feature |

When comparing circular chloroplast genomes, SSRs spanning the selected linearization boundary should be interpreted carefully.

# 9.3. Long Sequence Repeats

## 9.3.1. Definition

**Long sequence repeats**, abbreviated as **LSRs**, are relatively long sequence segments that occur at two or more positions within a genome.

Unlike SSRs, their copies do not need to be located directly adjacent to one another.

Long repeats may contribute to:

- Genome structural variation.
- Recombination.
- Insertions or deletions.
- Rearrangements.
- Genome instability.
- Evolutionary divergence.

## 9.3.2. Repeat Categories

The repeat categories most commonly reported in chloroplast genome studies include:

| Repeat type | Description |
|---|---|
| Forward repeat | Both copies occur in the same orientation |
| Palindromic repeat | The second copy corresponds to the reverse complement of the first |
| Reverse repeat | The second copy occurs in reverse order |
| Complementary repeat | The second copy is complementary to the first |

Vmatch directly supports searches for direct and palindromic matches. Additional classifications may require alternative search settings or post-processing.

## 9.3.3. LSR Identification Using Vmatch

Long sequence repeats are detected using **Vmatch v2.3.1**.

Vmatch performs large-scale sequence matching and allows the user to specify:

- Minimum repeat length.
- Match orientation.
- Number of allowed mismatches.
- Gap conditions.
- Sequence-similarity criteria.

## 9.3.4. Build a Vmatch Index

```bash
mkvtree \
  -db <path/to/genome.fasta> \
  -dna \
  -pl \
  -allout \
  -indexname <path/to/output_index>
```

### Command Explanation

| Option | Description |
|---|---|
| `-db` | Input FASTA sequence |
| `-dna` | Specifies DNA input |
| `-pl` | Stores sequence-position information |
| `-allout` | Generates all required index components |
| `-indexname` | Prefix for the generated index |

## 9.3.5. Identify Repeated Sequences

```bash
vmatch \
  -l 20 \
  -d \
  -p \
  <path/to/input_index> \
  > <path/to/output.txt>
```

### Command Explanation

| Option | Description |
|---|---|
| `-l 20` | Sets the minimum repeat length to 20 bp |
| `-d` | Reports direct matches |
| `-p` | Reports palindromic matches |
| `<path/to/input_index>` | Prefix of the Vmatch index |
| `>` | Redirects output to a text file |

The minimum repeat-length threshold should be selected consistently across all genomes being compared.

## 9.3.6. Recommended LSR Output Fields

| Field | Description |
|---|---|
| Repeat ID | Unique repeat identifier |
| Repeat type | Forward, palindromic, reverse, or complementary |
| Repeat length | Length of the repeated segment |
| Copy 1 start | Start coordinate of the first copy |
| Copy 1 end | End coordinate of the first copy |
| Copy 2 start | Start coordinate of the second copy |
| Copy 2 end | End coordinate of the second copy |
| Region | LSC, SSC, IRa, or IRb |
| Annotation | CDS, intron, or intergenic region |
| Sequence identity | Similarity between repeat copies |

## 9.4. Comparison Between SSRs and LSRs

| Characteristic | SSR | LSR |
|---|---|---|
| Full name | Simple sequence repeat | Long sequence repeat |
| Alternative name | Microsatellite | Long repeat |
| Basic structure | Short motif repeated in tandem | Longer sequence occurring at multiple positions |
| Typical motif or repeat length | Motifs of 1–6 bp | Commonly ≥20 or ≥30 bp, depending on the study |
| Copy arrangement | Adjacent tandem copies | Copies may be separated |
| Main software | MISA | Vmatch |
| Common application | Marker development and population genetics | Structural and evolutionary genome analysis |

# 10. Comparative Chloroplast Genome Visualization Using mVISTA

## 10.1. Overview

mVISTA is used to compare complete chloroplast genome sequences and visualize sequence conservation relative to an annotated reference genome.

The analysis can help identify:

- Conserved coding regions.
- Variable intergenic regions.
- Divergent introns.
- Potential mutation hotspots.
- Differences among LSC, SSC, and IR regions.
- Regions suitable for molecular-marker development.

One chloroplast genome should be selected as the reference sequence. The reference should preferably have a complete and manually reviewed annotation.

## 10.2. Prepare mVISTA Input Files

A custom Python script, `MvistaAnnotation.py`, is used to convert an annotated GenBank file into:

1. A nucleotide FASTA file.
2. An mVISTA-compatible annotation file.

```bash
python MvistaAnnotation.py \
  -i <path/to/input.gbf> \
  -a <path/to/output.mvista.ann> \
  -f <path/to/output.fasta>
```

Example:

```bash
python MvistaAnnotation.py \
  -i BG/Annotation/CPGAVAS2/BG.gbf \
  -a mVISTA/BG.mvista.ann \
  -f mVISTA/BG.fasta
```

## 10.3. Command Explanation

| Argument | Description |
|---|---|
| `python` | Runs the script using the Python interpreter |
| `MvistaAnnotation.py` | Custom GenBank-to-mVISTA conversion script |
| `-i` | Input annotated GenBank file |
| `-a` | Output mVISTA annotation file |
| `-f` | Output FASTA file |

## 10.4. Input GenBank Requirements

The GenBank file should contain:

- A complete chloroplast genome sequence.
- A valid `FEATURES` section.
- Gene coordinates.
- Strand information.
- Gene names.
- CDS features.
- tRNA features.
- rRNA features.
- Intron and exon coordinates where applicable.

The exact feature types processed depend on the implementation of `MvistaAnnotation.py`.

## 10.5. Expected Outputs

The command generates:

```text
mVISTA/
├── BG.fasta
└── BG.mvista.ann
```

### FASTA Output

The FASTA file contains the complete nucleotide sequence:

```fasta
>BG
ATGCGTATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA
TCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
```

### mVISTA Annotation Output

A simplified annotation structure may appear as:

```text
> 1000 2500 rbcL
1000 2500 exon

< 5000 6200 ndhF
5000 6200 exon
```

In this format:

- `>` represents a feature on the forward strand.
- `<` represents a feature on the reverse strand.
- The first two numbers indicate start and end coordinates.
- The final field contains the gene name.
- Subsequent lines define exon or untranslated-region coordinates.

The exact output format depends on the custom script.

## 10.6. Create the Output Directory

The output directory must exist before running the script.

```bash
mkdir -p mVISTA
```

## 10.7. Submit Files to mVISTA

The generated FASTA and annotation files may be submitted through the mVISTA interface.

Before submission:

- Ensure that all genomes are in the same orientation.
- Standardize genome starting coordinates.
- Verify the SSC orientation.
- Confirm that sequence names are short and informative.
- Remove unexpected spaces or special characters from FASTA headers.
- Ensure that the reference annotation coordinates match the reference FASTA sequence.

## 10.8. Interpretation of mVISTA Results

mVISTA plots generally display:

- Reference gene annotations.
- Genome coordinates.
- Percentage sequence identity.
- Conserved coding regions.
- Conserved noncoding regions.
- Regions with reduced sequence similarity.

Highly variable regions may be investigated as candidates for:

- DNA barcodes.
- Species-specific markers.
- Population-genetic markers.
- Phylogenetic loci.

However, apparent divergence may also result from:

- Poor alignment.
- Different genome orientations.
- Annotation errors.
- Missing sequence data.
- Incorrect IR boundaries.
- Structural rearrangements.

# 11. Reproducibility and Quality Control

## 11.1. Software Versions

The software versions used in the analysis should be recorded explicitly.

| Software | Version used |
|---|---:|
| FastQC | 0.12.1 |
| MultiQC | 1.33 |
| fastp | 1.0.1 |
| GetOrganelle | 1.7.7.1 |
| Bandage | 0.9.0 |
| GeSeq | 2.03 |
| OGDRAW | 1.3.1 |
| MAFFT | 7.490 |
| MEGA | 12 |
| CodonW | 1.4.4 |
| KaKs_Calculator | 3.0 |
| DnaSP | Record the exact local version |
| MISA | 2.1 |
| Vmatch | 2.3.1 |

## 11.2. Recommended Reproducibility Records

For each analysis, retain:

- Command-line commands.
- Software versions.
- Conda environment files.
- Configuration files.
- Input checksums.
- Log files.
- Intermediate files.
- Final output tables.
- Manual corrections.
- Reference accession numbers.
- Analysis dates.

## 11.3. Suggested Directory Structure

```text
Chloroplast_Genome_Analysis/
├── README.md
├── Images/
│   └── Workflow.svg
├── Raw_Data/
├── QC/
│   ├── FastQC_raw/
│   ├── MultiQC_raw/
│   ├── FastQC_clean/
│   └── MultiQC_clean/
├── Clean_Reads/
├── Assembly/
├── Annotation/
│   ├── GeSeq/
│   └── CPGAVAS2/
├── Genome_Features/
├── Phylogeny/
├── RSCU/
├── Selection_Pressure/
├── Nucleotide_Diversity/
├── Repeats/
│   ├── SSR/
│   └── LSR/
├── mVISTA/
├── Scripts/
└── Results/
```

# 12. Recommended Methods Description

The following paragraph may be adapted for a manuscript:

Raw paired-end reads were assessed using FastQC v0.12.1, and quality-control results were summarized using MultiQC v1.33. Adapter removal, quality trimming, and read filtering were performed using fastp v1.0.1. Cleaned reads were assembled de novo using GetOrganelle v1.7.7.1 with the embryophyte plastid database. Assembly graphs were inspected using Bandage v0.9.0. Complete chloroplast genomes were annotated using GeSeq v2.03 or CPGAVAS2, and graphical genome maps were generated using OGDRAW v1.3.1. Structural and compositional genome characteristics were summarized from the annotated GenBank files using custom Python scripts. Complete chloroplast genome sequences were aligned using MAFFT v7.490, and phylogenetic trees were reconstructed in MEGA 12 with 1,000 bootstrap replicates. Relative synonymous codon usage was calculated using CodonW v1.4.4. Selection pressure was evaluated using Ka/Ks ratios generated with KaKs_Calculator v3.0 following codon-aware alignment. Nucleotide diversity was calculated using DnaSP. Simple sequence repeats were identified using MISA v2.1, whereas long sequence repeats were detected using Vmatch v2.3.1. Comparative chloroplast genome conservation was visualized using mVISTA.

# 13. References and Tool Resources

## Quality Control and Read Processing

- FastQC: <https://www.bioinformatics.babraham.ac.uk/projects
