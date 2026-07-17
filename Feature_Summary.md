# Chloroplast Genome Feature Extractor

A Python script for extracting comprehensive genomic features from annotated chloroplast genome GenBank files.

This tool was developed to automate the extraction of structural, compositional, and gene-related information from chloroplast genomes, facilitating comparative genomics and downstream evolutionary analyses.

---

## Features

The script automatically extracts the following information from one or multiple annotated chloroplast GenBank files.

### 1. Genome structure

- Genome length
- Large Single Copy (LSC) length
- Small Single Copy (SSC) length
- Inverted Repeat A (IRa) length
- Inverted Repeat B (IRb) length

---

### 2. Genome composition

- Overall GC content
- GC content of
  - LSC
  - SSC
  - IR regions
- Overall AT content

---

### 3. Coding sequence statistics

- Total CDS length
- CDS percentage of the whole genome
- Number of protein-coding genes
- Number of tRNA genes
- Number of rRNA genes

---

### 4. Codon position composition

The script calculates nucleotide composition for coding sequences at:

- First codon position
- Second codon position
- Third codon position

including

- AT%
- GC%

These statistics are commonly used in codon usage and molecular evolution analyses.

---

### 5. Intron-containing genes

Automatically detects genes containing introns and reports

- Gene name
- Number of exons
- Number of introns
- Length of each exon
- Length of each intron
- Genomic coordinates
- Genome region (LSC / SSC / IRa / IRb)

---

### 6. Gene annotation summary

For every annotated gene, the script reports

- Gene name
- Gene type
- Start position
- End position
- Strand
- Gene length
- Genome region

---

### 7. Region assignment

Each genomic feature is automatically assigned to

- LSC
- SSC
- IRa
- IRb

using the annotated repeat regions in the GenBank file.

---

## Input

Supported formats

```
.gb
.gbk
.gbf
.genbank
```

The script can process

- a single GenBank file
- multiple files
- an entire directory
- directories recursively

---

## Output

Depending on the selected options, the script generates several tab-delimited tables.

Typical outputs include

```
genome_summary.tsv
gene_summary.tsv
intron_table.tsv
gene_coordinates.tsv
```

These tables can be directly imported into

- Microsoft Excel
- LibreOffice Calc
- R
- Python (pandas)
- GraphPad Prism

---

## Installation

Python ≥ 3.9 is recommended.

Install Biopython

```bash
pip install biopython
```

---

## Usage

### Single genome

```bash
python chloroplast_genome_parser.py \
-i genome.gbk \
-o results/
```

### Multiple genomes

```bash
python chloroplast_genome_parser.py \
-i GenBank_folder \
-o results/
```

### Recursive search

```bash
python chloroplast_genome_parser.py \
-i GenBank_folder \
-o results/ \
-r
```

---

## Example directory

```
Project/

├── GenBank/
│   ├── Bidens_alba.gbk
│   ├── Bidens_pilosa.gbk
│   ├── Bidens_aurea.gbk
│
├── results/
│
└── chloroplast_genome_parser.py
```

---

## Coordinate system

The script reports genomic coordinates using the conventional biological coordinate system

- 1-based
- inclusive

instead of Python's 0-based indexing.

---

## Intron detection

Genes containing

```
join(...)
```

locations are automatically recognized as intron-containing genes.

For each gene, the script calculates

- exon coordinates
- exon lengths
- intron coordinates
- intron lengths

without requiring any manual annotation.

---

## Handling duplicated genes

Genes duplicated within the inverted repeat (IR) regions are retained as separate genomic features.

This allows accurate reporting of

- duplicated CDS
- duplicated tRNA
- duplicated rRNA

commonly found in chloroplast genomes.

---

## Typical applications

This script is useful for

- Chloroplast genome characterization
- Comparative plastome analysis
- Genome annotation validation
- Molecular evolution studies
- Phylogenomics
- Chloroplast genome publications

---

## Advantages

Compared with manual extraction, this script

- avoids transcription errors
- provides reproducible results
- processes dozens or hundreds of genomes automatically
- generates publication-ready tables
- greatly reduces analysis time

---

## Citation

If you use this script in your research, please cite your associated publication describing the chloroplast genome analyses.

---

## License

MIT License

---

## Author

Developed for chloroplast genome comparative analyses and annotation statistics.
