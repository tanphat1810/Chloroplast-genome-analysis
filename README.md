# Chloroplast Genome Analysis

Quy trình phân tích tổng quát bộ gen lục lạp.

![Chloroplast genome analysis workflow](Images/Workflow.svg)

## 1. Data preprocessing
Raw paired-end FASTQ files are assessed for sequencing quality using **FastQC v0.12.1**, and quality reports are summarized with **MultiQC v1.33**. The reads are then processed using **Fastp** to trim low-quality bases, remove adapter sequences, and filter low-quality reads, producing high-quality paired-end FASTQ files for subsequent analyses.
### 1.1. Quality control using FastQC

```bash
fastqc \
  <path/to/input.fastq.gz> \
  -o <path/to/fastqc_output_directory> \
  -t 16
```
Với dữ liệu pair_end
```bash
fastqc \
  <path/to/input_R1.fastq.gz> \
  <path/to/input_R2.fastq.gz> \
  -o <path/to/fastqc_output_directory> \
  -t 16
```
### 1.2. Generate a summary report using MultiQC

```bash
multiqc \
  <path/to/fastqc_output_directory> \
  -o <path/to/multiqc_output_directory>
```
### 1.3. Read filtering and trimming using fastp
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
## 2. Genome assembly
## 2. Chloroplast genome assembly

Các reads paired-end đã được làm sạch, bao gồm file forward và reverse, được sử dụng để lắp ráp de novo bộ gen lục lạp bằng **GetOrganelle v1.7.7.1**. Quá trình lắp ráp được thực hiện với cơ sở dữ liệu `embplant_pt`, được thiết kế cho bộ gen lục lạp của thực vật.
Sau quá trình lắp ráp, GetOrganelle tạo ra các file kết quả, trong đó hai định dạng chính được sử dụng trong quy trình này là file trình tự **FASTA** và file đồ thị lắp ráp **GFA**.
Khi quá trình lắp ráp hoàn chỉnh, file FASTA có tên chứa `complete` thường đại diện cho trình tự bộ gen lục lạp đã được lắp ráp hoàn chỉnh và được sử dụng làm dữ liệu đầu vào cho các bước phân tích tiếp theo. File GFA được sử dụng để trực quan hóa và kiểm tra cấu trúc đồ thị lắp ráp bằng **Bandage v0.9.0**.
### Khởi động môi trường Conda
```bash
conda activate getorg
```
### Lắp ráp de novo bộ gen lục lạp bằng GetOrganelle
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
3. Genome annotation
File fasta sau khi được lắp ráp hoàn chỉnh sẽ được chú giải bằng công cụ trực tuyến Geseq (V2.03) hoặc có thể thay thế bằng công cụ trực tuyến CPGAVAS2.
File kết quả chứa thông tin chú giải bộ gen lục lạp là file genbank sẽ được trực quan hóa thành bản đồ gen bằng công cụ OGDRAW (V1.3.1) 
