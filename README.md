# Chloroplast Genome Analysis

Quy trình phân tích tổng quát bộ gen lục lạp.

![Chloroplast genome analysis workflow](Images/Workflow.svg)

## 1. Data preprocessing

Dữ liệu thô từ máy giải trình tự ở định dạng FASTQ được kiểm tra chất lượng bằng **FastQC v0.12.1**. Sau đó, **MultiQC v1.33** được sử dụng để tổng hợp các báo cáo chất lượng.

**fastp** được sử dụng để cắt bỏ các nucleotide chất lượng thấp ở đầu reads, loại bỏ reads chất lượng thấp và adapter nhằm làm sạch dữ liệu. Kết quả đầu ra là các file FASTQ paired-end đã được làm sạch.

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
###1.2. Generate a summary report using MultiQC

```bash
multiqc \
  <path/to/fastqc_output_directory> \
  -o <path/to/multiqc_output_directory>
```
###1.3. Read filtering and trimming using fastp
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
##2. Genome assembly
Sau đó các file Forward và Reverd sẽ được lắp ráp de novol thành bộ gen hoàng chỉnh bằng công cụ Get organelle (V1.7.7.1) với bộ dữ liệu embplant_pt (dành cho lục lạp). Kết quả được 2 file chính: file fasta hoàn chỉnh và file gfa. Thông thường kết quả fasta sẽ có dạng complete nếu lắp ráp hoàng chỉnh, đây là file kết quả chính được dùng để phân tích chính. Đồng thời file gfa sẽ được dùng để trực quan hóa dữ liệu kết quả lắp ráp bằng công cụ Bandage (V0.9.0).
3. Genome annotation
File fasta sau khi được lắp ráp hoàn chỉnh sẽ được chú giải bằng công cụ trực tuyến Geseq (V2.03) hoặc có thể thay thế bằng công cụ trực tuyến CPGAVAS2.
File kết quả chứa thông tin chú giải bộ gen lục lạp là file genbank sẽ được trực quan hóa thành bản đồ gen bằng công cụ OGDRAW (V1.3.1) 
