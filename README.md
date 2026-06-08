# In-House Bioinformatics Pipeline
---

## 🧬 In-House Bioinformatics Pipeline

Pipeline bioinformatika sederhana untuk pembelajaran analisis data sekuensing berbasis **Illumina** dan **Oxford Nanopore Technology (ONT)**.

Dirancang untuk membantu pengguna memahami alur kerja bioinformatika mulai dari quality control, read mapping, pengolahan file BAM, hingga variant calling. Pipeline ini cocok digunakan dalam kegiatan pelatihan, workshop, maupun pembelajaran mandiri di bidang genomik dan metagenomik.

> ⚠️ **Disclaimer**
>
> Pipeline ini hanya ditujukan untuk tujuan pendidikan dan penelitian. Tidak diperuntukkan untuk analisis klinis, diagnostik, maupun pelaporan resmi.

---
  <br>

## 🔧 Persyaratan Aplikasi
Sebelum menjalankan pipeline, pastikan Anda telah menginstal software berikut:

| Software | Fungsi |
|-----------|---------|
| **Fastp** | Quality control, trimming, filtering reads |
| **BWA-MEM** | Mapping reads Illumina ke reference genome |
| **Minimap2** | Mapping reads ONT/PacBio ke reference genome |
| **Samtools** | Manipulasi file SAM/BAM |
| **Bcftools** | Variant calling dan manipulasi file VCF/BCF |
<br>

---
### Installation References

#### Fastp
https://github.com/OpenGene/fastp

#### BWA-MEM
https://github.com/lh3/bwa

#### Minimap2
https://github.com/lh3/minimap2

#### Samtools
https://github.com/samtools/samtools

#### Bcftools
https://github.com/samtools/bcftools

---

## 📂 Pipeline Structure

```text
project/
│
├── Installation.sh
├── Illumina.sh
├── ONT.sh
├── config.txt
│
├── reference/
│   └── genome.fasta
│
├── bed/
│   └── target.bed
│
└── fastq/
    └── sample.fastq.gz
```

---


## 🚀 Running the Pipeline

### Illumina Workflow

Untuk data hasil sekuensing Illumina:

```bash
bash Illumina.sh config.txt
```


### Oxford Nanopore Technology (ONT) Workflow

Untuk data hasil sekuensing ONT:

```bash
bash ONT.sh config.txt
```

---

