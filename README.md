# In-House Bioinformatics Pipeline


## 🧬 In-House Bioinformatics Pipeline

Pipeline bioinformatika sederhana untuk pembelajaran analisis data sekuensing berbasis **Illumina** dan **Oxford Nanopore Technology (ONT)**.

Dirancang untuk membantu pengguna memahami alur kerja bioinformatika mulai dari quality control, read mapping, pengolahan file BAM, hingga variant calling. Pipeline ini cocok digunakan dalam kegiatan pelatihan, workshop, maupun pembelajaran mandiri di bidang genomik dan metagenomik.

> ⚠️ **Disclaimer**
>
> Pipeline ini hanya ditujukan untuk tujuan pendidikan dan penelitian. Tidak diperuntukkan untuk analisis klinis, diagnostik, maupun pelaporan resmi.


  <br>

---

## 🔧 Tools Reference
Sebelum menjalankan pipeline, pastikan Anda telah menginstal software berikut:

| Software | Fungsi | Link |
|-----------|---------|---------|
| **Fastp** | Quality control, trimming, filtering reads | https://github.com/OpenGene/fastp |
| **BWA-MEM** | Mapping reads Illumina ke reference genome | https://github.com/lh3/bwa |
| **Minimap2** | Mapping reads ONT/PacBio ke reference genome | https://github.com/lh3/minimap2 |
| **Samtools** | Manipulasi file SAM/BAM | https://github.com/samtools/samtools |
| **Bcftools** | Variant calling dan consensus untuk ONT file | https://github.com/samtools/bcftools |
| **iVar** | Membuat consensus untuk Illumina file | https://github.com/andersen-lab/ivar |

<br>


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

