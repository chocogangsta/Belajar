In-House Bioinformatics Pipeline

Dikembangkan untuk keperluan pembelajaran analisis bioinformatika dasar

Pipeline ini disusun untuk memberikan pengalaman praktis dalam menganalisis data sekuensing, khususnya shotgun metagenomics dan Oxford Nanopore Technology (ONT).
Cocok digunakan oleh peneliti, mahasiswa, atau teknisi laboratorium yang ingin memahami alur kerja analisis bioinformatika secara menyeluruh.

  ⚠️ Catatan Penting:
      > Pipeline ini hanya ditujukan untuk keperluan pembelajaran dan bukan untuk analisis resmi atau diagnosis klinis.

🔧 Persyaratan Aplikasi
Sebelum menjalankan pipeline, pastikan Anda telah menginstal software berikut:

Aplikasi 

1. Fastp (quality control, trimming, filtering, dll)
   > https://github.com/OpenGene/fastp

2. BWA	mem (Mapping to Reference) for Illumina 
   > https://github.com/lh3/bwa

3. Minimap2 (Mapping to Reference) for Oxford Nanopore Technology / Pacbio
   > https://github.com/lh3/minimap2

4. Samtools	(Manipulate SAM/BAM File from Mapping)
   > https://github.com/samtools/samtools

5. Bcftools	(Varian calling dan manipulating VCF/BCF file) 
   > https://github.com/samtools/bcftools


🚀 Cara Menjalankan Pipeline
       > bash Installation.sh <

File >config.txt< adalah file konfigurasi utama. Anda harus menyesuaikan isinya sesuai dengan:

       1. Lokasi FastQ File        <
       2. Lokasi Reference Genome  <
       3. Lokasi bed File          <

🚀 Cara Menjalankan Pipeline
       1. For Illumina
          > bash Illumina.sh config.txt <
       2. For ONT
          > bash ONT.sh config.txt <

