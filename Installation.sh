#!/bin/bash

# Unduh installer Miniconda
sudo wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /opt/miniconda-installer.sh
source ~/miniconda3/bin/activate
conda init --all

# Jalankan installer
sudo bash /opt/miniconda-installer.sh

# Sumber .bashrc agar perubahan diterapkan
source ~/.bashrc

# Install BWA
sudo apt install bwa

# Install Samtools
sudo apt install samtools


echo "Instalasi selesai!"
