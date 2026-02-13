# 

## Téléchargez les fichiers suivants pour les exercices d’aujourd’hui

### Description des fichiers

| nom du fichier | type de fichier |  description du fichier |
| - | - | - | - |
| reference_mt.fasta | FASTA | génome de référence mitochondrial humain |
| mt_multiSNP_2000_R1.fastq | FASTQ | 2000 lectures mitochondriales humaines (lectures avant) |
| mt_multiSNP_2000_R2.fastq | FASTQ | 2000 lectures mitochondriales humaines (lectures arrière)|


1.  Chargez les modules pour **bwa**, **samtools** et **bcftools**
    
2.  Créez tous les fichiers d’index pour le génome de référence
    
3.  Alignez les lectures sur le génome de référence
    
4.  Compressez, triez et indexez le fichier SAM
    
5.  Filtrez le fichier BAM pour retirer les lectures ayant un score d’alignement inférieur à 20
    
6.  Appelez les variants (créez un VCF) à partir du fichier BAM filtré
    
7.  Filtrez le VCF pour retirer les SNP ayant un score de qualité inférieur à 30
    
8.  Quel SNP a été retiré du VCF lors du processus de filtrage à l’étape 7 ?
