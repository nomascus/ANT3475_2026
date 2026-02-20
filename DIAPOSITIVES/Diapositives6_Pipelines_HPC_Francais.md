
# ANT3475 Cours 6 : Pipelines et HPC

Voici les fichiers à télécharger pour le cours d’aujourd’hui. Veuillez commencer par créer un nouveau répertoire comme suit.

```
mkdir ~/scratch/6_HPC
cd ~/scratch/6_HPC

wget https://raw.githubusercontent.com/nomascus/ANT3475_2026/refs/heads/main/FILES/6_HPC/reference_mt.fasta
wget https://raw.githubusercontent.com/nomascus/ANT3475_2026/refs/heads/main/FILES/6_HPC/mt_multiSNP_2000_R1.fastq 
wget https://raw.githubusercontent.com/nomascus/ANT3475_2026/refs/heads/main/FILES/6_HPC/mt_multiSNP_2000_R2.fastq
wget https://raw.githubusercontent.com/nomascus/ANT3475_2026/refs/heads/main/FILES/6_HPC/reads_sample1_F.fastq
wget https://raw.githubusercontent.com/nomascus/ANT3475_2026/refs/heads/main/FILES/6_HPC/reads_sample1_R.fastq
wget https://raw.githubusercontent.com/nomascus/ANT3475_2026/refs/heads/main/FILES/6_HPC/metadataHPC.txt
wget https://raw.githubusercontent.com/nomascus/ANT3475_2026/refs/heads/main/FILES/6_HPC/metadata_batch.txt
wget https://raw.githubusercontent.com/nomascus/ANT3475_2026/refs/heads/main/FILES/6_HPC/mtDNA_mappingVCF.sh
wget https://raw.githubusercontent.com/nomascus/ANT3475_2026/refs/heads/main/FILES/6_HPC/mtDNA_mappingVCF2.sh
wget https://raw.githubusercontent.com/nomascus/ANT3475_2026/refs/heads/main/FILES/6_HPC/mtDNA_mappingVCF2.sbatch
```

----------

# Tout assembler dans un pipeline

Nous avons fait beaucoup de bioinformatique la semaine dernière lorsque nous avons généré ces fichiers VCF ! Cependant, il n’est pas très efficace de taper toutes ces lignes de texte une par une, surtout si vous voulez exécuter cela sur des dizaines ou des centaines d’échantillons. C’est pourquoi nous créons des scripts appelés **pipelines** qui enchaînent tous ces programmes afin que nous puissions exécuter l’ensemble du processus en une seule ligne de commande.

Un pipeline est un long fichier texte que nous créons dans vi et qui se termine par **.sh** si nous l’exécutons en ligne de commande ou **.sbatch** si nous l’exécutons à distance (nous y reviendrons plus tard).

Voici un exemple de pipeline pour tout ce que nous avons fait la semaine dernière.

`vi mtDNA_mappingVCF.sh`

```
#!/bin/bash

############################
# 0. Load modules
############################

module load bwa
module load samtools
module load bcftools

############################
# 1. Define input files.
############################

REF=reference_mt.fasta
R1=reads_sample1_F.fastq
R2=reads_sample1_R.fastq
SAMPLE='Sample1'

############################
# 2. Index reference files 
############################

echo "Indexing reference..."

samtools faidx $REF
bwa index $REF

############################
# 3. Align reads to reference
############################

echo "Aligning reads..."

bwa mem $REF $R1 $R2 > ${SAMPLE}.sam

############################
# 4. Convert SAM → sorted BAM
############################

echo "Converting to BAM..."

samtools view -bS ${SAMPLE}.sam | \
samtools sort -o ${SAMPLE}.sorted.bam

############################
# 5. Index BAM
############################

samtools index ${SAMPLE}.sorted.bam

############################
# 6. Basic BAM filtering
# (remove low mapping quality reads)
############################

echo "Filtering BAM..."

samtools view -b -q 20 ${SAMPLE}.sorted.bam > ${SAMPLE}.filtered.bam
samtools index ${SAMPLE}.filtered.bam

############################
# 7. Call variants
############################

echo "Calling variants..."

bcftools mpileup -Ou -f $REF ${SAMPLE}.filtered.bam | \
bcftools call -mv -Ov -o ${SAMPLE}.vcf

############################
# 8. Filter VCF
# (keep QUAL > 30)
############################

echo "Filtering VCF..."

bcftools filter -i 'QUAL>30' ${SAMPLE}.vcf -o ${SAMPLE}.filtered.vcf

echo "Pipeline complete!"
```

----------

1.  La première ligne est toujours **#!/bin/bash**. Il s’agit d’une ligne spéciale qui indique à l’ordinateur qu’il s’agit d’un programme bash. Ne l’oubliez pas !
    
2.  Pour toutes les autres lignes, un **#** signifie que la ligne est un commentaire. Les lignes de commentaire **ne sont pas exécutées** par l’ordinateur (il existe quelques exceptions comme #! et #SBATCH). Vous pouvez les utiliser pour prendre des notes, rendre votre code plus lisible et le diviser en blocs visibles.
    
3.  Dans le bloc 0, nous chargeons les trois programmes : samtools, bwa et bcftools.
    
4.  Dans le bloc 1, nous indiquons au programme quels fichiers d’entrée utiliser. Ici, nous déclarons quatre variables : une pour le genome de référence, une pour les lectures avant, une pour les lectures arrière et une pour le nom de l’échantillon. Cela est utile, car si nous voulons exécuter le pipeline sur un autre ensemble d’échantillons, de références ou de lectures, nous n’avons qu’à modifier les noms de ces fichiers une seule fois.
  
5.  Dans les blocs restants, nous appelons simplement les mêmes commandes que la semaine dernière pour aligner les reads et traiter le VCF, sauf que nous utilisons maintenant les noms de variables déclarés ci-dessus. Nous avons également ajouté quelques lignes `echo` pour afficher de l’information à l’écran afin de voir quelle étape est en cours d’exécution.
    

Pour exécuter tout le pipeline qui aligne les reads et appelle les VCF pour cet ensemble de données de séquençage, il suffit de taper :

`bash mtDNA_mappingVCF.sh`

et tout le pipeline s’exécute !

----------

## Partie HPC 2

Vous vous souvenez lorsque nous avons discuté des variables environnementales comme $1, $2, $3 ? Je vous avais dit de ne pas les utiliser, car elles ont une fonction spéciale. Si nous définissons une variable à $1 dans le pipeline, elle prend le premier élément tapé après le programme et le stocke dans $1. De même, si nous définissons une variable à $2, elle prend le deuxième élément tapé après le programme et le stocke dans $2, et ainsi de suite.

Que se passe-t-il si nous modifions le bloc 1 de notre pipeline comme suit :

`vi mtDNA_mappingVCF2.sh`

```
############################
# 1. Define input files.
############################

REF=$1 
R1=$2
R2=$3
SAMPLE=$4

```

Lorsque nous appelons le programme `bash mtDNA_mappingVCF2.sh`, il génère plusieurs erreurs, car il s’attend à recevoir quatre informations supplémentaires (dans ce cas la référence, deux fichiers fastq et un nom d’échantillon) dont il a besoin pour s’exécuter.

Donc, si nous exécutons :

```
bash mtDNA_mappingVCF2.sh reference_mt.fasta reads_sample1_F.fastq reads_sample1_R.fastq Sample1
```

bash stockera ces quatre chaînes de texte dans $1, $2, $3 et $4 pour utilisation dans notre script. Pour rendre cela plus clair, nous pouvons ensuite stocker ces variables dans des variables plus explicites comme suit :
```
REF=$1 
R1=$2
R2=$3
SAMPLE=$4
```

| Environmental variable |Stored information |Variable name in pipeline
| - | - | - | - |
| $1 |reference_mt.fasta | $REF|
| $2 |reads_sample1_F.fastq | $R1|
| $3 |reads_sample1_R.fastq | $R2|
| $4 |Sample1 | $SAMPLE|

Cela donne beaucoup de flexibilité au pipeline, car nous pouvons simplement modifier l’information sur la ligne de commande et exécuter le même pipeline ! Cela devient très utile si nous voulons exécuter le pipeline sur plusieurs échantillons ou ensembles de reads en même temps. Nous pouvons faire cela avec une boucle **while** et un fichier de métadonnées.

```
[orkin@login1 6_Unix]$ cat metadataHPC.txt
reference_mt.fasta  reads_sample1_F.fastq  reads_sample1_R.fastq  Sample1
reference_mt.fasta  mt_multiSNP_2000_R1.fastq  mt_multiSNP_2000_R2.fastq  Sample1
```

```
while read line; do bash mtDNA_mappingVCF2.sh $line; done <metadataHPC.txt
```

Et voilà ! Avec cette seule ligne de code, nous sommes passés de fastq à VCF filtré pour deux échantillons !

----------

## Qu’est-ce qu’un HPC

Jusqu’à présent, tout notre travail s’est déroulé sur l’ordinateur infonuagique en ligne de commande. Bien que l’interface en ligne de commande que nous apprenons soit différente de l’interface graphique à laquelle vous êtes habitués sous Windows ou Mac, vous vous demandez probablement pourquoi tout cela est nécessaire. L’une des principales raisons d’utiliser la ligne de commande est que nous pouvons exécuter de nombreuses tâches en même temps. Contrairement à un ordinateur personnel, nous pouvons planifier des dizaines, des centaines, voire des milliers de tâches simultanément avec seulement quelques lignes de code.

Un cluster de calcul haute performance (HPC) est un ensemble de nombreux ordinateurs connectés entre eux et gérés de façon à exécuter efficacement de gros calculs. Au lieu d’exécuter les analyses sur votre ordinateur portable, vous soumettez des tâches au cluster, qui les exécute sur des machines distantes puissantes.

Un cluster typique comporte trois parties :

1.  **Le nœud de connexion (login node)**. C’est là où nous avons travaillé jusqu’à présent. C’est l’endroit où vous vous connectez via ssh, modifiez des scripts, déplacez des fichiers et soumettez des tâches avec `sbatch`. Nous ne l’utilisons pas pour des calculs lourds, car c’est un espace partagé.
    
2.  **Les nœuds de calcul (compute nodes)**. Ce sont les ordinateurs puissants auxquels nous envoyons nos scripts pour exécuter notre code. Vous ne vous connectez jamais directement à ces nœuds.
    
3.  **Le planificateur de tâches (SLURM job scheduler)**. SLURM est un programme qui gère toutes les tâches envoyées aux nœuds de calcul par les utilisateurs.     

Pensez à un restaurant très occupé. Le rôle de l’hôte (SLURM) est d’asseoir une longue file de clients (les tâches en attente dans la file) à des tables vides de différentes tailles (les nœuds de calcul) afin qu’ils puissent manger.

----------

## Envoyer une tâche (job) au cluster

Pour envoyer notre tâche aux nœuds de calcul, nous devons ajouter les lignes suivantes au début de notre script (juste sous #!/bin/bash) :

```
#SBATCH --job-name=Joe_vcf_pipeline
#SBATCH --cpus-per-task=1  
#SBATCH --mem=1G  
#SBATCH --time=00:10:00  
#SBATCH --output=%x_%j.out  
#SBATCH --error=%x_%j.err
```
Analysons cela. Les lignes **#SBATCH** indiquent au cluster les ressources dont votre tâche a besoin avant son exécution.

`#SBATCH --job-name=bwa_pipeline`  
Cela donne à votre tâche un nom lisible par un humain, ce qui rend les fichiers de sortie plus faciles à interpréter. Vous pouvez le changer pour ce que vous voulez.

`#SBATCH --cpus-per-task=1`  
Demande le nombre de cœurs CPU pour la tâche. Plus il y a de CPU, plus la tâche s’exécute rapidement.

`#SBATCH --mem=1G`  
Demande la quantité de RAM pour la tâche. Trop peu de RAM et la tâche plantera, mais trop de RAM et elle restera en attente dans la file (queue) avant de démarrer.

`#SBATCH --time=00:10:00`  
Définit le temps maximal demandé pour l’exécution de la tâche. Les demandes de temps plus longues restent plus longtemps en attente dans la file, mais si votre tâche n’est pas terminée avant le temps demandé, elle sera arrêtée. Le format est heures:minutes:secondes.

`#SBATCH --output=%x_%j.out`  
Spécifie un nom de fichier pour la sortie standard (ce qui serait normalement affiché à l’écran). %x est un code pour le nom de la tâche, et %j pour l’identifiant de la tâche (job ID), un numéro unique que l’ordinateur assigne à chaque tâche. Par exemple, cela peut créer un fichier appelé Joe_bwa_pipeline_12.out.

`#SBATCH --error=%x_%j.err`  
Spécifie un nom de fichier pour la sortie d’erreur standard (les messages d’erreur normalement affichés à l’écran). Cela peut créer un fichier appelé Joe_bwa_pipeline_12.err.
 
**Dans ce pipeline :** cela ressemble à ceci :

```vi  mtDNA_mappingVCF2.sbatch```
```
#!/bin/bash  
#SBATCH --job-name=Joe_vcf_pipeline  
#SBATCH --cpus-per-task=1  
#SBATCH --mem=1G  
#SBATCH --time=00:10:00  
#SBATCH --output=%x_%j.out  
#SBATCH --error=%x_%j.err

############################
# 0. Load modules
############################

module load bwa
module load samtools
module load bcftools

############################
# 1. Define input files.
############################

REF=$1
R1=$2
R2=$3
SAMPLE=$4

############################
# 2. Index reference files # don't need to run
############################

# echo "Indexing reference..."

# samtools faidx $REF 
# bwa index $REF

############################
# 3. Align reads to reference
############################

echo "Aligning reads..."

bwa mem $REF $R1 $R2 > ${SAMPLE}.sam

############################
# 4. Convert SAM → sorted BAM
############################

echo "Converting to BAM..."

samtools view -bS ${SAMPLE}.sam | \
samtools sort -o ${SAMPLE}.sorted.bam

############################
# 5. Index BAM
############################

samtools index ${SAMPLE}.sorted.bam

############################
# 6. Basic BAM filtering
# (remove low mapping quality reads)
############################

echo "Filtering BAM..."

samtools view -b -q 20 ${SAMPLE}.sorted.bam > ${SAMPLE}.filtered.bam

samtools index ${SAMPLE}.filtered.bam

############################
# 7. Call variants
############################

echo "Calling variants..."

bcftools mpileup -Ou -f $REF ${SAMPLE}.filtered.bam | \
bcftools call -mv -Ov -o ${SAMPLE}.vcf

############################
# 8. Filter VCF
# (keep QUAL > 30)
############################

echo "Filtering VCF..."

bcftools filter -i 'QUAL>30' ${SAMPLE}.vcf -o ${SAMPLE}.filtered.vcf

echo "Pipeline complete!"
```

----------

## Exécution sur le cluster

Pour envoyer la tâche aux nœuds de calcul, nous utilisons la commande **sbatch** au lieu de bash, et nous pouvons enregistrer le fichier du pipeline avec une extension différente afin de pouvoir l’identifier. Une fois la tâche soumise, le cluster nous indique :

```
[orkin@login1 TMP]$ sbatch mtDNA_mappingVCF2.sbatch
Submitted batch job 23
```

Nous pouvons ensuite vérifier l’état de la tâche avec la commande **squeue**.

```
[orkin@login1 TMP]$ squeue
JOBID USER  ACCOUNT NAME  ST  TIME_LEFT NODES CPUS TRES_PER_N MIN_MEM NODELIST (REASON)
23  orkin def-sponsor0 bwa_pipeline  R 0:58 1  1  N/A  1G nodecpupool1 (None)
```

Cela nous donne des informations sur chaque tâche  dans la file d’attente. Les éléments les plus importants sont : qui exécute la tâche (USER), le nom de la tâche (NAME), son statut (ST) et le temps restant (TIME_LEFT). Notez les statuts suivants : **PD** = en attente dans la file ; **CF** = en cours d’initialisation ; **R** = en cours d’exécution.

Si le cluster est resté inactif pendant un certain temps, vous pourriez voir votre tâche bloquée en état CF pendant quelques minutes. Une fois la première tâche lancée, les suivantes s’exécuteront rapidement. Si vos tâches restent bloquées et que vous souhaitez toutes les retirer de la file, vous pouvez taper **scancel** suivi de votre identifiant utilisateur. Par exemple : `scancel orkin`

Notez également que nous pouvons exécuter le pipeline sur un très grand nombre d’échantillons en même temps en utilisant une boucle **while** avec un fichier de métadonnées.

Chaque ligne du fichier de métadonnées correspond à une liste séparée par des tabulations contenant les éléments nécessaires ($1 $2 $3 $4) pour exécuter le pipeline.

```
[orkin@login1 TMP]$ cat metadata_batch.txt
reference_mt.fasta reads_sample1_F.fastq  reads_sample1_R.fastq  Sample1
reference_mt.fasta reads_sample1_F.fastq  reads_sample1_R.fastq  Sample2
reference_mt.fasta reads_sample1_F.fastq  reads_sample1_R.fastq  Sample3
reference_mt.fasta reads_sample1_F.fastq  reads_sample1_R.fastq  Sample4
reference_mt.fasta  mt_multiSNP_2000_R1.fastq  mt_multiSNP_2000_R2.fastq  Sample5
reference_mt.fasta  mt_multiSNP_2000_R1.fastq  mt_multiSNP_2000_R2.fastq  Sample6
reference_mt.fasta  mt_multiSNP_2000_R1.fastq  mt_multiSNP_2000_R2.fastq  Sample7
reference_mt.fasta  mt_multiSNP_2000_R1.fastq  mt_multiSNP_2000_R2.fastq  Sample8
```

```
while read line; do sbatch mtDNA_mappingVCF2.sbatch $line; done < metadata_batch.txt
```


Et voilà, nous avons exécuté le pipeline 8 fois sur différents échantillons ! 

**Remarque** : il ne s'agit pas de 8 ensembles de lectures différents, car il s'agit uniquement d'un exemple. Dans la réalité, les paires de fichiers fastq seraient différentes pour chaque échantillon et auraient leurs propres noms uniques.
