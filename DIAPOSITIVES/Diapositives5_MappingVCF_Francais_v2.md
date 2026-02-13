# Alignement des reads et appel de variants

Fichiers à télécharger pour le cours d’aujourd’hui

```
wget https://raw.githubusercontent.com/nomascus/ANT3475_2026/refs/heads/main/FILES/Mapping_VCF/reference.fasta
wget https://raw.githubusercontent.com/nomascus/ANT3475_2026/refs/heads/main/FILES/Mapping_VCF/reads_sample1_F.fastq
wget https://raw.githubusercontent.com/nomascus/ANT3475_2026/refs/heads/main/FILES/Mapping_VCF/reads_sample1_R.fastq
```

Imaginez que vous venez de télécharger des données de séquençage d’ADN provenant de 200 échantillons humains et que votre objectif est d’identifier comment les populations sont subdivisées selon la géographie. Malheureusement, on ne peut pas commencer directement à travailler sur les questions biologiques, parce que tout ce que la machine de séquençage nous donne, c’est un énorme mélange de **courts fragments d’ADN, appelés "reads" (lectures de séquençage))**, dans un type de fichier appelé FASTQ (notez que c’est différent du FASTA). 

Afin d'analyser les données, nous devons générer un fichier VCF qui ne contient que les positions nucléotidiques variables dans notre échantillon. En d'autres termes, nous devons passer des lectures de séquençage brutes aux SNP. C'est l'objectif d'aujourd'hui.

### Que sont les fichiers FASTQ ?

Les fichiers FASTQ contiennent quatre lignes de texte pour chaque lecture :

Exemple :

```
@read1_Forward
ACTGATCGTACGATCGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
+
IIIIIIIIIIHHHHHHHHHHFFFFFFFFFFDDDDDDDDDD!!!!!!!!!!
@read2_Forward
GCTAGCTAGCATCGATCGATCGTAGCTAGCATCGATCGTAGCTAGCATCG
+
IIIIIIHHHHHHFFFFFFDDDDDD!!!!!!!!!IIIIHHHHFFFFDDDDD


```

Ce sont les données brutes de la machine de séquençage. Elles contiennent :

1.  Un en-tête qui commence avec de l’information sur le read qu’on regarde. Les en-têtes FASTQ commencent toujours par le caractère @.
    
2.  Une ligne de séquence nucléotidique
    
3.  Un signe plus
    
4.  Une ligne de caractères correspondant au score de qualité pour chaque nucléotide. Autrement dit, à quel point la machine de séquençage est confiante que le premier nucléotide du read est A, le deuxième est C, etc. Ce score de qualité nous permet de retirer les reads peu fiables de notre analyse.
    

_Remarque_ : Le séquençage de l’ADN se fait habituellement dans les directions avant et arrière le long du brin d’ADN ; donc, on obtient souvent deux fichiers FASTQ par échantillon, un fichier avant et un fichier arrière.

Avant de pouvoir commencer à étudier les questions biologiques, nous devons :

1.  Déterminer où les lectures se placent dans le génome
    
2.  Nettoyer les mauvaises données
    
3.  Trouver les différences par rapport au génome de référence pour chaque individu
    

# 1. Alignement (FASTQ → SAM)

Avant de faire quoi que ce soit, nous devons charger trois programmes qui sont stockés sur le serveur. Pour faire ça, on utilise la commande **module load**. Si ces programmes ne fonctionnent pas plus tard, vérifiez d’abord que vous les avez bien chargés !

```
module load samtools
module load bwa
module load bcftools
```

Le processus d’alignement prend un fichier FASTQ et assigne chaque read à la position la plus probable dans un génome de référence. Avec assez de reads, on peut couvrir tout le génome.

Imaginez que vous avez une copie déchiquetée d’un document. C’est un peu comme prendre chaque morceau de papier déchiqueté et l’aligner avec une autre copie du document pour trouver où il appartient.

![Hands-on: Mapping / Mapping / Sequence analysis](https://training.galaxyproject.org/training-material/topics/sequence-analysis/images/mapping/mapping.png)

La différence ici est qu’on utilise un génome de référence, qui est simplement un fichier FASTA avec tous les nucléotides alignés sur tous les chromosomes de l’organisme qu’on veut étudier. Par exemple, le génome humain fait environ 3 milliards de nucléotides, et on doit trouver où chaque read s’insère.

Comme un fichier de génome de référence peut être très gros, on doit créer un index qui aide les autres programmes à retrouver les reads alignés plus rapidement. C’est comparable à l’index à la fin d’un livre de recettes qui vous aide à trouver la page d’une recette.

```
samtools faidx reference.fasta 
bwa index reference.fasta 
```

Ici, on utilise les commandes **samtools faidx** et **bwa index** pour créer les fichiers d’index. Les fichiers d’index seront créés automatiquement et vont se terminer par .amb, .ann, .bwt, .fai, .pac et .sa.

Pour aligner les reads, on utilise un programme appelé **bwa mem**.

```
bwa mem -t 1 reference.fasta reads_sample1_F.fastq reads_sample1_R.fastq > mappedReads_sample1.sam
```

On lance bwa mem, puis on lui donne trois fichiers :

1.  Le génome de référence
    
2.  Les lectures avant (souvent compressés en .gz)
    
3.  Les lectures arrière (souvent compressés en .gz)
    

et ensuite on sauvegarde le résultat (on redirige la sortie de bwa mem) dans un nouveau fichier, appelé un fichier SAM.

L’option **-t 1** indique qu’on veut utiliser un seul processeur sur l’ordinateur. Pour un très gros fichier, il faudrait augmenter ce nombre, mais ici c’est un petit exemple.

### Qu’est-ce qu’un fichier SAM ?

SAM (Sequence Alignment/Map) est un fichier texte qui contient :

-   Où chaque read s’aligne dans le génome
    
-   La qualité d’alignement (MAPQ)
    
-   Des détails sur l’alignement
    

Voici un exemple de fichier SAM.

```
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:1000

read1	0	chr5	1540	60	50M	*	0	0	ACTGATCGTACGATCGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII!!!!!!!!!#####
read2	0	chr1	292	15	50M	*	0	0	GCTAGCTAGCATCGATCGATCGTAGCTAGCATCGATCGTAGCTAGCATCG	#####!!!!!!!!!IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
```

On n’a pas besoin d’entrer dans les détails, mais essentiellement ça veut dire :

Le lecture 1 devrait commencer à la position 1540 sur le chromosome 5, et le programme est confiant qu’il appartient à cet endroit (MAPQ = 60),

Le lecture 2 devrait commencer à la position 292 sur le chromosome 1, et le programme n’est pas très confiant qu’il appartient à cet endroit (MAPQ = 15),

Dans un vrai fichier SAM, il y aura des millions de lectures qui couvrent chaque position nucléotidique plusieurs fois.

# 2. Compresser, trier et indexer le fichier SAM

Maintenant qu’on a notre fichier d’alignement, on doit l’organiser et le compresser, parce que beaucoup de programmes “en aval” demandent que les lectures soient ordonnés par chromosome et par position nucléotidique. Le fichier SAM qu’on a est énorme et les lectures sont dans un ordre aléatoire. Pour faire ça, on utilise un ensemble de programmes appelé **samtools**.

```
samtools view -bS mappedReads_sample1.sam | samtools sort -o mappedReads_sample1.sorted.bam
```

Décomposons ça :

1.  **samtools view** ouvre le fichier SAM. **-bS** indique que l’entrée est un SAM et que la sortie sera une version compressée appelée un fichier **BAM**.
    
2.  On “pipe” la sortie vers **samtools sort**, qui trie le fichier BAM par chromosome et par position nucléotidique.
    
3.  **-o** indique que le nom du fichier de sortie sera mappedReads_sample1.sorted.bam.
    

Il reste une petite étape de préparation. Comme le fichier BAM peut être très gros, on doit créer un fichier d’index qui aide les autres programmes à retrouver les lectures alignés plus rapidement.

```
samtools index mappedReads_sample1.sorted.bam
```

Ici, on utilise la commande **samtools index** pour créer le fichier d’index. Le fichier d’index sera créé automatiquement et va se terminer par .bai

```
mappedReads_sample1.sorted.bam.bai 
```

# 3. Filtrer le fichier BAM

Rappelez-vous qu’il y a un score de qualité d’alignement (MAPQ) pour chaque lecture dans notre fichier BAM. L’ordinateur fait de son mieux pour trouver une position dans le génome pour chaque read, mais parfois il aligne un lecutre au mauvais endroit, ou à plusieurs endroits du génome. Cela peut arriver pour plusieurs raisons : ADN bactérien contaminant dans l’échantillon, ADN dégradé, génome de référence de mauvaise qualité, ou erreurs de séquençage. On doit retirer les lectures mal alignés, pour ne garder que des séquences qui viennent du bon organisme et qui sont alignées au bon endroit.

```
samtools view -b -q 20 mappedReads_sample1.sorted.bam > mappedReads_sample1.sorted.filtered.bam
```

Ici encore, on utilise **samtools view** pour ouvrir le fichier BAM avec l’option -b, puis on ajoute **-q 20**, ce qui retire les lecture avec un score MAPQ inférieur à 20. Cette étape améliore beaucoup la fiabilité des données. Ensuite, on sauvegarde le résultat dans un nouveau fichier avec le symbole de redirection **>**.

# 4. Appel de variants (BAM → VCF)

Maintenant qu’on a un fichier BAM propre, c’est le moment de créer le **VCF**. Rappelez-vous que c’est le **Variant Call File**, qui nous donne une liste de toutes les positions variables dans le génome, ainsi que le génotype de chaque individu dans notre échantillon.

### Appeler des variants avec bcftools

```
bcftools mpileup -Ou -f reference.fasta mappedReads_sample1.sorted.filtered.bam |
bcftools call -mv -Ov -o HumanSamples.vcf
```

**bcftools** est une autre suite de programmes pour travailler avec des données génomiques. On va utiliser deux de ses programmes.

1.  `mpileup` regarde chaque position du génome et vérifie combien de lectures couvrent cette position. Par exemple :
    

![](https://github.com/nomascus/ANT3475_2026/blob/main/FILES/mpileup.png?raw=true)

**-Ou** veut dire que le type de sortie est un VCF non compressé

**-f** spécifie le fichier du génome de référence

2.  `call` identifie les positions où l’échantillon est différent du génome de référence.
    

![](https://github.com/nomascus/ANT3475_2026/blob/main/FILES/mpileupCall2.png?raw=true)

Ici, vous pouvez voir qu’il y a deux SNP potentiels : un A/G possible à la position 10 et un A/C possible à la position 18.

Avec **bcftools call**, on prend la sortie de **bcftools mpileup** du pipe, on garde seulement les positions variables, et on convertit le format en VCF.

Voici ce que veulent dire les autres options de bcftools call :

**-mv** utiliser l’appel de variants standard et ne garder que les positions variables  
**-Ov** le format de sortie est VCF  
**-o** le nom du fichier de sortie

Le résultat est un **fichier VCF** qui liste seulement les positions variables qui diffèrent du génome de référence chez un ou plusieurs individus. Voici un exemple.

```
##fileformat=VCFv4.2
##source=demo_pipeline
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1 sample2
chr1	105	.	A	G	50	PASS	.	GT	0/1	0/1
chr1	220	.	C	T	12	LowQual	.	GT	0/0	0/1
chr1	350	.	G	A	99	PASS	.	GT	1/1 1/0
```

Comme on l’a vu, un VCF contient les lignes suivantes :

Lignes d’en-tête (métadonnées) qui commencent par ## :

```
##fileformat=VCFv4.2
##source=demo_pipeline
```

Une ligne d’en-tête des colonnes qui commence par # :  
`#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT sample1 sample2`

et des lignes pour chaque SNP :

```
chr1	105	.	A	G	50	PASS	.	GT	0/1	0/1
chr1	220	.	C	T	12	LowQual	.	GT	0/0	0/1
chr1	350	.	G	A	99	PASS	.	GT	1/1 1/0
```

# 5. Filtrer le VCF

Le VCF qu’on a créé inclut tous les SNP potentiels, même ceux pour lesquels le programme n’est pas très confiant. Par exemple, si une position de référence est A et qu’elle est couverte par 9 lectures (7 A et 2 T), c’est difficile d’être confiant que c’est un vrai SNP. Dans ce cas, la confiance dans l’appel du SNP peut être faible.

Par prudence, on devrait retirer ces SNP de faible confiance.

```
bcftools filter -i 'QUAL>30' HumanSamples.vcf -o HumanSamples.filtered.vcf
```

Ici, on utilise un autre programme de bcftools appelé **bcftools filter** pour garder seulement les SNP avec un score QUAL d’au moins 30.

L'option **-i** signifie « information » et nous permet d'écrire une chaîne complexe de paramètres de filtrage. Pour l'instant, nous allons simplement lui indiquer les scores de qualité, mais n'oubliez pas d'inclure les guillemets.

----------

## Résumé des types de fichiers et des programmes

| Type de fichier| Description | Programme à utiliser|
| - | - | - |
| FASTQ | Reads de séquençage| bwa mem |
| SAM | Fichier d’alignement | samtools view, samtools sort |
| BAM | Fichier d’alignement compressé | samtools filter, bcftools mpileup, bcftools call |
| BAI | Fichier d’index BAM | samtools index |
| .amb, .ann, .bwt, .fai, .pac, and .sa | Fichiers d’index du génome de référence | samtools index, bwa index|
| VCF | Fichier d’appel de variants| bcftools filter|
