# Boucles Bash

Pour accéder à la leçon d’aujourd’hui, téléchargez le fichier html appelé 4_Unix4_2026.html depuis StudiUM. Il devrait s’ouvrir avec votre navigateur web.

Commençons par créer un nouveau dossier dans votre cluster pour la leçon de cette semaine. Après vous être connecté, nous commencerons par votre dossier personnel. Descendez ensuite d’un répertoire jusqu’au dossier scratch. Créez ensuite un nouveau dossier appelé 4_Unix4, et déplacez-vous dans ce dossier.

```
cd ~/scratch
mkdir 4_Unix4 
cd 4_Unix4
```
Jusqu'à présent, tout le code que nous avons utilisé a fonctionné ligne par ligne. Vous tapez quelque chose dans le terminal et cela fonctionne sur un seul fichier ou effectue une seule opération. Une des choses les plus puissantes à propos d'un langage de programmation comme bash est qu'il vous permet d'écrire des commandes qui agissent sur plusieurs fichiers l'un après l'autre sans entrée supplémentaire. Ce type de commande est appelé une **boucle**. Il existe quelques types de boucles différentes qui peuvent être écrites en bash, mais nous allons nous concentrer sur deux d'entre elles, la boucle **for** et la boucle **while**.

## Boucles for

Disons que nous avions une liste de noms : Charlotte, Joe, Jimmy, Rose, Sarah et que nous voulions les afficher à l'écran avec le mot "Bonjour" devant chaque nom sur une ligne séparée. Nous pourrions le faire en tapant les commandes suivantes une par une comme suit :

```
orkin@ANT3475:~$ echo "Bonjour Charlotte"
Bonjour Charlotte
orkin@ANT3475:~$ echo "Bonjour Joe"
Bonjour Joe
orkin@ANT3475:~$ echo "Bonjour Jimmy"
Bonjour Jimmy
orkin@ANT3475:~$ echo "Bonjour Rose"
Bonjour Rose
orkin@ANT3475:~$ echo "Bonjour Sarah"
Bonjour Sarah
```

C'est assez lent et inefficace cependant. Nous pourrions améliorer les choses un peu en ajoutant des points-virgules et des caractères de nouvelle ligne entre chaque ligne. Un point-virgule indique à bash que la commande est terminée. Cela vous permet d'écrire plusieurs commandes sur une seule ligne à exécuter par bash l'une après l'autre (commandes sérielles). Par exemple :

```
orkin@ANT3475:~$ echo "Bonjour Charlotte"; echo "Bonjour Joe"; echo "Bonjour Jimmy"; echo "Bonjour Rose"; echo "Bonjour Sarah";
Bonjour Charlotte
Bonjour Joe
Bonjour Jimmy
Bonjour Rose
Bonjour Sarah
```

C'est mieux, mais nous écrivons toujours cinq commandes séparées avec beaucoup de texte répété. Ne serait-il pas génial de simplifier cela pour dire à bash d'imprimer "Bonjour" devant une liste de noms ? C'est là que nous pouvons utiliser une boucle for, qui a la structure suivante :

**for une nouvelleVariable dans une liste ; do une commande à la $newVariable ; done**

|Action | Code|
|---|--|
|type de boucle|for |
|variable | prenom|
|emplacement|in|
|liste|Charlotte Joe Jimmy Rose Sarah|
|début|do;|
|commande| echo $prenom;|
|fin|done;|

```
orkin@ANT3475:~$ for prenom in Charlotte Joe Jimmy Rose Sarah; do echo "Bonjour $prenom"; done
Bonjour Charlotte
Bonjour Joe
Bonjour Jimmy
Bonjour Rose
Bonjour Sarah
```

Notez ici l'importance des points-virgules, qui servent à séparer les trois parties de la boucle

```for prenom in Charlotte Joe Jimmy Rose Sarah;```
```do "echo $prenom";```
```done```

Faisons quelque chose de un peu plus complexe. Supposons que nous ayons un dossier avec 50 fichiers dedans que nous devons renommer. Exécutez le script shell [makeFiles.sh](https://raw.githubusercontent.com/nomascus/ANT3814/main/FILES/makeFiles.sh) pour créer 50 fichiers vides dans un répertoire appelé "FILES", puis répertoriez le contenu du dossier

```
orkin@ANT3475:~$ bash makeFiles.sh
orkin@ANT3475:~$ ls FILES/
file1.txt  file14.txt file19.txt file23.txt file28.txt file32.txt file37.txt file41.txt file46.txt file50.txt
file10.txt file15.txt file2.txt  file24.txt file29.txt file33.txt file38.txt file42.txt file47.txt file6.txt
file11.txt file16.txt file20.txt file25.txt file3.txt  file34.txt file39.txt file43.txt file48.txt file7.txt
file12.txt file17.txt file21.txt file26.txt file30.txt file35.txt file4.txt  file44.txt file49.txt file8.txt
file13.txt file18.txt file22.txt file27.txt file31.txt file36.txt file40.txt file45.txt file5.txt  file9.txt
```

Remarquez que les fichiers 1-9 ne sont pas listés en premier. C'est parce que bash, comme de nombreux systèmes d'exploitation, classe les caractères selon l'ordre ascii :

```blanc ! " # $ % & ' ( ) * + , - . /0 1 2 3 4 5 6 7 8 9 : ; < = > ? @ A B C D E F G H I J K L M N O P Q R S T U V W X Y Z[ \\] ˆ_ a b c d e f g h i j k l m n o p q r s t u v w x y z { } ~\```

Dans ce système de classement, le caractère point vient avant les chiffres, donc le fichier1.txt est classé avant le fichier10.txt et le fichier4.txt ne vient qu'après le fichier39.txt.

Rendons les noms de fichiers plus faciles à classer en ajoutant un zéro aux noms des 9 premiers fichiers. Nous pouvons le faire avec une boucle for qui liste les nombres de 1 à 9 :

```
for number in 1 2 3 4 5 6 7 8 9; do mv file${number}.txt file0${number}.txt ; done
```

Attention, cependant ! Créer une boucle qui renomme vos fichiers est dangereux, car vous pouvez accidentellement mélanger tous vos fichiers ou les effacer en cas d'erreur dans votre code. **Testez toujours votre boucle en ajoutant echo** à la commande d'abord pour voir si le code est correct.

 ```
orkin@ANT3475:~/FILES$ for number in 1 2 3 4 5 6 7 8 9; do echo "mv file${number}.txt file0${number}.txt" ; done
mv file1.txt file01.txt
mv file2.txt file02.txt
mv file3.txt file03.txt
mv file4.txt file04.txt
mv file5.txt file05.txt
mv file6.txt file06.txt
mv file7.txt file07.txt
mv file8.txt file08.txt
mv file9.txt file09.txt
```

Ici, nous utilisons la boucle pour imprimer d'abord les commandes à l'écran. Une fois que tout semble bon, vous pouvez supprimer la commande echo et exécuter le code pour renommer les fichiers.

```
orkin@ANT3475:~/FILES$ ls 
file01.txt file06.txt file11.txt file16.txt file21.txt file26.txt file31.txt file36.txt file41.txt file46.txt
file02.txt file07.txt file12.txt file17.txt file22.txt file27.txt file32.txt file37.txt file42.txt file47.txt
file03.txt file08.txt file13.txt file18.txt file23.txt file28.txt file33.txt file38.txt file43.txt file48.txt
file04.txt file09.txt file14.txt file19.txt file24.txt file29.txt file34.txt file39.txt file44.txt file49.txt
file05.txt file10.txt file15.txt file20.txt file25.txt file30.txt file35.txt file40.txt file45.txt file50.txt
```

Si vous avez beaucoup de nombres dans votre boucle for, vous n'avez pas besoin de les écrire tous. Vous pouvez utiliser {1..9} à la place de 1 2 3 4 5 6 7 8 9. Les deux nombres à l'intérieur des accolades séparés par deux points indiquent une plage de nombres pour la boucle

```
for number in {1..9}; do echo "mv file${number}.txt file0${number}.txt" ; done
```

## Boucles while

Un des types de boucles les plus puissants est la boucle while, car elle vous permet de dire à bash de faire quelque chose tant qu'une condition est vraie. Le plus souvent, cela est utilisé pour lire un fichier et exécuter une commande en fonction des informations contenues dans chaque ligne. 

```
while condition est vraie; do 
    code à exécuter; 
done
```

Nous pouvons le faire à la ligne de commande en une seule ligne de texte. Disons que nous voulons lire un fichier et imprimer les lignes avec des informations devant chaque ligne

```
orkin@ANT3475:~/FILES$ lineNumber=1; while read line; do echo "Line $counter $line"; ((lineNumber++)) ; done < apol.fasta
Line 1 >seq1
Line 2 GAATTCAAGTTCTTGTGCGCACACAAATCCAATAAAAACTATTGTGCACACAGACGCGAC
Line 3 TTCGCGGTCTCGCTTGTTCTTGTTGTATTCGTATTTTCATTTCTCGTTCTGTTTCTACTT
Line 4 AACAATGTGGTGATAATATAAAAAATAAAGCAATTCAAAAGTGTATGACTTAATTAATGA
Line 5 GCGATTTTTTTTTTGAAATCAAATTTTTGGAACATTTTTTTTAAATTCAAATTTTGGCGA
Line 6 AAATTCAATATCGGTTCTACTATCCATAATATAATTCATCAGGAATACATCTTCAAAGGC
Line 7 AAACGGTGACAACAAAATTCAGGCAATTCAGGCAAATACCGAATGACCAGCTTGGTTATC
Line 8 AATTCTAGAATTTGTTTTTTGGTTTTTATTTATCATTGTAAATAAGACAAACATTTGTTC
Line 9 CTAGTAAAGAATGTAACACCAGAAGTCACGTAAAATGGTGTCCCCATTGTTTAAACGGTT
Line 10 GTTGGGACCAATGGAGTTCGTGGTAACAGTACATCTTTCCCCTTGAATTTGCCATTCAAA
Line 11 ATTTGCGGTGGAATACCTAACAAATCCAGTGAATTTAAGAATTGCGATGGGTAATTGACA
Line 12 TGAATTCCAAGGTCAAATGCTAAGAGATAGTTTAATTTATGTTTGAGACAATCAATTCCC
Line 13 CAATTTTTCTAAGACTTCAATCAATCTCTTAGAATCCGCCTCTGGAGGTGCACTCAGCCG
Line 14 CACGTCGGGCTCACCAAATATGTTGGGGTTGTCGGTGAACTCGAATAGAAATTATTGTCG
Line 15 CCTCCATCTTCATGGCCGTGAAATCGGCTCGCTGACGGGCTTCTCGCGCTGGATTTTTTC
Line 16 ACTATTTTTGAATACATCATTAACGCAATATATATATATATATATTTAT
```

Ici, nous donnons à une variable appelée ```lineNumber``` une valeur de 1. Ensuite, nous commençons la boucle while, disant à bash de lire les lignes dans un fichier. Nous lui dirons où se trouve le fichier plus tard, mais il est important de savoir que **while** et **read** sont des commandes interprétées par bash, **line** est juste un nom de variable que nous choisissons pour représenter chaque ligne de métadonnées alors qu'elle passe à travers la boucle que nous écrivons. Nous pourrions utiliser n'importe quel nom que nous voulons, mais ```line``` a du sens. Ensuite, nous utilisons ```echo``` parce que nous voulons imprimer à l'écran le mot 'Ligne' suivi de la valeur du compteur et ensuite la ligne du fichier. ((lineNumber++)) est un astuce pour ajouter 1 à la variable lineNumber à chaque fois que la boucle s'exécute. Enfin, nous fermons la boucle avec ```done``` et lui donnons le fichier à lire (apol.fasta). N'oubliez pas que le signe inférieur à ```<``` indique que nous envoyons des informations du fichier (que nous avons enregistré comme une variable à la ligne 2) à droite dans le programme à gauche


Essayons un exemple plus difficile, mais plus réaliste. Disons que nous avons des lectures de séquençage de cinq échantillons que nous voulons mapper sur différents génomes de référence en utilisant le programme **bwa**. Nous avons la feuille de données de métadonnées suivante avec une liste

 d'informations sur chaque échantillon dans le format suivant:

|Espèce|Individu|Génome de référence
|-|-|-| 
|Lemur_catta|LC_01|mLemcat1.fa
|Lemur_catta|LC_12|mLemcat1.fa
|Pan_paniscus|PP_09|Ptro3.fa
|Homo_sapiens|HS_01|hg38.fa
|Homo_sapiens|HS_22|hg38.fa

Si nous voulons exécuter bwa, nous devons lui indiquer les informations suivantes
 
```
bwa mem forward_sequencing_reads.fasta reverse_sequencing_reads.fasta reference_genome.fa -o outputFileName.sam
```

Pour la première ligne de notre table de métadonnées, cela ressemblerait à ceci:

``` bwa mem Lemur_catta_LC_01_F.fasta Lemur_catta_LC_01_R.fasta mLemcat1.fa -o Lemur_catta_LC_01.sam```

Mais ce serait pénible d'écrire toutes ces lignes de texte, surtout si vous travaillez avec des centaines d'échantillons. Nous pouvons écrire un **script** en utilisant une boucle while pour générer le code pour un nombre infini de tâches bwa!

```
#!/bin/bash

filename='bwa_metadata.txt'

while read line; do

   ID=$(echo $line | awk '{print $1_$2}')
   ref=$(echo $line | awk '{print $3}')

   echo "bwa mem ${ID}_F.fasta ${ID}_R.fasta $ref -o ${ID}.sam"

done < $filename
```

Ce script est un fichier texte que nous sauvegardons sous le nom de bwa_runBatch.sh. Il contient la boucle while complète et quelques instructions supplémentaires pour indiquer à bash comment l'exécuter. C'est un peu compliqué, alors décomposons-le ligne par ligne.

**Ligne 1** 
> #!/bin/bash

Cette ligne spéciale en haut du script indique à bash dans quel langage le script est écrit. Cette ligne sera un peu différente pour d'autres langages comme python ou perl.

**Ligne 2**
>filename='bwa_metadata.txt'

Nous sauvegardons le nom du fichier de métadonnées que nous voulons utiliser dans une variable appelée $filename. 

**Ligne 3**
> while read line; do

Cette ligne commence la boucle while, indiquant à bash de lire les lignes dans un fichier. Nous lui dirons où se trouve le fichier plus tard. Il est important de savoir que **while** et **read** sont des commandes interprétées par bash, **line** est juste un nom de variable que nous choisissons pour représenter chaque ligne de métadonnées lorsqu'elle passe à travers le programme que nous écrivons. Nous pourrions utiliser n'importe quel nom, mais "line" a du sens.  

**Ligne 4**
>ID=$(echo $line | awk '{print $1"_"$2}')

Maintenant que la première ligne de métadonnées a été lue par bash, nous pouvons commencer à faire des choses avec. Ici, nous utilisons awk pour séparer les deux premières colonnes (Lemur_catta et LC_01) et les fusionner en une seule chaîne séparée par un trait de soulignement (Lemur_catta_LC_01). Ensuite, nous le sauvegardons dans une variable appelée ```ID``` pour le stocker et l'utiliser plus tard dans le script.

**Ligne 5**
>  ref=$(echo $line | awk '{print $3}')

Comme dans la ligne 4, nous sauvegardons une autre partie de la ligne de métadonnées. Cette fois, il s'agit de la colonne 3 (mLemcat1.fa), qui est le nom du génome de référence. Nous le sauvegardons dans une variable appelée ```ref``` pour l'utiliser plus tard dans le script. Maintenant que la première ligne de métadonnées a été analysée et stockée dans des variables, nous pouvons utiliser les variables pour exécuter bwa. 

**Ligne 6**
>echo "bwa mem ${ID}_F.fasta ${ID}_R.fasta $ref -o ${ID}.sam"

Tout d'abord, nous exécutons une commande ```echo```, car nous voulons afficher les résultats de tout ce qui se trouve entre guillemets à l'écran. Nous commençons par appeler le programme bwa mem, puis nous lui fournissons quatre informations : les fichiers de lecture de séquençage avant et arrière, le génome de référence et le nom du fichier de sortie. Nous insérons les variables que nous avons créées avec awk ($ID et $ref) dans la chaîne de texte.

${ID}_F.fasta -> Lemur_catta_LC_01_F.fasta
${ID}_R.fasta -> Lemur_catta_LC_01_R.fasta
$ref -> mLemcat1.fa 
${ID}.sam -> Lemur_catta_LC_01.sam

**Ligne 7**
>done < $filename

Enfin, nous clôturons la boucle avec ```done```. C'est également à ce moment que nous donnons au fichier de métadonnées comme entrée à la boucle. Rappelez-vous que le signe inférieur ``` <``` indique que nous envoyons des informations du fichier (que nous avons sauvegardé sous forme de variable à la ligne 2) sur la droite dans le programme sur la gauche.

Alors exécutons notre script!

```
orkin@ANT3475:~/FILES$ bash bwa_batchRun.sh
bwa mem Lemur_catta_LC_01_F.fasta Lemur_catta_LC_01_R.fasta mLemcat1.fa -o Lemur_catta_LC_01.sam
bwa mem Lemur_catta_LC_12_F.fasta Lemur_catta_LC_12_R.fasta mLemcat1.fa -o Lemur_catta_LC_12.sam
bwa mem Pan_paniscus_PP_09_F.fasta Pan_paniscus_PP_09_R.fasta Ptro3.fa -o Pan_paniscus_PP_09.sam
bwa mem Homo_sapiensHS_01_F.fasta Homo_sapiens_HS_01_R.fasta hg38.fa -o Homo_sapiens_HS_01.sam
bwa mem Homo_sapiens_HS_22_F.fasta Homo_sapiens_HS_22_R.fasta hg38.fa -o Homo_sapiens_HS_22.sam
```
Regardez ça!


# Éditeurs de texte

[L'introduction suivante à VI est tirée de CSHL Programming For Biology et a été initialement rédigée par Sofia Robb]

Il est souvent nécessaire de créer et d'écrire dans un fichier tout en utilisant le terminal. Cela rend essentiel l'utilisation d'un éditeur de texte de terminal. Il existe de nombreux éditeurs de texte. Certains de nos préférés sont Emacs et vim. Nous allons vous initier à un éditeur de texte simple appelé `vi`.

## Introduction à vi

Qu'est-ce que **vi** ?

> **vi** est un éditeur de texte en ligne de commande. vi est inclus dans chaque installation Linux. Vous n'avez jamais besoin de l'installer.

Qu'est-ce qu'un éditeur de texte en ligne de commande ?

> Un éditeur de texte en ligne de commande est un éditeur de texte que vous utilisez à partir de la ligne de commande. Dans la plupart des éditeurs de texte en ligne de commande, ne vous attendez pas à pouvoir cliquer pour naviguer. Vous devrez vous déplacer avec des frappes de touches du clavier. Les deux éditeurs de texte les plus populaires sont **vi** et **emacs**. Vous pouvez utiliser l'un ou l'autre, mais nous commencerons par **vi** car les frappes de touches sont moins complexes que dans **emacs**.

Pourquoi devrais-je me soucier des éditeurs de texte en ligne de commande ?

> Si vous êtes connecté à une machine distante, un éditeur de texte en ligne de commande est le moyen le plus rapide, le plus facile et le plus efficace d'écrire des fichiers texte.

### Pour commencer avec vi

__Ouvrir un fichier__

Sur la ligne de commande, tapez `vi` suivi d'un nom de fichier.

```bash
orkin@ANT3475:~$ vi <fichier>
```
Essayons :

```bash
orkin@ANT3475:~$ vi premier_fichier_vi.txt
```
Vous verrez ceci dans votre terminal.

```
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
"premier_fichier_vi.txt" [New File]
```
Remarquez le nom du fichier en bas.

Si vous **n'incluez pas** de nom de fichier, vous verrez quelque chose de similaire à ceci :


```bash
~
~
~
~
~                                                VIM - Vi IMproved
~
~                                                version 8.0.1283
~                                            by Bram Moolenaar et al.
~                                   Vim is open source and freely distributable
~
~                                          Become a registered Vim user!
~                                 type  :help register<Enter>   for information
~
~                                 type  :q<Enter>               to exit
~                                 type  :help<Enter>  or  <F1>  for on-line help
~                                 type  :help version8<Enter>   for version info
~
~
~

```

Lisez le message et tapez `:q<Enter>` pour **q**uitter ou sortir.

__vi a deux modes.__

1. **Mode Insertion**

2. **Mode Commande**

**Mode Insertion** est pour taper le contenu de votre fichier. Toutes les frappes de touches seront interprétées comme des caractères que vous souhaitez voir dans votre fichier.

**Mode Commande** est pour utiliser des commandes. Toutes les frappes de touches seront interprétées comme des commandes et ***non*** comme faisant partie de votre fichier. Les commandes courantes sont pour supprimer, copier, rechercher, remplacer et enregistrer.

### Créer, Écrire et Enregistrer un Fichier

__Créer__

À partir de la ligne de commande, ouvrez un nouveau fichier en tapant

`vi premier_fichier_vi.txt` 

__Écrire__

Commencez à taper du contenu. Pour cela, nous devons entrer dans le **Mode Insertion**. 

Pour cela, tapez `i`.

Votre session vi ressemblera maintenant à ceci :

```bash
~
~
~
~
~
~
-- INSERT (paste) --
```

Remarquez l'`INSERT` en bas de l'écran.

Commencez à taper le contenu de votre fichier. N'oubliez pas que toutes les frappes de touches sont celles que vous voulez voir dans votre fichier et que votre souris ne fonctionnera pas.

__Enregistrer__

Maintenant que le fichier contient du contenu, passons en **Mode Commande** pour enregistrer notre fichier.

1. Appuyez sur la touche `<ESC>` pour entrer dans le **Mode Commande**.
2. tapez `:w` (deux points suivis d'un w) pour **Enregistrer (Write)**

Si vous voulez taper encore du contenu, entrez en **Mode Insertion** (`i`). 

Si, au lieu de cela, vous voulez quitter, puisque vous êtes déjà en **Mode Commande**, vous pouvez utiliser les touches de sortie `:q`

### Activités courantes et commandes vi

Passez en **Mode Commande** pour toutes les commandes. Si vous n'êtes pas sûr d'être en **Mode Commande**, appuyez simplement sur la touche `<esc>`. Cela n'affectera pas si vous êtes déjà en **Mode Commande**

__Enregistrer et Quitter__

N'oubliez pas de passer en **Mode Commande** avec la touche `<esc>`

|frappe de touche | signification | 
|-----------|---------|  
|`:wq` | Enregistrer (**w**rite) et **q**uitter |   
|`:q!` | **q**uitter sans Enregistrer !!! |  
|`:w`  | Enregistrer (**w**rite) uniquement |  

 
La plupart des commandes dans vi sont exécutées dès que vous appuyez sur une séquence de touches. Toute commande commençant par deux points ( : ) nécessite que vous appuyiez sur `<enter>` pour compléter la commande.

__Se déplacer__

N'oubliez pas de passer en **Mode Commande** avec la touche `<esc>`

|frappe de touche | signification | 
|-----------|---------|  
|Touches fléchées  | déplacer le curseur |  
|`j`, `k`, `h`, `l` | déplacer le curseur vers le bas, vers le haut, vers la gauche et vers la droite (similaire aux touches fléchées) |  
|`0` (zéro) | déplacer le curseur au début de la ligne actuelle |  


|`^` (accent circonflexe) | déplacer le curseur au début de la ligne actuelle |  
|`$`  |  déplacer le curseur à la fin de la ligne actuelle |  
|`:n` | aller à la n-ième ligne du fichier |  
|`nG` | aller à la n-ième ligne (par exemple 5G pour aller à la 5e ligne) |  
|`G` | aller à la dernière ligne |  
|`w` | aller au début du mot suivant |  
|`nw`| avancer de n mots (par exemple, 2w pour avancer de deux mots) |  
|`b` | aller au début du mot précédent |  
|`nb` | reculer de n mots |  


__Supprimer du contenu__

N'oubliez pas de passer en **Mode Commande** avec la touche `<esc>`

|frappe de touche | signification | 
|-----------|---------|  
|`x`        | supprimer un seul caractère  |   
|`nx`       | supprimer n caractères (par exemple, 5x pour supprimer cinq caractères)  |  
|`dd`       | supprimer la ligne actuelle  |  
|`dn`       | d suivi d'une commande de déplacement. Supprimez jusqu'à où la commande de déplacement vous aurait mené. (par exemple, d5w signifie supprimer 5 mots)|  

__Annuler__

N'oubliez pas de passer en **Mode Commande** avec la touche `<esc>`

|frappe de touche | signification | 
|-----------|---------|  
| `u` | Annuler la dernière action (vous pouvez continuer à appuyer sur u pour continuer à annuler)  |  
| `U` | (Remarque: majuscule)** - Annuler toutes les modifications de la ligne actuelle  |  

### Autres astuces utiles

|frappe de touche | signification | 
|-----------|---------|  
| `:set number` | afficher les numéros de ligne | 
| `:set nonumber` | désactiver les numéros de ligne |  
| `:/[texte de recherche]` | rechercher [ce texte] dans votre fichier |  
| `n` | aller à l'occurrence suivante de votre recherche |  

[Unix4 ProblemSet](https://raw.githubusercontent.com/nomascus/ANT3814/main/PROBLEM_SETS/Unix4_problemSet_Francais.md)
