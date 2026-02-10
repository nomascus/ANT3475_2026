# Unix2: Problème défini

  

Utilisez la ligne de commande pour répondre à toutes les questions suivantes. Pour chaque question, copiez la ligne de code que vous avez utilisée pour répondre à la question et enregistrez-la dans un fichier texte.

  

1) Affichez votre UTILISATEUR et HOSTNAME à l'écran

  

2) Utilisez une variable d'environnement pour afficher un nombre aléatoire à l'écran

  

3) Enregistrez les variables suivantes en bash : votre nom, votre groupe préféré, le titre de votre livre préféré, une liste de 3 films, le nombre 7, le nombre 3.141592654

  

4) Quel type de données est chaque variable ?

  

5) Utilisez les variables pour écrire une phrase qui liste votre nom, votre groupe préféré, le titre de votre livre préféré et les 3 films

  

6) Multipliez ensemble les deux variables numériques

  

7) Ajoutez un quatrième film à la variable qui est votre liste de 3 films en utilisant la méthode $( )

  

8) Ajoutez un cinquième film à la variable qui est votre liste de 3 films en utilisant la méthode ${ }

  

Téléchargez la séquence génique du gène de l'insuline humaine (INS) à partir d'ici :

https://raw.githubusercontent.com/nomascus/ANT3814/main/FILES/INS_gene.fasta

  

9) Utilisez grep -v pour extraire la ligne d'en-tête et enregistrez-la en tant que variable appelée $headerINS

  

10) Utilisez grep pour extraire les lignes de séquence et enregistrez-les en tant que variable appelée $sequenceINS.

  

11) Affichez $sequenceINS à l'écran. Pourquoi y a-t-il des espaces entre chaque ligne de séquence ?

  

12) Utilisez sed pour supprimer les espaces entre les lignes de séquence. Notez ici que vous devez utiliser le caractère spécial \s à la place d'un espace blanc dans sed.

  

13)

Téléchargez le fichier vcf suivant: https://raw.githubusercontent.com/nomascus/ANT3814/main/FILES/chr22.1kg.phase3.v5a.trimmed.vcf.gz

Décompressez le fichier à l'aide de la commande suivante : gunzip chr22.1kg.phase3.v5a.trimmed.vcf.gz

Ensuite, le fichier doit se terminer par .vcf au lieu de .vcf.gz

  

Les génotypes dans ce VCF sont indiqués comme suit :

  

Allèle de référence homozygote :  0|0

Allèle alternatif homozygote :  1|1

Hétérozygote :  0|1 ou 0|1

  

13) Regardez le fichier. Combien y a-t-il de lignes, et combien de lignes sont des lignes d'en-tête par rapport aux lignes de données ? Combien d'individus y a-t-il dans le VCF ?

  

14) Utilisez l'approche cat <( ) < ( ) pour afficher les lignes d'en-tête et les 10 dernières lignes

  

15) Imprimez une version simplifiée du VCF qui inclut uniquement le chromosome, la position, l'allèle de référence, l'allèle alternatif et les génotypes des individus. Enregistrez ce fichier sous le nom de clean.vcf

  

BONUS : Utilizer grep -E pour créer une version de clean.vcf qui n'inclut que les lignes avec au moins 1 position qui est hétérozygote. Combien y a-t-il de ces positions ?

