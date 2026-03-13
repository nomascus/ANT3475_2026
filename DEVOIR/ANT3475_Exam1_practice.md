# Exam Practice

**Lors de l'examen, vous pouvez utiliser vos diapositives de cours et vos devoirs, mais vous NE POUVEZ PAS utiliser Internet ou l'IA. Si je vous vois ouvrir un autre programme ou navigateur Internet, vous obtiendrez un zéro à l'examen.**


1) Écrivez la commande pour déplacer un fichier nommé `results.txt` du répertoire courant vers un répertoire nommé `output` situé dans le même répertoire.

2) Quelle est la différence entre les guillemets simples et les guillemets doubles en bash ?

3) Comment concaténez-vous les valeurs de deux variables, `var1` et `var2`, et stockez-vous le résultat dans une nouvelle variable `var3` ?

4) Écrivez une expression régulière pour correspondre à une adresse e-mail.

5) Écrivez une commande Perl pour substituer toutes les instances de ATG par TAC dans la ligne de séquence suivante ATGCGAAGGATGAAG.

6) Quel est le résultat du script suivant ? Les 5 premières lignes du fichier `metadata.txt` sont :

primates  
carnivores  
cétacés  
félidés  
canidés
```
#!/bin/bash

while read line; do
   alignment=$(echo $line | awk '{print $1}')
   echo "iqtree2 -s ${alignment}.aligned.gb.fasta -m MFP -b 1000"
done < metadata.txt
```

1) mv results.txt output

2) Guillemets simples `' '` : tout est pris littéralement.
Guillemets doubles `" "` : les variables (`$var`), substitutions (`$(...)`) et certains échappements (`\`) sont interprétés.

3) var3=${var1}${var2}

4)  \^\w+@\w+\\.\w+$

5) echo ATGCGAAGGATGAAG | perl -lpe 's/ATG/TAC/g'

6) Pour chaque ligne du fichier metadata.txt, le script prend le texte de la première colonne, l'enregistre dans une variable appelée alignment, puis l'insère dans la commande  commençant par iqtree2 et l'affiche à l'écran.
