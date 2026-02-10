# Révision rapide des bases Unix et ensemble de problèmes

commandes utiles
=========================

| commande         | description                              |
| --------------- | ---------------------------------------- |
| `ls`            | affiche le contenu du répertoire                  |
| `cd`            | changer de répertoire                         |
| `mkdir`         | créer un répertoire                         |
| `rm`            | supprimer des fichiers et répertoires. Faites attention, il est facile de supprimer plus que ce que vous voulez. |
| `head`          | affiche les premières lignes dans la fenêtre du terminal |
| `tail`          | affiche les dernières lignes dans la fenêtre du terminal |
| `sort`          | trie les lignes                          |
| `uniq`          | affiche les lignes uniques                  |
| `grep`          | trouve les lignes contenant un motif  |
| `wc`            | compte le nombre de lignes, de caractères et de mots |
| `mv`            | déplace des fichiers                               |
| `cp`            | copie des fichiers                               |
| `date`          | renvoie la date et l'heure actuelles        |
| `pwd`           | renvoie le nom du répertoire de travail            |
| `ssh`           | connexion à distance                             |
| `scp`           | copie sécurisée à distance                       |
| `~`             | raccourci pour votre répertoire personnel         |
| `man <commande>` | page de manuel pour la commande, par exemple `man ls` pour obtenir la page de manuel pour `ls` |





### Essayez ces exemples


Allons dans un répertoire contenant de nombreux fichiers et listons ces fichiers

```bash
cd /bin/
ls
```

__Quelle est la différence entre ces deux commandes?__

Essayez-les toutes les deux!!
```
ls -l
ls -lt
```


__Pipes__ 

Vous pouvez chaîner plusieurs commandes avec un pipe `|`, de sorte que la sortie standard de la première commande soit 'piped' dans l'entrée standard de la deuxième commande.

Essayez!!
```
ls -lt | head
```



__Points-virgules__

Vous pouvez chaîner plusieurs commandes en mettant un point-virgule `;` après chaque commande. Ici, les commandes seront exécutées séquentiellement, mais toute sortie ne sera pas transmise d'une commande à l'autre.

Essayez!!
```
date ; sleep 2 ; date
```
> Si vous voulez en savoir plus sur `sleep`, tapez `man sleep`


__Télécharger un fichier__. 

Changez de répertoire vers votre répertoire personnel. Vous avez probablement les autorisations pour écrire dans votre répertoire personnel. Maintenant, utilisez `wget` ou `curl` pour télécharger des fichiers. Sur certains systèmes, seul l'un de ces deux peut être disponible

```
cd ~

curl -O https://raw.githubusercontent.com/nomascus/ANT3814/main/PROBLEM_SETS/cuffdiff.txt
```
Notez que '-O' est la lettre O, pas le chiffre zéro 0

__Rediriger STDOUT__  
Vous pouvez rediriger la sortie d'une commande dans un fichier. 
```
cd ~
grep Chr7 cuffdiff.txt > fav_chr_cuffdiff.txt
```



__Ajouter STDOUT à la fin d'un fichier qui existe déjà__

Vous pouvez ajouter la sortie d'une commande à un fichier

```
grep Chr9 cuffdiff.txt >> fav_chr_cuffdiff.txt
```



__Rediriger STDERR__  

Vous pouvez rediriger STDERR vers un fichier.

Revoyons ce qu'est réellement STDERR.
```
cat blablabla.txt
```
> Le fichier blablabla.txt n'existe pas, donc nous obtenons `cat: blablabla.txt: No such file or directory` imprimé dans le terminal. Ce message est étiqueté par le système d'exploitation comme un message d'erreur ou STDERR. 

STDERR est un type de sortie étiquetée que nous pouvons rediriger

```
cat blablabla.txt 2> errors.txt
```
> Nous pouvons rediriger les messages d'erreur, alias STDERR, vers un nouveau fichier appelé comme nous le voulons


Que se passe-t-il lorsque vous essayez de rediriger STDOUT?
```
cat blablabla.txt > errors.txt
```
> `cat: blablabla.txt: No such file or directory` est toujours imprimé à l'écran car nous redirigeons seulement STDOUT vers notre fichier. Il n'y a pas de STDOUT dans ce cas et notre fichier sera vide. Comment le vérifieriez-vous?




__Rediriger STDOUT et STDERR__

Vous pouvez rediriger à la fois STDOUT et STDERR vers **deux fichiers séparés** en une seule commande.

```
# affichez-le d'abord dans le terminal
cat fav_chr_cuffdiff.txt blablabla.file

# rediriger vers deux fichiers, STDOUT vers out.txt, STDERR vers err.txt 
cat fav_chr_cuffdiff.txt blablabla.file 1> out.txt 2> err.txt

# cela fait la même chose, voyez-vous la différence?
cat fav_chr_cuffdiff.txt blablabla.file > out.txt 2> err.txt

```
> Examinez le contenu de `out.txt` et `err.txt`

Vous pouvez également rediriger à la fois STDOUT et STDERR vers **le même** fichier.
```
cat fav_chr_cuffdiff.txt blablabla.file &> all_out_err.txt 
```
> Consultez ce qui se trouve dans `all_out_err.txt`

Ensemble de problèmes
===========

1. Connectez-vous à votre machine. 

2. Quel est le chemin complet de votre répertoire personnel?

3. Remontez d'un répertoire?
    - Combien de fichiers contient-il?
    - Combien de répertoires?

4. Créez un répertoire appelé `ANT3814_problemsets` dans votre répertoire personnel.

5. Accédez à ce nouveau répertoire appelé problemsets. Vérifiez que vous êtes dans le bon répertoire en utilisant `pwd`.

6. Utilisez `wget` pour copier <https://raw.githubusercontent.com/nomascus/ANT3814/main/PROBLEM_SETS/sequences.nt.fa> depuis le web dans votre répertoire

 problemsets. Si `wget` n'est pas disponible sur votre système, utilisez `curl -O` comme alternative. (Notez que '-O' est la lettre O, pas le chiffre zéro 0)

7. Sans utiliser d'éditeur de texte, utilisez des commandes Unix pour trouver ces informations pour le fichier `sequences.nt.fa`.
  Ce fichier peut être trouvé ici <https://raw.githubusercontent.com/nomascus/ANT3814/main/PROBLEM_SETS/sequences.nt.fa>
      - Combien de lignes contient ce fichier?   
      - Combien de caractères?    (Astuce: consultez les options de wc)
      - Quelle est la première ligne de ce fichier?    (Astuce: lisez la page de manuel de head)
      - Quelles sont les 3 dernières lignes?    (Astuce: lisez la page de manuel de tail)
      - Combien de séquences y a-t-il dans le fichier?    (Astuce: utilisez grep) (Remarque : Le début d'une séquence est indiqué par un caractère `>`.)    

8. Renommez `sequences.nt.fa` en `cancer_genes.fasta`. (Astuce: lisez la page de manuel pour mv)

9. Copiez ce fichier distant, cuffdiff.txt, dans votre répertoire problemset. Voici l'URL que vous pouvez utiliser : <https://raw.githubusercontent.com/nomascus/ANT3814/main/PROBLEM_SETS/cuffdiff.txt>

Utilisez `wget` pour copier <https://raw.githubusercontent.com/nomascus/ANT3814/main/PROBLEM_SETS/cuffdiff.txt> depuis le web dans votre répertoire problemsets. Si `wget` n'est pas disponible sur votre système, utilisez `curl -O` comme alternative. (Notez que '-O' est la lettre O, pas le chiffre zéro 0)




10. Faites ce qui suit pour `cuffdiff.txt`. Les descriptions de chaque colonne du fichier sont dans le tableau ci-dessous.
    - Regardez les premières lignes du fichier
    - Triez le fichier par changement de log 'log2(fold_change)', du plus élevé au plus bas, et enregistrez dans un nouveau fichier dans votre répertoire appelé sorted.cuffdiff.out
    - Triez le fichier (log du changement le plus élevé au plus bas), puis affichez uniquement les 100 premières lignes. Enregistrez dans un fichier appelé `top100.sorted.cuffdiff.out`.
    - Triez le fichier par changement de log, affichez les 100 premiers, affichez uniquement la première colonne. Il s'agira d'une liste des 100 gènes avec le plus grand changement d'expression. Assurez-vous que votre liste est triée par nom de gène et est unique. Enregistrez cette liste dans un fichier appelé `differentially.expressed.genes.txt`.



__Format du fichier Cuffdiff__

| Numéro de colonne | Nom de la colonne       | Exemple              | Description                              |
| ------------- | ----------------- | -------------------- | ---------------------------------------- |
| 1             | Identifiant testé         | XLOC_000001          | Un identifiant unique décrivant le transcrit, le gène, le transcrit principal ou le CDS testé |
| 2             | Identifiant testé         | XLOC_000001          | Un identifiant unique décrivant le transcrit, le gène, le transcrit principal ou le CDS testé |
| 3             | gène              | Lypla1               | Le nom du gène ou les identifiants du gène qui sont testés |
| 4             | locus             | chr1:4797771-4835363 | Coordonnées génomiques pour faciliter la navigation vers les gènes ou transcrits testés. |
| 5             | échantillon 1          | Foie                | Libellé (ou numéro si aucun libellé n'est fourni) du premier échantillon testé |
| 6             | échantillon 2          | Cerveau                | Libellé (ou numéro si aucun libellé n'est fourni) du deuxième échantillon testé |
| 7             | Statut du test       | NOTEST               | Peut être OK (test réussi), NOTEST (pas assez d'alignements pour le test), LOWDATA (trop complexe ou séquencé de manière trop superficielle), HIDATA (trop de fragments dans le locus), ou FAIL,  lorsque une matrice de covariance mal conditionnée ou une autre exception numérique empêche le test. |
| 8             | FPKMx             | 8.01089              | FPKM du gène dans l'échantillon x             |
| 9             | FPKMy             | 8.551545             | FPKM du gène dans l'échantillon y             |
| 10            | log2(FPKMy/FPKMx) | 0.06531              | Le log (base 2) du changement y/x  |
| 11            | test stat         | 0.860902             | La valeur de la statistique de test utilisée pour calculer la signification du changement observé dans FPKM |
| 12            | p value           | 0.389292             | La valeur p non corrigée de la statistique de test |
| 13            | q value           | 0.985216             | La valeur p corrigée de la FDR de la statistique de test |
| 14            | significatif       | no                   | Peut être "oui" ou "non", en fonction de la valeur de p qui est supérieure à la FDR après correction de Benjamini-Hochberg pour les tests multiples |
