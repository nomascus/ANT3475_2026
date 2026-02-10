# Unix 3

Pour accéder à la leçon d’aujourd’hui, téléchargez le fichier html appelé 3_Unix3_2026.html depuis StudiUM. Il devrait s’ouvrir avec votre navigateur web.

Commençons par créer un nouveau dossier dans votre cluster pour la leçon de cette semaine. Après vous être connecté, nous commencerons par votre dossier personnel. Descendez ensuite d’un répertoire jusqu’au dossier scratch. Créez ensuite un nouveau dossier appelé 3_Unix3, et déplacez-vous dans ce dossier.

```
cd ~/scratch
mkdir 3_Unix3 
cd 3_Unix3
```
Copiez les fichiers file1.txt, file2.txt, file3.txt et file4.txt du répertoire que vous avez créé la dernière fois (2_Unix_2) dans votre nouveau répertoire 3_Unix_3.

Après avoir fait cela, téléchargez les fichiers suivants pour cette leçon en copiant, collant et exécutant chaque ligne de ce code.

```
wget https://raw.githubusercontent.com/nomascus/ANT3814/refs/heads/main/FILES/FOXP2.fasta
wget https://raw.githubusercontent.com/nomascus/ANT3814/refs/heads/main/FILES/FOXP2_long.fasta 
wget https://raw.githubusercontent.com/nomascus/ANT3814/main/FILES/apol.fasta
wget https://raw.githubusercontent.com/nomascus/ANT3814/main/FILES/unix3.fasta
wget https://raw.githubusercontent.com/nomascus/ANT3814/main/FILES/nobody.txt

```
Site web utile pour tester les expressions régulières
https://regex101.com/


## Expressions régulières

Jusqu'à présent, nous avons introduit quelques façons de rechercher et de remplacer du texte en ligne de commande. La première était l'utilisation du caractère générique **\***, qui en ligne de commande signifie "un ou plusieurs de n'importe quel caractère". Par exemple,
```
[orkin@login1 3_Unix3]$ ls *.txt
file1.txt
file2.txt
file3.txt
file4.txt
```

Dans ce cas, nous demandons à ```ls``` de lister tous les fichiers se terminant par .txt. Nous pourrions faire de même  de cette manière

```
[orkin@login1 3_Unix3]$ ls fi*txt
file1.txt
file2.txt
file3.txt
file4.txt
```
Notez qu'il n'y a pas de différence dans le résultat si nous demandons n'importe quel fichier se terminant par '.txt' ou n'importe quel fichier commençant par 'fi'  et se terminant par 'xt'. Cependant, si nous recherchons des fichiers se terminant par '.tx', nous n'en trouverons aucun

```
[orkin@login1 ~]$ ls *tx
ls: *.tx: Aucun fichier ou dossier de ce type
```

Nous avons également introduit ```grep``` pour trouver des lignes correspondantes dans un fichier et le langage `sed`  qui nous permet d'éditer du texte en continu en utilisant des caractères de recherche et de remplacement

```
[orkin@login1 ~]$ echo "This is ANT3745" | sed 's/74/47/'
This is ANT3814
```

Cependant, il existe un moyen beaucoup plus puissant de manipuler du texte en ligne de commande et dans des fichiers en utilisant **des expressions régulières**, qui est un langage de correspondance de motifs qui utilise des caractères spécifiques ayant des significations secondaires. Nous avons déjà appris quelques-uns de ces caractères, comme **\n**, qui signifie nouvelle ligne et **\t** qui signifie tabulation. Dans ce cas, le backslash ```\``` est appelé le "caractère d'échappement", car il indique que le langage ne doit pas interpréter les deux caractères littéralement, mais plutôt que mettre un backslash devant le caractère bascule vers la signification secondaire du caractère (ou il échappe à la signification standard). Ces combinaisons de caractères sont appelées 'atomes', car elles peuvent être combinées en formes plus complexes. Tous les atomes d'expression régulière n'utilisent pas le caractère d'échappement, mais beaucoup le font.

### Atomes d'expressions régulières

|Type|Atome|Description|Exemple
|-|-|-|-|
Caractère|[ ]| l'un des caractères entre crochets| /[AtG]/
Caractère|[^]| PAS un caractère entre crochets|/[^AtG]/
Caractère|[a-z]|un lettre minuscule ordinaire|/GAT[a-z]ACA/
Caractère|[A-Z]|un lettre majuscule ordinaire|/GAT[A-Z]ACA/
Caractère|\d| un seul chiffre (digit) [0-9]| /file\d\.txt/
Caractère|\D| PAS un chiffre (digit) [^0-9]|/\D\D\D\D\d.txt/
Caractère|\w| un seul caractère de mot (word) [A-Za-z0-9_]|/file3\.\wxt/
Caractère|\W| PAS un caractère de mot (word) [^A-Za-z0-9_]|file\W\Wtxt
Caractère|\s| un caractère d'espace (space) blanc [ \n\t\r]| /Chr1\s536/
Caractère|\S| PAS d'espace (space) blanc [^ \n\t\r]| //Chr1\s\S36//
Caractère|.| n'importe quel caractère unique| /GA.TT/
Caractère|\\.| un seul point| /GA\\.TT/
Ancre|^| l'atome correspond au début de la chaîne de texte|/^>Fasta header/
Ancre|$| l'atome correspond à la fin de la chaîne de texte |/ATGACA$/ 

#### Exemples

Rappelez-vous que lors de la recherche, nous utiliserons la structure des deux barres obliques
 ```/ce que je cherche/```

 ```/p...e/```  correspondrait aux chaînes 'pomme', 'pTd9e' et 'pomme de terre ' .  
 
```/p[omme][omme][omme]e```  correspondrait à 'pomme', mais aussi à des chaînes comme 'pemme', 'pmmme' et 'pooom'

```/\w\w\w\d\d\d/``` correspondrait à ANT3475 et ANT2440

```/\d\d\d-d\d\d\d/``` correspondrait à un numéro de téléphone comme 555-1212. Il correspondrait également à un modèle comme celui-ci 5555-1212, mais pas à un modèle comme celui-ci 514-555-1212

```/^\d\d\d-d\d\d\d/``` correspondrait à 555-1212 ainsi qu'à 555-121212
 
```/^\d\d\d-d\d\d\d$/``` ne correspondrait qu'à un numéro de téléphone à sept chiffres comme 555-1212 (trois chiffres suivis d'un tiret et de quatre chiffres)

### Quantifieurs d'expressions régulières

Plutôt que de taper les mêmes caractères encore et encore, nous pouvons utiliser un groupe d'expressions régulières appelé **quantifieurs** qui indiquent le nombre de répétitions pour un caractère

|Type|Atome|Description|
|-|-|-|
Quantifieur|?| l'atome correspond zéro fois ou exactement une fois
Quantifieur|*| l'atome correspond zéro ou plusieurs fois
Quantifieur|+| l'atome correspond une ou plusieurs fois
Quantifieur|{3}| l'atome correspond à 3 fois
Quantifieur|{2,5}| l'atome entre 2 et 5 fois inclusivement
Quantifieur|{4,}| l'atome correspond à 4 fois ou plus 

```/pomm?e/``` correspondrait à pomme et pome
```/pommS?e/``` correspondrait à pomme ou pommSe
```/p.+e/```  correspondrait aux chaînes 'pomme', 'pluie' 'pTdse' et 'pomme de terre ' 
```/^p.+e$/```  correspondrait aux chaînes 'pomme', 'pluie' 'pTdse' et 'pomme de terre ' 
```/p.*e/```  correspondrait aux chaînes 'pomme', 'pluie', pTd9e', 'pe' et 'pomme de terre'   
``` /^\d{3}-\d{4}$/``` ne correspondrait qu'à un numéro à sept chiffres comme 555-1212 (trois chiffres suivis d'un tiret et de quatre chiffres)
 
 #### Exemple de sites de méthylation

Essayons quelque chose de plus compliqué et réaliste. Supposons que nous voulions rechercher un site de méthylation à l'intérieur d'une séquence nucléotidique. Pour ce faire, la séquence doit répondre aux critères suivants :

Commence par G ou A
suivi de C
suivi d'un seul nucléotide ou de rien
suivi de G

Cela pourrait correspondre à l'un des éléments suivants

GCAG
GCTG
GCGG
GCCG
GCG
ACAG
ACTG
ACGG
ACCG
ACG

Rechercher chacun de ces éléments individuellement serait très lent, mais avec des expressions régulières, nous pouvons rechercher n'importe lequel d'entre eux en une seule ligne de texte!

```/[GA]C[ATCG]?G/'```


### Soit... ou

Parfois, nous voulons rechercher l'un des deux motifs complexes. Nous pouvons le faire avec un caractère de pipe ```|``` et un ensemble de parenthèses ```( )```. Regardons quelques lignes de la chanson enfantine, "Alouette".

*alouette, gentille alouette
alouette, je te plumerai
je te plumerai la tête
je te plumerai le bec
je te plumerai les yeux
je te plumerai le cou
je te plumerai les ailes
je te plumerai les pattes
je te plumerai la queue
je te plumerai le dos*

Supposons que nous voulions rechercher les lignes qui contiennent 'la queue' ou 'le bec'.

```/je te plumerai (la queue|le bec)/```

Cela correspond à toute chaîne qui :
- contient les mots je te plumerai séparés par un espace simple
- est suivie de la queue ou du bec séparé par un seul espace

Que se passerait-il si nous excluions l'énoncé (soit|soit) ?
```/je te plumerai/```

Et si nous ajoutions un point d'interrogation à la fin ?

```/je te plumerai (la queue|le bec)?/```

Et si nous ajoutions un espace après la queue ?

```/je te plumerai (la queue |le bec)?/```

## Expressions régulières en ligne de commande avec perl

Essayons quelques-unes de ces expressions. Plutôt que d'utiliser ```sed```, qui est un peu limité dans sa fonctionnalité, nous pouvons utiliser un autre langage plus robuste appelé **perl**. Perl est un langage de codage extrêmement puissant qui était plus courant, mais qui est principalement remplacé par **python**. Néanmoins, perl reste le meilleur langage pour les expressions régulières et peut être facilement appelé depuis la ligne de commande. Si nous voulons transmettre du texte à perl pour utiliser des expressions régulières, tapez ```perl -lpe ' '```. Les options '-lpe' indiquent à perl que la commande entrante doit être exécutée ( e ) ligne par ligne ( l ) et imprimée ( p ). Attention cependant ! L'option ```-lpe``` imprimera toutes les lignes du fichier (même celles qui ne correspondent pas), en ne changeant que celles que vous lui indiquez. Une manière plus sûre est d'indiquer à perl d'utiliser l'option ```-lne``` qui n'imprime pas automatiquement toutes les lignes. Si vous voulez n'imprimer que les lignes qui correspondent aux conditions, vous pouvez dire à perl d'"imprimer si", en utilisant la structure suivante :

Pour la recherche :  ```perl -lne 'print if //' fichier.txt```
Pour la substitution :  ```perl -lne 'print if s///' fichier.txt```

Au lieu de mettre le fichier à la fin, vous pouvez également utiliser cat pour afficher le fichier ou echo pour transmettre le texte à perl de cette manière

Pour la correspondance :  ```cat fichier.txt | perl -lne 'print if //' ```
Pour la substitution :  ```echo 'voici un peu de texte' |perl -lne 'print if s///' ```
 
Commençons par une chaîne de nucléotides et recherchons quelques motifs. Tout comme sed, perl utilise des barres obliques pour séparer le texte à rechercher et à remplacer. 
 
```
[orkin@login1 ~]$ echo 'Ce texte est dans le mauvais ordre' | perl -lpe 's/mauvais/bon/ '
Ce texte est dans le bon ordre 
```

### Correspondances multiples

Parfois, nous avons plus d'une instance de texte correspondant dans une chaîne

```pomme raisin framboise pomme orange```

mais une expression régulière ne fonctionnera que sur la première correspondance qu'elle trouve

```
[orkin@login1 ~]$ echo 'pomme raisin framboise pomme orange' | perl -lpe 's/pomme/ananas/'
ananas raisin framboise pomme orange
```
si vous voulez agir sur toutes les correspondances dans une ligne, vous devez faire une recherche **globale** en ajoutant 'g' à la fin de l'expression. 

```
[orkin@login1 ~]$ echo 'pomme raisin framboise pomme orange' | perl -lpe 's/pomme/ananas/g'
ananas raisin framboise ananas orange
```
C'est particulièrement important lorsqu'on travaille avec un texte très répétitif, comme une séquence nucléotidique. Supposons que vous vouliez changer une séquence d'ADN en ARN en remplaçant tous les nucléotides de Thymine par de l'Uracile.

```
[orkin@login1 ~]$ echo 'ATCGTCGTATGCGTATGCTGATA' | perl -lpe 's/T/U/'
AUCGTCGTATGCGTATGCTGATA
```

Comme vous pouvez le voir, perl ne traite que le premier Thymine, ce qui n'est pas très utile. Essayez maintenant avec une recherche globale

```
[orkin@login1 ~]$ echo 'ATCGTCGTATGCGTATGCTGATA' | perl -lpe 's/T/U/g'
AUCGUCGUAUGCGUAUGCUGAUA
```

### Sous-patterns

Les expressions régulières permettent également de capturer des sous-patterns dans le texte et de les déplacer en divisant le texte en parties entre parenthèses. Par exemple, si nous voulions changer

```Ce ordre est dans le bon texte``` en ```Ce texte est dans le bon ordre```

nous commençons par diviser en quatre parties, qui correspondent à des variables scalaires numérotées en perl

|Texte|Variable| expression régulière
|-|-|-|
|Ce|$1|^(\w+\s)| 
|ordre|$2|(ordre)|
|est dans le bon|$3|(.+)
|texte|$4|(texte)$|

Ensuite, nous utilisons une substitution ``` s/// ``` pour réorganiser les variables
```
echo 'Ce ordre est dans le bon texte' | perl -lpe 's/^(\w+ )(ordre)(.*)(texte)/$1$4$3$2/'
Ce texte est dans le bon ordre
```

Ou supposons que nous ayons un fichier fasta avec un en-tête très désordonné et que nous voulions le réorganiser et le simplifier

```
>NC_19939459.0 Homo_sapiens_A87w889 FOXP2 ISOFORM_X1.9982
ATGTCGATGCTGACTGATGCTACTAGCTGGAGGGATGTGATATTACGGCTAGTTGCGCGCCCTGTGCATGAT
```

Changeons l'en-tête en ```>FOXP2 Homo_sapiens```

```
[orkin@login1 ~]$ perl -lpe 's/.*(Homo_sapiens).*(FOXP2).*/>$2 $1/' FOXP2.fasta
>FOXP2 Homo_sapiens
ATGTCGATGCTGACTGATGCTACTAGCTGGAGGGATGTGATATTACGGCTAGTTGCGCGCCCTGTGCATGAT
```

Et si nous avions un fichier fasta avec des séquences de plusieurs espèces dedans?

```
>NC_19939459.0 Homo_sapiens_A87w889 FOXP2 ISOFORM_X1.9982
ATGTCGTTGCTGACTGATGCTACTAGCTGGAGGGATGTGATATTACGGCTAGTTGCGCGCCCTGTGCATGAT
>NC_2943459.0 Pan_troglodytes_A4 FOXP2
ATGTCGATGCTGACTGATGCTACTAGCTGGAGGGATGTGATATTACGGCTAGTTGCGCGCCCTGTGCATGAT
>NC_5559.0 Gorilla_gorilla_K49 FOXP2 ISOFORM_x2s
TTGTCGATGCTGACTGATGCTACTAGCTGGAGGGATGTGATATTACGGCTAGTTGCGCGCCCTGTGCATGAT
>KA_ds39459.0 Nomascus_concolor_A87w889 FOXP2 isoform 1
CTGTGGATGCCGACTGATGCTACTAGCTGGAGGGATGTGATATTACGGCTAGTTGCGCGCCCTGTGCATGAT
```
Nous devrions utiliser une expression régulière plus générique
```
[orkin@login1 ~]$ perl -lpe 's/>.* (\w+_\w+) (FOXP2) .*/>$2 $1/' FOXP2_long.fasta
>FOXP2 Homo_sapiens_A87w889
ATGTCGTTGCTGACTGATGCTACTAGCTGGAGGGATGTGATATTACGGCTAGTTGCGCGCCCTGTGCATGAT
>FOXP2 Pan_troglodytes_A4
ATGTCGATGCTGACTGATGCTACTAGCTGGAGGGATGTGATATTACGGCTAGTTGCGCGCCCTGTGCATGAT
>FOXP2 Gorilla_gorilla_K49
TTGTCGATGCTGACTGATGCTACTAGCTGGAGGGATGTGATATTACGGCTAGTTGCGCGCCCTGTGCATGAT
>FOXP2 Nomascus_concolor_A87w889
CTGTGGATGCCGACTGATGCTACTAGCTGGAGGGATGTGATATTACGGCTAGTTGCGCGCCCTGTGCATGAT
```

Cela nous rapproche assez, mais que se passe-t-il si nous voulons supprimer le texte après le trait de soulignement supplémentaire?

```
[orkin@login1 ~]$ perl -lpe 's/>.* (\w+_\w+)_\w+ (FOXP2) .*/>$2 $1/' FOXP2_long.fasta
>FOXP2 Homo_sapiens
ATGTCGTTGCTGACTGATGCTACTAGCTGGAGGGATGTGATATTACGGCTAGTTGCGCGCCCTGTGCATGAT
>FOXP2 Pan_troglodytes
ATGTCGATGCTGACTGATGCTACTAGCTGGAGGGATGTGATATTACGGCTAGTTGCGCGCCCTGTGCATGAT
>FOXP2 Gorilla_gorilla
TTGTCGATGCTGACTGATGCTACTAGCTGGAGGGATGTGATATTACGGCTAGTTGCGCGCCCTGTGCATGAT
>FOXP2 Nomascus_concolor
CTGTGGATGCCGACTGATGCTACTAGCTGGAGGGATGTGATATTACGGCTAGTTGCGCGCCCTGTGCATGAT
```


