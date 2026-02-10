
# ANT3475 : Bioinformatique évolutive humaine

## Lecture 1 : La ligne de commande Unix

Le matériel de cours d’aujourd’hui est modifié et traduit en français à partir de : 

Programming for Biology [programmingforbiology.org](http://programmingforbiology.org) et a été initialement rédigé par Simon Prochnik et Sofia Robb. Une grande partie de ce contenu est disponible gratuitement en ligne. Il est également proposé sous forme d’un excellent [cours intensif de 2 semaines]([https://meetings.cshl.edu/courses.aspx?course=c-info&year=24]) au Cold Spring Harbor Laboratory pour apprendre Python dans le cadre de la bio-informatique.

# Grande idée

## Pourquoi ?

Pourquoi est-il important pour **les anthropologues biologistes** d’apprendre à programmer ?

Vous avez peut-être déjà une idée de la réponse à cette question puisque vous êtes ici.

Nous croyons fermement que savoir programmer est tout aussi essentiel que savoir réaliser une électrophorèse sur gel ou configurer une réaction PCR. Les données que nous obtenons aujourd’hui à partir d’une seule expérience peuvent être considérables. Souvent, ces données doivent être reformattées, filtrées et analysées d’une manière spécifique. La programmation vous permet de réaliser ces tâches de manière **efficace** et **reproductible**.

Quels sont nos conseils pour réussir dans ce cours de programmation ?

1. Pratiquez, pratiquez, pratiquez. Passez autant de temps que possible à coder réellement.
2. Écrivez seulement une ou deux lignes de code, puis testez-les. Si vous écrivez trop de lignes d’un coup, il sera plus difficile de déboguer en cas d’erreur.
3. Les erreurs ne sont pas des échecs. Chaque erreur est une occasion d’apprendre. Chaque erreur que vous corrigez est un grand succès. Corriger des erreurs vous aidera à consolider vos connaissances.
4. Ne passez pas trop de temps à essayer de résoudre un problème tout seul. Bien qu’il soit formateur de tenter de le résoudre par vous-même, il peut être frustrant de rester bloqué pendant longtemps. Je suis là pour vous aider, alors n’hésitez pas à demander.
5. Les cours magistraux sont importants, mais la pratique est plus importante encore.
6. Notre objectif principal est de vous apprendre progressivement à résoudre les problèmes par vous-même.

---

# Unix (linux, mac OS)

### Aperçu d’Unix

#### Qu’est-ce que la ligne de commande ?

La plupart des outils de bio-informatique sont conçus pour fonctionner dans la ligne de commande, sans interface graphique. Bien souvent, un outil en ligne de commande est plus polyvalent qu’un outil graphique, car vous pouvez aisément combiner plusieurs outils en scripts automatisés qui accomplissent des tâches personnalisées de manière séquentielle, sans intervention humaine. Presque tout le travail en bio-informatique se fait sur un système d’exploitation de type UNIX comme mac OSX ou linux.

Derrière l’interface graphique (GUI) de Mac OSX se cache le système d’exploitation (OS). Il est basé sur une version d’Unix. Linux est assez similaire et constitue également un système d’exploitation très répandu en bio-informatique, avec des variantes comme Red Hat, Ubuntu et autres.

La ligne de commande vous donne accès au cœur du système d’exploitation et facilite l’écriture de logiciels et de scripts personnalisés.

Dans ce cours, nous allons utiliser Linux et R pour écrire des scripts et les exécuter depuis la ligne de commande.

Le système Unix est complexe et il y a toujours plus à apprendre. Nous n’aborderons ici que l’essentiel pour démarrer.

### Les bases

#### Se connecter à votre ordinateur

Presque tout notre travail dans ce cours sera effectué sur un ordinateur virtuel dans le cloud. Cependant, nous devons y accéder d’une manière particulière. Si vous utilisez un ordinateur Mac ou Linux, c’est très simple : il vous suffit d’ouvrir un programme appelé **Terminal**. Sous Windows, les choses sont un peu plus compliquées, car ce système n’est pas basé sur Unix. Nous utiliserons un programme appelé **PowerShell**. Une version devrait déjà être installée, mais vous souhaiterez peut-être passer à la version 7.

##### Installation de Powershell
1. Ouvrez le **Microsoft Store**  
2. Recherchez **PowerShell**  
3. Installez  

#### Créer un compte sur notre serveur linux dans le cloud

1. Rendez-vous sur ce site pour créer un compte. Veuillez choisir un nom d’utilisateur court, qui soit simplement votre prénom. Assurez-vous de retenir votre mot de passe !  
   https://jupyter.ant3475.calculquebec.cloud/  
2. Attendez que je vous donne l’accès au serveur  

#### Faire apparaître la ligne de commande

Pour faire apparaître la ligne de commande sur Windows, ouvrez **PowerShell**. Cela ouvrira une fenêtre semblable à :

![enter image description here](https://github.com/nomascus/ANT3814/blob/main/FILES/windowsTerminal.png?raw=true)

Pour faire apparaître la ligne de commande sur Mac, ouvrez **Terminal** depuis le Launchpad. Cela ouvrira une fenêtre semblable à :

![enter image description here](https://github.com/nomascus/ANT3814/blob/main/FILES/macTerminal.png?raw=true)

Pour faire apparaître la ligne de commande sur Linux, utilisez le **Finder** pour naviguer jusqu’à _Applications->Utilities_ et double-cliquez sur l’application **Terminal**. Cela ouvrira une fenêtre similaire à :

![enter image description here](https://github.com/nomascus/ANT3814/blob/main/FILES/UbuntuTerminal.png?raw=true)

Vous pouvez ouvrir plusieurs fenêtres de Terminal ou PowerShell en même temps. Cela s’avère souvent pratique. Une façon courante de travailler est d’ouvrir votre éditeur de code dans une fenêtre, d’exécuter votre code dans une deuxième fenêtre et d’observer vos données et vos fichiers dans une troisième.

Je vous conseille de placer l’icône de Terminal dans la barre de raccourcis en bas de votre écran, puisque vous l’utiliserez fréquemment.

#### Se connecter à notre serveur dans le cloud

Pour accéder au serveur, vous devez taper la commande **ssh**, qui signifie “secure shell”. Elle ouvre un tunnel chiffré vers l’ordinateur distant dans le cloud.


```
ssh userid@ant3475.calculquebec.cloud
```
Par exemple, pour mon nom d’utilisateur :

```
ssh orkin@ant3475.calculquebec.cloud
```

Appuyez ensuite sur la touche Enter.  

Si le terminal vous demande de confirmer une empreinte de sécurité, tapez **yes** et appuyez sur Enter. Vous ne devriez avoir à accorder l'autorisation que la première fois.

Ensuite, entrez le mot de passe que vous avez défini en créant votre compte, puis appuyez sur enter. Vous ne verrez pas de caractères s’afficher pendant que vous saisissez votre mot de passe, ce qui est une mesure de sécurité. Votre écran de terminal devrait ressembler à ceci :

![enter image description here](https://github.com/nomascus/ANT3814/blob/main/FILES/cloudTerminal.png?raw=true)

Vous êtes maintenant connecté à un cluster de serveurs Linux dans le cloud. Bravo !

#### D’accord, je suis connecté… et maintenant ?

La fenêtre de terminal exécute un **shell** appelé « bash ». Le shell est une boucle qui :
1. Affiche un invite de commande (prompt)  
2. Lit une ligne de texte que vous tapez au clavier  
3. Décompose cette ligne en une ou plusieurs commandes  
4. Exécute ces commandes (qui envoient généralement une sortie à l’écran)  
5. Retourne à l’étape 1  

Il existe plusieurs shells aux noms étranges comme **bash**, **sh**, **csh**, **tcsh**, **ksh**, **zsh**. La partie « sh » signifie shell. Chaque shell a sa propre syntaxe et ses particularités. **bash** est standard dans les systèmes Linux, mais Mac OS est récemment passé de bash à **zsh**, bien qu’ils soient très similaires.

#### L’invite de commande (prompt)

La plupart des tâches en bio-informatique se font en lançant des logiciels dans un shell, donc il est utile de prendre un peu de temps pour bien maîtriser le shell.

Voici un exemple d’invite de commande :

```
[orkin@login1] ~$
```

En voici un autre :

```
$ 1:12PM>
```

Et encore un autre :


```
[CharlesDarwin@narval3 projects]$
```

Ce que vous voyez dépend de la configuration faite par l’administrateur du système ou de vos propres paramètres (mais nous n’aborderons pas ce dernier point ici).

L’invite de commande vous indique que le shell est prêt à recevoir vos instructions. La plupart des commandes s’exécutent presque instantanément. Pour une commande plus longue, l’invite ne réapparaîtra que lorsqu’elle sera terminée et que le système sera prêt à accepter votre prochaine requête.

#### Exécution de commandes

Vous tapez une commande, puis appuyez sur **Enter** pour l’exécuter. Si la commande produit une sortie, elle s’affiche à l’écran. Exemple :


```
[orkin@login1] ~$ ls 
[orkin@login1] ~$ projects scratch
```

Ici, la commande est `ls`, qui affiche la liste des fichiers et répertoires dans le répertoire courant (nous verrons cela sous peu). Sous la liste, l’invite de commande réapparaît.

Certaines commandes prennent plus de temps à s’exécuter. Tant que la commande n’est pas terminée, l’invite de commande ne réapparaît pas. Vous pouvez continuer à travailler en parallèle, soit en ouvrant un nouveau shell (dans le menu File de Terminal), soit en exécutant la commande en arrière-plan avec un **&** :


```
[orkin@login1] ~$ long_running_application &
[orkin@login1] ~$
```
La commande s’exécutera en arrière-plan. Si elle produit une sortie, celle-ci s’affichera dans la fenêtre du terminal. Vous pouvez vouloir enregistrer cette sortie dans un fichier (on parlera plus tard de la redirection).

#### Édition de la ligne de commande

La plupart des shells offrent la possibilité d’éditer la ligne de commande. Tant que vous n’avez pas appuyé sur **Enter**, vous pouvez revenir sur la ligne et la modifier grâce au clavier. Voici quelques raccourcis utiles :

- **Backspace** : Supprime le caractère précédent et recule d’un espace.  
- **Flèches gauche/droite** : Déplacent le curseur d’un caractère à gauche ou à droite.  
- **control-a (^a)** : Déplace le curseur au début de la ligne. (Astuce : A est la première lettre de l’alphabet)  
- **control-e (^e)** : Déplace le curseur à la fin de la ligne. (Astuce : E pour End)  
- **control-d (^d)** : Supprime le caractère situé sous le curseur. D = Delete.  
- **control-k (^k)** : Supprime de la position du curseur jusqu’à la fin de la ligne. k = kill. (Le texte est placé dans un “kill buffer”, équivalent à un « couper »)  
- **control-y (^y)** : Colle le contenu du kill buffer à la position du curseur. y = yank. (Équivalent à « coller »)  
- **Flèches haut/bas** : Défilent dans l’historique des commandes. Vous pouvez ainsi réexécuter une commande précédente et éventuellement la modifier avant de l’exécuter.  
- **control-c (^c)** :  Envoyez un signal d'interruption à la commande pour l'arrêter. Pensez à C pour « cancel » (annuler).
- 
Il existe aussi quelques commandes de shell utiles :

- `history` : Affiche l’historique des commandes récentes, numérotées.  
- `!<number>` : Réexécute une ancienne commande en se basant sur son numéro (affiché par `history`).  
- `!!` : Réexécute la dernière commande.  

**bash** propose l’auto-complétion et la correction orthographique. Si vous tapez une partie d’une commande et appuyez sur la touche **Tab**, il vous affichera toutes les complétions possibles. Exemple :


```
[orkin@login1] ~$ le< less lessecho lessfile lesskey lesspipe let lexgrog
[orkin@login1] ~$ le
```
 Si vous appuyez sur **Tab** après avoir partiellement saisi une commande, mais avant **Enter**, **bash** vous proposera aussi la liste des fichiers existants, car de nombreuses commandes opèrent sur des fichiers.

#### Wildcards (caractères génériques)

Vous pouvez utiliser des caractères génériques pour désigner des fichiers. `*` remplace zéro ou plusieurs caractères, `?` remplace n’importe quel caractère unique. Par exemple, pour lister tous les fichiers ayant l’extension “.txt”, exécutez `ls` avec la forme générique `*.txt` :


```
[orkin@login1] ~$ ls *.txt 
data.txt readme.txt results.txt transitions.txt
[orkin@login1] ~$
```

Il existe également des formes plus avancées de motifs génériques que vous pouvez découvrir dans la page de manuel de **bash**. Par exemple, pour faire correspondre les fichiers commençant par “f” ou “g” et se terminant par “.txt”, vous pouvez mettre les caractères entre crochets `[fg]`. Exemple :


```
[orkin@login1] ~$ ls [rd]*.txt data.txt readme.txt results.txt 
[orkin@login1] ~$
```
Vous pouvez aussi indiquer des intervalles de caractères comme `[a-e]` ou `[1-5]`.

### Les répertoires et leur organisation

Un ordinateur est équipé d’un espace de stockage pour les scripts, données, images, le système d’exploitation, etc. Historiquement, on a eu des disquettes, puis des disques durs, et maintenant des SSD (Solid State Drives). Parlons de l’organisation de ce stockage afin de vous aider à localiser vos travaux. Les répertoires (ou dossiers) se créent à l’intérieur d’autres répertoires. Un répertoire est spécial : le **répertoire racine** (root directory), noté `/`. Les fichiers associés sont placés à l’intérieur de répertoires pour rester organisés, créant ainsi une structure en arborescence. On peut la représenter comme un arbre à l’envers.

**Exemple de structure de répertoires Linux depuis la racine :**

![enter image description here](https://github.com/nomascus/ANT3814/blob/main/FILES/directory_structure.png?raw=true)

#### Home sweet home (votre répertoire personnel `~`)

Lorsque vous vous connectez, vous vous trouvez dans votre répertoire personnel (home directory), généralement dans `/home` (sur mac OS, c’est dans `/Users`). Si votre identifiant est `jdesjardins`, votre répertoire s’appelle alors `/home/jdesjardins`. On appelle cela un chemin d’accès ou **path**. Unix est rempli d’abréviations pour vous éviter de trop taper. Vous pouvez ainsi raccourcir ce chemin en `~`. Le caractère `~` s’appelle “tilde” et se trouve en haut à gauche de votre clavier.

Dans votre répertoire personnel, vous pouvez sauvegarder, supprimer, ouvrir et déplacer des fichiers et d’autres répertoires. En général, vous pouvez lire, mais pas modifier, les fichiers situés en dehors de ce répertoire.

Pour voir le contenu de votre répertoire personnel, utilisez la commande `ls` :

![enter image description here](https://github.com/nomascus/ANT3814/blob/main/FILES/ls-home.png?raw=true)

Vous pouvez modifier le comportement d’une commande avec des options, par exemple `ls -aF` (l’option `a` pour tout lister, même les fichiers cachés qui commencent par un `.` et l’option `F` pour “decorated listing”). Cela ajoute par exemple un `/` à la fin des noms de répertoires :

![enter image description here](https://github.com/nomascus/ANT3814/blob/main/FILES/ls-aF.png?raw=true)

Ne supprimez pas les fichiers cachés. Beaucoup sont cruciaux pour la configuration de vos commandes et autres programmes. Par exemple, le fichier `.bashrc` renferme la configuration pour le shell **bash**. Vous pourrez l’éditer plus tard (une fois que vous saurez ce que vous faites) pour changer l’invite de commande ou le chemin de recherche des commandes.

#### Changer de répertoire

Vous pouvez vous déplacer d’un répertoire à l’autre grâce à la commande `cd` (change directory). Indiquez le nom du répertoire où vous souhaitez aller, ou laissez vide pour retourner au répertoire personnel. Utilisez `pwd` (print working directory) pour voir où vous êtes (ou fiez-vous à l’invite de commande, si configurée) :


```
[orkin@login1] $ cd projects
[orkin@login1 projects]$ pwd /home/orkin/projects
[orkin@login1 projects]$ cd
[orkin@login1 ~] $ cd /
[orkin@login1 /] $ ls -aF
```

![enter image description here](https://github.com/nomascus/ANT3814/blob/main/FILES/ls-root.png?raw=true)

Chaque répertoire contient deux répertoires cachés, `.` et `..`. `.` désigne toujours le répertoire courant, `..` son répertoire parent. Cela vous permet par exemple de remonter d’un niveau :


```
[orkin@login1 scratch]$ cd ..
```
ou de faire des déplacements plus complexes :


```
[orkin@login1 scratch] cd ../../orkin/scratch
```
Cette commande remonte de deux niveaux, puis descend dans un répertoire appelé `orkin` à l’intérieur d’un répertoire appelé `scratch`.

Si vous êtes perdu, la commande `pwd` vous indique le chemin complet de votre répertoire actuel :


```
[orkin@login1 scratch]$ pwd
/home/orkin/scratch
```

#### Chemins absolus et chemins relatifs

Vous pouvez spécifier un chemin de deux façons : un chemin absolu commence toujours par `/` et est lu à partir de la racine, tandis qu’un chemin relatif commence par n’importe quel autre caractère et est lu à partir de votre position courante. Voici deux façons d’accéder au répertoire `projects` dans le répertoire personnel d’`orkin`. En partant de `scratch`, on peut remonter d’un niveau, puis entrer dans `projects` (chemin relatif), ou bien partir de la racine `/`, puis `home`, `orkin` et `projects` (chemin absolu).


```
[orkin@login1 scratch] cd ../project  
[orkin@login1 project] $ cd
[orkin@login1 /] $ cd /home/orkin/project 
[orkin@login1 project] $
```

### Commandes Unix essentielles

À l’exception de certaines commandes intégrées dans le shell (par ex. `history`, `for`, `if`), toutes les commandes Unix sont de véritables programmes exécutables. Lorsque vous tapez le nom d’une commande, le shell cherche dans tous les répertoires listés dans la variable d’environnement PATH un exécutable de ce nom. S’il le trouve, il exécute la commande. Sinon, il vous indique qu’il ne l’a pas trouvée.

La plupart des commandes résident dans `/bin`, `/sbin`, `/usr/bin`, `/usr/sbin/` ou `/usr/local/bin`.

#### Obtenir des informations sur les commandes

La commande `man` vous donne un résumé de la commande. Par exemple :


```
[orkin@login1] $ man wc Formatting page, please wait... WC(1) WC(1)

NAME wc - print the number of bytes, words, and lines in files ...
```

#### Arguments et options en ligne de commande

De nombreuses commandes acceptent des arguments, souvent le nom de fichiers sur lesquels opérer. La plupart acceptent également des **options** (appelées “switches”) pour ajuster leur comportement. Certaines commandes reconnaissent les options “courtes” (-c, -l, etc.), d’autres reconnaissent des options “longues” (-- chars, --lines, etc.).  

Le programme `wc` (word count) est un exemple : il reconnaît `-c`, `-w` et `-l` (ou `--chars`, `--words`, `--lines`) afin de compter les caractères, les mots et les lignes d’un fichier texte. Les deux exemples ci-dessous comptent le nombre de caractères et de lignes du fichier `/var/log/README` :


```
[orkin@login1 ~] $ wc -c -l /var/log/README 25 150 /var/log/README
[orkin@login1 ~] $ wc --chars --lines /var/log/messages 25 150 /var/log/README
```

Vous pouvez même regrouper plusieurs options courtes :


```
[orkin@login1 ~] $ wc -cl 
/var/log/README 25 150 /var/log/README
```

Beaucoup de commandes affichent un résumé rapide lorsque vous les appelez avec `-h` ou `--help`.

#### Espaces, tabulations et caractères de fin de ligne

Le shell utilise l’espace pour séparer les arguments. Si vous voulez inclure un espace (ou autre caractère blanc, tabulation, nouvelle ligne, etc.) dans un argument, placez ce texte entre guillemets simples. Exemple :


```
mail -s 'An important message' 'Bob Ghost [bob@ghost.org](mailto:bob@ghost.org)'
```
Ici, `-s` signifie “subject” pour le courriel. Comme le sujet contient des espaces, nous plaçons cette portion entre guillemets. De même pour le nom et l’adresse, qui contiennent aussi des espaces.

Certains caractères non-imprimables ont des **codes d’échappement** :

| Code d’échappement | Description                                                  |
| ------------------ | ------------------------------------------------------------ |
| \n                 | Nouvelle ligne                                              |
| \t                 | Tabulation                                                  |
| \r                 | Retour chariot                                              |
| \a                 | Caractère de cloche (ding !)                                (https://www.rapidtables.com/code/text/ascii-table.html) |

### Commandes utiles

Voici quelques commandes couramment utilisées. Utilisez `man` pour en savoir plus. Certaines pourraient être utiles pour les exercices ;-)

#### Manipuler les répertoires

| Commande  | Description                                                                 |
| --------- | --------------------------------------------------------------------------- |
| `ls`      | Affiche le contenu d’un répertoire. Souvent utilisé sous la forme `ls -F`, `ls -l` ou `ls -a`. |
| `mv`      | Renomme ou déplace un fichier ou un répertoire.                            |
| `cp`      | Copie un fichier. `cp -r` (récursif) pour copier des répertoires.          |
| `rm`      | Supprime (efface) un fichier. Exemple : `rm -rf olddata.tsv`               |
| `mkdir`   | Crée un répertoire.                                                        |
| `ln`      | Crée un lien symbolique ou dur. `ln -s` crée un lien symbolique.           |
| `chmod`   | Modifie les permissions d’un fichier ou d’un répertoire.                   |

#### Lire des fichiers

| Commande | Description                                                                                 |
| -------- | ------------------------------------------------------------------------------------------- |
| `cat`    | Programme pour concaténer. Souvent employé pour voir le contenu d’un fichier dans le terminal. |
| `less`   | Fait défiler un fichier page par page. Très utile pour les gros fichiers inaccessibles à un éditeur de texte. |
| `head`   | Affiche les premières lignes d’un fichier. `head -3` affiche par exemple les trois premières lignes. |
| `tail`   | Affiche la fin d’un fichier. `tail -f` permet de suivre les nouvelles lignes ajoutées en temps réel. |

#### Analyser et traiter des fichiers

| Commande         | Description                                                                      |
| ---------------- | -------------------------------------------------------------------------------- |
| `wc`             | Compte le nombre de mots, de lignes et/ou de caractères dans un ou plusieurs fichiers. |
| `sort`           | Trie les lignes d’un fichier par ordre alphabétique ou numérique (`-g` ou `-n`), éventuellement en ordre inverse (`-r`). |
| `uniq`           | Supprime les lignes dupliquées dans un fichier.                                  |
| `cut`            | Supprime des colonnes dans chaque ligne d’un fichier ou de plusieurs fichiers.    |
| `grep`           | Filtre un fichier pour n’afficher que les lignes correspondant à un motif. Peut s’inverser (`-v`) pour les lignes ne correspondant pas au motif. |
| `gzip` (`gunzip`) | Compresse (ou décompresse) un fichier.                                          |
| `tar`            | Archive ou désarchive un répertoire entier dans un fichier unique.                |

#### Éditer des fichiers

| Commande  | Description                                                           |
| --------- | --------------------------------------------------------------------- |
| `tr`      | Remplace un caractère par un autre. Peut aussi supprimer des caractères. |
| `vi`      | Lance l’éditeur de texte vi (pratique pour les experts).             |
| `echo`    | Affiche du texte à l’écran. Exemple : `echo 'Hello World!'`          |

#### Se connecter à d’autres ordinateurs

| Commande                | Description                                                                 |
| ----------------------- | --------------------------------------------------------------------------- |
| `ssh`                  | Permet de se connecter de façon sécurisée (chiffrée) à une machine distante. |
| `scp`                  | Copie (sécurisée) de fichiers vers/depuis des machines distantes.            |

#### Standard I/O et redirection

Les commandes Unix envoient leurs résultats vers la sortie standard (l’écran) et reçoivent leurs entrées depuis l’entrée standard (clavier).

Chaque programme Unix démarre avec trois connexions au monde extérieur, appelées **flux** (streams) :

| Nom du flux        | Description                                                                                      |
| ------------------ | ------------------------------------------------------------------------------------------------ |
| standard input     | Correspond au clavier (ce que vous tapez).                                                       |
| standard output    | Correspond à l’écran (ce que le programme affiche normalement).                                  |
| standard error     | Correspond aussi à l’écran, mais séparé pour afficher les messages d’erreur.                     |

Le mot “initialement” implique qu’on peut en fait **rediriger** ces flux vers un fichier, un autre périphérique ou un autre programme. C’est un mécanisme très pratique.

#### Exemple simple

Le programme `wc` compte les lignes, les caractères et les mots envoyés à son entrée standard. On peut l’utiliser de manière interactive :


```

[orkin@login1 ~] $ wc Alouette, gentille alouette,  
Alouette, je te plumerai
Je te plumerai la tete 
Je te plumerai la tete 
Alouette, alouette 
A a a ah ^d 
6 23 126

```

Dans cet exemple, j’ai lancé la commande `wc`. Elle attend que je tape un texte (ici quelques vers). Lorsque j’ai fini, j’appuie sur **Control-D** pour indiquer la fin (End-of-File). `wc` m’affiche alors le nombre de lignes (6), de mots (23) et de caractères (126) du texte.

On veut souvent compter les lignes d’un gros fichier et sauvegarder le résultat ailleurs. D’abord, on envoie le résultat de `wc` vers stdout :

```
[orkin@login1 FILES] wc big_file.fasta
82 103 5024 big_file.fasta
```

Puis on **redirige** cette sortie vers un fichier avec `>` :


```
[orkin@login1 files] wc big_file.fasta > count.txt
```

Si vous regardez le contenu de **count.txt** avec `cat`, vous verrez les données que vous venez d’enregistrer :


```
[orkin@login1 FILES] $ cat count.txt
82 103 5024
```

#### Caractères spéciaux de redirection

Voici la liste de base pour **bash** :

| Commande de redirection | Description                                                                                          |
| ----------------------- | ---------------------------------------------------------------------------------------------------- |
| `< myfile.txt`          | Redirige le contenu du fichier vers l’entrée standard.                                              |
| `> myfile.txt`          | Redirige la sortie standard vers un fichier.                                                         |
| `>> logfile.txt`        | Ajoute la sortie standard à la fin du fichier.                                                      |
| `1> myfile.txt`         | Redirige uniquement la sortie standard vers le fichier (même chose que `>`).                         |
| `2> myfile.txt`         | Redirige uniquement la sortie d’erreur.                                                              |
| `&> myfile.txt`         | Redirige à la fois la sortie standard et la sortie d’erreur vers le fichier.                        |

On peut les combiner. Par exemple :


```
[orkin@login1 FILES] $ grep root /etc/passwd > search.out 2> search.err
```
Ici, on recherche les lignes contenant “root” dans le fichier `/etc/passwd`, on envoie les résultats dans **search.out** et les erreurs éventuelles dans **search.err**.

#### Filtres, noms de fichiers et entrée standard

De nombreuses commandes Unix agissent comme des **filtres** : elles lisent des données à partir d’un fichier ou de l’entrée standard, les transforment et écrivent le résultat sur la sortie standard. Le plus souvent, si elles reçoivent un nom de fichier en argument, elles l’utiliseront comme entrée ; sinon, elles liront l’entrée standard. Exemple, ces deux commandes sont équivalentes :


```
[orkin@login1 FILES] grep 'Colobus' big_file.fasta 
[orkin@login FILES1] grep 'Colobus' < big_file.fasta
```

La première recherche “Colobus” directement dans le fichier **big_file.fasta**, la seconde lit l’entrée standard redirigée depuis **big_file.fasta**.

Parfois, vous souhaitez qu’un filtre agisse à la fois sur plusieurs fichiers et sur l’entrée standard. Beaucoup de commandes acceptent `-` comme synonyme de l’entrée standard :


```
[orkin@login1 FILES] grep 'Colobus' big_file.fasta bigger_file.fasta -
```
Ceci recherche “Colobus” dans les fichiers **big_file.fasta**, **bigger_file.fasta** et enfin dans l’entrée standard (ici le clavier, puisque non redirigée).

#### Standard I/O et pipelines

Le grand avantage du shell Unix est de pouvoir enchaîner plusieurs commandes avec des **pipes**. Exemple :


```
[orkin@login1 FILES] grep 'Colobus' big_file.fasta | wc -l 
3
```
Ici, `grep` recherche les lignes contenant “Colobus” dans **big_file.fasta**, et `wc -l` compte le nombre de lignes qui lui sont passées. Le symbole `|` (pipe) relie la sortie standard de `grep` à l’entrée standard de `wc`.

#### Autres idiomes de pipe

Les pipes sont très puissants. Voici des exemples courants :

**Compter le nombre de lignes où un motif n’apparaît pas**  

```
[orkin@login1 FILES] grep -v "Colobus" big_file.fasta | wc -l 
79
```

**Extraire les lignes uniques dans un fichier**  

```
[orkin@login1 FILES] sort long_file.txt | uniq > unique.out
```

**Concaténer plusieurs listes et supprimer les doublons**  

```
[orkin@login1 FILES] cat file1 file2 file3 file4 | sort | uniq
```

**Compter les lignes uniques**  

```
[orkin@login1 FILES]$ sort long_file.txt | uniq | wc -l
```
Pour compter le nombre d’occurrences de chaque ligne dans plusieurs fichiers :  

```
[orkin@login1 FILES]$ sort file1.txt file2.txt file3.txt file4.txt |uniq -c
```

**Faire un listing de répertoire très long et le défiler**  

```
[orkin@login1 FILES]$ less bigger_file.fasta 
```
Pour des lignes très longues, ajoutez `-S` :


```
[orkin@login1 FILES]$ less -S bigger_file.fasta
```

### Permissions

Souvent, plusieurs utilisateurs partagent le même ordinateur. Afin d’éviter qu’ils ne se gênent mutuellement, chaque fichier ou répertoire possède un ensemble de permissions (lire, écrire, exécuter) pour le propriétaire, le groupe et les autres. Pour les visualiser :


```
[orkin@login1 scratch]$ ls -l
total 8
drwxr-x---. 2 orkin orkin 4096 Jan 15 19:54 FILES
drwxr-x---. 2 orkin orkin  6 Jan 15 16:49 MyCode
drwxr-x---. 2 orkin orkin  6 Jan 15 19:14 New_dossier
drwxr-x---. 3 orkin orkin 19 Jan 15 16:48 RawData
-rw-r-----. 1 orkin orkin 33 Jan 15 16:57 count.txt
-rw-r-----. 1 orkin orkin  0 Jan 15 15:16 data.txt
-rw-r-----. 1 orkin orkin  0 Jan 15 15:16 readme.txt
-rw-r-----. 1 orkin orkin  0 Jan 15 15:16 results.txt
drwxr-x---. 2 orkin orkin 32 Jan 15 14:10 tmp
-rw-r-----. 1 orkin orkin  0 Jan 15 15:16 transitions.txt
```

Le premier caractère indique le type (`d` pour un répertoire, `l` pour un lien, `-` pour un fichier). Ensuite, il y a trois groupes de trois caractères, indiquant les permissions pour l’utilisateur, le groupe et les autres (par exemple, `rwx`). Par exemple :


```
Usr Grp Autr d rwx rwx rwx
```

Les utilisateurs appartiennent à des groupes, et “autres” désigne tous les autres utilisateurs. 

La commande `chmod` (change mode) sert à modifier les permissions. Vous pouvez ajouter `+` ou retirer `-` des permissions pour l’utilisateur `u`, le groupe `g`, les autres `o` ou tout le monde `a`.

Pour rendre **print.py** exécutable comme un programme :


```
[orkin@login1 FILES] $ ls -l print.py 
-rw-r--r-- 1 orkin staff 24 Jan 9 04:22 print.py 
[orkin@login1 FILES] $ chmod a+x print.py
[orkin@login1 FILES] $ ls -l print.py 
-rwxr-xr-x 1 orkin staff 24 Jan 9 04:22 print.py
[orkin@login1 FILES] $
```

Pour empêcher les autres de lire ou d’exécuter le script :


```
[orkin@login1 FILES] $ chmod o-rx print.py
[orkin@login1 FILES] $ ls -l print.py 
-rwxr-x--- 1 dbrown staff 23 Oct 16 13:37 print.py
```


<div style="page-break-after: always;"></div>

### Plus d’Unix

Voici quelques commandes plus avancées d’Unix qui peuvent s’avérer très utiles. Quand vous aurez un peu de temps, vous pourrez les explorer davantage (entre parenthèses, nous indiquons les pages du *Linux Phrasebook Second Edition* de Scott Granneman ou bien des liens vers des tutoriels en ligne) :

- `awk` (Linux Phrasebook p.194-198) ([tutoriel en ligne](https://www.tutorialspoint.com/awk/index.htm))  
- `sed` (Linux Phrasebook p.188-194) ([tutoriel en ligne](https://www.tutorialspoint.com/sed/index.htm))  
- `perl` one-liners ([tutoriel en ligne](https://catonmat.net/introduction-to-perl-one-liners))  
- `for` loops ([tutoriel en ligne](https://www.tutorialspoint.com/unix/for-loop.htm))

---
### [Lien vers le problème Unix 1](https://github.com/nomascus/ANT3814/blob/main/PROBLEM_SETS/Unix1_problemSet.md)
```
