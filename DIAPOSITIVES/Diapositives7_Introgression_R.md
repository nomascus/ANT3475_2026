# R Studio

R est un langage de programmation largement utilisé en informatique statistique et biologique. Il est similaire à Unix, mais la syntaxe du code est un peu différente. Nous n’allons pas passer beaucoup de temps à apprendre R, mais nous l’utiliserons pour quelques applications dans ce cours. Commençons par charger l’interface du studio R.

Allez d’abord sur le site web suivant et connectez-vous avec le login et le mot de passe de votre compte cloud ANT3814 :

[https://jupyter.ant3475.calculquebec.cloud/hub/login](https://jupyter.ant3814.calculquebec.cloud/hub/login)

   
1.  Dans la fenêtre Launcher sur la droit, cliquez sur l’icône RStudio (un cercle bleu). Cela devrait ouvrir un terminal blanc avec plusieurs fenêtres, qui est l’interface de R Studio.
    
2.  Nous devons maintenant installer quelques logiciels dans R. Comme c’est la première fois que nous utilisons R, nous devrons installer quelques gros logiciels, mais vous n’aurez à le faire que la première fois ! Ensuite, les choses se dérouleront assez rapidement.
    

Nous devons installer trois programmes :  **admixr**  **rmarkdown**  et  **ggplot2**  qui sont les principaux programmes que nous utiliserons aujourd’hui.

```
install.packages('admixr'); install.packages('rmarkdown'); install.packages('ggplot2')
```

Pendant ce temps, commençons la leçon d’aujourd’hui.

## Flux génétique, hybridation, mélange génétique et introgression

**Arbre phylogénétique des espèces**  : Un arbre phylogénétique qui montre le schéma général de la façon dont les espèces partagent une histoire de population ancestrale commune.

**Arbre génétique**  : Les gènes individuels (et leurs allèles) peuvent avoir une histoire différente au sein des espèces à mesure qu’elles se séparent les unes des autres. Nous pouvons examiner les phylogénies des gènes individuels, que l’on appelle « arbres génétiques ».

**Tri incomplet des lignées**  : Lorsque le polymorphisme génétique ne se regroupe pas au sein d’une seule espèce. Il se regroupe dans l’ancêtre commun de plusieurs espèces.

**Flux génétique**  : Lorsque des individus migrent d’une population à une autre, ils peuvent introduire des allèles nouveaux ou de faible fréquence dans la population secondaire. Cela peut changer la fréquence des allèles de la population secondaire.

**Admixture génomique**  : Des populations (ou espèces) précédemment isolées se croisent (s’hybrident) ce qui entraîne un mélange génétique des deux populations d’origine.

**Introgression génomique**  : Incorporation d’allèles d’une population (ou espèce) dans le pool génétique d’une deuxième population divergente.

Ce sont des concepts très similaires. Généralement, l’hybridation, le mélange génétique et l’introgression font référence au flux génétique entre espèces ou des populations plus éloignées.

[Diapositives ABBA BABA](https://github.com/nomascus/ANT3814/blob/main/DIAPOSITIVES/ABBA-BABA_Francais.pdf)

## Statistique D

REMARQUE : cette section est modifiée à partir du tutoriel admixr ([https://bodkan.net/admixr/articles/01-tutorial.html](https://bodkan.net/admixr/articles/01-tutorial.html))

Supposons que nous soyons intéressés par la question suivante :  _« Quelles populations humaines aujourd’hui montrent des preuves de mélange avec les Néandertaliens ? »_

Par exemple, si nous voulons tester si certains humains modernes d’Europe se sont mélangés avec des Néandertaliens, ce qui augmenterait leur affinité génétique avec ce groupe archaïque par rapport aux humains modernes d’Afrique (dont les ancêtres n’ont jamais rencontré les Néandertaliens), nous aurions besoin d’un groupe de trois populations/espèces et d’un groupe externe :

(Français, Yoruba, Néandertalien, Chimpanzé)

Nous pouvons appeler ces quatre groupes les populations W, X, Y et Z.

Une façon de voir cela est d’utiliser la statistique D suivante :

### Sites BABA - # Sites ABBA  
***
### Sites BABA + # Sites ABBA

Les statistiques 𝐷 sont basées sur la comparaison des proportions de motifs de sites BABA et ABBA observés dans les données. Un écart significatif de 𝐷 par rapport à zéro indique un excès de partage d’allèles entre la première et la troisième population (𝐷 positif), ou un excès de partage d’allèles entre la deuxième et la troisième population (𝐷 négatif). Si nous obtenons 𝐷 qui n’est pas significativement différent de 0, cela suggère que les première et deuxième populations forment un clade et ne diffèrent pas dans le taux de partage d’allèles avec la troisième population (c’est l’hypothèse nulle contre laquelle les données sont comparées).

## Devoir_7

Nous devons télécharger quelques fichiers pour travailler. Cliquez sur l’autre onglet de votre navigateur, qui devrait ressembler au terminal ssh, mais en blanc et non en noir. L’invite de commande sera différente, mais fonctionnera comme un terminal normal.

Tout d’abord, nous devons créer un nouveau dossier pour le cours de cette semaine (7_introgression), alors  **créons ce dossier à l’emplacement suivant, en remplaçant USERID par votre nom d’utilisateur**.

```
cd ~/scratch
mkdir 7_introgression
```

téléchargeons et extrayons quelques fichiers :

```
cd 7_introgression
wget https://raw.githubusercontent.com/nomascus/ANT3814/refs/heads/main/PROBLEM_SETS/Introgression_problemSet_francais_2025.rmd 
wget --no-check-certificate https://bioinf.eva.mpg.de/admixr/snps.tar.gz
tar -xvf snps.tar.gz

```

Ensuite, nous devons installer une version de admixtools sur la ligne de commande.

```
git clone https://github.com/DReichLab/AdmixTools.git
cd AdmixTools
module load openblas/0.3.24
module load gsl/2.7
cd src
make clobber; make all; make install
cd ../bin
```

Une fois admixtools installé, nous devons aider R à trouver son emplacement. Veillez à  **changer USERID par votre identifiant dans les commandes suivantes**.

```
echo "export PATH=\"$PATH:/home/USERID/scratch/7_introgression/AdmixTools/bin\"" >> ~/.bashrc
source ~/.bashrc
cp /cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v4/Compiler/gcccore/r/4.4.0/lib64/R/etc/Renviron ~/.Renviron
echo "PATH=$PATH" >> ~/.Renviron

```

Passez maintenant à l’onglet Rstudio de votre navigateur. Allez dans la fenêtre en bas à droite de l’écran et cliquez sur  _Files_. Vous obtiendrez ainsi une liste des dossiers et des fichiers de votre répertoire personnel (~) sur l’ordinateur dématérialisé. Naviguez jusqu’au nouveau dossier (7_introgression) en cliquant sur les icônes de dossier comme vous le feriez sur votre ordinateur personnel. Cliquez sur le fichier de l’ensemble de problèmes que nous venons de télécharger sur l’ordinateur cloud (Introgression_problemSet_francais_2025.rmd).

Cela devrait ouvrir un fichier R markdown dans la fenêtre en haut à gauche, qui vous permettra de cliquer sur le code R. Les fichiers R markdown comportent des blocs de code séparés en gris et un espace pour les commentaires en blanc.

Passons au deuxième bloc de code, qui se trouve sous « > Charger admixr et les données ». Vous pouvez exécuter chaque ligne de code individuellement, en allant sur cette ligne et en cliquant (sur pc) sur ctrl-enter ou (sur mac) sur cmd-enter. Remarquez qu’il y a probablement du code qui s’exécute dans la fenêtre de la console ci-dessous. Il se peut également qu’une petite fenêtre contenant une figure ou un tableau apparaisse sous le bloc de code une fois que le code a fini de s’exécuter. Vous pouvez également copier et coller le code dans la fenêtre de la console et appuyer sur Entrée pour l’exécuter.

En général, dans R, quelques lignes de code peuvent être traitées comme une seule ligne de code si elles se trouvent entre des parenthèses (comme les lignes 19 à 22) ou des signes plus (comme les lignes 65 à 68). Dans ce cas, il vous suffit de cliquer sur ctrl-enter ou cmd-enter à la ligne 19 ou 65 pour exécuter le groupe entier. Vous pouvez également cliquer sur la flèche verte en haut du bloc de code pour exécuter toutes les lignes de code du bloc. Cependant, il est préférable d’exécuter le code ligne par ligne, afin de voir ce qu’il fait.

Ok, cela devrait être suffisant pour vous permettre de commencer à travailler sur les problèmes de cette semaine. Les questions sont écrites dans le fichier rmarkdown, mais vous devez coller vos réponses dans le quiz (devoir_7) sur StudiUM.
