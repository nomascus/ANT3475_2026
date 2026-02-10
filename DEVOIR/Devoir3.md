ANT3475 : Série d'exercices 3

1.  Dans le fichier  [nobody.txt](https://raw.githubusercontent.com/nomascus/ANT3814/main/FILES/nobody.txt), trouvez toutes les occurrences de « Nobody » et affichez la ligne. Commencez par utiliser  `grep`, puis utilisez une expression régulière Perl. Veillez à n'imprimer QUE les lignes contenant « Nobody ».
    
2.  Dans le fichier  [nobody.txt](https://raw.githubusercontent.com/nomascus/ANT3814/main/FILES/nobody.txt), utilisez Perl pour remplacer toutes les occurrences de « Nobody » par votre nom préféré et écrivez un fichier de sortie avec le nom de cette personne (par exemple, Michael.txt). .
    
3.  À l'aide d'expressions régulières, recherchez toutes les lignes d'en-tête FASTA dans  [unix3.fasta](https://raw.githubusercontent.com/nomascus/ANT3814/main/FILES/unix3.fasta). Notez que le format d'un en-tête dans un fichier FASTA est une ligne qui commence par un symbole supérieur à et qui est suivie d'un texte (par exemple  `>seqName description`  où seqName est le nom ou l'identifiant de la séquence. L'identifiant ne peut pas contenir d'espaces. La description qui suit peut contenir des espaces.
    
4.  Si une ligne correspond au format d'un en-tête FASTA, extrayez l'identifiant de séquence (seqID), qui correspond à tout ce qui précède le premier espace, et la description (seqDescription), qui correspond à tout ce qui suit le premier espace, à l'aide de sous-motifs (groupes). Astuce : vérifiez comment rechercher les espaces et les caractères « non espaces ». 
    
    - Imprimez les informations de séquence dans ce format :  `id:seqID desc:seqDescription`
5.  L'enzyme ApoI a un site de coupure d'enzyme de restriction : RAATTY où R et Y sont des nucléotides dégénérés. L'enzyme coupera la séquence après le premier R. Recherchez les codes IUPAC sur Google pour identifier les possibilités de nucléotides pour R et Y. Écrivez une expression régulière pour trouver et imprimer toutes les lignes où l'enzyme peut couper la séquence dans la séquence suivante  [apoI.fasta](https://raw.githubusercontent.com/nomascus/ANT3814/main/FILES/apol.fasta) .
    
6.  GTTGCCTGAAATGGCGGAACCTTGAAA a une longueur de 27 nucléotides, ce qui signifie qu'il devrait être composé de 9 codons. Utilisez des expressions régulières pour le diviser en 9 codons séparés d'abord par des tabulations, puis par des nouvelles lignes. Indice : vous devrez utiliser un quantificateur, un sous-motif et une recherche globale.

