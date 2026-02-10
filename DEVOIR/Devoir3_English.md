ANT3475: Problem set 3

1.  In the file  [nobody.txt](https://raw.githubusercontent.com/nomascus/ANT3814/main/FILES/nobody.txt), find every occurrence of "Nobody" and print out the line. First do this with  `grep`, then do it with a perl regular expression. Make sure to print ONLY the lines with "Nobody"
    
2.  In the file  [nobody.txt](https://raw.githubusercontent.com/nomascus/ANT3814/main/FILES/nobody.txt), use perl to substitute every occurrence of 'Nobody' with your favorite name and write an output file with that person's name (ex. Michael.txt).
    
3.  Using regular expressions, find all the FASTA header lines in  [unix3.fasta](https://raw.githubusercontent.com/nomascus/ANT3814/main/FILES/unix3.fasta). Note that the format for a header in a FASTA file is a line that starts with a greater than symbol and is followed by some text (e.g.  `>seqName description`  where seqName is the sequence name or identifier. The identifier cannot have spaces in it. The description that follows it can have spaces.)
    
4.  If a line matches the format of a FASTA header, extract the sequence identifier (seqID), which is everything before the first space character and description (seqDescription), which is everything after the first space character) using sub patterns (groups). Hint: check how to search for space characters and "not space" characters 
    
    -   Print sequence information in this format:  `id:seqID desc:seqDescription`
5.  The enzyme ApoI has a restriction enzyme cut site: RAATTY where R and Y are degenerate nucleotides. the enzyme will cut the sequence after the first R. Search google for the IUPAC codes to identify the nucleotide possibilities for the R and Y. Write a regular expression to find and print all the lines where the enzyme be a cut the sequence in the following sequence  [apoI.fasta](https://raw.githubusercontent.com/nomascus/ANT3814/main/FILES/apol.fasta).
    
6.  GTTGCCTGAAATGGCGGAACCTTGAAA is 27 nucleotide long, which means that it should be composed of 9 codons. Use regular expressions to split this into 9 codons separated first by tabs, then by new lines. Hint: you will need to use a quantifier and a subpattern and a global search.
