#!/bin/bash

filename='bwa_metadata.txt'

while read line; do

   ID=$(echo $line | awk '{print $1"_"$2}')
   ref=$(echo $line | awk '{print $3}')

   echo "bwa mem ${ID}_F.fasta ${ID}_R.fasta $ref -o ${ID}.sam"

done < $filename
