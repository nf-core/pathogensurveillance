#!/bin/bash

ref="Phytophthora ramorum"  # get the choosed Reference Sequence

CODEX=$(echo "$ref" | awk '{print toupper(substr($1,1,3))toupper(substr($2,1,2))}')

grep '^>' *.fna > buscos.fof # list of names to modify

while IFS= read -r fastafile; do
	#echo "Processing file $fastafile"

	# Extracting file name
	fasta=$(echo $fastafile | awk -F ':' '{print $1}')
	#header=$(echo $fastafile | awk -F ':' '{print $2":"$3}')
	header=$(grep '>' $fasta)
	
	# Genearting a unique IDentifier
	NUM=$(echo $fasta | cut -d 'a' -f1)
	echo "NUM $NUM"
	ID="$CODEX$NUM"
	echo "Identifier: $ID"

        # Replace the headers in the FASTA file
	echo "Processing file fasta $fasta with header $header"
        sed -i "s/$header/>$ID \| [$ref]/" "$fasta"
done < buscos.fof
