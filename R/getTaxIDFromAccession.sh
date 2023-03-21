#Run with Command cat myIDs.txt | xargs ./getTaxIDFromAccession.sh > myListOfTaxa.tsv
for ACC in $@
do
    #echo "${ACC}"
    curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${ACC}&rettype=fasta&retmode=xml" |\
    grep TSeq_defline |\
    cut -d '>' -f 2 |\
    cut -d '<' -f 1 |\
    tr -d "\n"
    echo "\t${ACC}"
    

    sleep 0.1
    #break
done
