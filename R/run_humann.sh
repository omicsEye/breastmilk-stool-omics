#!/bin/bash

# run: bash run_humann.sh
input='/mys3bucket/input'
output='/mys3bucket/output'
i=0
pair='_R2_'
for file in $input/*.fastq.gz
do
        t=$(basename $file .fastq.gz)
        if [ -f $output/${t}/${t}_pathcoverage.tsv ]
        then
                continue
        fi
        if [[ $t == *$pair* ]]; then
                continue
        fi
        echo ${file}
        echo ${t}
        # humann -i /mys3bucket/input/${file} -o /mys3bucket/output/${t} --nucleotide-database /opt/biobakery_workflows_databases/humann/chocophlan/ --protein-database /opt/biobakery_workflows_databases/humann/uniref/  --search-mode uniref90 --memory-use maximum --thread 16  --bypass-translated-search --resume &
        #disown
        sleep 2
        i=$((i + 1))
        if [ $i -ge 16 ]; then
                i=0
                sleep 2
        fi


done
echo "Done running humann"
