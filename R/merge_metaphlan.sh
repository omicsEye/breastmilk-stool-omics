#!/bin/sh

#  merge_metaphlan.sh
#  
#
#  Created by Ali Rahnavard on 5/6/21.
#  
python3 ./utils/merge_metaphlan_tables.py  ~/Documents/INOVA/bugs_list/*.tsv>  ~/Documents/INOVA/species_abundance_table.txt
