#preprocessing

tr "\n" "\t" < SSU.fasta | sed -e 's/>/\n>/g' | sed '/^\n*$/d' > SSU.temp

tr "\n" "\t" < full_ITS.fasta | sed -e 's/>/\n>/g' | sed '/^\n*$/d' > full_ITS.temp

tr "\n" "\t" < LSU.fasta | sed -e 's/>/\n>/g' | sed '/^\n*$/d' > LSU.temp

#matching SSU + ITS

gawk 'BEGIN {FS=OFS="\t"} FNR==NR{a[$1]=$2;next} ($1 in a) {print $1,a[$1],$2}' SSU.temp full_ITS.temp > SSU-ITS.temp

#matching SSU+ITS with LSU

gawk 'BEGIN {FS=OFS="\t"} FNR==NR{a[$1]=$2$3;next} ($1 in a) {print $1,a[$1],$2}' SSU-ITS.temp LSU.temp > SSU-ITS-LSU.temp

#Formatting output; Final output = SSU-ITS-LSU.fasta

sed -e 's/\t/\n/' < SSU-ITS-LSU.temp | sed -e 's/\t//' > SSU-ITS-LSU.fasta

rm *.temp

#check how many seqs it has

grep -c “>" SSU-ITS-LSU.fasta
