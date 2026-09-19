# in the directory containing hydrolytic enzyme amino acid sequences

for file in *.fasta

do
withpath="${file}"
filename=${withpath##*/}
base="${filename%.fasta}"
echo "${base}"

signalp6 --fastafile ./"${base}".fasta --organism other --output_dir ./"${base}" --format txt --mode fast

done
