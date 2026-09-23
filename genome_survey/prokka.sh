# Split genomes by Kingdom 

for file in *.fna
do
withpath="${file}"
filename=${withpath##*/}
base="${filename%.fna}"
echo "${base}"
prokka --cpus 20 --evalue 1e-50 --kingdom Bacteria --outdir ./prokkawithpfam_e50/"${base}" --prefix "${base}" --locustag "${base}" ./"${base}".fna
done

for file in *.fna
do
withpath="${file}"
filename=${withpath##*/}
base="${filename%.fna}"
echo "${base}"
prokka --cpus 20 --evalue 1e-50 --kingdom Archaea --outdir ./prokkawithpfam_e50/"${base}" --prefix "${base}" --locustag "${base}" ./"${base}".fna
done
