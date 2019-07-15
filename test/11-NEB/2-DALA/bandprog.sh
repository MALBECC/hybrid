name=DALA
rm salida.dat col1 col2
grep  Energy-band $name.out| grep -v eV |awk '{print $2}' > col1
grep  Energy-band $name.out| grep -v eV |awk '{print $3}' > col2
paste col1 col2 > Band-progresion.dat
rm col1 col2

