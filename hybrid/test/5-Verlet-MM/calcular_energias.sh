#Copia name_original.out a un archivo temporal
#Extra energia potencial (segun qm,mm,qmmm), cinética y total
#

name_original=GCG
name=XX12131234
tipo=mm # qm, mm, qmmm

cp $name_original.out $name.out

if [ $tipo == qmmm ] 
then
grep "Etots:" $name.out|awk '{print $2}' > potencial
elif [ $tipo == mm ]
then
grep "Esolv:" $name.out|awk '{print $2}' > potencial
elif [ $tipo == qm ]
then
grep "E lio:" $name.out|awk '{print $3}' > potencial
fi

grep "hybrid: Kinetic Energy (eV):" $name.out|awk '{print $5}' > cinetica
grep "hybrid: Total Energy + Kinetic (eV):" $name.out|awk '{print $7}' > total
grep "Begin move" $name.out|awk '{print $4}' > paso


paste paso potencial cinetica total > energias_todas
sed  -i '1i #paso    potencial      cinetica         total'  energias_todas.dat
paste paso total > energia_total.dat

rm $name.out potencial cinetica total total
#primervalor_potencial=`head -1 potencial_temp`
#primervalor_cinetica=`head -1 cinetica_temp`
#primervalor_total=`head -1 total_temp`

#echo "awk '{print $1 - primervalor_potencial}' potencial_temp > potencial" >> awk.sh
#echo "awk '{print $1 - primervalor_cinetica}' cinetica_temp > cinetica " >> awk.sh
#echo "awk '{print $1 - primervalor_total}' total_temp > total " >> awk.sh

#sed -i s/primervalor_potencial/$primervalor_potencial/g awk.sh
#sed -i s/primervalor_cinetica/$primervalor_cinetica/g awk.sh
#sed -i s/primervalor_total/$primervalor_total/g awk.sh

#chmod +x awk.sh 
#./awk.sh 

#rm *_temp

