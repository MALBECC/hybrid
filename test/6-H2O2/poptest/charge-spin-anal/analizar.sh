grep -A 7 "        ║   MULLIKEN POPULATION ANALYSIS  ║" mulliken >mullonly
grep -A 7 "        ║     SPIN POPULATION ANALYSIS    ║" mulliken >spinonly

grep "        ║    1   ║     8     ║" mullonly |awk '{print $6}' >mull1
grep "        ║    2   ║     8     ║" mullonly |awk '{print $6}' >mull2
grep "        ║    3   ║     1     ║" mullonly |awk '{print $6}' >mull3
grep "        ║    4   ║     1     ║" mullonly |awk '{print $6}' >mull4

grep "        ║    1   ║     8     ║" spinonly |awk '{print $6}' >spin1
grep "        ║    2   ║     8     ║" spinonly |awk '{print $6}' >spin2
grep "        ║    3   ║     1     ║" spinonly |awk '{print $6}' >spin3
grep "        ║    4   ║     1     ║" spinonly |awk '{print $6}' >spin4

NUMOFLINES=$(cat mull1 | wc -l )

#for i in `seq 1 $NUMOFLINES`;do
#  echo $i >> step
#done

awk '{print $1}' cut5.rce > step

paste step mull1 mull2 mull3 mull4 spin1 spin2 spin3 spin4 > todo.dat

rm mullonly spinonly step mull1 mull2 mull3 mull4 spin1 spin2 spin3 spin4

