#!/bin/bash
#rm *.out 2> captured_stderr
rm *.csv 2> captured_stderr
maxmodellen="1000"
lenmdlconc="15"
kmax="20"
Tcutoff="20"
quadpts="20"
echo $kmax          1>  "Data/ramprem.csv"
echo $maxmodellen   1>> "Data/ramprem.csv"
echo $lenmdlconc    1>> "Data/ramprem.csv"
echo $Tcutoff       1>> "Data/ramprem.csv"
echo $quadpts       1>> "Data/ramprem.csv"
./compile.sh  || exit 1

for basisset in "631G" "R31G" "631pG" "R31pG" 
do
for moleculename in "C" "C2H6" "CO2" "NH3"
do

cd Data
cp Molecules/$moleculename.inp .
cp -rf BasisSets/$basisset .
echo $moleculename 1>"moleculename.csv"
echo $basisset 1>"basisset.csv"
python createinput.py 1>captured_stderr
python quick.py  1> captured_stderr 2> captured_stderr
qchem allgaus.inp > allgaus.out
python parseoutput.py 
rm $moleculename.inp 
cd ..

for lenmdlconc in "15"
do
echo $kmax          1>  "Data/ramprem.csv"
echo $maxmodellen   1>> "Data/ramprem.csv"
echo $lenmdlconc    1>> "Data/ramprem.csv"
echo $Tcutoff       1>> "Data/ramprem.csv"
echo $quadpts       1>> "Data/ramprem.csv"

echo 'STARTING RampItUp'
./RampItUp.exe 1> HFramp.out
echo 'FINISHING RampItUp'
echo $moleculename
echo $basisset
grep 'Hartrampenergy' HFramp.out
grep -A 3 'Hartrampenergy' HFramp.out
grep 'Ramp - All-Gaussian atomisation energy' HFramp.out
done
done
rm -rf Data/$basisset
done

rm Data/*.out
rm Data/captured_stderr
rm Data/*.csv
rm Data/allgaus*
rm captured_stderr
