#!/bin/bash

rm Libraries/*
rm ObjectFiles/*
for dir in "SCF" "ModelRG" "Multipoles" "gggg" "RIgg" "Initial_Setup" "One_Elec_Int" 
do
    cp DEBUG.h $dir/.
    cd $dir
rm *.f*~
for file in *.f*
do
  echo $file
  gfortran -O3 -cpp -pg -g -c $file || exit 1
done
for file in *.o
do
  echo $file
  ar r $dir.a $file || exit 1
  mv $file ../.
done
mv $dir.a ../.
rm DEBUG.h
cd ..

done
gfortran -O3 -g -c -cpp master.f90 || exit 1
gfortran -O3 -g -c -cpp twoelec.f90 || exit 1
gfortran -O3 -g -c -cpp twoelecinside.f90 || exit 1
gfortran -O3 -g -c -cpp SRgenRIsp.f90
gfortran -O3 -g -c -cpp SRgenRIss.f90
gfortran -O3 -g -c -cpp SRgenRIpp.f90
gfortran -O3 -g -c -cpp readQchem.f90 

mv Libraries/* .
mv ObjectFiles/* .

gfortran -O3 -g -o RampItUp.exe master.o twoelec.o twoelecinside.o SRgenRIsp.o SRgenRIss.o SRgenRIpp.o SCF.a One_Elec_Int.a Multipoles.a Initial_Setup.a RIgg.a ModelRG.a gggg.a -framework Accelerate readQchem.o 

mv *.a Libraries/.
mv *.o ObjectFiles/.
#cp RampItUp.exe ../Exe/.


