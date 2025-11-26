#!/bin/sh

if [ -e OUTCAR ] ; then
  rm OUTCAR
fi
echo "# total energy taken from dmol3 *.outmol file" >> OUTCAR
echo "# unit is converted from Ha to eV" >> OUTCAR
n=`grep "opt==" *.outmol | tail -1 | awk '{print $3}'`
m=`echo $n*27.2113863|bc`
echo "  free  energy   TOTEN  =       "$m" eV" >> OUTCAR
