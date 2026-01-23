export mspath=/home/zyy/Accelrys/MaterialsStudio8.0
export PATH=$mspath/etc/CASTEP/bin:$PATH
export PATH=$mspath/etc/DMol3/bin:$PATH
export PATH=$mspath/etc/Scripting/bin:$PATH

inputfile=cluster 
RunDMol3.sh -np 6 $inputfile
