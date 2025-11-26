#!/bin/bash
#submit all the vasp jobs
for x in `ls POSCAR*`
do
    ind=`echo $x | sed 's/POSCAR//'`
    rm -rf dmol3_$ind
    mkdir dmol3_$ind
    cp cluster.input ms8.sh outmol2OUTCAR.sh dmol3_$ind
    cp cluster.car$ind dmol3_$ind/cluster.car
    cp cluster.car$ind dmol3_$ind/cluster.car.bak
    cd dmol3_$ind
 #   sbatch vasp.pbs |tail -1|awk '{print $NF}' > job_ID.dat
#PBS
    oldpwd=`pwd`
    sh ms8.sh
    sleep 1s
    sh outmol2OUTCAR.sh
    sleep 1s
    cd ..
done

#now check whether all the calculations finish

#stat=0
#while [ $stat -eq 0 ]; do
#    sleep 60s
#    stat=1
#    /opt/gridview/pbs/dispatcher/bin/qstat | grep xggong | awk '{printf("%6d\n", $1)}' >back
#    n=$(wc -l back | awk '{printf("%d\n", $1)}')
#    nu=$(($n+1))
#    for x in `ls -d vasp_*/`
#    do
#        jobid=$(head -1 $x/job_ID.dat|awk '{printf("%6d\n", $1)}')
#        for((l=1;l<$nu;l++)); do
#            jobone=$(sed -n "${l}p" back | awk '{printf("%6d\n", $1)}')
#            if [[ $jobid == $jobone ]]; then
#                stat=0
#            fi
#        done
#    done
#    rm back
#done
#
#rm -rf POSCAR_*
