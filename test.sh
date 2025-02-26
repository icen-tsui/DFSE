#!/bin/sh
#SBATCH -J HfO2-2d-12-atoms
#SBATCH -N 1
#SBATCH -n 48
#SBATCH -t 168:0:0
#SBATCH -p wholenodeQ
#SBATCH --output=%j.out
#SBATCH --error=%j.err

export LD_LIBRARY_PATH=/public/software/compiler/intel-compiler/2021.3.0/mkl/lib/intel64/:$LD_LIBRARY_PATH
#module load mpi/intelmpi/2021.3.0
#module load apps/vasp/intelmpi/5.4.4
cd $SLURM_SUBMIT_DIR
#mpirun -np 48 /public/software/apps/vasp/intelmpi/5.4.4/bin/vasp_std >out.vasp
source /public/home/win0812/intel/oneapi/setvars.sh --force
export PATH=/public/home/win0812/intel/oneapi/mkl:$PATH

export PATH=/public/home/win0812/gcc/gcc-12.1.0/objdir/bin:/public/home/win0812/gcc/gcc-12.1.0/objdir/lib64:$PATH
export LD_LIBRARY_PATH=/public/home/win0812/gcc/gcc-12.1.0/objdir/lib/:/public/home/win0812/gcc/gcc-12.1.0/objdir/lib64/:$LD_LIBRARY_PATH

ulimit -s unlimited
chmod 777 de.x 
conda activate  ferroelectric_search
# Output files
results_file="./results.dat"
log_file="./run.log"
num_parallel=48    # Number of parallel processes
IM2ODE_runtime=300
Perturbation_times=2
MAX_LOOPS=5 #for sturcture relaxzation
MAX_TIME1=7200 #for one time vasp max running time in second when relax the reference phase
MAX_TIME2=12800 #for  one time vasp max running time in second when relax the ferroelectric candidates
#please change VASP_COMMAND for your environment
VASP_COMMAND="mpirun -np 48 /public/home/win0812/software/vasp.5.4.4/bin/vasp_std"
# Function to generate random numbers within a specified range
generate_random_integer() {
    local min=$1
    local max=$2
    # 使用awk生成随机整数
    awk -v seed="$RANDOM" -v min="$min" -v max="$max" 'BEGIN {
        srand(seed)
        print int(min + rand() * (max - min + 1))
    }'
}

generate_random_decimal() {
    local min=$1
    local max=$2
    awk -v seed="$RANDOM" -v min="$min" -v max="$max" 'BEGIN {
        srand(seed)
        printf "%.1f\n", min + rand() * (max - min)
    }'
}


# Ranges for the variables
DIS_min=1.5
DIS_max=4.0
Volumn_min=240.0
Volumn_max=300.0
AREA_min=24.0
AREA_max=27.0
LAYER_HEIGHT_min=5.0
LAYER_HEIGHT_max=7.0

# Maximum runtime for each simulation in seconds
max_runtime=120  # Example: 300 seconds

# Clear log file and results file
> $log_file
> $results_file

# Function to generate configuration and run the simulation
run_simulation() {
    local i=$1
    local work_dir="work_dir_$i"
    local config_file="$work_dir/de.in"
    local results_dir="$work_dir/results"
    local results_individual="$results_dir/de_ini_1"
    
    # Create a unique work directory for this simulation
    mkdir -p "$results_dir"

    # Generate random values for the variables

    DIS1_1=$(generate_random_decimal $DIS_min $DIS_max)
    DIS1_2=$(generate_random_decimal $DIS_min $DIS_max)
    DIS1_3=$(generate_random_decimal $DIS_min $DIS_max)
    DIS2_1=$(generate_random_decimal $DIS_min $DIS_max)
    DIS2_2=$(generate_random_decimal $DIS_min $DIS_max)
    DIS2_3=$(generate_random_decimal $DIS_min $DIS_max)
    DIS3_1=$(generate_random_decimal $DIS_min $DIS_max)
    DIS3_2=$(generate_random_decimal $DIS_min $DIS_max)
    DIS3_3=$(generate_random_decimal $DIS_min $DIS_max)
    Volumn=$(generate_random_integer  $Volumn_min  $Volumn_max)
    AREA=$(generate_random_integer $AREA_min $AREA_max)
    LAYER_HEIGHT=$(generate_random_integer $LAYER_HEIGHT_min $LAYER_HEIGHT_max)
    # Create the configuration file
    cat > "$config_file" <<EOL
SystemName=Hf4O8
NumberOfSpecies=2
NumberOfElements=4 8
NameOfElements=Hf  O 
DistanceOfAtom=
DIS1=$DIS1_1 $DIS1_2 #$DIS1_3
DIS2=$DIS2_1 $DIS2_2 #$DIS2_3
#DIS3=$DIS3_1 $DIS3_2 $DIS3_3
Population=10
MaxStep=1
#Volumn=$Volumn
#De_ratio=0.6
SelectiveDynamics=F
#Symmetry=T
#spg_front=123
#spg_rear=142
#Pickup=T         #########################F
#Pickup_step=10    #########################no
#Mystruct=2
#choice of ab initio software (default: VASP)
#PWMAT=F
#The unit of pressure is GPa
#PSTRESS = 15

#Multi-Objective=T
#hardness=F
#rcut=2.5
#ionicity=0.0
#ESflag=T
#ES_mod=1
#ES_Eg=0.4
#ES_opt=1.0

#HSE=T
#HSE_population=5
#LDA_population=15
#energy_cut=-7.35
#gap_cut=5.0
#LDA_ES_Eg=0.70
#LDA_Es_opt=0.85
#HSE_ES_Eg=1.50


Q2D=T
vacuum_layer=15
Area=$AREA
Layer_height=$LAYER_HEIGHT

#cluster=T
#model_ball=1.0
#init_radius=1.5
#cluster_ctr_x=0.5
#cluster_ctr_y=0.5
#cluster_ctr_z=0.5
#model_shell=0.0
#shell_radius_out=2.0
#shell_radius_in=1.8
#shell_ctr_x=0.50
#shell_ctr_y=0.50
#shell_ctr_z=0.50
#model_plate=0.0
#plate_radius=2.8
#plate_height=0.2
#plate_ctr_x=0.5
#plate_ctr_y=0.5
#plate_ctr_z=0.32
#cluster_substrate=T
#fix_atom=T
#SubstrateElements=4 0 0

#for surface model, substrate is needed, otherwise you should switch to Q2D
#SelectiveDynamics=T is strongly recommended
#model_surface=T
#surface_height=5.5

#ribbon is built upon surface model
#a directin is periodic, b is perpendicular to the surface
#model_ribbon=T
#ribbon_b_min=4.1850
#ribbon_b_max=6.2925

#for grain boundary generation
#model_gb=T
#the following three are cartisian coordinations
#bottom_height=15.053
#gb_height=3.893
#top_height=16.204
#transverse freedom for the top layer (direct)
#transverse_a=1.0
#transverse_b=1.0


#fix_lat=T
#fix_a=4.26000
#fix_b=2.46000
#fix_c=30.0000            #  53.0350
#fix_alpha=90
#fix_beta=90
#fix_gama=90
#find_defect=T

#dimer_mode=T
#In dimer_mode, NumberOfElements should be even
#only available in bulk system
#dimer_dis=1.2

#polarization_3d = .TRUE.
#crystal_system_3d = 2 3 4 5 6 7 
#crystal_system_count_3d = 6

polarization_2d = .TRUE.
crystal_system_2d = 1 2 3 4 5 6  
crystal_system_count_2d = 6

FQFE=.FALSE.
COEFFICIEN=0.06


EOL

    #echo "Iteration $i: Updated DIS1 to $DIS1_1 $DIS1_2, DIS2 to $DIS2_1 $DIS2_2, Area to $AREA, Layer_height to $LAYER_HEIGHT"
    echo "Iteration $i: Updated DIS1 to $DIS1_1 $DIS1_2, DIS2 to $DIS2_1 $DIS2_2, Volumn to $Volumn"

    # Log the current values
    #echo "Iteration $i: DIS1=$DIS1_1 $DIS1_2, DIS2=$DIS2_1 $DIS2_2,  Area=$AREA, Layer_height=$LAYER_HEIGHT" >> "$log_file"
    echo "Iteration $i: Updated DIS1 to $DIS1_1 $DIS1_2, DIS2 to $DIS2_1 $DIS2_2, Volumn to $Volumn" >> "$log_file"

    # Change to work directory and run ./de.x with a timeout
    (cd "$work_dir" && timeout $max_runtime ../de.x)

    # Check the exit status of ./de.x
    if [ $? -eq 124 ]; then
        echo "Iteration $i: ./de.x exceeded the maximum runtime of $max_runtime seconds and was terminated." >> "$log_file"
    else
        echo "Iteration $i: ./de.x completed successfully." >> "$log_file"
    fi

    # Check if de_ini_1 file was updated and copy it to the results directory
    if [ -f "$results_individual" ]; then
        cat "$results_individual" >> "$results_file"
    else
        echo "Iteration $i: No results found for this iteration." >> "$log_file"
    fi

    sleep 10
}

export -f run_simulation generate_random_decimal generate_random_integer
export DIS_min DIS_max AREA_min AREA_max Volumn_max Volumn_min LAYER_HEIGHT_min LAYER_HEIGHT_max max_runtime log_file results_file

seq 1 $IM2ODE_runtime | xargs -n 1 -P $num_parallel -I {} bash -c 'run_simulation "$@"' _ {}
python3 screening_POSCAR.py
# Delete work directories 
for i in $(seq 1 $IM2ODE_runtime); do
    if [ -d "work_dir_$i" ]; then
        rm -rf "work_dir_$i"
    fi
done
sleep 10

cd ./results
cp ../INCAR_* ./
cp ../POTCAR ./
cp ../*.py ./
judge_string="reached required accuracy - stopping structural energy minimisation"
main_dir=$(pwd)

python3 get_POSCAR.py
rm de_ini_1
for i  in POSCAR*;
do
cat $i >> de_ini_1
done

sleep 20
rm POSCAR*

python3 ferroelectric_search.py 'generate_refer_phase'
# Cyclic optimisation reference phase, self-consistency, berryphase
for i in refer-POSCAR* ; do
    echo $i >> "$main_dir/refer_polar_info.dat"
    echo $i >> "$main_dir/main.log"
    cp INCAR_* $i/
    cp POTCAR $i/
    cp POSCAR POSCAR-origin
    cd $i/
    vaspkit -task 102 -kpr 0.04
    cp INCAR_1 INCAR
    START_TIME=$(date +%s)
    timeout $MAX_TIME1 $VASP_COMMAND | tee -a relax.log
    END_TIME=$(date +%s)
    
    ELAPSED_SECONDS=$((END_TIME - START_TIME))
    ELAPSED_MINUTES=$((ELAPSED_SECONDS / 60))
    ELAPSED_HOURS=$((ELAPSED_MINUTES / 60))

    if [ $ELAPSED_HOURS -gt 0 ]; then
        echo "runtime_$count $ELAPSED_HOURS hours $((ELAPSED_MINUTES % 60)) minutes" >> relax.log
    else
        echo "runtime_$count $ELAPSED_MINUTES minutes" >> relax.log
    fi

    if [ $ECLAPSE_TIME -ge $MAX_TIME ]; then
         count=$MAX_LOOPS
    fi
    found=0
    count=0
     
    while [ $found -eq 0 ] && [ $count -lt $MAX_LOOPS ];do
       ((count++))
       if grep -q "$judge_string" "relax.log";then
	   found=1
	   echo "relax has done" >> "$main_dir/main.log"
           cp CONTCAR POSCAR
           vaspkit -task 102 -kpr 0.04
           cp INCAR_2 INCAR
           $VASP_COMMAND
	   result=$(python3 "$main_dir/ferroelectric_search.py" ./EIGENVAL)

           if [ "$result" = "not_metal" ]; then
               cp INCAR_3 INCAR
	       $VASP_COMMAND
               grep 'Total electronic dipole moment' OUTCAR | awk '{print $6, $7, $8}' >> "$main_dir/refer_polar_info.dat"
               grep 'Ionic dipole moment' OUTCAR | awk '{print $5, $6, $7}' >> "$main_dir/refer_polar_info.dat"
           else
               echo "Sorry, seems this phase is metal, will set the berry phase value as 0" >> "$main_dir/refer_polar_info.dat"
           fi
	else
           cp CONTCAR POSCAR
           vaspkit -task 102 -kpr 0.04
           cp INCAR_1 INCAR
	   $VASP_COMMAND | tee -a  relax.log
        fi
    done
   if [ $found -eq 0 ];then
	   echo "Sorry, seems this phase is not a good one" >> "$main_dir/refer_polar_info.dat"
	   echo "relax failed" >> "$main_dir/main.log"
   fi
   phonopy --symmetry  --tolerance=1.0e-3  -c CONTCAR >> sym.log
   cd ../
done

python3 ferroelectric_search.py 'generate_ferroelectric_phase' $Perturbation_times

# Cyclic optimisation of ferroelectric phase, self-consistency, berryphase
for i in POSCAR-* ; do
    echo $i >> "$main_dir/target_polar_info.dat"
    echo $i >> "$main_dir/main.log"
    echo $i >> "$main_dir/target_polar_energy.dat"
    cp INCAR_* $i/
    cp POTCAR $i/
    cd $i/
    cp POSCAR POSCAR-origin
    vaspkit -task 102 -kpr 0.04
    cp INCAR_4 INCAR
    START_TIME=$(date +%s)
    timeout $MAX_TIME2 $VASP_COMMAND | tee -a relax.log
    END_TIME=$(date +%s)

    ELAPSED_SECONDS=$((END_TIME - START_TIME))
    ELAPSED_MINUTES=$((ELAPSED_SECONDS / 60))
    ELAPSED_HOURS=$((ELAPSED_MINUTES / 60))

    if [ $ELAPSED_HOURS -gt 0 ]; then
        echo "runtime_$count $ELAPSED_HOURS hours $((ELAPSED_MINUTES % 60)) minutes" >> relax.log
    else
        echo "runtime_$count $ELAPSED_MINUTES minutes" >> relax.log
    fi

    found=0
    count=0

    while [ $found -eq 0 ] && [ $count -lt $MAX_LOOPS ];do
       ((count++))
       if grep -q "$judge_string" "relax.log";then
	   found=1
	   echo "relax has done" >> "$main_dir/main.log"
           cp CONTCAR POSCAR
           vaspkit -task 102 -kpr 0.04
           cp INCAR_2 INCAR
           $VASP_COMMAND
	   result=$(python3 "$main_dir/ferroelectric_search.py" ./EIGENVAL)

           if [ "$result" = "not_metal" ]; then
               cp INCAR_3 INCAR
	       $VASP_COMMAND
               grep 'Total electronic dipole moment' OUTCAR | awk '{print $6, $7, $8}' >> "$main_dir/target_polar_info.dat"
               grep 'Ionic dipole moment' OUTCAR | awk '{print $5, $6, $7}' >> "$main_dir/target_polar_info.dat"
           else
               echo "Sorry, seems this phase is metal, will set the berry phase value as 0" >> "$main_dir/target_polar_info.dat"
           fi
	else
           cp CONTCAR POSCAR
           vaspkit -task 102 -kpr 0.04
           cp INCAR_4 INCAR
	   $VASP_COMMAND | tee -a  relax.log
        fi
    done
   if [ $found -eq 0 ];then
	   echo "Sorry, seems this phase is not a good one" >> "$main_dir/target_polar_info.dat"
	   echo "relax failed" >> "$main_dir/main.log"
   fi
   atoms_numbers_info=$(sed -n '7p' POSCAR)
   atoms_numbers_sum=$(echo $atoms_numbers_info | awk '{sum=0; for(i=1;i<=NF;i++) sum+=$i; print sum}')
   energy=$(grep 'free energy' OUTCAR | awk 'END{print $5}')
   energy_per_atom=$(echo "scale=6; $energy / $atoms_numbers_sum" | bc)
   echo "$energy_per_atom" >> "$main_dir/target_polar_energy.dat"
   unset atoms_numbers_info; unset atoms_numbers_sum; unset energy; unset energy_per_atom
   cd ../
done
# sort out the refer_phase and ferroelectric phase
python3 "$main_dir/ferroelectric_search.py" 'process_data'


