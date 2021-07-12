#!/bin/sh
#PBS -N solv_screen
#PBS -q default
#PBS -l nodes=1:ppn=24,walltime=100:00:00
#PBS -m e

cd ${PBS_O_WORKDIR}
#Determine number of isoconfigurational trajectories (=n1*n2)
n1=25
n2=4

#Equilibrate system

sed -e "s/RSEED1/${RANDOM}/" -e "s/RSEED2/${RANDOM}/" -e "s/RSEED3/${RANDOM}/" -e "s/RSEED4/${RANDOM}/" -e "s/RSEED5/${RANDOM}/" "equilibration.lmp" > "testlmps"
//home/shared/lammps-images/lmp_ubuntu_simple -p 1  -i "testlmps" &
wait

mv master* data/
rm testlmps

#Run isoconfigurational trajectories

F=0.001
for i in `seq 1 $n2`
do
	for j in `seq 1 $n1`
	do
		k=$((($i-1) * $n1 + $j))
		rseed=${RANDOM}
		sleep 0.1
		sed -e "s/RSEED1/${rseed}/" -e "s/FLAG/${k}off/" -e "s/TRAJOUT/dumpoff-${k}/" -e "s/VVAL/0.0/" "in2.images" > "suboff-${k}"
		//home/shared/lammps-images/lmp_ubuntu_simple -p 1  -i "suboff-${k}" &
		sleep 0.1
		sed -e "s/RSEED1/${rseed}/" -e "s/FLAG/${k}on/" -e "s/TRAJOUT/dumpon-${k}/" -e "s/VVAL/${F}/" "in2.images" > "subon-${k}"
		//home/shared/lammps-images/lmp_ubuntu_simple -p 1  -i "subon-${k}" &
	done
	wait
done 
wait

#Recalculate energies with correct potentials

for i in `seq 1 $n2`
do
	for j in `seq 1 $n1`
	do
		k=$((($i-1) * $n1 + $j))
		sleep 0.1
		sed -e "s/FLAG/${k}off/" -e "s/TRAJOUT/fieldoff-${k}/" -e "s/INNAME/dumpoff-${k}/" "in3.images" > "suboff-${k}"
		//home/shared/lammps-images/lmp_ubuntu_simple -p 1 -i "suboff-${k}" &
		sleep 0.1
		sed -e "s/FLAG/${k}on/" -e "s/TRAJOUT/fieldon-${k}/" -e "s/INNAME/dumpon-${k}/" "in3.images" > "subon-${k}"
		//home/shared/lammps-images/lmp_ubuntu_simple -p 1 -i "subon-${k}" &
	done
	wait
done
wait

mv fieldo* data/
mv dumpo* data/

rm log.lammps*
rm suboff*
rm subon*
rm screen.0

