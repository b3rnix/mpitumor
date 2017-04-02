if [[ -z "$1" ]]; then
	echo 'Debe especificar el numero de procesos'
	exit
fi

mpirun -ppn $1 ./main_mpi
