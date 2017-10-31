#include "mpi_util.h"

void MpiUtil::allgatherv_char_array(char* local_input, int local_input_size, char** output, int* output_size, MPI_Comm comm) {
	int nproc;
	MPI_Comm_size(comm, &nproc);
	int local_sizes[nproc];
	int displacements[nproc];
	MPI_Allgather(&local_input_size, 1, MPI_INT, local_sizes, 1, MPI_INT, comm);
	*output_size = 0;
	for (int i = 0; i < nproc; i++) {
		displacements[i] = *output_size;
		*output_size += local_sizes[i];
	}
	*output = new char[*output_size];
	MPI_Allgatherv(local_input, local_input_size, MPI_CHAR, *output, local_sizes, displacements, MPI_CHAR, comm);
	//&(*output)[0]
}

void MpiUtil::broadcast_char_array(char** array, int* array_size, int root, MPI_Comm comm) {
	int me;
	MPI_Comm_rank(comm, &me);
	MPI_Bcast(array_size, 1, MPI_INT, root, comm);
	if (me != root) {
		*array = new char[*array_size];
	}
	MPI_Bcast(*array, *array_size, MPI_CHAR, root, comm);
}

void MpiUtil::gatherv_char_array(char* local_input, int local_input_size, char** output, int* output_size, int root, MPI_Comm comm) {
	int nproc, me;
	MPI_Comm_size(comm, &nproc);
	MPI_Comm_rank(comm, &me);
	int local_sizes[nproc];
	int displacements[nproc];
	MPI_Gather(&local_input_size, 1, MPI_INT, local_sizes, 1, MPI_INT, root, comm);
	*output_size = 0;
	if (me == root) {
		for (int i = 0; i < nproc; i++) {
			displacements[i] = *output_size;
			*output_size += local_sizes[i];
		}
	}
	*output = new char[*output_size];
	MPI_Gatherv(local_input, local_input_size, MPI_CHAR, *output, local_sizes, displacements, MPI_CHAR, root, comm);
}
