#pragma once
#include "mpi.h"
#include <vector>

class MpiUtil {
public:
	template <typename Container>
		void static allgatherv_string_container(Container* input, Container* output, MPI_Comm comm = MPI_COMM_WORLD) {
			int local_input_size = 0, output_size = 0;
			std::vector<char> local_input_vector = concat_strings(input, &local_input_size);
			char local_input_char_array[local_input_size];
			memcpy(local_input_char_array, &local_input_vector[0], local_input_size);
			char* output_char_array = NULL;
			allgatherv_char_array(local_input_char_array, local_input_size, &output_char_array, &output_size, comm);
			std::vector<char> output_vector(output_char_array, output_char_array + output_size);
			extract_strings(output_vector.begin(), output_vector.end(), output);
			if (output_char_array) delete[] output_char_array;
	}

	template <typename Container>
		void static gatherv_string_container(Container* input, int root, Container* output, MPI_Comm comm = MPI_COMM_WORLD) {
			int local_input_size = 0, output_size = 0;
			std::vector<char> local_input_vector = concat_strings(input, &local_input_size);
			char local_input_char_array[local_input_size];
			memcpy(local_input_char_array, &local_input_vector[0], local_input_size);
			char* output_char_array = NULL;
			gatherv_char_array(local_input_char_array, local_input_size, &output_char_array, &output_size, root, comm);
			int me;
			MPI_Comm_rank(MPI_COMM_WORLD, &me);
			if (me == root) {
				std::vector<char> output_vector(output_char_array, output_char_array + output_size);
				extract_strings(output_vector.begin(), output_vector.end(), output);
			}
			if (output_char_array) delete[] output_char_array;
	}

	template <typename Container>
		void static broadcast_string_container(Container* input_output, int root, MPI_Comm comm = MPI_COMM_WORLD) {
			int me;
			MPI_Comm_rank(comm, &me);
			int total_size = 0;
			char* char_array = NULL;
			if (me == root) {
				std::vector<char> input_char_vector = concat_strings(input_output, &total_size);
				char_array = new char[total_size];
				memcpy(char_array, &input_char_vector[0], total_size);
			}
			//std::cout << "total size: " << total_size << " me:" << me << std::endl;
			//MPI_Barrier(MPI_COMM_WORLD);
			broadcast_char_array(&char_array, &total_size, 0, comm);
			//std::cout << "total size: " << total_size << " me:" << me << std::endl;
			//MPI_Barrier(MPI_COMM_WORLD);
			//if (me == 0) std::cout << "waiting for broadcast" << std::endl;
			if (me != 0) {
				//std::cout << "me: " << me << " received: " << char_array << std::endl;
				std::vector<char> output_char_vector(char_array, char_array + total_size);
				//std::cout << "extracting strings" << std::endl;
				extract_strings(output_char_vector.begin(), output_char_vector.end(), input_output);
				//std::cout << "strings extracted" << std::endl;
				//input_output -> shrink_to_fit();
			}
			//std::cout << " deleteing array " << me << std::endl;
        		if (char_array) delete[] char_array;
			//std::cout << " array deleted " << me << std::endl;
	}

	template<class Iter, typename Container>
		void static extract_strings(Iter input_begin, Iter input_end, Container* output) {
			Iter current = input_begin;
			while (current != input_end) {
				if (*current == '\0') {
					std::string str(input_begin, current);
					output -> insert(output -> begin(), str);
					input_begin = current + 1;
				}
				current++;
			}
	}

	template<typename Container>
		std::vector<char> static concat_strings(Container* strs, int* total_size) {
			std::vector<char> char_vector;
			*total_size = 0;
			for (auto iter = strs -> begin(); iter != strs -> end(); iter++) {
				char_vector.insert(char_vector.end(), iter -> begin(), iter -> end());
				char_vector.push_back('\0');
				*total_size += iter -> length() + 1;
			}
			return char_vector;
	}
private:
	void static broadcast_char_array(char** array, int* array_size, int root, MPI_Comm comm);
	void static gatherv_char_array(char* local_input, int local_input_size, char** output, int* output_size, int root, MPI_Comm comm);
	void static allgatherv_char_array(char* local_input, int local_input_size, char** output, int* output_size, MPI_Comm comm);
};
