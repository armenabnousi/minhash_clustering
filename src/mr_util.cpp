#include "mr_util.h"

void MrToolbox::gather_kmv(MapReduce* mr, int procs, bool convert, int* num_kmvs, int* num_kvs) {
	int n_kmvs = 0, n_kvs = 0;
	n_kvs = mr -> reduce(convert_kmv_to_kv, &n_kmvs);
	mr -> gather(procs);
	if (convert) n_kmvs = mr -> convert();
	else MPI_Allreduce(&n_kmvs, num_kmvs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	if (num_kmvs) {
		*num_kmvs = n_kmvs;
		*num_kvs = n_kvs;
	}
}

void MrToolbox::convert_kmv_to_kv(char* key, int keybytes, char* values, int nvalues, int* valuebytes, KeyValue* kv, void* int_ptr) {
	int* num_kmv = (int*) int_ptr;
	*num_kmv = *num_kmv + 1;
	for (int i = 0, offset = 0; i < nvalues; offset += valuebytes[i], i++) {
		kv -> add(key, keybytes, values + offset, valuebytes[i]);
	}
}

void MrToolbox::collect_keys(MapReduce* mr, std::vector<char>* collected_keys) {
	int me, nprocs;
	MPI_Comm_rank(MPI_COMM_WORLD, &me);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	CollectionContainer local_container(0);
	//if (me == 0) std::cout << "collecting keys: me: " << me << std::endl;
	mr -> scan(collect_local_keys, &local_container);
	//if (me == 0) std::cout << "local size is: " << local_container.local_size << std::endl;
	char local_array[local_container.local_size];
	prepare_array(&local_container, local_array);
	//if (me == 0) std::cout << "array prepared to be sent" << std::endl;
	int all_sizes[nprocs];
	MPI_Gather(&(local_container.local_size), 1, MPI_INT, all_sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);
	//if (me == 0) std::cout << "local sizes received for gatherv" << std::endl;
	int total_size = 0;
	char* collected_keys_array = NULL;
	int displacements[nprocs];
	if (me == 0) {
		for (int i = 0; i < nprocs; i++) {
			displacements[i] = total_size;
			total_size += all_sizes[i];	
		}
	}
	collected_keys_array = new char[total_size];
	//if (me == 0) std::cout << "sending all keys " << std::endl;
	MPI_Gatherv(local_array, local_container.local_size, MPI_CHAR, collected_keys_array, all_sizes, 
			displacements, MPI_CHAR, 0, MPI_COMM_WORLD);
	//if (me == 0) std::cout << "all keys received; inserting in vector" << std::endl;
	collected_keys -> insert(collected_keys -> begin(), collected_keys_array, collected_keys_array + total_size);
	if(collected_keys_array) delete[] collected_keys_array;
	collected_keys_array = NULL;
}

void MrToolbox::collect_local_keys(char* key, int keybytes, char* values, int nvalues, int* valuebytes, void* toolbox_ptr) {
	CollectionContainer* toolbox = (CollectionContainer*) toolbox_ptr;
	std::vector<char> key_vector(key, key + keybytes);
	toolbox -> local_size += keybytes;
	(toolbox -> local_list).push_back(key_vector);
}

void MrToolbox::prepare_array(CollectionContainer* container, char keys_array[]) {
	int offset = 0;
	//std::cout << "local container size is: " << (container -> local_list).size() << std::endl;
	int me;
	MPI_Comm_rank(MPI_COMM_WORLD, &me);
	//if (me == 0) std::cout << "proc0 num keys is: " << (container -> local_list).size() << std::endl;
	while (!((container -> local_list).empty())) {
		std::vector<char> current = (container -> local_list).back();
		//if (me == 0) {
		//	std::cout << "current key in proc 0 is: ";
		//	for (int i = 0; i < current.size(); i++) {
		//		std::cout << current[i];
		//	}
		//	std::cout << std::endl;
		//}
		std::copy(current.begin(), current.end(), keys_array + offset);
		//std::cout << "current copied as: " << keys_array + offset;
		offset += current.size();
		(container -> local_list).pop_back();
	}
	//if (me == 0) std::cout << "all copied " << std::endl;
}

void MrToolbox::remove_redundant(MapReduce* mr) {
	mr -> reduce(redundant_agglomerator, NULL);
	mr -> convert();
	mr -> reduce(redundant_remover, NULL);
	mr -> convert();
}

void MrToolbox::remove_redundant_local(MapReduce* mr, bool is_string) {
	if (is_string) {
		remove_redundant_string_manual(mr);
	} else {
		mr -> reduce(redundant_agglomerator, NULL);
		mr -> convert();
		int kvs = mr -> reduce(redundant_remover, NULL);
		mr -> convert();
	}
	//std::cout << "kvs after removing local redundancies: " << kvs << std::endl;
}

void MrToolbox::remove_redundant_string_manual(MapReduce* mr) {
	mr -> reduce(emit_unique_strings, NULL);
	mr -> convert();
}

void MrToolbox::emit_unique_strings(char* key, int keybytes, char* values, int nvalues, int* valuebytes, KeyValue* kv, void* ptr) {
	std::unordered_set<std::string> value_set;
	for (int i = 0, offset = 0; i < nvalues; offset += valuebytes[i], i++) {
		value_set.insert((std::string) (values + offset));
	}
	for (auto iter = value_set.begin(); iter != value_set.end(); iter++) {
		int valuesize = iter -> length() + 1;
		char value[valuesize];
		strcpy(value, iter -> c_str());
		kv -> add(key, keybytes, value, valuesize);
	}
}

void MrToolbox::open_kmv(MapReduce* mr) {
	mr -> reduce(emit_kvs, NULL);
}

void MrToolbox::emit_kvs(char* key, int keybytes, char* values, int nvalues, int* valuebytes, KeyValue* kv, void* ptr) {
	for (int i = 0, offset = 0; i < nvalues; offset += valuebytes[i], i++) {
		kv -> add(key, keybytes, values + offset, valuebytes[i]);
	}
}

void MrToolbox::redundant_agglomerator(char* key, int keybytes, char* values, int nvalues, int* valuebytes, KeyValue* kv, void* ptr) {
        int val_offset = 0;
        for (int val_index = 0; val_index < nvalues; val_index++)
        {
                char combinedkv[keybytes + 1 + valuebytes[val_index]];
                memcpy(combinedkv, key, keybytes);
                combinedkv[keybytes] = SEPARATOR;
                memcpy(combinedkv + keybytes + 1, values + val_offset, valuebytes[val_index]);
                kv -> add(combinedkv, keybytes + 1 + valuebytes[val_index], NULL, 0);
                val_offset += (valuebytes[val_index]);
        }
}

void MrToolbox::redundant_remover(char* key, int keybytes, char* values, int nvalues, int* valuebytes, KeyValue* kv, void* ptr) {
        int separator_pos = -1;
        for (int i = 0; i < keybytes; i++)
        {
                if ( key[i] == SEPARATOR)
                {
                        separator_pos = i;
                        break;
                }
        }
        char splitted_key[separator_pos];
        char splitted_val[keybytes - separator_pos - 1];
        memcpy(splitted_key, key, separator_pos);
        memcpy(splitted_val, key + separator_pos + 1, keybytes - separator_pos - 1);
        kv -> add(splitted_key, separator_pos, splitted_val, keybytes - separator_pos - 1);
}
/*
void MrToolbox::set_mr_params(MapReduce* mr, int pagesize, int timer, int verbosity, int outofcore) {
	mr -> memsize = pagesize;
	mr -> verbosity = verbosity;
	mr -> timer = timer;
	mr -> outofcore = outofcore;
}
*/
void MrToolbox::set_mr_params(MapReduce* mr, int pagesize) {
	//if (pagesize > 1024) pagesize = 1024;
	mr -> verbosity = 0;
	mr -> memsize = pagesize;
	mr -> timer = 0;
	mr -> outofcore = -1;
}


