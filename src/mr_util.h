#pragma once
#include "mpi.h"
#include "mapreduce.h"
#include "keyvalue.h"
#include <list>
#include <vector>
#include <unordered_set>
using namespace MAPREDUCE_NS;

#define SEPARATOR '!'

typedef struct CollectionContainer {
	int local_size;
	std::list< std::vector<char> > local_list;
	CollectionContainer(int initial_size) : local_size(initial_size) { }
} CollectionContainer;

class MrToolbox {

public:
	void static gather_kmv(MapReduce* mr, int procs = 1, bool convert = true, int* num_kmvs = NULL, int* num_kvs = NULL);
	void static collect_keys(MapReduce* mr, std::vector<char>* collection);
	void static remove_redundant(MapReduce* mr);
	void static remove_redundant_local(MapReduce* mr, bool is_string);
	void static remove_redundant_string_manual(MapReduce* mr);
	void static open_kmv(MapReduce* mr);
	void static emit_kvs(char* key, int keybytes, char* values, int nvalues, int* valuebytes, KeyValue* kv, void* ptr);
	void static collect_local_keys(char* key, int keybytes, char* values, int nvalues, int* valuebytes, void* toolbox_ptr);
	void static prepare_array(CollectionContainer* container, char keys_array[]);
	void static redundant_agglomerator(char* key, int keybytes, char* values, int nvalues, int* valuebytes, KeyValue* kv, void* ptr);
	void static redundant_remover(char* key, int keybytes, char* values, int nvalues, int* valuebytes, KeyValue* kv, void* ptr);
	void static emit_unique_strings(char* key, int keybytes, char* values, int nvalues, int* valuebytes, KeyValue* kv, void* ptr);
	//void static set_mr_params(MapReduce* mr, int pagesize, int timer, int verbosity, int outofcore);
	void static set_mr_params(MapReduce* mr, int pagesize);
private:
	void static convert_kmv_to_kv(char* key, int keybytes, char* values, int nvalues, int* valuebytes, KeyValue* kv, void* int_ptr);
};
