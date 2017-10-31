#pragma once
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <fstream>

class ShinglingToolbox {
public:
        int* hash_randoms = NULL;
        int Prime_num;
        int num_hash;
	int max_hash;
        int shingle_size;
        int kmer_len;
	int latest_hashed = 0;
	bool multishingle;
	bool last_iter = false;
	bool use_persistent_hashing = true;
	std::unordered_map<std::string, std::string>* local_proteins;
	std::unordered_map<std::string, std::vector<std::vector<unsigned int>>*> persistent_hash_function_min_values;
	std::unordered_map<std::string, std::vector<std::vector<unsigned int>>*> persistent_vertex_hashes;
	std::unordered_map<std::string, std::unordered_set<std::string>> persistent_graph;
	std::unordered_map<std::string, unsigned long long> vertex_hash;
	ShinglingToolbox() : hash_randoms(NULL), num_hash(0) { }
	ShinglingToolbox(int prime_num, int num_hash, int shingle_size, int kmer_len, char* random_filename, int max_hash, 
				bool use_persistent_hashing = true);
	~ShinglingToolbox();
	void set_hash_randoms(int* rand_array, int max_hash, char* filename);
};
