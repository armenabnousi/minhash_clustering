#pragma once
#include "mpi.h"
#include "mapreduce.h"
#include "keyvalue.h"
#include <string>
#include "limits.h"
#include "sstream"
#include <fstream>
#include "iostream"
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <math.h>
#include <cstdio> //for deleting existing cluster files - remove()
#include <algorithm> //for reversing std::vector
#include <sys/stat.h>
#include "fasta_util.h"
#include "message_value.h"
#include "union_find.h"
#include "mpi_util.h"
#include "mr_util.h"
#include "shingling_toolbox.h"
using namespace MAPREDUCE_NS;

#define ALPHABET_SIZE 20
#define NUMBERED_NAMES 1        //WILL ACCEPT UP TO 1000,000 SEQUENCES, EACH WITH MAXIMUM 998 DOMAINS.
#define SEPARATOR '!'
#define HTM_CC_LABEL_INDICATOR 'l'

typedef struct Clusters {
	std::vector< std::set <std::string> >* clusters;
	std::vector< std::string >* cluster_labels;
	Clusters(std::vector< std::set< std::string> >* clusters, std::vector< std::string >* labels) : clusters(clusters), 
													cluster_labels(labels) { }
} Clusters;

typedef struct DistributeUtil {
	Clusters clusters;
	bool single_proc; //container that should be distributed is only on one processor or all
	DistributeUtil(Clusters clusters, bool single_proc) : clusters(clusters), single_proc(single_proc) { }	
} DistributeUtil;

typedef struct ClusterLimits {
	int min_cluster_size_limit;
	int largest_cluster_size;
	ClusterLimits(int min_limit) : min_cluster_size_limit(min_limit), largest_cluster_size(0) { }
} ClusterLimits;

typedef struct MpiGatherParameters {
	int local_clusters_size;
	int global_clusters_size;
	int* proc_clusters_sizes = NULL;
	int num_processors;
	int root;
	int me;
	~MpiGatherParameters() {
		if (proc_clusters_sizes) delete[] proc_clusters_sizes;
	}
} MpiGatherParameters;

typedef struct MpiGatherToolbox {
	MpiGatherParameters* gather_params;
	std::vector< std::set<std::string> >* clusters;
	std::vector< std::string >* cluster_labels;
	int largest_cluster_size;
	MpiGatherToolbox(MpiGatherParameters* gather_params, std::vector<std::set<std::string>>* clusters, 
				std::vector<std::string>* cluster_labels) : 
				gather_params(gather_params), clusters(clusters), cluster_labels(cluster_labels) { }
} MpiGatherToolbox;

class ShinglingClusterGen {
	ShinglingToolbox shingling_params;
	int shingling_iterations;
	int connected_component_iterations;
	std::string output_filename;
	std::string connected_component_algorithm;
	static int me;
	static int log_proc;
	static int log_scale;
	//static int latest_hashed;
	static std::unordered_map<std::string, int> seq_owner_proc;
	//static std::unordered_map<std::string, std::unordered_set<std::string>> persistent_graph;
public:
	//static std::unordered_map<std::string, std::vector<std::vector<unsigned int>>*> persistent_hash_function_min_values;
	//static std::unordered_map<std::string, std::vector<std::vector<unsigned int>>*> persistent_vertex_hashes;
	ShinglingClusterGen(int shingling_iters, int cc_iters, 
						std::string output_prefix, std::string cc_alg, int proc_me, int log_scale, int log_proc,
						int prime_num, int num_hash, int shingle_size, int kmer_len, char* random_filename, 
						bool fixed_randoms, int max_hash, bool use_persistent_hashing);
	ShinglingClusterGen(int shingling_iters, int cc_iters, 
						std::string output_prefix, std::string cc_alg, int prime_num, int num_hash, 
						int shingle_size, int kmer_len, char* random_filename, bool fixed_randoms, int max_hash, 
						bool use_persistent_hashing);
	ShinglingClusterGen() { }
	~ShinglingClusterGen();
	void update_hashnum(int hashnum);
	void cluster(MapReduce* mr, std::string input_fname, int nmap, int sepdelta, 
					std::vector< std::set<std::string> >* clusters, int* largest_cluster_size, bool rename, 
					bool collect_in_root, MpiGatherParameters gather_params, int min_cluster_limit, 
					std::unordered_map<std::string, std::string>* local_sequences, 
					std::unordered_set<std::string> clusters_to_remove,
					double* kmer_shingling_time = NULL, double* vertex_shingling_time = NULL, 
					double* connected_comp_time = NULL, double time_breakdowns[5] = NULL, 
					std::string finalized_seq_list_filename = "");
	void mrmpi_cluster(MapReduce* mr, std::string input_fname, int nmap, int sepdelta, 
				std::unordered_map<std::string, std::string>* local_sequences, 
				MpiGatherParameters gather_params, bool collect_in_root,
				std::vector< std::set<std::string> >* final_clusters,
				double* kmer_shingling_time, double* vertex_shingling_time, double* connected_comp_time,
				double time_breakdowns[5], std::string finalized_seq_list_filename);
	std::string set_output_filename(ShinglingToolbox* params, std::string prefix, int shingling_iters, int cc_iters);
	static void read_proteins(int, char *, int, KeyValue *, void *hash_parameters);
	static void establish_sequence_ownership(char* name, int name_length, char* sequence, int seq_length, void* local_seqs_ptr);
	//static void populate_seq_owners_proc(std::unordered_map<std::string, std::string>* local_seqs);
	static void generate_local_seqs_kmer_shingles(int itask, KeyValue* kv, void* shingling_params);
	static void generate_kmer_min_shingles(KeyValue* kv, std::string seq_name, std::string sequence, 
						ShinglingToolbox* shingling_parameters, int itask, bool use_persistent_hashing);
	static void initialize_array_by_max_uint(unsigned int*** array, int rows, int cols);
	static unsigned long long compute_kmer_int_value(std::string kmer, std::string previous_kmer, 
							unsigned long long previous_kmer_value, int kmer_len);
	static int compute_hash_value(unsigned long long kmer_int_value, int random_num1, int random_num2, int Prime_num);
	static void update_hash_func_min_values(unsigned int kmer_hashed_value, unsigned int* hash_function_min_values, int shingle_size);
	static void update_hash_func_min_values(unsigned int hashed_value, std::vector<std::vector<unsigned int>>* hash_table, 
										int row, int max_col);
	static void emit_min_shingle_kvs(unsigned int** hash_function_min_values, KeyValue* kv, int num_hash, int shingle_size, 
						char* name, int name_length);
	static void emit_min_shingle_kvs(std::vector<std::vector<unsigned int>> hash_function_min_values, KeyValue* kv, 
                                                int num_hash, int shingle_size, char* name, int name_length);
	static void emit_min_shingle_kvs(std::vector<std::vector<unsigned int>>* hash_table, KeyValue* kv, int nrows, int ncols,
                                                        char* name, int name_length);
	static void delete_2d_array(unsigned int** array, int rows, int cols);
	static void draw_graph(char* key, int keybytes, char* values, int nvalues, int* valuebytes, KeyValue *kv, void* ptr);
	static void redundant_agglomerator(char* key, int keybytes, char* values, int nvalues, int* valuebytes, KeyValue* kv, void* ptr);
	static void redundant_remover(char* key, int keybytes, char* values, int nvalues, int* valuebytes, KeyValue* kv, void* ptr);
	static void generate_vertex_neighbor_min_shingles(char* sequence_name, int name_length, char* neighbors_list, int number_of_neighbors, 
							int* neighbor_name_lengths, KeyValue* kv, void* shingling_parameters);
	static unsigned long long compute_numbered_protname_int_value(std::string name);
	static void update_persistent_vertex_hashes(unsigned int neighbor_name_hash_value, int hash,
                                unordered_map<std::string, std::vector<std::vector<unsigned int>>*>* persistent_vertex_hashes,
                                std::string current_seq, int ncols);
	static void label_propagate(char* key, int keybytes, char* values, int nvalues, int* valuebytes, KeyValue* kv, void* ptr);
	static void extract_messages_info(char* messages, int* message_lengths, int num_messages, MessageValue* current_label, 
					MessageValue* min_label, std::vector<MessageValue>* neighbors, int* neighbor_count);
	static void label_finalize(char* key, int keybytes, char* values, int nvalues, int* valuebytes, KeyValue* kv, void* ptr);
	static void rename_sequences(char* cluster_label, int label_length, char* seq_names, int num_members, int* seq_name_lengths, 
				KeyValue* kv, void* ptr);
	static void remove_small_clusters(char* cluster_label, int label_length, char* member_names, int num_members,
						int* member_name_lengths, KeyValue* kv, void* min_limit);
	static std::vector< std::set<std::string> > collect_clusters(MapReduce* mr, MpiGatherParameters gather_params, 
									int* largest_cluster_size, std::vector<std::string>* cluster_labels);
	static void local_cluster_agglomerator(char* cluster_name, int cluster_name_length, char* members, 
							int num_members, int* member_name_lengths, KeyValue* kv, void* gather_parameters);
	static void send_local_clusters(char* empty_key, int keybytes, char* agglomerated_clusters, int nvalues_one, 
							int* agg_clusters_size, KeyValue* kv, void* gather_toolbox_ptr);
	static void extract_clusters_from_array(char* clusters_array, int array_size, MpiGatherToolbox* gather_toolbox);
	inline bool file_exists (const std::string& name);
	std::vector< std::set<std::string> > read_shingling_clusters(std::string filename, int* largest_cluster_size, bool rename, 
									int min_cluster_limit, std::vector<std::string>* cluster_labels,
									std::unordered_set<std::string> clusters_to_remove);
	static void distribute_clusters(int itask, KeyValue* kv, void* clusters_ptr);
	std::vector< std::string > split_string(std::string str, char delim);
	static void add_label_to_cluster(char* key, int keylength, char* values, int nvalues, int* valuebytes, KeyValue* kv, void* ptr);
	void perform_label_propagate(MapReduce* mr, int nnodes, int nedges, int connected_component_iteration, double* connected_comp_time); 
	void perform_union_find(MapReduce* mr, std::vector<std::set<std::string>>* final_clusters, bool collect_in_root,
					MpiGatherParameters gather_params);
	void perform_hash_to_min_connected_component(MapReduce* mr, int nnodes, int nedges, int connected_component_iterations, 
						double* connected_comp_time);
	void static hash_to_min_hasher(char* node, int node_length, char* members, int nmembers, int* member_lengths,
					KeyValue* kv, void* label_unchanged);
	void static find_label_and_min_neighbor(std::unordered_set<std::string>* neighbors, std::string* min_neighbor, std::string* label);
	void static save_completed_components(char* node, int node_length, char* members, int nmembers, int* member_lengths, 
						KeyValue* kv, void* ptr);
	int static send_to_owner_proc(char* key, int keybytes);
	void static add_special_characters(char* cluster_label, int label_size, char* members, int num_members, int* member_lengths,
						KeyValue* kv, void* ptr);
	//static void clean_up_static_members();
	void remove_finalized_sequences(std::unordered_map<std::string, std::string>* local_sequences, 
						std::string finalized_seq_list_filename);
	static void print_details(char* key, int key_length, char* values, int nvalues, int* valuebytes, void* ptr) {
		if (nvalues > 1) {
			int me; MPI_Comm_rank(MPI_COMM_WORLD, &me);
			std::cerr << me << "\tprinting\t" << key << "\t" << nvalues << std::endl;
		}
	}
};
