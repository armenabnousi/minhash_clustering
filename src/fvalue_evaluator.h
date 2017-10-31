#pragma once
#include "mpi.h"
#include "mapreduce.h"
#include "keyvalue.h"
#include "nmi_messages.h"
#include "message_value.h"
#include "mr_util.h"
#include <vector>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <math.h>
#include <algorithm>
#include <float.h>
#include <iostream>
#include "sstream"
#include <fstream>
using namespace MAPREDUCE_NS;

#define SIZE_LABEL '$'

typedef struct NmiToolbox {
	int seq_count;
	double H_X;
	double H_Y;
	double H_X_cond_Y;
	double H_Y_cond_X;
	bool working_on_X;
	int member_count_X;
	int member_count_Y;
	KeyValue* kv; //this will be a pointer to the kv for intersection_mr
	NmiToolbox(int seq_count) : seq_count(seq_count), H_X(0), H_Y(0), H_X_cond_Y(0), H_Y_cond_X(0), working_on_X(true), 
					member_count_X(0), member_count_Y(0) { }
} NmiToolbox;

typedef struct ClusterReadToolbox {
	std::vector< std::unordered_set <std::string> >* clusters;
	std::vector< std::string >* cluster_labels;
	bool unlabeled;
	int* size_limit;
	ClusterReadToolbox(std::vector< std::unordered_set<std::string> >* clusters, bool unlabeled, 
				int* size_limit = NULL, std::vector< std::string >* labels = NULL) : 
				clusters(clusters), unlabeled(unlabeled), cluster_labels(labels), size_limit(size_limit) { }
} ClusterReadToolbox;

typedef struct ClusterSizeInspectionUtil {
	int min_allowed_size;
	int* largest_cluster_size = NULL;
	ClusterSizeInspectionUtil(int min_allowed_size, int* largest) : min_allowed_size(min_allowed_size), largest_cluster_size(largest) { }
} ClusterSizeInspectionUtil;

class ClusterEvaluator {
	int seq_count;
	int log_scale;
	int log_proc;
	int me;
	static int largest_shingling_cluster_size;
	static int largest_pfam_cluster_size;
public:
	ClusterEvaluator(int num_seq) : seq_count(num_seq) { }
	ClusterEvaluator(int num_seq, int log_scale, int log_proc, int me) : seq_count(num_seq), log_scale(log_scale), 
										log_proc(log_proc), me(me) { }
	void compute_nmi(std::vector< std::set<std::string> > test_clusters, 
				std::vector< std::set<std::string> > base_clusters, 
				double* nmi_lfk, double* nmi_max, int print_match_stats_cluster_size_limit);
	void compute_nmi(MapReduce* test_mr, MapReduce* base_mr, int root_processor, double* nmi_lfk, double* nmi_max);
	void compute_nmi(MapReduce* test_mr, MapReduce* base_mr, int root_processor, double* nmi_lfk, double* nmi_max, bool rename, 						int min_cluster_size_limit, int* largest_cluster_size);
	static void compute_local_cluster_H_X(char* label, int label_length, char* members, int num_members,
                                                                int* member_lengths, void* local_nmi_metrics);
	static void emit_cluster_members_as_key(char* label, int label_length, char* members, int num_members,
								int* member_lengths, void* intersection_mr_ptr);
	static void generate_all_pairs(char* null_key, int keybytes, char* all_labels, int clusters_count,
								int* message_lengths, KeyValue* kv, void* ptr);
	static void couple_intersecting_clusters(char* common_member, int member_length, char* labels, int num_labels,
								int* label_lengths, KeyValue* kv, void* ptr);
	static void compute_pair_entropy(char* paired_labels, int key_length, char* paired_cluster_sizes, int intersection_size,
								int* valuebytes, KeyValue* kv, void* ptr);
	static void compute_single_cluster_entropy(char* cluster_label, int label_length, char* pairwise_entropies, 
					int paired_clusters_count, int* size_of_double, void* ptr);
	static void compute_clusterings_conditional_entropy(char* clustering_symbol, int one, char* single_cluster_entropies, 
								int num_clusters, int* size_of_double, KeyValue* kv, void* root_ptr);
	void read_pfam_clusters(std::string filename, int* largest_cluster_size, int min_cluster_size_limit, MapReduce* pfam_mr, int nmap, bool domained);
	void read_shingling_clusters(std::string filename, int* largest_cluster_size, bool rename, int min_cluster_size, 
					MapReduce* shingling_mr, int nmap, int largest_clusters_count, std::string clustering_type);
	void compute_nmi(std::vector<std::unordered_set<std::string>>* test_clusters, 
					std::vector<std::unordered_set<std::string>>* base_clusters, double* nmi_lfk, double* nmi_max, 
					int min_cluster_size, int largest_clusters_count, int* largest_test_cluster_size,
					int* largest_base_cluster_size);
	void static distribute_clusters(int itask, KeyValue* kv, void* clusters_ptr);
	void read_pclust_clusters(std::string filename, int* largest_cluster_size, int min_cluster_limit,
					MapReduce* pclust_mr, int nmap, bool rename);
	void read_pclust_clusters_on_all_procs(std::string filename, int min_cluster_limit, int* largest_cluster_size,
						std::vector< std::unordered_set<std::string> >* clusters, bool rename);
private:
	void compute_entropies(std::vector< std::set<std::string> > X, 
				std::vector< std::set<std::string> > Y, double* H_X, double* H_X_cond_Y, 
				int print_match_stats_cluster_size_limit);
	void compute_H_X(std::vector< std::set<std::string> > X, double* H_X);
	void compute_H_X_cond_Y(std::vector< std::set<std::string> > X, std::vector< std::set<std::string> > Y, double* H_X_cond_Y,
					int print_match_stats_cluster_size_limit);
	static double compute_h(int w, int n);
	static double compute_H_star_Xi_cond_Yj(int a, int b, int c, int d, int seq_count, double* H_Xi_cond_Yj);
	void aggregate_local_metrics(NmiToolbox* nmi_toolbox, int root, int me);
	void compute_nmi_from_entropies(double H_X, double H_Y, double H_X_cond_Y, double H_Y_cond_X, double* nmi_lfk, double* nmi_max);
	void static rename_sequences(char* cluster_label, int label_length, char* seq_names, int num_members,
                                                int* seq_name_lengths, KeyValue* kv, void* ptr);
	void static remove_small_clusters(char* cluster_label, int label_length, char* member_names, int num_members,
                                                int* member_name_lengths, KeyValue* kv, void* cluster_size_limit);
	std::vector< std::unordered_set<std::string> > read_pfam_clusters_on_all_proc(std::string filename, int* largest_cluster_size, int min_cluster_size_limit, bool wholeseq);
	std::vector< std::string > split_string(std::string str, char delim);
	//void static distribute_clusters(int itask, KeyValue* kv, void* clusters_ptr);
	std::vector< std::unordered_set<std::string> > read_shingling_clusters_on_all_proc(std::string filename, int* largest_cluster_size, 
											bool rename, int min_cluster_limit, 
											std::vector<std::string>* cluster_labels);
	void find_largest_clusters(std::vector<std::unordered_set<std::string>>* shingling_clusters,
					std::vector<int>* selected_clusters_indices, int count);
	void static send_intersection_sizes(char*, int, char*, int, int*, void*);
	void static compute_cluster_entropies(char*, int, char*, int, int*, void*);
	void convert_clusters_to_mr(MapReduce* mr, std::vector<std::unordered_set<std::string>>* clusters, int min_cluster_size, 
					int largest_clusters_count, int* largest_cluster_size);
	void find_largest_clusters(std::vector<std::unordered_set<std::string>>* shingling_clusters, int count, 
						std::vector< std::unordered_set< std::string >>* selected_clusters);
};
