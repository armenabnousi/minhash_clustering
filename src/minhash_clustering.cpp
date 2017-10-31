#include "mpi.h"
#include "mapreduce.h"
#include "keyvalue.h"
#include "shingling_cluster_generator_optimized.h"
#include <string.h>
#include "iostream"
#include <vector>
#include <set>
#include <list>
#include <math.h>
#include <algorithm> //for set_intersection(), set_difference() and distance()
#include "fvalue_evaluator.h"
#include "fasta_util.h"
//#include "merge_evaluator.h"
#include "graph_formater.h"
#include <iomanip>
#include <unistd.h>
#include <stdio.h>
using namespace std;
using namespace MAPREDUCE_NS;

#define LOG_SCALE 1
#define LOG_PROC 0
#define SEPARATOR '!'		//USED IN REMOVING REDUNDANT KV's
#define ROOT_PROCESSOR 0
#define PRIME_NUM 1280000003

void parse_input_args(char** args, int* nmap, int* pagesize, int* sepdelta, std::string* input_fname, std::string* output_prefix, 
			int* seq_count, int* kmer_len, int* shingle_size, int* max_num_hash, int* shingling_iterations, 
			int* connected_component_iterations, std::string* nmi_output_filename, int* min_cluster_size, bool* rename,
			bool* collect_in_root, int* print_match_stats_cluster_size_limit, std::string* connected_component_algorithm,
			double* sim_score_threshold, double* sim_length_threshold, string* similarity_function, string* finalized_directory,
			std::string* community_detection_method, std::string* setname, int* start_hash_num, std::string* metis_prefix,
			std::string* fvalue_filename, double* fvalue_threshold, int* hash_difference, std::string* grappolo_command,
			std::string* grappolo_script);
void print_nmi_values(std::string nmi_output_filename, int largest_test_cluster_size, double nmi_max, double nmi_lfk, 
                        int test_num_hash, int base_num_hash, double test_cluster_time, double base_cluster_time, double nmi_time, 
			double merge_time, double previous_nmi, int clusters_to_remove, double sim_score_threshold, 
			double sim_length_threshold, int seqs_to_align, int sequences_to_remove_count, int remaining_sequences_count,
			double test_kmer_shingle_time, double base_kmer_shingle_time, double test_vertex_shingle_time, 
			double base_vertex_shingle_time, double test_cc_time, double base_cc_time, double time_breakdowns[5]);
unordered_map<std::string, std::string> read_all_proteins(std::string filename);
std::string generate_metis_filename(std::string metis_prefix, std::string fasta_filename, int num_hash, int shingling_iterations, int kmer_len);
void print_sizes(char* key, int keybytes, char* values, int nvalues, int* valubytes, void* ptr) {
	std::cout << key << "\t" << nvalues << std::endl;
}

int me;

/* ---------------------------------------------------------------------- */

int main(int narg, char **args)
{
	//std::cout << narg << " " << args[22] << std::endl;
	if (narg != 31) {
		printf("ERROR: Usage profiler.exe (number of maps) (page size) (separation delta) (input file name) (output prefix) (number of sequences) (kmer length) (shingle size) (max number of hash functions) (number of shingling iterations) (number of test connected component iterations) (nmi result filename) (min cluster size) (rename) (collect clusters in root) (min size of cluster that we want its match stats be printed) (connected component algorithm) (score threshold) (length threshold) (similarity function) (finalized clusters directory) (community detection method) (setname) (start hash num) (metis prefix) (f1-value filename) (fvalue threshold) (hash difference) (grappolo command prefix) (grappolo script name)");
		return(1);
	}
	int clusters_to_remove_from_step = 0;
	int num_seqs_to_align = 5;
	bool allow_repeats = false;
	std::string finalized_filename_prefix;
	//std::cout << finalized_filename_prefix << std::endl;
	double test_cluster_time = 0, test_kmer_shingle_time = 0, test_cc_time = 0, base_cluster_time = 0, base_kmer_shingle_time = 0, 
		base_cc_time = 0, nmi_time = 0,	merge_time = 0, t_start = 0, t_end = 0, test_vertex_shingle_time = 0, 
		base_vertex_shingle_time = 0;
	double sim_score_threshold = 0, sim_length_threshold = 0;
	double nmi_lfk = 0, nmi_max = 0;
	MPI_Init(&narg, &args);
	int nmap = 0, pagesize = 0, sepdelta = 0, seq_count = 0, kmer_len = 0, shingle_size = 0, test_num_hash = 0, base_num_hash = 0, shingling_iterations = 0, 
		connected_component_iterations = 0, min_cluster_limit = 0, max_num_hash = 0;
	bool rename = true, collect_in_root = false;
	int print_match_stats_cluster_size_limit = 0;
	std::string input_fname = "", output_prefix = "", nmi_output_filename = "", connected_component_algorithm = "";
	std::string similarity_function = "", finalized_directory = "";
	std::string community_detection_method = "";
	std::string setname, metis_prefix, fvalue_filename, grappolo_command, grappolo_script;
	int start_hash_num = 0, hash_difference = 0;
	double fvalue_threshold;
	parse_input_args(args, &nmap, &pagesize, &sepdelta, &input_fname, &output_prefix, &seq_count, &kmer_len, &shingle_size, 
				&max_num_hash, &shingling_iterations, &connected_component_iterations, 
				&nmi_output_filename, &min_cluster_limit, &rename, &collect_in_root, &print_match_stats_cluster_size_limit,
				&connected_component_algorithm, &sim_score_threshold, &sim_length_threshold, &similarity_function,
				&finalized_directory, &community_detection_method, &setname, &start_hash_num, &metis_prefix, 
				&fvalue_filename, &fvalue_threshold, &hash_difference, &grappolo_command, &grappolo_script);
	test_num_hash = start_hash_num;
	base_num_hash = test_num_hash - hash_difference;
	if (allow_repeats) {
		finalized_filename_prefix = finalized_directory + "/finalized_clusters_freerepeat";
	} else {
		finalized_filename_prefix = finalized_directory + "/finalized_clusters_norepeat";
	}

	MPI_Comm_rank(MPI_COMM_WORLD, &me);
	int num_processors = 0;
	MPI_Comm_size(MPI_COMM_WORLD, &num_processors);
	//setting parameters for cluster generator instantiation
	int max_hash_possible = 2000;
	bool use_persistent_hashing = true;
	//ShinglingParameters test_shingling_params(PRIME_NUM, test_num_hash, shingle_size, kmer_len, "random_hash_abs", 2000, 
	//						use_persistent_hashing);
	//ShinglingParameters base_shingling_params(PRIME_NUM, base_num_hash, shingle_size, kmer_len, "random_hash_abs", 2000,
	//						use_persistent_hashing);

	MapReduce *test_mr = new MapReduce(MPI_COMM_WORLD);
        MrToolbox::set_mr_params(test_mr, pagesize);
        MapReduce *base_mr = new MapReduce(MPI_COMM_WORLD);
        MrToolbox::set_mr_params(base_mr, pagesize);
	MpiGatherParameters gather_params;
        gather_params.num_processors = num_processors;
        gather_params.root = ROOT_PROCESSOR;
        gather_params.me = me;

	//std::unordered_map<std::string, std::string> all_proteins = read_all_proteins(input_fname);
        int largest_test_cluster_size = 0, largest_base_cluster_size = 0;
        std::vector< std::set<std::string> > test_clusters;
        std::vector< std::set<std::string> > base_clusters;
	std::unordered_set<std::string> clusters_to_remove;
	std::unordered_map<std::string, std::string> local_sequences;

	double max_test_cluster_time = 0, max_base_cluster_time = 0, max_nmi_time = 0, max_merge_time = 0, max_test_kmer_shingle_time = 0, 
		max_base_kmer_shingle_time = 0, max_test_vertex_shingle_time = 0, max_base_vertex_shingle_time = 0, max_test_cc_time = 0, max_base_cc_time = 0;
	int remaining_sequences_count = 0, sequences_to_remove_count = 0;
	double previous_nmi = 0.0;
	bool rerun_base = true, increment_hash = true;
	int loop_count = 0;
	//std::cout << "starting loop" << std::endl;
	//std::cout << "calling the constructor" << std::endl;
	ShinglingClusterGen cluster_gen_test(shingling_iterations, 
						connected_component_iterations, output_prefix, connected_component_algorithm, 
						me, LOG_SCALE, LOG_PROC, PRIME_NUM, test_num_hash, shingle_size, kmer_len, 
						"random_hash_abs", 2000, use_persistent_hashing);
	ShinglingClusterGen cluster_gen_base(shingling_iterations, 
						connected_component_iterations, output_prefix, connected_component_algorithm, 
						me, LOG_SCALE, LOG_PROC, PRIME_NUM, base_num_hash, shingle_size, kmer_len, 
						"random_hash_abs", 2000, use_persistent_hashing);
	//std::cout << "constructor completed" << std::endl;
	double fvalue = 0;
	for (test_num_hash, base_num_hash; test_num_hash < max_num_hash && fvalue < fvalue_threshold; test_num_hash++, base_num_hash++, loop_count++) {
		if (me == 0) remove(fvalue_filename.c_str());
		cluster_gen_test.update_hashnum(test_num_hash);
		cluster_gen_base.update_hashnum(base_num_hash);
		//test_mr -> timer = 1;
		test_cluster_time = 0, base_cluster_time = 0, nmi_time = 0, merge_time = 0;
		test_kmer_shingle_time = 0, base_kmer_shingle_time = 0, test_vertex_shingle_time = 0, base_vertex_shingle_time = 0;
		test_cc_time = 0, base_cc_time = 0;
		//std::string clfname = "/home/armen.abnousi/shingle/xv_checkhash" + std::to_string(test_num_hash) + "_doubled";
		//char tc_filename[200];
		//strcpy(tc_filename, clfname.c_str());
		MPI_Barrier(MPI_COMM_WORLD);
		

		//TEST SET GRAPH GENERATION
		double test_time_breakdowns[5] = {0, 0, 0, 0, 0};
		if (LOG_SCALE > 0 && (LOG_PROC == me || LOG_PROC == -1)) std::cout << "clustering test set " << test_num_hash << std::endl;
		t_start = MPI_Wtime();
		cluster_gen_test.cluster(test_mr, input_fname, nmap, sepdelta, &test_clusters, &largest_test_cluster_size, rename,
						collect_in_root, gather_params, min_cluster_limit,
						&local_sequences, clusters_to_remove, &test_kmer_shingle_time, &test_vertex_shingle_time,
						&test_cc_time, test_time_breakdowns);
		t_end = MPI_Wtime();
		test_cluster_time = t_end - t_start;
		
		t_start = MPI_Wtime();
		std::string fasta_filename = input_fname.substr(input_fname.rfind('/') + 1);
		std::string test_metis_filename = generate_metis_filename(metis_prefix, fasta_filename, test_num_hash, shingling_iterations, kmer_len);
		if (community_detection_method != "no_com_detection") {
			GraphFormater metis_generator;
			metis_generator.generate_metis_from_mapreduce_graph(test_mr, test_metis_filename, community_detection_method);
		}
		t_end = MPI_Wtime();
		double test_metis_generation_time = t_end - t_start;
		t_start = MPI_Wtime();
		MPI_Barrier(MPI_COMM_WORLD);

		//BASE SET GRAPH GENERATION
		double base_time_breakdowns[5] = {0, 0, 0, 0, 0};
		std::string base_metis_filename = generate_metis_filename(metis_prefix, fasta_filename, base_num_hash, shingling_iterations, kmer_len);
                if (LOG_SCALE > 0 && (LOG_PROC == me || LOG_PROC == -1)) std::cout << "clustering base set " << base_num_hash << std::endl;
                t_start = MPI_Wtime();
		if (!std::ifstream(base_metis_filename)) {
                	cluster_gen_base.cluster(base_mr, input_fname, nmap, sepdelta, &base_clusters, &largest_base_cluster_size, rename,
                                                collect_in_root, gather_params, min_cluster_limit,
                                                &local_sequences, clusters_to_remove, &base_kmer_shingle_time, &base_vertex_shingle_time,
                                                &base_cc_time, base_time_breakdowns);
                	t_end = MPI_Wtime();
                	base_cluster_time = t_end - t_start;

                	t_start = MPI_Wtime();
                	if (community_detection_method != "no_com_detection") {
                	        GraphFormater metis_generator;
                	        metis_generator.generate_metis_from_mapreduce_graph(base_mr, base_metis_filename, community_detection_method);
                	}
		} else {
			base_cluster_time = 0;
		}
                t_end = MPI_Wtime();
                double base_metis_generation_time = t_end - t_start;
 //               t_start = MPI_Wtime();
                MPI_Barrier(MPI_COMM_WORLD);

		//GATHERING RUNTIME INFORMATION
		double max_test_kmer_hashing, max_test_kmer_shingle_collating, max_test_kmer_graph, 
			max_test_costum_hash_collating, max_test_redundant_removal, max_test_metis_time, max_base_cluster_time;
		MPI_Reduce(&test_kmer_shingle_time, &max_test_kmer_shingle_time, 1, MPI_DOUBLE, MPI_MAX, ROOT_PROCESSOR, MPI_COMM_WORLD);
		MPI_Reduce(&test_vertex_shingle_time, &max_test_vertex_shingle_time, 1, MPI_DOUBLE, MPI_MAX, ROOT_PROCESSOR, MPI_COMM_WORLD);
		MPI_Reduce(&test_cluster_time, &max_test_cluster_time, 1, MPI_DOUBLE, MPI_MAX, ROOT_PROCESSOR, MPI_COMM_WORLD);
		MPI_Reduce(&test_time_breakdowns[0], &max_test_kmer_hashing, 1, MPI_DOUBLE, MPI_MAX, ROOT_PROCESSOR, MPI_COMM_WORLD);
		MPI_Reduce(&test_time_breakdowns[1], &max_test_kmer_shingle_collating, 1, MPI_DOUBLE, MPI_MAX, ROOT_PROCESSOR, MPI_COMM_WORLD);
		MPI_Reduce(&test_time_breakdowns[2], &max_test_kmer_graph, 1, MPI_DOUBLE, MPI_MAX, ROOT_PROCESSOR, MPI_COMM_WORLD);
		MPI_Reduce(&test_time_breakdowns[3], &max_test_costum_hash_collating, 1, MPI_DOUBLE, MPI_MAX, ROOT_PROCESSOR, MPI_COMM_WORLD);
		MPI_Reduce(&test_time_breakdowns[4], &max_test_redundant_removal, 1, MPI_DOUBLE, MPI_MAX, ROOT_PROCESSOR, MPI_COMM_WORLD);
		MPI_Reduce(&test_metis_generation_time, &max_test_metis_time, 1, MPI_DOUBLE, MPI_MAX, ROOT_PROCESSOR, MPI_COMM_WORLD);
		MPI_Reduce(&base_cluster_time, &max_base_cluster_time, 1, MPI_DOUBLE, MPI_MAX, ROOT_PROCESSOR, MPI_COMM_WORLD);
		if (me == 0) {
			std::cerr << "kmer_hashing: " << max_test_kmer_hashing << "\nkmer_collating: " << max_test_kmer_shingle_collating << "\nkmer_graph: " << max_test_kmer_graph << "\ncustom_hash_collating: " << max_test_costum_hash_collating << "\nredun_removal: " << max_test_redundant_removal << "\nmetis_generation: " << max_test_metis_time << "\n***kmer_shingling_time: " << max_test_kmer_shingle_time << "\n***vertex_shingling_time: " << max_test_vertex_shingle_time << "\n***cc_time: " << max_test_cc_time << "\n***test_cluster_time: " << max_test_cluster_time << "\n**base_cluster_time: " << max_base_cluster_time << std::endl;
		}

		delete test_mr;
		test_mr = new MapReduce(MPI_COMM_WORLD);
		MrToolbox::set_mr_params(test_mr, pagesize);
		delete base_mr;
		base_mr = new MapReduce(MPI_COMM_WORLD);
		MrToolbox::set_mr_params(base_mr, pagesize);
		if (false) {
		//////////////////////////////////////////
		if (me == 0) {
			t_start = MPI_Wtime();
			std::string grappolo_args = " \"" + std::to_string(test_num_hash) +
                                " " + test_metis_filename + " " + base_metis_filename + " " + fvalue_filename +
                                "\" ";
			std::string command = grappolo_command +  grappolo_args + grappolo_script + "; exit";
			//std::string command = "cd /home/armen.abnousi/shingle; qsub -F \"" + std::to_string(test_num_hash) + 
			//	" " + test_metis_filename + " " + base_metis_filename + " " + fvalue_filename + 
			//	"\" script_cluster_fval" + std::to_string(nmap) + ".sh; exit;";
			//std::cout << "calling the next job script using command\n\t" << command << std::endl;	
			system(command.c_str());
                }

		//if (me == 0) std::cout << "waiting for the other jobs to complete: " << fvalue_filename << std::endl;
		while(!std::ifstream(fvalue_filename)) {
			sleep(10);
			//if (me == 0) std::cout << "waiting for the other jobs to complete inside" << std::endl;
		}
		//std::cout << "jobs completed, reading the new fvalue" << std::endl;
		//read the fvalue
		ifstream infile(fvalue_filename);
		if (infile.good()) {
			std::string line;
			std::getline(infile, line);
			fvalue = atof(line.c_str());
		} else {
			std::cout << "cannot read fvalue file" << std::endl;
			exit(1);
		}
		infile.close();
		t_end = MPI_Wtime();
		double wait_time = 0;
		if (me == 0) {
			wait_time = t_end - t_start;
			std::cout << "\nwait time: " << wait_time << std::endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);
		}
		//////////////////////////////////////////
		if (me == 0) std::cout << "iteration completed...\tfvalue = " << fvalue << "\n******************\n" << std::endl;
	}
	delete test_mr;
	delete base_mr;
	MPI_Finalize();
	return(0);
}

void parse_input_args(char** args, int* nmap, int* pagesize, int* sepdelta, std::string* input_fname, std::string* output_prefix, 
                        int* seq_count, int* kmer_len, int* shingle_size, int* max_num_hash, int* shingling_iterations, 
			int* connected_component_iterations, std::string* nmi_output_filename, int* min_cluster_size, bool* rename,
			bool* collect_in_root, int* print_match_stats_cluster_size_limit, std::string* connected_component_algorithm,
			double* sim_score_threshold, double* sim_length_threshold, std::string* similarity_function, 
			std::string* finalized_directory, std::string* community_detection_method, std::string* setname, 
			int* start_hash_num, std::string* metis_prefix, std::string* fvalue_filename, double* fvalue_threshold,
			int* hash_difference, std::string* grappolo_command, std::string* grappolo_script) {
	*nmap = atoi(args[1]);
	*pagesize = atoi(args[2]);
	*sepdelta = atoi(args[3]);
	*input_fname = std::string(args[4]);
	*output_prefix = std::string(args[5]);
	*seq_count = atoi(args[6]);
	*kmer_len = atoi(args[7]);
	*shingle_size = atoi(args[8]);
	*max_num_hash = atoi(args[9]);
	*shingling_iterations = atoi(args[10]);
	*connected_component_iterations = atoi(args[11]);
	*nmi_output_filename = std::string(args[12]);
	int min_size_limit = atoi(args[13]);
	*min_cluster_size = min_size_limit > 2 ? min_size_limit : 2;
	*rename = atoi(args[14]);
	*collect_in_root = atoi(args[15]);
	*print_match_stats_cluster_size_limit = atoi(args[16]);
	*connected_component_algorithm = args[17];
	*sim_score_threshold = atof(args[18]);
	*sim_length_threshold = atof(args[19]);
	*similarity_function = args[20];
	*finalized_directory = args[21];
	*community_detection_method = args[22];
	*setname = args[23];
	*start_hash_num = atoi(args[24]);
	*metis_prefix = args[25];
	*fvalue_filename = args[26];
	*fvalue_threshold = atof(args[27]);
	*hash_difference = atoi(args[28]);
	*grappolo_command = args[29];
	*grappolo_script = args[30];
}

void print_nmi_values(std::string nmi_output_filename, int largest_test_cluster_size, double nmi_max, double nmi_lfk, 
			int test_num_hash, int base_num_hash, double test_cluster_time, double base_cluster_time, 
			double nmi_time, double merge_time, double previous_nmi, int clusters_to_remove, double sim_score_threshold,
			double sim_length_threshold, int seq_to_align, int sequences_to_remove_count, int remaining_sequences_count,
			double test_kmer_shingle_time, double base_kmer_shingle_time, double test_vertex_shingle_time,
			double base_vertex_shingle_time, double test_cc_time, double base_cc_time, double time_breakdowns[5]) {
	int me;
	MPI_Comm_rank(MPI_COMM_WORLD, &me);
	if (me == ROOT_PROCESSOR) {
                std::ofstream outfile;
                outfile.open(nmi_output_filename, std::ios_base::app);
		if (test_cluster_time < 0.001) test_cluster_time = 0;
		if (base_cluster_time < 0.001) base_cluster_time = 0;
		if (nmi_time < 0.001) nmi_time = 0;
		if (merge_time < 0.001) merge_time = 0;
		
		double total_time = test_cluster_time + base_cluster_time + nmi_time + merge_time;
		outfile.precision(5);
                outfile << test_num_hash << "\t" << base_num_hash << "\t" << test_kmer_shingle_time << "\t" << test_vertex_shingle_time <<
				"\t" << test_cc_time << "\t" << test_cluster_time << "\t" << base_kmer_shingle_time << "\t" <<
				base_vertex_shingle_time << "\t" << base_cc_time << "\t" << base_cluster_time << "\t" << 
				nmi_time << "\t" << merge_time << "\t" << total_time << "\t" << largest_test_cluster_size << "\t" << 
				clusters_to_remove << "\t" << sequences_to_remove_count << "\t" << remaining_sequences_count << "\t" <<
				sim_score_threshold << "\t" << sim_length_threshold << "\t" << seq_to_align << "\t" <<
				nmi_lfk << "\t" << previous_nmi << "\t" << nmi_max << "\t" << time_breakdowns[0] << "\t" 
				<< time_breakdowns[1] << "\t" << time_breakdowns[2] << "\t" << time_breakdowns[3] << "\t"
				<< time_breakdowns[4] << std::endl;
                outfile.close();
	}
}

unordered_map<std::string, std::string> read_all_proteins(std::string filename) {
	FastaToolbox fasta_toolbox;
	std::list<Protein> prots = fasta_toolbox.read_fasta(&filename[0]);
	unordered_map<std::string, std::string> seq_map;
	for (auto iter = prots.begin(); iter != prots.end(); iter++) {
		seq_map.insert(std::make_pair(iter -> name, iter -> sequence));
	}
	return seq_map;
}

std::string generate_metis_filename(std::string metis_prefix, std::string fasta_filename, int num_hash, int shingling_iterations, int kmer_len) {
	std::string metis_filename = metis_prefix + fasta_filename + "_hash" + std::to_string(num_hash) + "_m" + 
					std::to_string(shingling_iterations) + "_kmer" + std::to_string(kmer_len);
	return metis_filename;
}
