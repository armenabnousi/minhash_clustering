#include "mpi.h"
#include "mapreduce.h"
#include "keyvalue.h"
#include <unistd.h>
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
			std::string* grappolo_script, bool* fixed_randoms, bool* run_grappolo);

void parse_input_args(char** args, int* nmap, int* pagesize, int* sepdelta, std::string* input_fname, int* seq_count, int* kmer_len, 
			int* shingle_size, int* max_num_hash, int* shingling_iterations,
                        int* min_cluster_size, int* start_hash_num, std::string* metis_prefix,
                        std::string* fvalue_filename, double* fvalue_threshold, int* hash_difference, std::string* grappolo_command,
                        std::string* grappolo_script, bool* fixed_randoms, bool* run_grappolo);

void parse_input_args(int narg, char** args, int* nmap, int* pagesize, int* sepdelta, std::string* input_fname, int* seq_count, int* kmer_len,
                        int* shingle_size, int* max_num_hash, int* shingling_iterations,
                        int* min_cluster_size, int* start_hash_num, std::string* metis_prefix,
                        std::string* fvalue_filename, double* fvalue_threshold, int* hash_difference, std::string* grappolo_command,
                        std::string* grappolo_script, bool* fixed_randoms, bool* run_grappolo);

std::unordered_map<std::string, std::string> read_all_proteins(std::string filename);
std::string generate_metis_filename(std::string metis_prefix, std::string fasta_filename, int num_hash, int shingling_iterations, int kmer_len);
void print_sizes(char* key, int keybytes, char* values, int nvalues, int* valubytes, void* ptr) {
	std::cout << key << "\t" << nvalues << std::endl;
}

int me;

/* ---------------------------------------------------------------------- */

int main(int narg, char **args)
{
	//std::cout << narg << " " << args[22] << std::endl;
	//if (narg != 32) {
	if (narg != 20) {
		printf("ERROR: Usage profiler.exe (number of maps) (page size) (separation delta) (input file name) (number of sequences) (kmer length) (shingle size) (max number of hash functions) (number of shingling iterations) (min cluster size) (start hash num) (metis prefix) (f1-value filename) (fvalue threshold) (hash difference) (grappolo command prefix) (grappolo script name) (use fixed randoms) (run grappolo)")
		//printf("ERROR: Usage profiler.exe (number of maps) (page size) (separation delta) (input file name) (output prefix) (number of sequences) (kmer length) (shingle size) (max number of hash functions) (number of shingling iterations) (nmi result filename) (min cluster size) (rename) (collect clusters in root) (min size of cluster that we want its match stats be printed) (connected component algorithm) (score threshold) (length threshold) (similarity function) (finalized clusters directory) (community detection method) (setname) (start hash num) (metis prefix) (f1-value filename) (fvalue threshold) (hash difference) (grappolo command prefix) (grappolo script name)");
		return(1);
	}
	
	double test_cluster_time = 0, test_kmer_shingle_time = 0, test_cc_time = 0, base_cluster_time = 0, base_kmer_shingle_time = 0, 
		base_cc_time = 0, nmi_time = 0,	merge_time = 0, t_start = 0, t_end = 0, test_vertex_shingle_time = 0, 
		base_vertex_shingle_time = 0;
	MPI_Init(&narg, &args);
	int num_processors = 0;
	MPI_Comm_size(MPI_COMM_WORLD, &num_processors);
	
	int nmap = num_processors, pagesize = 4096, sepdelta = 30000, seq_count = 0, kmer_len = 6, shingle_size = 2;
	int shingling_iterations = 2,  min_cluster_limit = 10, max_num_hash = 2000;
	bool rename = false, collect_in_root = false, fixed_randoms = false;
	std::string input_fname = "", output_prefix = "",  connected_component_algorithm = "nocc", community_detection_method = "usc_louvain";
	std::string metis_prefix, fvalue_filename, grappolo_command, grappolo_script;
	int start_hash_num = 41, hash_difference = 40;
	double fvalue_threshold = 0.9;
	parse_input_args(args, &nmap, &pagesize, &sepdelta, &input_fname, &seq_count, &kmer_len, &shingle_size, 
				&max_num_hash, &shingling_iterations, &min_cluster_limit, &start_hash_num, &metis_prefix, 
				&fvalue_filename, &fvalue_threshold, &hash_difference, &grappolo_command, &grappolo_script, &fixed_randoms);
	int test_num_hash = start_hash_num;
	int base_num_hash = test_num_hash - hash_difference;

	MPI_Comm_size(MPI_COMM_WORLD, &num_processors);
	//setting parameters for cluster generator instantiation
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
	int connected_component_iterations = -1;
	//std::cout << "starting loop" << std::endl;
	//std::cout << "calling the constructor" << std::endl;
	ShinglingClusterGen cluster_gen_test(shingling_iterations, 
						connected_component_iterations, output_prefix, connected_component_algorithm, 
						me, LOG_SCALE, LOG_PROC, PRIME_NUM, test_num_hash, shingle_size, kmer_len, 
						"../data/random_hash_abs", fixed_randoms, max_num_hash, use_persistent_hashing);
	ShinglingClusterGen cluster_gen_base(shingling_iterations, 
						connected_component_iterations, output_prefix, connected_component_algorithm, 
						me, LOG_SCALE, LOG_PROC, PRIME_NUM, base_num_hash, shingle_size, kmer_len, 
						"../data/random_hash_abs", fixed_randoms, max_num_hash, use_persistent_hashing);
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
		if (run_grappolo) {
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

void parse_input_args(int narg, char** args, int* nmap, int* pagesize, int* sepdelta, std::string* input_fname, int* seq_count, int* kmer_len,
                        int* shingle_size, int* max_num_hash, int* shingling_iterations,
                        int* min_cluster_size, int* start_hash_num, std::string* metis_prefix,
                        std::string* fvalue_filename, double* fvalue_threshold, int* hash_difference, std::string* grappolo_command,
                        std::string* grappolo_script, bool* fixed_randoms) {
	while ((c = getopt (narg, argv, "m:p:q:f:n:k:s:e:r:l:h:o:v:t:d:g:c:x:u:")) != -1)	{
		switch(c): {
		case 'm': 
			*nmap = atoi(optarg);
			break;
		case 'p': 
			*pagesize = atoi(optarg);
			break;
		case 'q':
			*sepdelta = atoi(optarg);
			break;
		case 'f':
			*input_fname = optarg;
			break;
		case 'n':
			*seq_count = atoi(optarg);
			break;
		case 'k': //equivalent to w in (w,c)-sketch
			*kmer_len = atoi(optarg);
			break;
		case 's': //equivalent to c in (w,c)-sketch
			*shingle_size = atoi(optarg);
			break;
		case 'e':
			*max_num_hash = atoi(optarg);
			break;
		case 'r':
			*shingling_iterations = atoi(optarg);
			break;
		case 'l':
			*min_cluster_size = (atoi(optarg) > 2) ? atoi(optarg) : 2;
			break;
		case 'h': 
			*start_hash_num = atoi(optarg);
			break;
		case 'o':
			*metis_prefix = optarg;
			break;
		case 'v':
			*fvalue_filename = optarg;
			break;
		case 't':
			*fvalue_threshold = atof(optarg);
			break;
		case 'd':
			*hash_difference = atoi(optarg);
			break;
		case 'g':
			*grappolo_command = optarg;
			break;
		case 'c':
			*grappolo_script = optarg;
			break;
		case 'x':
			*fixed_randoms = (atoi(optarg) == 0) ? false : true;
			break;
		case 'u':
			*run_grappolo = (atoi(optarg) == 0) ? false: true;
		case '?':
			cout << "Error: " << "option -" << optopt << " not recognized or requires an argument!\n";
			return 1;
		default:
			abort();	
		}
	}
}

void parse_input_args(char** args, int* nmap, int* pagesize, int* sepdelta, std::string* input_fname, int* seq_count, int* kmer_len, 
			int* shingle_size, int* max_num_hash, int* shingling_iterations,
                        int* min_cluster_size, int* start_hash_num, std::string* metis_prefix,
                        std::string* fvalue_filename, double* fvalue_threshold, int* hash_difference, std::string* grappolo_command,
                        std::string* grappolo_script, bool* fixed_randoms) {
	*nmap = atoi(args[1]);
        *pagesize = atoi(args[2]);
        *sepdelta = atoi(args[3]);
        *input_fname = std::string(args[4]);
        *seq_count = atoi(args[5]);
        *kmer_len = atoi(args[6]);
        *shingle_size = atoi(args[7]);
        *max_num_hash = atoi(args[8]);
        *shingling_iterations = atoi(args[9]);
        int min_size_limit = atoi(args[10]);
        *min_cluster_size = min_size_limit > 2 ? min_size_limit : 2;
        *start_hash_num = atoi(args[11]);
        *metis_prefix = args[12];
        *fvalue_filename = args[13];
        *fvalue_threshold = atof(args[14]);
        *hash_difference = atoi(args[15]);
        *grappolo_command = args[16];
        *grappolo_script = args[17];
        *fixed_randoms = (atoi(args[18]) == 0) ? false : true;
	*run_grappolo = (atoi(args[19]) == 0) ? false : true;
}

void parse_input_args(char** args, int* nmap, int* pagesize, int* sepdelta, std::string* input_fname, std::string* output_prefix, 
                        int* seq_count, int* kmer_len, int* shingle_size, int* max_num_hash, int* shingling_iterations, 
			std::string* nmi_output_filename, int* min_cluster_size, bool* rename,
			bool* collect_in_root, int* print_match_stats_cluster_size_limit, std::string* connected_component_algorithm,
			double* sim_score_threshold, double* sim_length_threshold, std::string* similarity_function, 
			std::string* finalized_directory, std::string* community_detection_method, std::string* setname, 
			int* start_hash_num, std::string* metis_prefix, std::string* fvalue_filename, double* fvalue_threshold,
			int* hash_difference, std::string* grappolo_command, std::string* grappolo_script, bool* fixed_randoms) {
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
	*nmi_output_filename = std::string(args[11]);
	int min_size_limit = atoi(args[12]);
	*min_cluster_size = min_size_limit > 2 ? min_size_limit : 2;
	*rename = atoi(args[13]);
	*collect_in_root = atoi(args[14]);
	*print_match_stats_cluster_size_limit = atoi(args[15]);
	*connected_component_algorithm = args[16];
	*sim_score_threshold = atof(args[17]);
	*sim_length_threshold = atof(args[18]);
	*similarity_function = args[19];
	*finalized_directory = args[20];
	*community_detection_method = args[21];
	*setname = args[22];
	*start_hash_num = atoi(args[23]);
	*metis_prefix = args[24];
	*fvalue_filename = args[25];
	*fvalue_threshold = atof(args[26]);
	*hash_difference = atoi(args[27]);
	*grappolo_command = args[28];
	*grappolo_script = args[29];
	*fixed_randoms = atoi(args[30]) == 0 ? false : true;
	*run_grappolo = atoi(args[31]) == 0 ? false : true
}

std::unordered_map<std::string, std::string> read_all_proteins(std::string filename) {
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
