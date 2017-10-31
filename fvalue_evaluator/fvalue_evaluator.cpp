#include "mpi.h"
#include "mapreduce.h"
#include "keyvalue.h"
#include <string.h>
#include "iostream"
//#include "metrics/cluster_evaluator.h"
#include "../src/fvalue_evaluator.h"
#include "../src/mr_util.h"
using namespace std;
using namespace MAPREDUCE_NS;

#define LOG_SCALE 2
#define LOG_PROC -1
#define ROOT_PROCESSOR 0

#define PFAM_COMPARISON true

void parse_input_arguments(char** args, int* nmap, int* pagesize, std::string* shingling_fname, std::string* pfam_fname, int* seq_count, 
				int* min_cluster_limit,	int* print_match_stats_cluster_size_limit, int* largest_cluster_count, int* hashes, 
				std::string* base_type, std::string* fvalue_filename, std::string* setname, bool* rename);

/* ---------------------------------------------------------------------- */

int main(int narg, char **args)
{
	std::cout << "inside main function" << std::endl;
	if (narg != 14) {
		printf("ERROR: Usage profiler.exe (number of maps) (page size) (shingling clusters fname) (pfam fname) (number of sequences) (min cluster size limit) (cluster size limit for printing matching cluster stats) (top largest shingling clusters count for comparison with pfam) (hashnum) (base type) (fvalue output) (setname) (rename)");
		return(0);
	}
	int nmap, pagesize, seq_count, min_cluster_limit, print_match_stats_cluster_size_limit, largest_cluster_count, hashes;
	std::string shingling_fname, pfam_fname, fvalue_output_fname, setname;
	std::string base_type;
	bool rename;
	parse_input_arguments(args, &nmap, &pagesize, &shingling_fname, &pfam_fname, &seq_count, &min_cluster_limit, 
				&print_match_stats_cluster_size_limit, &largest_cluster_count, &hashes, &base_type, 
				&fvalue_output_fname, &setname, &rename);
	//std::cout << "hello" << std::endl;
	//std::cout << "largest cluster count initial: " << largest_cluster_count << std::endl;
	int me;
	MPI_Init(&narg, &args);
	MPI_Comm_rank(MPI_COMM_WORLD, &me);
	
	MapReduce *shingling_mr = new MapReduce(MPI_COMM_WORLD);
	MrToolbox::set_mr_params(shingling_mr, pagesize);
	MapReduce *pfam_mr = new MapReduce(MPI_COMM_WORLD);
	MrToolbox::set_mr_params(pfam_mr, pagesize);

	//READING THE PFAM CLUSTERING
	int pfam_largest_cluster_size, shingling_largest_cluster_size;
	
	ClusterEvaluator clust_eval(seq_count);
	//if (me == 0) std::cout << "reading shingling clusters" << std::endl;
	clust_eval.read_shingling_clusters(shingling_fname, &shingling_largest_cluster_size, rename, min_cluster_limit,
                                                shingling_mr, nmap, largest_cluster_count, "as_test");
        //clust_eval.read_pclust_clusters(shingling_fname, &shingling_largest_cluster_size, min_cluster_limit, shingling_mr, nmap, rename);
	//if (me == 0) std::cout << "reading pfam clusters" << std::endl;
	if (base_type == "shingling") {
	//	std::cout << "reading shingling base clusters from: " << pfam_fname << std::endl;
		clust_eval.read_shingling_clusters(pfam_fname, &pfam_largest_cluster_size, rename, min_cluster_limit,
							pfam_mr, nmap, largest_cluster_count, "as_base");
	} else if (base_type == "pfam") {
	//	std::cout << "reading pfam base clusters from: " << pfam_fname << std::endl;
		clust_eval.read_pfam_clusters(pfam_fname, &pfam_largest_cluster_size, min_cluster_limit, pfam_mr, nmap, rename);
	} else if (base_type == "pclust") {
	//	std::cout << "reading pclust base clusters from: " << pfam_fname << std::endl;
		clust_eval.read_pclust_clusters(pfam_fname, &pfam_largest_cluster_size, min_cluster_limit, pfam_mr, nmap, rename);
	}
	double nmi_lfk, nmi_max;
	//if (me == 0) std::cout << "computing nmi for largest base: " << pfam_largest_cluster_size << " largest_test" << shingling_largest_cluster_size << std::endl;
	
	//pfam_mr -> print("/home/armen.abnousi/shingle/datasets/large70/large70_pfammed_mrmpiformated", 0, -1, 1, 5, 5);
	clust_eval.compute_nmi(shingling_mr, pfam_mr, ROOT_PROCESSOR, &nmi_lfk, &nmi_max);

	if (me == 0) {
	//	std::cout << "set\thash\tbase\tbase_largest\tsh_largest\tsh_f1\tavg_f1\n";
	//	std::cout << setname << "\t" << hashes << "\t" << base_type << "\t" << pfam_largest_cluster_size << "\t" 
	//			<< shingling_largest_cluster_size << "\t" << nmi_lfk << "\t" << nmi_max << std::endl;
	//	std::ofstream outfile;
	//	outfile.open(fvalue_output_fname, std::ios_base::app);
	//	outfile << base_type << "\t" << setname << "\t" << hashes << "\t" << pfam_largest_cluster_size << "\t"
	//		<< shingling_largest_cluster_size << "\t" << nmi_lfk << "\t" << nmi_max << "\n";
	//	outfile.close();
		std::ofstream outfile;
		outfile.open(fvalue_output_fname);
		outfile << nmi_max;
		outfile.close();
	}
		
	delete pfam_mr;
	delete shingling_mr;
	MPI_Finalize();

	return(0);
}

void parse_input_arguments(char** args, int* nmap, int* pagesize, std::string* shingling_fname, std::string* pfam_fname, 
				int* seq_count, int* min_cluster_limit, int* print_match_stats_cluster_size_limit, 
				int* largest_cluster_count, int* hashes, std::string* base_type, std::string* fvalue_output, 
				std::string* setname, bool* rename) {
	*nmap = atoi(args[1]);
	*pagesize = atoi(args[2]);
	*shingling_fname = args[3];
	*pfam_fname = (args[4]);
	*seq_count = atoi(args[5]);
	int min_size_limit = atoi(args[6]);
	*min_cluster_limit = min_size_limit > 2 ? min_size_limit : 2;
	*print_match_stats_cluster_size_limit = atoi(args[7]);
	*largest_cluster_count = atoi(args[8]);
	*hashes = atoi(args[9]);
	*base_type = args[10];
	//std::cout << "arg10 is: " << args[10] << " which makes base_type:" << base_type << std::endl;
	*fvalue_output = args[11];
	*setname = args[12];
	if (atoi(args[13]) == 0) {
		*rename = false;
	} else {
		*rename = true;
	}
}
