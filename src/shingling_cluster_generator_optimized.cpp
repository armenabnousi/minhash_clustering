#include "shingling_cluster_generator_optimized.h"

ShinglingClusterGen::ShinglingClusterGen(int shingling_iters, int cc_iters, 
						std::string output_prefix, std::string cc_alg, int proc_me, int log_scale, int log_proc,
						int prime_num, int num_hash, int shingle_size, int kmer_len, char* random_filename, 
						int max_hash, bool use_persistent_hashing) :
						shingling_iterations(shingling_iters),
                                                connected_component_iterations(cc_iters), connected_component_algorithm(cc_alg),
						shingling_params(ShinglingToolbox(prime_num, num_hash, 
											shingle_size, kmer_len, random_filename, 
											max_hash, use_persistent_hashing)) {
	//std::cout << "first constructor called" << std::endl;
	
	me = me;
	log_scale = log_scale;
	log_proc = log_proc;
	output_filename = set_output_filename(&shingling_params, output_prefix, shingling_iters, cc_iters);
	//std::cout << " now the random: " << std::endl;
	//std::cout << shingling_params.hash_randoms[0] << " " << shingling_params.hash_randoms[1] << std::endl;
	//shingling_params = ShinglingToolbox(prime_num, num_hash, shingle_size, kmer_len, random_filename, max_hash, use_persistent_hashing);
	//shingling_params -> persistent_hash_function_min_values = &(ShinglingClusterGen::persistent_hash_function_min_values);
	//shingling_params -> persistent_vertex_hashes = &(ShinglingClusterGen::persistent_vertex_hashes);
	//shingling_params -> persistent_graph = &(persistent_graph);
	//std::cout << "current sizes: kmertable:" << persistent_hash_function_min_values.size() << " " << persistent_vertex_hashes.size() << std::endl;
	//if (persistent_hash_function_min_values.size() > 0) std::cout << "and cols for kmers: " << persistent_hash_function_min_values.begin() -> second -> size() << " " << persistent_hash_function_min_values.begin() -> second -> at(0).size();
//	ShinglingClusterGen(shingling_iters, cc_iters, output_prefix, cc_alg, prime_num, num_hash, shingle_size, kmer_len, 
//				random_filename, max_hash, use_persistent_hashing);
}

ShinglingClusterGen::ShinglingClusterGen(int shingling_iters, int cc_iters, 
						std::string output_prefix, std::string cc_alg, int prime_num, int num_hash, 
						int shingle_size, int kmer_len, char* random_filename, int max_hash, 
						bool use_persistent_hashing) : shingling_iterations(shingling_iters),
						connected_component_iterations(cc_iters), connected_component_algorithm(cc_alg),
						shingling_params(ShinglingToolbox(prime_num, num_hash,
                                                                                        shingle_size, kmer_len, random_filename,
                                                                                        max_hash, use_persistent_hashing)) {
	//std::cout << "second constructor called" << std::endl;
	//ShinglingToolbox shingling_params(prime_num, num_hash, shingle_size, kmer_len, random_filename, max_hash, use_persistent_hashing);
	output_filename = set_output_filename(&shingling_params, output_prefix, shingling_iters, cc_iters);
	//std::cout << "construction completed" << std::endl;
}

ShinglingClusterGen::~ShinglingClusterGen() { }

void ShinglingClusterGen::update_hashnum(int hashnum) {
	shingling_params.num_hash = hashnum;
}

/*
void ShinglingClusterGen::clean_up_static_members() {
	int me; MPI_Comm_rank(MPI_COMM_WORLD, &me);
	for (auto iter = persistent_hash_function_min_values.begin(); iter != persistent_hash_function_min_values.end(); iter++) {
		//std::cerr << me << " size is: " << persistent_hash_function_min_values -> size() << std::endl;
		delete iter -> second;
		//std::cerr << me << " deleted value, size is: " << persistent_hash_function_min_values -> size() << std::endl;
		//persistent_hash_function_min_values -> erase(iter);
		//std::cerr << me << " deleted key, size is: " << persistent_hash_function_min_values -> size() << std::endl;
	}
	persistent_hash_function_min_values.clear();

	//delete persistent_hash_function_min_values;
}
*/

void ShinglingClusterGen::cluster(MapReduce* mr, std::string input_fname, int nmap, int sepdelta, 
					std::vector< std::set<std::string> >* clusters, int* largest_cluster_size, bool rename, 
					bool collect_in_root, MpiGatherParameters gather_params, int min_cluster_limit, 
					std::unordered_map<std::string, std::string>* local_sequences, 
					std::unordered_set<std::string> clusters_to_remove,
					double* kmer_shingling_time, double* vertex_shingling_time, double* connected_comp_time,
					double time_breakdowns[5], std::string finalized_seq_list_filename) {
	//min_cluster_limit works only if rename==true
	//std::cout << "clustering called" << std::endl;
	std::vector< std::set<std::string> > temp_clusters;
/*	if (file_exists(output_filename)) {
		//if (log_scale > 0 && (log_proc == me || log_proc == -1)) std::cerr << "reading file" << std::endl;
		std::vector<std::string> cluster_labels;
		temp_clusters = read_shingling_clusters(output_filename, largest_cluster_size, rename, min_cluster_limit, &cluster_labels,
							clusters_to_remove);
		int me;
		MPI_Comm_rank(MPI_COMM_WORLD, &me);
		if (me == 0) std::cerr << "largest_cluster_size_is: " << *largest_cluster_size << std::endl;
		if (!collect_in_root) {
			//if (log_scale > 1 && (log_proc == me || log_proc == -1)) std::cerr << "generating mr for " << 
			//								temp_clusters.size() << " clusters" << std::endl;
			int nprocs;
			MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
			Clusters clusters(&temp_clusters, &cluster_labels);
			bool single_proc_owns_container = false;
			DistributeUtil distribute_util(clusters, single_proc_owns_container);
			//if (me == 0) std::cerr << "distributing clusters" << std::endl;
			mr -> map(nprocs, distribute_clusters, &distribute_util);
		//	exit(0);
			mr -> convert();
			//if (me == 0) std::cerr << "clusters distributed" << std::endl;
			//mr -> print("clusterscheck", 0, -1, 1, 5, 5);
		}
	} else {
*/		//if (shingling_params -> use_persistent_hashing) {
		//	shingling_params -> latest_hashed = latest_hashed;
		//}
		mrmpi_cluster(mr, input_fname, nmap, sepdelta, local_sequences, gather_params, collect_in_root, &temp_clusters, 
				kmer_shingling_time, vertex_shingling_time, connected_comp_time, time_breakdowns, 
				finalized_seq_list_filename);
		//if (shingling_params -> use_persistent_hashing) {
		//	latest_hashed = std::max(shingling_params -> num_hash, latest_hashed);
		//}
		//mr -> print("checkingccs", 0, -1, 1, 5, 5);
		//exit(0);
		if (rename) {
			if (log_scale > 1 && (log_proc == me || log_proc == -1)) std::cerr << "renaming sequences" << std::endl;
			mr -> reduce(rename_sequences, NULL);
			mr -> collate(NULL);
			bool is_string = true;
			MrToolbox::remove_redundant_local(mr, is_string);
			//mr -> convert();
			//mr -> print("renamed_clusts", 0, -1, 1, 5, 5);
			ClusterLimits limits(min_cluster_limit);
			*largest_cluster_size = 0;
			mr -> reduce(remove_small_clusters, &limits);
			*largest_cluster_size = limits.largest_cluster_size;
		//std::cerr << "hereitsallgeneratedkv " << collect_in_root << std::endl;
			int clusters_count = mr -> convert();
		 	int local_largest_cluster_sizes[gather_params.num_processors];
			MPI_Gather((char*) largest_cluster_size , 1, MPI_INT, local_largest_cluster_sizes, 1, MPI_INT, 
						gather_params.root, MPI_COMM_WORLD);
			int me;
			MPI_Comm_rank(MPI_COMM_WORLD, &me);
			if (me == gather_params.root) {
				for (int i = 0; i < gather_params.num_processors; i++) {
					*largest_cluster_size = *largest_cluster_size < local_largest_cluster_sizes[i] ? 
									local_largest_cluster_sizes[i] : *largest_cluster_size;
				}
				//std::cerr << "largestclusterfoundinthisround: " << *largest_cluster_size << std::endl;
			}
			//mr -> print("mrclusterscheck", 0, -1, 1, 5, 5);
			//std::cerr << "clusterscountathtispoint: " << clusters_count << std::endl;
		}
		if (collect_in_root && connected_component_algorithm != "union_find" ) {
			if (log_scale > 1 && (log_proc == me || log_proc == -1)) std::cerr << "collecting" << std::endl;
			std::vector<std::string> cluster_labels;
			temp_clusters = collect_clusters(mr, gather_params, largest_cluster_size, &cluster_labels);
		}
//	}
	clusters -> assign (temp_clusters.begin(), temp_clusters.end());
	shingling_params.latest_hashed = shingling_params.num_hash;
}

void ShinglingClusterGen::mrmpi_cluster(MapReduce* mr, std::string input_fname, int nmap, int sepdelta, 
					std::unordered_map<std::string, std::string>* local_sequences, 
					MpiGatherParameters gather_params, bool collect_in_root, 
					std::vector< std::set<std::string> >* final_clusters,
					double* kmer_shingling_time, double* vertex_shingling_time, double* connected_comp_time,
					double kmer_shingling_time_breakdown[5], std::string finalized_seq_list_filename) {
	//std::cout << "inside mrmpi" << std::endl;
	int me; MPI_Comm_rank(MPI_COMM_WORLD, &me);
	//shingling_params -> persistent_hash_function_min_values = persistent_hash_function_min_values;
	bool multishingle = shingling_iterations > 0 ? true : false;
	shingling_params.multishingle = multishingle;
	//if (me == 0) std::cerr << "hashrandoms: me:" << me << " hash_num:" << shingling_params -> num_hash << " h1:" << shingling_params -> hash_randoms[0] << " " << shingling_params -> hash_randoms[1] << std::endl;
	//if (me == 43) std::cerr << "hashrandoms: me:" << me << " hash_num:" << shingling_params -> num_hash << " h1:" << shingling_params -> hash_randoms[0] << " " << shingling_params -> hash_randoms[1] << std::endl;
	double tstart = MPI_Wtime();
	char** dirname = new char*[1];
        dirname[0] = new char[input_fname.length()+1];
        strcpy(&dirname[0][0], input_fname.c_str());
	int nproteins;
	//std::cerr << "entering for first time" << std::endl;
	if (shingling_params.latest_hashed == 0 || !(shingling_params.use_persistent_hashing)) { 
		//if (me == 0) std::cerr << "reading proteins anew" << std::endl;
		nproteins = mr -> map(nmap, 1, dirname, 0, 0, "\n>", sepdelta, &ShinglingClusterGen::read_proteins, local_sequences);
		//if (me == 0) std::cout << "md1: proteins read" << std::endl;
		//if (me < 5) std::cerr << "me : " << me << " at the beginning: " << local_sequences -> size() << std::endl;
		local_sequences -> clear();
		mr -> aggregate(NULL);
		//if (me == 0) std::cout << "md2 aggregated" << std::endl;
		mr -> scan(establish_sequence_ownership, local_sequences);
		//if (me == 0) std::cout << "md3 ownership established" << std::endl;
		if (finalized_seq_list_filename != "") {
			remove_finalized_sequences(local_sequences, finalized_seq_list_filename);
		}
		
		//populate_seq_owners_proc(local_sequences);
		//if (!(shingling_params -> use_persistent_hashing)) latest_hashed = 1;
	}
	int local_proteins_count = local_sequences -> size();
	int global_proteins_count;
	if (me == 0) std::cout << "starting regular computation..." << std::endl;
	MPI_Reduce(&local_proteins_count, &global_proteins_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	//if (me == 0) std::cerr << "last hashed: " << shingling_params.latest_hashed << std::endl;
	//if (me == 0) std::cerr << "proteins in the dataset: " << global_proteins_count << std::endl;
	//std::cerr << "firstifpassed" << std::endl;
	//if (latest_hashed < 9) {
	//	if (me < 5) std::cerr << "me : " << me << " next round:" << local_sequences -> size() << std::endl;
	//}
	//std::cerr << "map is done" << std::endl;
	//if (me == 0) std::cerr << "check dist: sequences mapped" << std::endl;
	//MPI_Barrier(MPI_COMM_WORLD);
	//mr -> kv_stats(2);
	delete[] dirname[0];
        delete[] dirname;
	if (log_scale > 1 && (log_proc == me || log_proc == -1)) std::cerr << "num proteins: " << nproteins << std::endl;
	//changed from aggregate(NULL)
	//if (me == 0) std::cerr << "check dist: aggregated (shouldn't change much)" << std::endl;
        //MPI_Barrier(MPI_COMM_WORLD);
        //mr -> kv_stats(2);
        shingling_params.local_proteins = local_sequences;
	double tstartbd = MPI_Wtime();
	//if (me == 0) std::cerr << "map to generate local seq kmer shingles" << std::endl;
	int emitted_shingles = mr -> map(nmap, generate_local_seqs_kmer_shingles, &shingling_params, 0);
//mr -> sort_keys(6);
	if (me == 0) std::cerr << "kmer shingles generated" << std::endl;
//std::string filename = "xv_kmer_doubled_" + std::to_string(shingling_params.num_hash);
//mr -> print(&filename[0], 0, -1, 1, 0, 5);
	double tendbd = MPI_Wtime();
	kmer_shingling_time_breakdown[0] = tendbd - tstartbd;
	//latest_hashed = shingling_params -> num_hash;
	//shingling_params -> latest_hashed = latest_hashed;
	tstartbd = MPI_Wtime();
	//if (me == 0) std::cerr << "check dist: min shingles generated - second map" << std::endl;
        //MPI_Barrier(MPI_COMM_WORLD);
        //mr -> kv_stats(2);
        //double tend = MPI_Wtime();
	if (log_scale > 1 && (log_proc == me || log_proc == -1)) std::cerr << "emitted shingles: " << emitted_shingles << std::endl;
	mr -> collate(NULL);
	if (log_scale > 1 && (log_proc == me || log_proc == -1)) std::cerr << "drawing initial graph " << std::endl;
	tendbd = MPI_Wtime();
	kmer_shingling_time_breakdown[1] = tendbd - tstartbd;
	tstartbd = MPI_Wtime();
	//double tend = MPI_Wtime();
	//if (me == 0) std::cerr << "check dist: collated by min shingles" << std::endl;
	//if (me == 0) std::cerr << " writing down the nvalues: " << std::endl;
	//MPI_Barrier(MPI_COMM_WORLD);
	//mr -> scan(print_details, NULL);
	//if (shingling_params -> num_hash == 5) mr -> print("nvals", 0, -1, 1, 5, 0);
	//MPI_Barrier(MPI_COMM_WORLD);
	//mr -> kmv_stats(2);
	if (me == 0) std::cerr << "drawing graph now... " << std::endl;
	int nedges = mr -> reduce(draw_graph, NULL);
	if (me == 0) std::cerr << "graph drawn" << std::endl;
	
//mr -> sort_keys(5);
//filename = "xv_reduce_draw_doubled_" + std::to_string(shingling_params.num_hash);
//mr -> print(&filename[0], 0, -1, 1, 5, 5);
	tendbd = MPI_Wtime();
	kmer_shingling_time_breakdown[2] = tendbd - tstartbd;
	tstartbd = MPI_Wtime();
	//int ncount = mr -> collate(send_to_owner_proc);
	int ncount = mr -> collate(NULL);
	//if (me == 0) std::cout << "initial node count: " << ncount << std::endl;
	tendbd = MPI_Wtime();
	kmer_shingling_time_breakdown[3] = tendbd - tstartbd;
	tstartbd = MPI_Wtime();
	int nnodes = 0;
	nedges = 0;
//mr -> sort_multivalues(5);
//filename = "xv_before_rrs_doubled_" + std::to_string(shingling_params.num_hash);
//mr -> print(&filename[0], 0, -1, 1, 5, 5);
	MrToolbox::remove_redundant_string_manual(mr);
//filename = "xv_after_rrs_doubled_" + std::to_string(shingling_params.num_hash);
//mr -> print(&filename[0], 0, -1, 1, 5, 5);
	/*
	mr** -> reduce(redundant_agglomerator, NULL);
	mr -> convert();
	nedges = mr -> reduce(redundant_remover, NULL);
	int nnodes = mr -> convert();
	*/
	tendbd = MPI_Wtime();
	kmer_shingling_time_breakdown[4] = tendbd - tstartbd;
	double tend = MPI_Wtime();
	if (kmer_shingling_time) {
	        *kmer_shingling_time = tend - tstart;
	}
	if (me == 0) std::cerr << "entering the iterations loop" << std::endl;
	tstart = MPI_Wtime();
	//mr -> print("interfiles_initgraph", 0, -1, 1, 5, 5);
	if (log_scale > 2 && (log_proc == me || log_proc == -1)) std::cerr << "redundant edges removed" << std::endl;
	if (log_scale > 1 && (log_proc == me || log_proc == -1)) std::cerr << "edges: " << nedges << "\t" << "nodes: " 
										<< nnodes << std::endl;
	if (log_scale > 1 && (log_proc == me || log_proc == -1)) std::cerr << "multishingle set to: " << multishingle << 
							" with sh iterations asked for: " << shingling_iterations << std::endl;
	for (int iteration = 1; iteration < shingling_iterations; iteration++) {
		multishingle = true;
		if (log_scale > 1 && (log_proc == me || log_proc == -1)) std::cerr << "running shingling iteration: " 
										<< iteration << std::endl;
		if (iteration == shingling_iterations - 1) {
			shingling_params.last_iter = true;
		}
//std::string filename = "/home/armen.abnousi/shingle/kmer_graph" + std::to_string(shingling_params -> num_hash) + "_optimized";
//mr -> print(&filename[0], 0, -1, 1, 5, 5);
		//if (me == 0) std::cerr << "generating vertex shingles" << std::endl;
		mr -> reduce(generate_vertex_neighbor_min_shingles, &shingling_params);
		if (me == 0) std::cerr << "vertex shingles generated" << std::endl;
		if (log_scale > 1 && (log_proc == me || log_proc == -1)) std::cerr << "shingles generated" << std::endl;
		mr -> collate(NULL);
		if (log_scale > 1 && (log_proc == me || log_proc == -1)) std::cerr << "collate completed" << std::endl;
		int edges = mr -> reduce(draw_graph, NULL);
		if (log_scale > 1 && (log_proc == me || log_proc == -1)) std::cout << "edges: " << edges << std::endl;
		//if (me == 0) std::cerr << "removong redundant..." << std::endl;
		bool is_string = true;
		mr -> convert();
		MrToolbox::remove_redundant_local(mr, is_string);
		MrToolbox::open_kmv(mr);
		//MrToolbox::remove_redundant(mr);
		if (log_scale > 1 && (log_proc == me || log_proc == -1)) std::cerr << "graph drawn" << std::endl;
		//mr -> collate(send_to_owner_proc);
		if (me == 0) std::cerr << "final collate " << std::endl;
		int ncount = mr -> collate(NULL);
		if (log_scale > 1 && (log_proc == me || log_proc == -1)) std::cerr << "collate completed" << std::endl;
		//mr -> print("interfiles_m2graph", 0, -1, 1, 5, 5);
	}
	if (me == 0) std::cerr << "finalizing stuff" << std::endl;
	//latest_hashed = shingling_params -> num_hash;
	//shingling_params -> latest_hashed = latest_hashed;
	shingling_params.last_iter = false;
	if (multishingle) {
		if (log_scale > 1 && (log_proc == me || log_proc == -1)) std::cout << "removing redundant edges" << std::endl;
		bool is_string = true;
		MrToolbox::remove_redundant_local(mr, is_string);
		//mr -> convert();
		//NODES AND EDGES ARE NOT COMPUTED ANYMORE!!!!
		if (log_scale > 1 && (log_proc == me || log_proc == -1)) std::cout << "done" << std::endl;
		/*
  		mr -> reduce(redundant_agglomerator, NULL);
		if (log_scale > 1 && (log_proc == me || log_proc == -1)) std::cout << "reduce done" << std::endl;
		mr -> convert();
		if (log_scale > 1 && (log_proc == me || log_proc == -1)) std::cout << "converted" << std::endl;
		nedges = mr -> reduce(redundant_remover, NULL);
		if (log_scale > 1 && (log_proc == me || log_proc == -1)) std::cout << "reduce2 donr" << std::endl;
		nnodes = mr -> convert();
		if (log_scale > 1 && (log_proc == me || log_proc == -1)) std::cout << "converted2" << std::endl;
		*/
		//mr -> print("interfiles_redun_removed_graph", 0, -1, 1, 5, 5);
	}
	tend = MPI_Wtime();
        if (vertex_shingling_time) {
		*vertex_shingling_time = tend - tstart;
	}
	//mr -> print("/home/armen.abnousi/shingle/nocc", 0, -1, 1, 5, 5);
	//exit(0);
	tstart = MPI_Wtime();
	if (connected_component_algorithm == "label_propagation") {
		//NODES AND EDGES ARE NOT COMPUTED ANYMORE!!!!!
		perform_label_propagate(mr, nnodes, nedges, connected_component_iterations, connected_comp_time);
	} else if (connected_component_algorithm == "union_find") {
		perform_union_find(mr, final_clusters, collect_in_root, gather_params);
	} else if (connected_component_algorithm == "hash_to_min") {
		//NODES AND EDGES ARE NOT COMPUTED ANYMORE!!!!!
		perform_hash_to_min_connected_component(mr, nnodes, nedges, connected_component_iterations, connected_comp_time);
	} else if (connected_component_algorithm == "nocc") {
		if (log_scale > 1 && (log_proc == me || log_proc == -1)) std::cout << "doing final touches... adding label to cluster" << std::endl;
		if (me == 0) std::cerr << "in nocc" << std::endl;
		mr -> reduce(add_label_to_cluster, NULL);
		mr -> convert();
		mr -> reduce(add_special_characters, NULL);
//mr -> sort_keys(5);
		int final_nodes_count = mr -> convert();
//mr -> sort_multivalues(5);
		//if(me == 0) std::cout << "final nodes count: " << final_nodes_count << std::endl;
	}
	/*
	mr -> reduce(redundant_agglomerator, NULL);
	mr -> convert();
	mr -> reduce(redundant_remover, NULL);
	mr -> convert();
	*/
	bool is_string = true;
	MrToolbox::remove_redundant_local(mr, is_string);
	//mr -> convert();
	tend = MPI_Wtime();
	if (connected_comp_time) {
		*connected_comp_time = tend - tstart;
	}
	//mr -> print("interfiles_finalized_redun_removed", 0, -1, 1, 5, 5);
	if (output_filename != "") {
		char output_name[output_filename.length() + 1];
		strcpy(output_name, output_filename.c_str());
		remove(output_name);
		//ifstream file_check(output_name);
		//if (file_check) {
		//	file_check.close();
		//	if (log_scale > 1 && (log_proc == me || log_proc == -1)) std::cerr << " output printed" << std::endl;
		//	remove(output_name);
		//}
//		mr -> print(output_name, 0, -1, 1, 5, 5);
		if (log_scale > 1 && (log_proc == me || log_proc == -1)) std::cerr << " output printed" << std::endl;
		//exit(0);
	}
	if (log_scale > 1 && (log_proc == me || log_proc == -1)) std::cerr << "clusters computed." << std::endl;
}

std::string ShinglingClusterGen::set_output_filename(ShinglingToolbox* params, std::string prefix, int shingling_iters, int cc_iters) {
	std::string filename = prefix + "_k" + std::to_string((long long int) params -> kmer_len) + 
					"sh" + std::to_string((long long int) params -> shingle_size) + 
					"m" + std::to_string((long long int) shingling_iters) + 
					"cc" + std::to_string((long long int) cc_iters) + 
					"h" + std::to_string((long long int) params -> num_hash);
	return filename;
}

void ShinglingClusterGen::read_proteins(int itask, char* text, int size, KeyValue* kv, void *ptr) {
	FastaToolbox fasta_toolbox;
	//std::cerr << "reading proteins" << std::endl;
	std::unordered_map<std::string, std::string>* local_seqs = (std::unordered_map<std::string, std::string>*) ptr;
//	if (local_seqs -> size() == 0) {
		int me; //used in setting the seq_owner_proc
		MPI_Comm_rank(MPI_COMM_WORLD, &me);
		if (log_scale > 1 && (log_proc == me || log_proc == -1)) std::cerr << "proc: " << me << " reading proteins" << std::endl;
		if (log_scale > 3 && (log_proc == me || log_proc == -1)) std::cerr << "proc: " << me <<	
										"received text:\n" << text << std::endl;
		std::list<Protein> proteins = fasta_toolbox.read_fasta(text, size, itask);
		if (log_scale > 1 && (log_proc == me || log_proc == -1)) {
			auto protein = proteins.begin();
			std::cerr << protein -> name << std::endl;
			protein = proteins.end();
			protein--;
			std::cerr << protein -> name << std::endl;
			std::cerr << itask << "readseqs: " << proteins.size() << std::endl;
		}
		for (auto protein = proteins.begin(); protein != proteins.end(); protein++) {
			int name_length = (protein -> name).length() + 1;
			int sequence_length = (protein -> sequence).length() + 1;
			char sequence_name[name_length];
			char sequence[sequence_length];
			strcpy(sequence_name, (protein -> name).c_str());
			strcpy(sequence, (protein -> sequence).c_str());
			kv -> add(sequence_name, name_length, sequence, sequence_length);
			local_seqs -> insert(std::make_pair(protein -> name, protein -> sequence));
		}
/*		populate_seq_owners_proc(local_seqs);
		//generate_min_shingles(proteins, kmer_len, kv, hash_params);
	} else {
		int me; MPI_Comm_rank(MPI_COMM_WORLD, &me);
		//if (log_scale > -1 && (log_proc == me || log_proc == -1)) std::cerr << " hashing local sequences" << std::endl;
		//std::cerr << "mapping: " << local_seqs -> size() << "seqs by: " << me << std::endl;
		for (auto protein = local_seqs -> begin(); protein != local_seqs -> end(); protein++) {
			int name_length = (protein -> first).length() + 1;
			int sequence_length = (protein -> second).length() + 1;
			char sequence_name[name_length];
			char sequence[sequence_length];
			strcpy(sequence_name, (protein -> first).c_str());
			strcpy(sequence, (protein -> second).c_str());
		//	std::cerr << "adding: " << sequence_name << " len:" << name_length << " seq:" << sequence_length << std::endl;
			kv -> add(sequence_name, name_length, sequence, sequence_length);
		}
	}
*/
}

void ShinglingClusterGen::generate_local_seqs_kmer_shingles(int itask, KeyValue* kv, void* shingling_params) {
	ShinglingToolbox* shingling_parameters = (ShinglingToolbox*) shingling_params;
	bool use_persistent_hashing = shingling_parameters -> use_persistent_hashing;
	int me; MPI_Comm_rank(MPI_COMM_WORLD, &me);
	//if (shingling_parameters -> num_hash < 12) {
	//	std::cerr << "itask: " << itask << " me: " << me << " nhash: " << shingling_parameters -> num_hash << " seqs: " << shingling_parameters -> local_proteins -> size() << std::endl;
	//}
	for (auto iter = shingling_parameters -> local_proteins -> begin(); iter != shingling_parameters -> local_proteins -> end(); iter++) {
		generate_kmer_min_shingles(kv, iter -> first, iter -> second, shingling_parameters, itask, use_persistent_hashing);
	}
}
void ShinglingClusterGen::generate_kmer_min_shingles(KeyValue* kv, std::string seq_name, std::string sequence, 
							ShinglingToolbox* shingling_parameters, int itask, bool use_persistent_hashing) {
	int me; MPI_Comm_rank(MPI_COMM_WORLD, &me);

//	ShinglingToolbox* shingling_parameters = (ShinglingToolbox*) shingling_params;
	int name_length = seq_name.length() + 1;
	char name[name_length];
	strcpy(name, seq_name.c_str());
	std::vector<std::vector<unsigned int>>* hash_function_min_values;
	int new_rows_count = (shingling_parameters -> num_hash) - (shingling_parameters -> latest_hashed);
	int total_rows_count = shingling_parameters -> num_hash;
	int ncols = shingling_parameters -> shingle_size;
	//if (itask == 0) std::cerr << "generateing " << new_rows_count << " new rows" << std::endl;
	//std::cout << "here: " << new_rows_count << " " << total_rows_count << " " << ncols << " " << shingling_parameters -> latest_hashed << " " << shingling_parameters -> num_hash << std::endl;
//	std::cerr << "for seq: " << name << " initializing a local table of nrows: " << new_rows_count << std::endl;
	if (new_rows_count > 0) {
		std::vector<unsigned int> empty_vector;
		hash_function_min_values = new std::vector<std::vector<unsigned int>>;
		hash_function_min_values -> insert(hash_function_min_values -> begin(), new_rows_count, empty_vector);
//		std::cerr << " initializied: " << hash_function_min_values -> size() << " rows" << std::endl;
//		unsigned int** hash_function_min_values;
		//if (me == 0) std::cerr << "inside hashrandoms: me:" << me << " hash_num:" << shingling_parameters -> num_hash << " h1:" << shingling_parameters -> hash_randoms[0] << " " << shingling_parameters -> hash_randoms[1] << std::endl;
        	//if (me == 43) std::cerr << "inside hashrandoms: me:" << me << " hash_num:" << shingling_parameters -> num_hash << " h1:" << shingling_parameters -> hash_randoms[0] << " " << shingling_parameters -> hash_randoms[1] << std::endl;
	
//		std::vector<std::vector<unsigned int>> hash_function_min_values(shingling_parameters -> num_hash, 
//								std::vector<unsigned int>(shingling_parameters -> shingle_size, UINT_MAX));
		if (log_scale > 3 && (log_proc == itask || log_proc == -1)) std::cerr << "proc " << itask << ": initializing array" << std::endl;
//		initialize_array_by_max_uint(&hash_function_min_values, shingling_parameters -> num_hash, shingling_parameters -> shingle_size);
		if (log_scale > 3 && (log_proc == itask || log_proc == -1)) std::cerr << "proc " << itask << ": array initialized" << std::endl;
		//std::string sequence = (std::string) sequence_array;
		int last_kmer_pos = sequence.length() - (shingling_parameters -> kmer_len);
		std::string kmer;
		std::string previous_kmer = "";
		unsigned long long kmer_int_value = 0;
		unsigned int kmer_hashed_value;
		for (int i = 0; i <= last_kmer_pos; i++, previous_kmer = kmer) {
			kmer = sequence.substr(i, shingling_parameters -> kmer_len);
			kmer_int_value = compute_kmer_int_value(kmer, previous_kmer, kmer_int_value, shingling_parameters -> kmer_len);
			for (int hash = shingling_parameters -> latest_hashed; hash < shingling_parameters -> num_hash; hash++) {
				kmer_hashed_value = compute_hash_value(kmer_int_value, shingling_parameters -> hash_randoms[hash * 2],
							shingling_parameters -> hash_randoms[hash * 2 + 1], shingling_parameters -> Prime_num);
				//if (shingling_parameters -> num_hash == 5 && kmer_hashed_value == 0) 
				//				std::cerr << "hashedto: 0, by: " << kmer << " intval: " << kmer_int_value 
				//					<< " a: " << shingling_parameters -> hash_randoms[hash * 2] << " b: " 
				//					<< shingling_parameters -> hash_randoms[hash * 2 + 1] << " P: " 
				//					<< shingling_parameters -> Prime_num << " hash:" << hash << std::endl;
				if (log_scale > 4 && (log_proc == itask || log_proc == -1)) 
					std::cerr << "proc " << itask << ": hashval computed " << kmer_int_value 
					<< " " <<  kmer_hashed_value << " " << hash << " " << shingling_parameters -> latest_hashed << std::endl;
				//update_hash_func_min_values(kmer_hashed_value, &kmer_hash_function_min_values[hash][0], shingling_parameters -> shingle_size);
				int current_row = hash - shingling_parameters -> latest_hashed;
				update_hash_func_min_values(kmer_hashed_value, hash_function_min_values, current_row, ncols);
				if (log_scale > 4 && (log_proc == itask || log_proc == -1)) 
								std::cerr << "proc " << itask << ": array updated" << std::endl;
			}
		}
		if (log_scale > 3 && (log_proc == -1)) {
			std::cerr << "new table for: " << name << std::endl;
			for (int j = 0; j < new_rows_count; j++) {
        		        for (int i = 0; 
					i < shingling_parameters -> shingle_size && i < hash_function_min_values -> at(j).size(); i++) { 
        		                std::cerr << "\t" << hash_function_min_values -> at(j)[i];
       		        	}       
        	        	std::cerr << std::endl;
        		}
		}
//	if (!(shingling_parameters -> multishingle)) {
		//append hash_function_min_values to the persistent_hash_table
//		std::cerr << "appending new rows to the persistent_hash_table" << std::endl;
		if (use_persistent_hashing) {
			if (shingling_parameters -> persistent_hash_function_min_values.find(seq_name) ==
				shingling_parameters -> persistent_hash_function_min_values.end()) {
				shingling_parameters -> persistent_hash_function_min_values.insert(
                	                      std::make_pair(seq_name, new std::vector<std::vector<unsigned int>>));
			}
//			std::cerr << "the name: " << name << " didn't exist in the map, creating now to append" << std::endl;
		//	shingling_parameters -> persistent_hash_function_min_values -> insert(
		//			std::make_pair(seq_name, hash_function_min_values));
		//	delete hash_function_min_values;
		//} else {
//			std::cerr << "existed: " << name << ". just appending" << std::endl;
			shingling_parameters -> persistent_hash_function_min_values.at(seq_name) -> insert(
					shingling_parameters -> persistent_hash_function_min_values.at(seq_name) -> end(), 
					hash_function_min_values -> begin(), hash_function_min_values -> end());
			delete hash_function_min_values;
		//}
		}
	}
	//emit shingles from persistent_hash_table
//		std::cerr << "emitting from the persistent_table" << std::endl;
/*	if (shingling_parameters -> num_hash == 9) {
                std::cerr << "printing table for " << name << " for num_hash 9, total_rows_count " << total_rows_count 
				<< " sizes:" << persistent_hash_function_min_values -> at(name) -> size() << ", " 
				<< persistent_hash_function_min_values -> at(name) -> at(0).size() << std::endl;
		for (int i = 0; i < total_rows_count; i++) {
                       	for (int j = 0; j < ncols; j++) {
                               	std::cerr << "\t" << persistent_hash_function_min_values -> at(name) -> at(i)[j];
                       	}
	                std::cerr << std::endl;
		}
	}
*/	if (use_persistent_hashing) {
		emit_min_shingle_kvs(shingling_parameters -> persistent_hash_function_min_values.at(name), kv, total_rows_count, ncols, name, name_length);
		//emit_min_shingle_kvs(hash_function_min_values, kv, new_rows_count, ncols, name, name_length);
		//delete hash_function_min_values;
	} else {
		emit_min_shingle_kvs(hash_function_min_values, kv, total_rows_count, ncols, name, name_length);
		delete hash_function_min_values;
	}
/*	} else {
		//emit shingles from hash_function_min_values
//		std::cerr << "emitting from new table and deleting the table" << std::endl;
		if (shingling_parameters -> num_hash == 7) {
			std::cerr << "printing the last row for " << name << " for num_hash 7, new_rows_count " << new_rows_count << std::endl;
			for (int i = 0; i < hash_function_min_values -> at(0).size(); i++) {
				std::cerr << "\t" << hash_function_min_values -> at(0)[i];
			}
			std::cerr << std::endl;
		}
		emit_min_shingle_kvs(hash_function_min_values, kv, new_rows_count, ncols, name, name_length);
	
		delete hash_function_min_values;
	}
*/	if (log_scale > 2 && (log_proc == itask || log_proc == -1)) std::cerr << "proc " << itask << ": all shingles computed" << std::endl;
//	emit_min_shingle_kvs(hash_function_min_values, kv, shingling_parameters -> num_hash, shingling_parameters -> shingle_size, 
//					name, name_length);
	if (log_scale > 3 && (log_proc == itask || log_proc == -1)) std::cerr << "proc " << itask << ": kvs emitted" << std::endl;
//	delete_2d_array(hash_function_min_values, shingling_parameters -> num_hash, shingling_parameters -> shingle_size);
	if (log_scale > 3 && (log_proc == itask || log_proc == -1)) std::cerr << "proc " << itask << ": array destroyed" << std::endl;
}

void ShinglingClusterGen::initialize_array_by_max_uint(unsigned int*** array, int rows, int cols) {
	*array = new unsigned int*[rows];
	for (int i = 0; i < rows; i++) {
		(*array)[i] = new unsigned int[cols];
	}
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			(*array)[i][j] = UINT_MAX;
		}
	}
}

unsigned long long ShinglingClusterGen::compute_kmer_int_value(std::string kmer, std::string previous_kmer,
						unsigned long long previous_kmer_value, int kmer_len) {
	unsigned long long result;
	if (previous_kmer == "") {
		result = 0;
		for (int i = 0; i < kmer_len; i++) {
			result += (kmer[i] - 'A') * pow(ALPHABET_SIZE, kmer_len - 1 - i);
		}
	} else {
		result = previous_kmer_value;
		result -= (previous_kmer[0] - 'A') * pow(ALPHABET_SIZE, kmer_len - 1);
		result *= ALPHABET_SIZE;
		result += kmer[kmer_len - 1] - 'A';
	}
	return result;
}

int ShinglingClusterGen::compute_hash_value(unsigned long long value, int random_num1, int random_num2, int Prime_num) {
	int hash_value = ((value * random_num1) + random_num2) % Prime_num;
	return hash_value;
}

void ShinglingClusterGen::update_hash_func_min_values(unsigned int hashed_value, unsigned int* hash_function_min_values, int shingle_size) {
	for (int small_me = 0; small_me < shingle_size; small_me++) {
		if (hashed_value < hash_function_min_values[small_me]) {
			unsigned int temp = hashed_value;
			hashed_value = hash_function_min_values[small_me];
			hash_function_min_values[small_me] = temp;
		}
	}
}

void ShinglingClusterGen::update_hash_func_min_values(unsigned int hashed_value, std::vector<std::vector<unsigned int>>* hash_table, 
								int row, int max_col) {
	//row is the hash function's number
	//number of columns is the number of shingles
	//std::cerr << "updating by:" << hashed_value << std::endl;
	//std::cout << "hash table nrows: " << hash_table -> size() << std::endl;
	//std::cout << "currently processing row: " << row << std::endl;
	//std::cout << "hash table cols: " << hash_table -> at(row).size() << std::endl;
	//std::cout << "max cols allowed: " << max_col << std::endl;
	auto row_begin_iter = hash_table -> at(row).begin();
	bool row_is_full = false;
	int insert_pos = 0;
	//std::cerr << "printing size " << std::endl;
	//std::cerr << " size is: " <<  hash_table -> at(row).size() << std::endl;
	//std::cerr << "size printed" << std::endl;
	for (insert_pos = 0; insert_pos < hash_table -> at(row).size(); insert_pos++) {
		//std::cerr << "in loop" << std::endl;
		//std::cerr << "comparing with : " << (hash_table -> at(row))[insert_pos] << std::endl;
		if (hashed_value <= (hash_table -> at(row))[insert_pos]) {
			//hash_table -> insert(row_begin_iter + i, hashed_value);
			break;
		}
	}
	//std::cerr << "inserting in pos: " << insert_pos;
	//std::cerr << " current_size:"  << hash_table -> at(row).size() << std::endl;
	if (insert_pos < max_col) {
		if (hash_table -> at(row).size() + 1 > max_col) {
			hash_table -> at(row).pop_back();
		}
		hash_table -> at(row).insert(row_begin_iter + insert_pos, hashed_value);
	}
	//std::cerr << "updated: ";
	//for (insert_pos = 0; insert_pos < hash_table -> at(row).size(); insert_pos++) {
	//	std::cerr << hash_table -> at(row).at(insert_pos) << " ";
	//}
	//std::cerr << std::endl;
}

void ShinglingClusterGen::emit_min_shingle_kvs(std::vector<std::vector<unsigned int>>* hash_table, KeyValue* kv, int nrows, int ncols,
							char* name, int name_length) {
	if (hash_table -> size() < nrows) {
		std::cerr << "for name: " << name << " there will be a new edge. nrows: " << nrows << " " << hash_table -> size() << std::endl;
	}
	for (int i = 0; i < nrows; i++) {
		if (hash_table -> at(i).size() >= ncols) {
			std::ostringstream convert;
			convert << i;
			for (int j = 0; j < ncols; j++) {
				convert << "$" << hash_table -> at(i).at(j);
			}
			std::string shingle = convert.str();
			char concated_shingles[shingle.length() + 1];
			strcpy(concated_shingles, shingle.c_str());
			kv -> add(concated_shingles, shingle.length() + 1, name, name_length);
			if (log_scale > 4 && (log_proc == -1)) std::cerr << "emitted: " << concated_shingles << "\t" << name << std::endl;
		}
	}
}

void ShinglingClusterGen::emit_min_shingle_kvs(unsigned int** hash_function_min_values, KeyValue* kv, int num_hash, int shingle_size,
				char* name, int name_length) {
	bool incomplete = false;
	for (int hash = 0; hash < num_hash; hash++) {
		std::ostringstream convert;
		for (int shingle = 0; !incomplete && shingle < shingle_size; shingle++) {
			convert << "$" << hash_function_min_values[hash][shingle];
			if (hash_function_min_values[hash][shingle] == UINT_MAX) {
				incomplete = true;
			}
		}
		if (!incomplete) {
			std::string shingle = convert.str();
			char concated_shingles[shingle.length() + 1];
			strcpy(concated_shingles, shingle.c_str());
			char sequence_name[name_length];
			strcpy(sequence_name, name);
			//std::cerr << "sending: " << concated_shingles << " and " << sequence_name << std::endl;
			kv -> add(concated_shingles, shingle.length() + 1, sequence_name, name_length);
			if (log_scale > 4 && (log_proc == -1)) std::cerr << "emitting: " << concated_shingles << "\t" << name << std::endl;
		}
	}
}

void ShinglingClusterGen::emit_min_shingle_kvs(std::vector<std::vector<unsigned int>> hash_function_min_values, KeyValue* kv, 
						int num_hash, int shingle_size, char* name, int name_length) {
        bool incomplete = false;
        for (int hash = 0; hash < num_hash; hash++) {
                std::ostringstream convert;
                for (int shingle = 0; !incomplete && shingle < shingle_size; shingle++) {
                        convert << "$" << hash_function_min_values[hash][shingle];
                        if (hash_function_min_values[hash][shingle] == UINT_MAX) {
                                incomplete = true;
                        }
                }
                if (!incomplete) {
                        std::string shingle = convert.str();
                        char concated_shingles[shingle.length() + 1];
                        strcpy(concated_shingles, shingle.c_str());
                        char sequence_name[name_length];
                        strcpy(sequence_name, name);
                        //std::cerr << "sending: " << concated_shingles << " and " << sequence_name << std::endl;
                        kv -> add(concated_shingles, shingle.length() + 1, sequence_name, name_length);
                        if (log_scale > 4 && (log_proc == -1)) std::cerr << "emitting: " << concated_shingles << "\t" << name << std::endl;
		}
	}
}

void ShinglingClusterGen::delete_2d_array(unsigned int** array, int rows, int cols) {
	for (int i = 0; i < rows; i++) {
		delete[] array[i];
	}
	delete[] array;
}

void ShinglingClusterGen::draw_graph(char* key, int keybytes, char* values, int nvalues, int* valuebytes, KeyValue *kv, void* ptr) {
	int val1_offset = 0;
        int total = 0;

        for (int val1_index = 0; val1_index < (nvalues - 1); val1_index++)
        {
                int val2_index = val1_index + 1;
                int val2_offset = val1_offset + valuebytes[val1_index];
                for (val2_index; val2_index < nvalues; val2_index++)
                {
                        if ( strcmp(values + val1_offset, values + val2_offset) != 0 )
                        {
				total += (valuebytes[val1_index] + valuebytes[val2_index]) * 2;
			}
		}
	}	


	val1_offset = 0;
	for (int val1_index = 0; val1_index < (nvalues - 1); val1_index++)
	{
		int val2_index = val1_index + 1;
		int val2_offset = val1_offset + valuebytes[val1_index];
		for (val2_index; val2_index < nvalues; val2_index++)
		{
			if ( strcmp(values + val1_offset, values + val2_offset) != 0 )
			{
				char node1[valuebytes[val1_index]];
				char node2[valuebytes[val2_index]];
				memcpy(node1, values + val1_offset, valuebytes[val1_index]);
				memcpy(node2, values + val2_offset, valuebytes[val2_index]);
				//if (total > 1000000000) {
				//	if (valuebytes[val1_index] > 20 || valuebytes[val2_index] > 20) std::cerr << " too big of a valuebyte" << std::endl;
				//	std::cerr << "found suspicous spot: " << node1 << " " << node2 << " nvalues: " << nvalues << std::endl;
				//}
				kv -> add(node2, valuebytes[val2_index], node1, valuebytes[val1_index]);
				kv -> add(node1, valuebytes[val1_index], node2, valuebytes[val2_index]);
				total += (valuebytes[val1_index] + valuebytes[val2_index]) * 2;
//				kv -> add(values + val2_offset, valuebytes[val2_index],
//						values + val1_offset, valuebytes[val1_index]);
//				kv -> add(values + val1_offset, valuebytes[val1_index],
//						values + val2_offset, valuebytes[val2_index]);
			}
			val2_offset += valuebytes[val2_index];
		}
		val1_offset += valuebytes[val1_index];
	}
	//std::cerr << "total: " << total << std::endl;
}

void ShinglingClusterGen::redundant_agglomerator(char* key, int keybytes, char* values, int nvalues, int* valuebytes, KeyValue* kv, void* ptr) {
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

void ShinglingClusterGen::redundant_remover(char* key, int keybytes, char* values, int nvalues, int* valuebytes, KeyValue* kv, void* ptr) {
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

void ShinglingClusterGen::generate_vertex_neighbor_min_shingles(char* sequence_name, int name_length, char* neighbors_list, 
								int number_of_neighbors, int* neighbor_name_lengths, KeyValue* kv, 
								void* shingling_parameters) {
	ShinglingToolbox* shingling_params = (ShinglingToolbox*) shingling_parameters;
	std::vector<std::vector<unsigned int>>* hash_function_min_values = NULL;
	int new_rows_count = (shingling_params -> num_hash) - (shingling_params -> latest_hashed);
	int total_rows_count = shingling_params -> num_hash;
	//int nrows = shingling_params -> num_hash;
	int ncols = shingling_params -> shingle_size;
	std::string current_seq = (std::string) sequence_name;
	int me; MPI_Comm_rank(MPI_COMM_WORLD, &me);
	if (shingling_params -> num_hash == 21 && current_seq == "46625_1") {
		int old_n = 0;
		if(shingling_params -> persistent_graph.find(current_seq) != shingling_params -> persistent_graph.end()) {
			std::cout << "wait i know my old neighbors" << std::endl;
			old_n = shingling_params -> persistent_graph.at(current_seq).size();
		}
		std::cout << "old neighbors: " << old_n << " new neighbors: " << number_of_neighbors << std::endl;
	}
	//if (me == 0) std::cout << "latest_hashed: " << shingling_params -> latest_hashed << std::endl;
	//std::cout << "numbers: " << new_rows_count << " " << total_rows_count << " " << ncols << " " << shingling_params -> latest_hashed << " " << shingling_params -> num_hash << " " << number_of_neighbors << std::endl;
	//std::cerr << "for vertex " << sequence_name << " initializing a local table of nrows: " << new_rows_count << std::endl;
	//std::cout << "maybe entering if? " << new_rows_count << std::endl;	
	if (new_rows_count > 0) {
		int previously_computed = 0;
		//std::vector<unsigned int> empty_vector;
		//hash_function_min_values = new std::vector<std::vector<unsigned int>>;
		//hash_function_min_values -> insert(hash_function_min_values -> begin(), new_rows_count, empty_vector);
		//std::cerr << " initializied: " << hash_function_min_values -> size() << " rows" << std::endl;
		//unsigned int** hash_function_min_values;
		//initialize_array_by_max_uint(&hash_function_min_values, shingling_params -> num_hash, shingling_params -> shingle_size);
		unsigned int neighbor_name_hash_value;
		int neighbor_offset = 0;
		//std::cerr << "for node : " << sequence_name << std::endl;
		if (shingling_params -> persistent_graph.find(current_seq) == shingling_params -> persistent_graph.end()) {
			shingling_params -> persistent_graph.insert(std::make_pair(current_seq, std::unordered_set<std::string>()));
		}
		if (shingling_params -> persistent_vertex_hashes.find(current_seq) != shingling_params -> persistent_vertex_hashes.end()) {
				previously_computed = shingling_params -> persistent_vertex_hashes.at(current_seq) -> size();
		}
		for (int neighbor_index = 0; 
				neighbor_index < number_of_neighbors; 
				neighbor_offset += neighbor_name_lengths[neighbor_index], neighbor_index++) {
			std::string neighbor_name = (std::string) (neighbors_list + neighbor_offset);
			bool new_neighbor = true;
			if (shingling_params -> persistent_graph.at(current_seq).find(neighbor_name) == 
				shingling_params -> persistent_graph.at(current_seq).end()) {
				shingling_params -> persistent_graph.at(current_seq).insert(neighbor_name);
			} else {
				new_neighbor = false;
			}
			unsigned long long name_int_value;
			if (NUMBERED_NAMES == 1) {				
				if ((shingling_params -> vertex_hash).find(neighbor_name) == (shingling_params -> vertex_hash).end()) {
					name_int_value = compute_numbered_protname_int_value(neighbor_name);
					shingling_params -> vertex_hash.insert(std::make_pair(neighbor_name, name_int_value));
				} else {
					//std::cout << "read from  table" << std::endl;
					name_int_value = shingling_params -> vertex_hash.at(neighbor_name);
				}
				//std::cerr << "neighbor name: " << neighbors_list + neighbor_offset << ", its value: " << name_int_value;
			} else {
				std::cerr << "\n\nERROR: PROVIDE A HASH FUNCTION FOR ALPHABETICAL PROTEIN NAMES TO NUMBERS\n\n" << std::endl;
			}
			//for (int hash = 0; hash < shingling_params -> num_hash; hash++) {
			//int previously_computed = 0;
			//if (previously_computed > shingling_params -> latest_hashed) std::cout << "ERROR: prevcomputed cannot be larger than latest hashed" << std::endl;
			if (new_neighbor) {
				for (int hash = 0; hash < previously_computed; hash++) {
					//if (shingling_params -> persistent_vertex_hashes -> find(current_seq) ==
					//	shingling_params -> persistent_vertex_hashes -> end()) {
					//	std::vector<unsigned int> empty_vector;
					//	shingling_params -> persistent_vertex_hashes -> insert(std::make_pair(current_seq, new std::vector<std::vector<unsigned int>>(shingling_params -> latest_hashed, empty_vector)));
					//}
					neighbor_name_hash_value = compute_hash_value(name_int_value, shingling_params -> hash_randoms[hash * 2],
				                                        shingling_params -> hash_randoms[hash * 2 + 1], shingling_params -> Prime_num);
					//std::cerr << ", neighbor hash value: " << neighbor_name_hash_value << std::endl;
					if (hash < shingling_params -> latest_hashed) {
						//std::cout << "going into update earlier results " << hash << " " << neighbor_index << " " << previously_computed << std::endl;
						update_persistent_vertex_hashes(neighbor_name_hash_value, hash, &(shingling_params -> persistent_vertex_hashes), current_seq, ncols);
						//std::cout << "prev results updated" << std::endl;
					} /*else {
						//std::cout << "updating new table " << hash << " " << neighbor_index << std::endl;
						update_hash_func_min_values(neighbor_name_hash_value, hash_function_min_values, hash - shingling_params -> latest_hashed, ncols);
						//std::cout << "new table updated" << std::endl;
					}*/
				        //update_hash_func_min_values(neighbor_name_hash_value, &hash_function_min_values[hash][0], shingling_params -> shingle_size);
				}
			}
		}
		//std::cout << "prev computed: " << previously_computed << " latest hased: " << shingling_params -> latest_hashed << std::endl;
		//if (previously_computed > shingling_params -> latest_hashed) std::cout << "ERROR: prevcomputed cannot be larger than latest hashed" << std::endl;
		std::vector<unsigned int> empty_vector;
                hash_function_min_values = new std::vector<std::vector<unsigned int>>;
		new_rows_count = shingling_params -> num_hash - previously_computed;
                hash_function_min_values -> insert(hash_function_min_values -> begin(), new_rows_count, empty_vector);
		std::unordered_set<std::string>* neighbors_list = &(shingling_params -> persistent_graph.at(current_seq));
		//std::cout << "updating for " << neighbors_list -> size() << " new neighbors " << std::endl;
		for (auto neighb = neighbors_list -> begin(); neighb != neighbors_list -> end(); neighb++) {
			//if (shingling_params -> vertex_hash.find(*neighb) == shingling_params -> vertex_hash.end()) {
			//	std::cout << "couldn't find in the vertex hash " << std::endl;
			//	std::cout << "vertex_hash size: " << shingling_params -> vertex_hash.size() << std::endl;
			//} else { std::cout << "vertex hash found" << std::endl; }
			unsigned long long name_int_value = shingling_params -> vertex_hash.at(*neighb);
			for (int hash = previously_computed; hash < shingling_params -> num_hash; hash++) {
				neighbor_name_hash_value = compute_hash_value(name_int_value, shingling_params -> hash_randoms[hash * 2],
								shingling_params -> hash_randoms[hash * 2 + 1], shingling_params -> Prime_num);
				update_hash_func_min_values(neighbor_name_hash_value, hash_function_min_values, hash - previously_computed, ncols);
			}
		}
		//std::cout << "adding to the persistent table" << std::endl;
		if (shingling_params -> persistent_vertex_hashes.find(current_seq) ==
			shingling_params -> persistent_vertex_hashes.end()) {
			shingling_params -> persistent_vertex_hashes.insert(
                                      std::make_pair(current_seq, new std::vector<std::vector<unsigned int>>));
		}
		shingling_params -> persistent_vertex_hashes.at(current_seq) -> insert(
				shingling_params -> persistent_vertex_hashes.at(current_seq) -> end(), 
				hash_function_min_values -> begin(), hash_function_min_values -> end());
		if (hash_function_min_values) delete hash_function_min_values;
		//std::cout << "added to the persistent" << std::endl;
		//emit_min_shingle_kvs(hash_function_min_values, kv, shingling_params -> num_hash, shingling_params -> shingle_size, sequence_name, name_length);
		//delete_2d_array(hash_function_min_values, shingling_params -> num_hash, shingling_params -> shingle_size);
	/*	if (shingling_params -> last_iter) {
			//append hash_function_min_values to the persistent_hash_table
			//std::cerr << "vertex: appending new rows to the persistent_hash_table" << std::endl;
			if (shingling_params -> persistent_hash_function_min_values -> find(std::string(sequence_name)) ==
				shingling_params -> persistent_hash_function_min_values -> end()) {
			//	std::cerr << "the name: " << sequence_name << " didn't exist in the map, creating now to append" << std::endl;
				shingling_params -> persistent_hash_function_min_values -> insert(
						std::make_pair(std::string(sequence_name), hash_function_min_values));
			//	std::cerr << "the name " << sequence_name << " is now added to the map." << std::endl;
			} else {
			//	std::cerr << "existed: " << sequence_name << ". just appending and deleting" << std::endl;
				shingling_params -> persistent_hash_function_min_values -> at(std::string(sequence_name)) -> insert(
						shingling_params -> persistent_hash_function_min_values -> at(std::string(sequence_name)) -> end(), 
						hash_function_min_values -> begin(), hash_function_min_values -> end());
				delete hash_function_min_values;
			}
			//emit shingles from persistent_hash_table
	//		std::cerr << "emitting from the persistent_table" << std::endl;
	//		std::cerr << "printing table for " << sequence_name << " for num_hash " << total_rows_count << std::endl;
	//		for (int i = 0; i < persistent_hash_function_min_values -> at(sequence_name) -> size(); i++) {
	//			for (int j =0; j < persistent_hash_function_min_values -> at(sequence_name) -> at(i).size(); j++) {
	//				std::cerr << "\t" << persistent_hash_function_min_values -> at(sequence_name) -> at(i).at(j);
	//			}
	//			std::cerr << std::endl;
	//		}
	//		std::cerr << std::endl;
			emit_min_shingle_kvs(shingling_params -> persistent_hash_function_min_values -> at(sequence_name), kv, total_rows_count, ncols, sequence_name, name_length);
		} else {
	*/		//emit shingles from hash_function_min_values
			//std::cout << "coming from inside if" << std::endl;
	}
	//emit_min_shingle_kvs(hash_function_min_values, kv, nrows, ncols, sequence_name, name_length);
	//std::cout << "emitting new kvs" << std::endl;
	if (shingling_params -> num_hash == 21 && current_seq == "46625_1") {
		std::cout << "table in optimized: " << std::endl;
		for (auto row = shingling_params -> persistent_vertex_hashes.at(current_seq) -> begin();
			row != shingling_params -> persistent_vertex_hashes.at(current_seq) -> end(); row++) {
			for (auto col = row -> begin(); col != row -> end(); col ++) {
				std::cout << *col << " ";
			}
			std::cout << std::endl;
		}
	}
	if (shingling_params -> persistent_vertex_hashes.at(current_seq) -> size() != total_rows_count) std::cout << "ERROR: for " << total_rows_count << " hashes there are " << shingling_params -> persistent_vertex_hashes.at(current_seq) -> size() << " rows in the table. for seq:" << sequence_name << std::endl;
	emit_min_shingle_kvs(shingling_params -> persistent_vertex_hashes.at(current_seq), kv, total_rows_count, ncols, sequence_name, name_length);
	//delete hash_function_min_values;
//	}
}

unsigned long long ShinglingClusterGen::compute_numbered_protname_int_value(std::string name) {
	unsigned long long result = 0;
	int seperator_index = name.find('_');
	if (seperator_index == std::string::npos) {
		result = std::stoi(name);
	} else {
		result = std::stoi(name.substr(0, seperator_index));
		result += (1e7 * ( std::stoi( name.substr(seperator_index + 1) ) + 1));
	}
	return result;
}

void ShinglingClusterGen::update_persistent_vertex_hashes(unsigned int neighbor_name_hash_value, int hash, 
				unordered_map<std::string, std::vector<std::vector<unsigned int>>*>* persistent_vertex_hashes, 
				std::string current_seq, int ncols) {
	auto row_begin_iter = persistent_vertex_hashes -> at(current_seq) -> at(hash).begin();
        bool row_is_full = false;
        int insert_pos = 0;

	for (insert_pos = 0; insert_pos < persistent_vertex_hashes -> at(current_seq) -> at(hash).size(); insert_pos++) {
		if (neighbor_name_hash_value <= (persistent_vertex_hashes -> at(current_seq) -> at(hash))[insert_pos]) {
			break;
		}
	}

	if (insert_pos < ncols) {
                if (persistent_vertex_hashes -> at(current_seq) -> at(hash).size() + 1 > ncols) {
                        persistent_vertex_hashes -> at(current_seq) -> at(hash).pop_back();
                }
                persistent_vertex_hashes -> at(current_seq) -> at(hash).insert(row_begin_iter + insert_pos, neighbor_name_hash_value);
	}
}

void ShinglingClusterGen::label_propagate(char* node, int node_length, char* messages, int num_messages, int* message_lengths, KeyValue* kv, void* ptr) {
	if (log_scale > 2 && (log_proc == me || log_proc == -1)) {
		int off = 0;
		for (int i = 0; i < num_messages; off += message_lengths[i], i++) {
			std::cerr << messages + off << ", ";
		}
		std::cerr << '\n';
	}
	MessageValue current_label;
	MessageValue min_label(node);
	MessageValue current_node(node);
	std::vector<MessageValue> neighbors;
	int neighbors_count = 0;
	//if (log_scale > 1 && (log_proc == -1)) std::cerr << "extracting messages" << std::endl;
	extract_messages_info(messages, message_lengths, num_messages, &current_label, &min_label, &neighbors, &neighbors_count);
	current_node.set_type(NEIGHBOR_NODE);
	//if (log_scale > 1 && (log_proc == -1)) std::cerr << "messages extracted" << std::endl;
	//re-send the neighbors list to yourself:
	for (int neighbor = 0; neighbor < neighbors_count; neighbor++) {
		kv -> add(current_node.get_message_array(), current_node.get_length(), 
				neighbors[neighbor].get_message_array(), neighbors[neighbor].get_length());
		if (log_scale > 2 && (log_proc == me || log_proc == -1)) std::cerr << "repeating neighbor for " << current_node.get_message_array() << ": " << neighbors[neighbor].get_message_array() << "s" << neighbors[neighbor].get_length() << "\n";
	}
	//if (log_scale > 1 && (log_proc == -1)) std::cerr << "neighbors re-generated" << std::endl;

	//send the label in your cluster to yourself
	MessageValue label_to_send(CURRENT_LABEL, min_label.get_name());
	label_to_send.set_type(CURRENT_LABEL);
	kv -> add(current_node.get_message_array(), current_node.get_length(), 
			label_to_send.get_message_array(), label_to_send.get_length());
	if (log_scale > 2 && (log_proc == me || log_proc == -1)) std::cerr << "setting label for " << current_node.get_message_array() << "s" << current_node.get_length() << " to: " << label_to_send.get_message_array() << "s" << label_to_send.get_length() << "\n";
	//if (log_scale > 1 && (log_proc == -1)) std::cerr << "new label sent" << std::endl;

	//if label has changed inform the neighbors
	if ( (current_label.get_message_array() == NULL) || strcmp(min_label.get_name(), current_label.get_name()) != 0) {
		label_to_send.set_type(INCOMING_LABEL);
		for (int neighbor = 0; neighbor < neighbors_count; neighbor++) {
			kv -> add(neighbors[neighbor].get_message_array(), neighbors[neighbor].get_length(), 
					label_to_send.get_message_array(), label_to_send.get_length());
			if (log_scale > 2 && (log_proc == me || log_proc == -1)) std::cerr << "informing: " << neighbors[neighbor].get_message_array() << "s" << neighbors[neighbor].get_length() << " of new label " << label_to_send.get_message_array() << "s" << label_to_send.get_length() << "\n";
		}
	}
	else {
		if (log_scale > 2 && (log_proc == me || log_proc == -1)) std::cerr << "label not changed" << "\n";
	}
	if (log_scale > 2 && (log_proc == me || log_proc == -1)) std::cerr << std::endl;
}

void ShinglingClusterGen::extract_messages_info(char* messages, int* message_lengths, int num_messages, MessageValue* current_label,
				MessageValue* min_label, std::vector<MessageValue>* neighbors, int* neighbor_count) {
	int message_offset = 0;
	MessageValue current_read;
	for (int message_index = 0; message_index < num_messages; message_index++) {
		//if (log_scale > 1 && (log_proc == -1)) std::cerr << "creating object" << std::endl;
		current_read = MessageValue(&messages[message_offset]);
		//if (log_scale > 1 && (log_proc == -1)) std::cerr << "object created" << std::endl;
		if ( strcmp(current_read.get_name(), min_label -> get_name() ) < 0 ) {
			*min_label = current_read;
			//if (log_scale > 1 && (log_proc == -1)) std::cerr << "min updated" << std::endl;
		}
		if (current_read.get_type() == CURRENT_LABEL) {
			//if (log_scale > 1 && (log_proc == -1)) std::cerr << "label found" << std::endl;
			*current_label = current_read;
		} else if (current_read.get_type() == NEIGHBOR_NODE) {
			//if (log_scale > 1 && (log_proc == -1)) std::cerr << "adding to neighbors list" << std::endl;
			//neighbors -> at(*neighbor_count) = current_read;
			neighbors -> push_back(current_read);
			*neighbor_count = *neighbor_count + 1;
		}
		//if (log_scale > 1 && (log_proc == -1)) std::cerr << "going for next message" << std::endl;
		message_offset += message_lengths[message_index];
	}
	//neighbors -> shrink_to_fit();
}

void ShinglingClusterGen::label_finalize(char* key, int keybytes, char* values, int nvalues, int* valuebytes, KeyValue* kv, void* ptr) {
	int label_index;
	int label_offset = 0;
	for (int val_index = 0; val_index < nvalues; val_index++)
	{
		if (values[label_offset] == CURRENT_LABEL)
		{
			label_index = val_index;
			break;	
		}
		label_offset += (valuebytes[val_index]);
	}
	if (label_index == nvalues)
	{
		std::cerr << "ERROR: label not found!" <<'\n';
	}
	kv -> add(values + label_offset, valuebytes[label_index], key, keybytes);
}

void ShinglingClusterGen::rename_sequences(char* cluster_label, int label_length, char* seq_names, int num_members, 
						int* seq_name_lengths, KeyValue* kv, void* ptr) {
	for (int index = 0, offset = 0; index < num_members; offset += seq_name_lengths[index], index++) {
		char* domained_seq_name = strtok(seq_names + offset + 1, "_");
		kv -> add(cluster_label, label_length, domained_seq_name, strlen(domained_seq_name) + 1);
	}
}

void ShinglingClusterGen::remove_small_clusters(char* cluster_label, int label_length, char* member_names, int num_members, 
							int* member_name_lengths, KeyValue* kv, void* cluster_limits) {
	ClusterLimits* limits = (ClusterLimits*) cluster_limits;
	int min_cluster_size_limit = limits -> min_cluster_size_limit;
	limits -> largest_cluster_size = limits -> largest_cluster_size < num_members ? num_members : limits -> largest_cluster_size;
	if (num_members >= min_cluster_size_limit) {
		for (int index = 0, offset = 0; index < num_members; offset += member_name_lengths[index], index++) {
			kv -> add(cluster_label, label_length, member_names + offset, member_name_lengths[index]);
		}
	}
}

std::vector< std::set<std::string> > ShinglingClusterGen::collect_clusters(MapReduce* mr, MpiGatherParameters gather_params, 
								int* largest_cluster_size, std::vector<std::string>* cluster_labels) {
	gather_params.local_clusters_size = 0;
	std::vector< std::set<std::string> > clusters;
	std::vector<std::string> local_cluster_labels;
	MpiGatherToolbox gather_toolbox(&gather_params, &clusters, &local_cluster_labels);
	mr -> reduce(local_cluster_agglomerator, &gather_toolbox);
	mr -> convert();
	if (log_scale > 2 && (log_proc == me || log_proc == -1)) std::cerr << "my size: " << gather_params.local_clusters_size << std::endl;
	if (gather_params.me == gather_params.root) {
		gather_params.proc_clusters_sizes = new int[gather_params.num_processors];
	}
	MPI_Gather(&gather_params.local_clusters_size, 1, MPI_INT, gather_params.proc_clusters_sizes, 1, MPI_INT, gather_params.root, MPI_COMM_WORLD);
	if (gather_params.me == gather_params.root) {
		for (int proc = 0; proc < gather_params.num_processors; proc++) {
			gather_params.global_clusters_size += gather_params.proc_clusters_sizes[proc];
			if (log_scale > 2 && (log_proc == me || log_proc == -1)) std::cerr << "rcvd: " << gather_params.proc_clusters_sizes[proc] << std::endl;
		}
	}
	gather_toolbox.largest_cluster_size = 0;
	gather_toolbox.gather_params = &gather_params;
	mr -> reduce(send_local_clusters, &gather_toolbox);
	*largest_cluster_size = gather_toolbox.largest_cluster_size;
	MpiUtil::gatherv_string_container(&local_cluster_labels, gather_params.root, cluster_labels, MPI_COMM_WORLD);
	std::reverse(cluster_labels -> begin(), cluster_labels -> end());
	return clusters;
}

void ShinglingClusterGen::local_cluster_agglomerator(char* cluster_name, int cluster_name_length, char* members, int num_members, int* member_name_lengths, KeyValue* kv, void* gather_toolbox_ptr) {	
	MpiGatherToolbox* gather_toolbox = (MpiGatherToolbox*) gather_toolbox_ptr;
	int cluster_size = 0;
	for (int i = 0; i < num_members; i++) {
		cluster_size += member_name_lengths[i];
	}
	//cluster_size -= num_members;
	cluster_size++;
	gather_toolbox -> gather_params -> local_clusters_size += cluster_size;
	char cluster[cluster_size];
	cluster[cluster_size - 1] = '$';
	memcpy(&cluster[0], members, cluster_size - 1); 
	//int copy_offset = 0;
	//for (int member = 0, offset = 0; member < num_members; offset += member_name_lengths[member], member++) {
	//	copy_offset = offset - member;
	//	strcpy(cluster + copy_offset, members + offset + 1);
	//}
	kv -> add(NULL, 0, cluster, cluster_size);
	gather_toolbox -> cluster_labels -> push_back(std::string(cluster_name));
}

void ShinglingClusterGen::send_local_clusters(char* empty_key, int keybytes, char* agglomerated_clusters, int nvalues_one, int* agg_clusters_size, KeyValue* kv, void* gather_toolbox_ptr) {
	MpiGatherToolbox* gather_toolbox = (MpiGatherToolbox*) gather_toolbox_ptr;
	int root = gather_toolbox -> gather_params -> root;
	int me = gather_toolbox -> gather_params -> me;
	int local_size = gather_toolbox -> gather_params -> local_clusters_size;
	int nprocs = gather_toolbox -> gather_params -> num_processors;
	if (me == root) {
		if (log_scale > 2 && (log_proc == me || log_proc == -1)) std::cerr << "sending: " << agglomerated_clusters << " and " << (int) (agglomerated_clusters + strlen(agglomerated_clusters))[0] << " and " << (int)(agglomerated_clusters + strlen(agglomerated_clusters) + 1)[0] << std::endl;
		extract_clusters_from_array(agglomerated_clusters, local_size, gather_toolbox);
	}
	for (int proc = 0; proc < nprocs; proc++) {
		if (proc != root) {
			if (me == root) {
				int rcv_size = gather_toolbox -> gather_params -> proc_clusters_sizes[proc];
				char* rcv_buf = new char[rcv_size];
				MPI_Recv(rcv_buf, rcv_size / sizeof(char), MPI_CHAR, proc, proc, MPI_COMM_WORLD, NULL);
				extract_clusters_from_array(rcv_buf, rcv_size, gather_toolbox);
				delete[] rcv_buf;
			} else if (me == proc) {
				MPI_Send(agglomerated_clusters, local_size / (sizeof(char)), MPI_CHAR, root, proc, MPI_COMM_WORLD);
			}
		}
	}
	if (me == root) {
		if (log_scale > 2 && (log_proc == me || log_proc == -1)) std::cerr << "number of extracted clusters: " << gather_toolbox -> clusters -> size() << std::endl;
		if (log_scale > 2 && (log_proc == me || log_proc == -1)) {
			for (auto iter = gather_toolbox -> clusters -> at(gather_toolbox -> clusters -> size() - 1).begin(); 
					iter != gather_toolbox -> clusters -> at(gather_toolbox -> clusters -> size() - 1).end(); iter++) {
				std::cerr << *iter <<" \t";
			}
		}
	}
}

void ShinglingClusterGen::extract_clusters_from_array(char* clusters_array, int array_size, MpiGatherToolbox* gather_toolbox) {
	int offset = 0;
	std::set< std::string > cluster;
	while (offset < array_size) {
		if (clusters_array[offset] == '$') {
			if (cluster.size() > 0) {
				if (cluster.size() > gather_toolbox -> largest_cluster_size) {
					gather_toolbox -> largest_cluster_size = cluster.size();
				}
				gather_toolbox -> clusters -> push_back(cluster);
				cluster.clear();
				offset++;
			}
		} else {
			int member_length = strlen(clusters_array + offset) + 1;
			char temp[member_length];
			strcpy(temp, clusters_array + offset);
			cluster.insert( (std::string) temp );
			offset += member_length;
		}
	}
}

bool ShinglingClusterGen::file_exists (const std::string& name) {
	struct stat buffer;   
	return (stat (name.c_str(), &buffer) == 0); 
}

std::vector< std::set<std::string> > ShinglingClusterGen::read_shingling_clusters(std::string filename, int* largest_cluster_size, 
											bool rename, int min_cluster_limit,
											std::vector<std::string>* cluster_labels,
											std::unordered_set<std::string> clusters_to_remove) {
	int me;
	MPI_Comm_rank(MPI_COMM_WORLD, &me);
	*largest_cluster_size = 0;
	if (me == 0) std::cerr << "reading file: " << filename << std::endl;
	std::ifstream infile(filename);
	std::string line;
	std::vector< std::set<std::string> > clusters;
	int index_counter = 0;
	while (std::getline(infile, line)) {
		std::set<std::string> cluster;
		std::vector< std::string > numbered_members = split_string(line, ' ');
		int members_size = numbered_members.size();
		std::string cluster_label = numbered_members[10].substr(0, numbered_members[10].length() - 1);
		if (clusters_to_remove.find(cluster_label) == clusters_to_remove.end()) {
			for (int i = 12; i < members_size; i++) {
				std::string member;
				if (rename) {
					int length = numbered_members[i].find('_') - 1;
					std::string member = numbered_members[i].substr(1, length);
					cluster.insert(member);
				} else {
					cluster.insert(numbered_members[i]);
				}
			}
			int cluster_size = cluster.size();
			if (cluster_size == 5923) {
				//std::cerr << "this is the largest cluster: \n" << cluster_label << std::endl;
			}
			if (cluster_size >= min_cluster_limit) {
				*largest_cluster_size = cluster_size > *largest_cluster_size ? cluster_size : *largest_cluster_size;
				clusters.push_back(cluster);
				cluster_labels -> push_back(cluster_label);
				//std::cerr << "cluster label: " << numbered_members[10] << " in index " << index_counter << std::endl;
				index_counter++;
			}
		}
	}
	return clusters;
}

std::vector< std::string > ShinglingClusterGen::split_string(std::string str, char delim) {
	std::vector< std::string > tokens;
	std::stringstream ss(str);
	std::string token;
	while (std::getline(ss, token, delim)) {
		tokens.push_back(token);
	}
	return tokens;
}

void ShinglingClusterGen::distribute_clusters(int itask, KeyValue* kv, void* distribute_util_ptr) {
	DistributeUtil* util = (DistributeUtil*) distribute_util_ptr;
	Clusters clusters = util -> clusters;
	int me, nprocs;
	MPI_Comm_rank(MPI_COMM_WORLD, &me);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	int clusters_count = clusters.clusters -> size();
	//std::cerr << "all clusters count: " << clusters_count << std::endl;
	int proc_share, start_index, end_index;
	if (util -> single_proc) {
		proc_share = clusters_count;
		start_index = 0;
		end_index = clusters_count;
	} else {		
		proc_share = ceil( (double) clusters_count / nprocs );
		start_index = me * proc_share;
		end_index = std::min((me + 1) * proc_share, clusters_count);
	}
	//std::cerr << "proc" << me << " start: " << start_index << " end: " << end_index << " share: " << proc_share << std::endl;
	for (int i = start_index; i < end_index; i++) {
		std::string label_str = clusters.cluster_labels -> at(i);
		int label_length = label_str.length() + 1;
		char label[label_length];
		strcpy(label, label_str.c_str());
		for (auto iter = clusters.clusters -> at(i).begin(); iter != clusters.clusters -> at(i).end(); iter++) {
			std::string member_str = *iter;
			int member_length = member_str.length() + 1;
			char member[member_length];
			strcpy(member, member_str.c_str());
			//std::cerr << "adding kv now: " << label << " " << member << std::endl;
			kv -> add(label, label_length, member, member_length);
		}
	}
	//std::cerr << "distributed" << std::endl;
}

void ShinglingClusterGen::perform_label_propagate(MapReduce* mr, int nnodes, int nedges, int connected_component_iterations, 
						double* connected_comp_time) { //hash_min algorithm
	double tstart = MPI_Wtime();
	if (log_scale > 1 && (log_proc == me || log_proc == -1)) std::cerr << "starting connected component" << std::endl;
	int nmessages = mr -> reduce (label_propagate, NULL);
	//mr -> print("interfiles_label1", 0, -1, 1, 5, 5);
	//mr -> print ("test_res", 0, -1, 1, 5, 5);
	while (connected_component_iterations != 0 && nmessages > nedges + nnodes) {
		if (log_scale > 2 && (log_proc == me || log_proc == -1))std::cerr << "in connected component loop" << std::endl;
		//mr -> print ("test_res", 0, -1, 1, 5, 5);
		mr -> collate(NULL);
		nmessages = mr -> reduce(label_propagate, NULL);
		//mr -> print("interfiles_labelnew", 0, -1, 1, 5, 5);
		connected_component_iterations--;
	}
	double tend = MPI_Wtime();
	if (connected_comp_time) {
		*connected_comp_time = tend - tstart;
	}
	if (log_scale > 1 && (log_proc == me || log_proc == -1)) std::cerr << "cc computed. last collate... " << std::endl;
	mr -> collate(NULL);
	if (log_scale > 1 && (log_proc == me || log_proc == -1)) std::cerr << "finalizing labels" << std::endl;
	mr -> reduce(label_finalize, NULL);
	//mr -> print("interfiles_finalized", 0, -1, 1, 5, 5);
	mr -> convert(); //changed  from collate
}

void ShinglingClusterGen::perform_union_find(MapReduce* mr, std::vector<std::set<std::string>>* final_clusters, bool collect_in_root, 
						MpiGatherParameters gather_params) {
	int me; MPI_Comm_rank(MPI_COMM_WORLD, &me);
	if (log_scale > 1 && (log_proc == gather_params.root || log_proc == -1)) std::cerr << "union finding" << std::endl;
	//MPI_Barrier(MPI_COMM_WORLD);

//	mr -> reduce(add_label_to_cluster, NULL);
//	mr -> convert();
	int largest_cluster_size = 0;
	std::vector<std::string> cluster_labels;
	if (log_scale > 1 && (log_proc == gather_params.root || log_proc == -1)) std::cerr << "collecting clusters" << std::endl;
	//MPI_Barrier(MPI_COMM_WORLD);
	//mr -> print("edgescheck", 0, -1, 1, 5, 5);
	std::vector< std::set< std::string> > clusters_vector = collect_clusters(mr, gather_params, 
											&largest_cluster_size, &cluster_labels);
	if (me == gather_params.root) {
		for (int i = 0; i < clusters_vector.size(); i++) {
			clusters_vector[i].insert(cluster_labels[i]);
		}
	}
	//std::cerr << "THE UNION FIND METHOD IS NOT IMPLEMENTED COMPLETELY.\n";
	//std::cerr << "IT MESSES THE CLUSTER LABELS SO THAT THEY ARE NOT GOOD FOR FINDING BAD MERGERS ANYMORE.\n";
	//std::cerr << "EXITING THE PROGRAM..." << std::endl;
	//exit(0);
	if (log_scale > 1 && (log_proc == gather_params.root || log_proc == -1)) std::cerr << "running ufo" << std::endl;
	//MPI_Barrier(MPI_COMM_WORLD);
	unordered_map<std::string, std::vector<std::string>> clusters_map;
	if (me == gather_params.root) {
		UnionFind ufo;
		ufo.run_connected_component(&clusters_vector, &clusters_map);
		if (log_scale > 1 && (log_proc == gather_params.root || log_proc == -1)) std::cerr << "ufo completed: " << 
												clusters_map.size() << std::endl;
	//	std::cerr << "sizecheck: " << (clusters_map.begin()) -> first << " " << ((clusters_map.begin()) -> second).size() << std::endl;
	//	MPI_Barrier(MPI_COMM_WORLD);
		clusters_vector.clear();
		cluster_labels.clear();
	//	MPI_Barrier(MPI_COMM_WORLD);
		//if (me == 0) std::cerr << "distributing clusters" << std::endl;
	//	MPI_Barrier(MPI_COMM_WORLD);
	}
	//MPI_Barrier(MPI_COMM_WORLD);
	if (!collect_in_root) {
		if (log_scale > 1 && (log_proc == gather_params.root || log_proc == -1)) std::cerr << "distributing clusters" << std::endl;
		if (me == gather_params.root) {
			clusters_vector.reserve(clusters_map.size());
			cluster_labels.reserve(clusters_map.size());
			for (auto iter = clusters_map.begin(); iter != clusters_map.end(); iter++) {
				clusters_vector.push_back(std::set<std::string>(iter -> second.begin(), iter -> second.end()));
				cluster_labels.push_back(iter -> first);
			}
	//		std::cerr << clusters_vector[0].size() << " for sizecheck " << cluster_labels[0] << std::endl;
		}
		int nprocs;
                MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
		bool single_proc_owns_container = true;
		Clusters clusters(&clusters_vector, &cluster_labels);
		DistributeUtil distribute_util(clusters, single_proc_owns_container);
		mr -> map(nprocs, distribute_clusters, &distribute_util);
		mr -> convert();
	} else {
		if (log_scale > 1 && (log_proc == gather_params.root || log_proc == -1)) std::cerr << "returning the container" << std::endl;
		if (me == gather_params.root) {
			final_clusters -> insert(final_clusters -> begin(), clusters_vector.begin(), clusters_vector.end());
		}
	}
}

void ShinglingClusterGen::perform_hash_to_min_connected_component(MapReduce* mr, int nnodes, int nedges, int connected_component_iterations, 
								double* connected_comp_time) {
	double t_start = MPI_Wtime();
	mr -> reduce(add_label_to_cluster, NULL); //appends the labelname to the cluster (make a copy of the key in values)
	mr -> convert();
	bool completed = false;
	int loop_counter = 0;
	int me; MPI_Comm_rank(MPI_COMM_WORLD, &me);
	while (!completed && connected_component_iterations != 0) {
		bool local_completed = true;
		//char fnamearray[100];
		//std::string filename = "cc_check" + std::to_string(loop_counter);
		//strcpy(fnamearray, filename.c_str());
		//mr -> reduce(MrToolbox::redundant_agglomerator, NULL);
		//mr -> convert();
		//mr -> reduce(MrToolbox::redundant_remover, NULL);
		//mr -> convert();
		mr -> reduce(hash_to_min_hasher, &local_completed);
		mr -> collate(NULL);
		//mr -> print(fnamearray, 0, -1, 1, 5, 5);
		MPI_Allreduce(&local_completed, &completed, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD);
		//std::cerr << "me: " << me << " local:" << local_completed << " completed:" << completed << std::endl;
	//	if (me == 0) std::cerr << "overall bool:" << completed << std::endl;
		connected_component_iterations--;
		loop_counter++;
	}
	if (me == 0) std::cerr << "loop counter: " << loop_counter << std::endl;
	mr -> reduce(save_completed_components, NULL);
	mr -> convert(); //change from collate
	//mr -> reduce(MrToolbox::redundant_agglomerator, NULL);
	//mr -> convert();
	//mr -> reduce(MrToolbox::redundant_remover, NULL);
	//mr -> convert();
	double t_end = MPI_Wtime();
	*connected_comp_time = t_end - t_start;
}

void ShinglingClusterGen::add_label_to_cluster(char* key, int keylength, char* values, int nvalues, int* valuebytes, KeyValue* kv, void* ptr) {
	kv -> add(key, keylength, key, keylength);
	for (int i = 0, offset = 0; i < nvalues; offset += valuebytes[i], i++) {
		kv -> add(key, keylength, values + offset, valuebytes[i]);
	}
}

void ShinglingClusterGen::hash_to_min_hasher(char* node, int node_length, char* members, int nmembers, int* member_lengths, 
							KeyValue* kv, void* label_unchanged_ptr) {
	bool* label_unchanged = (bool*) label_unchanged_ptr;
	std::unordered_set<std::string> neighbors;
	int total_members_length = 0;
	for (int i = 0; i < nmembers; total_members_length += member_lengths[i], i++);
	std::vector<char> members_vector(members, members + total_members_length);
	MpiUtil::extract_strings(members_vector.begin(), members_vector.end(), &neighbors);
	std::string min_neighbor, current_label;
	find_label_and_min_neighbor(&neighbors, &min_neighbor, &current_label);
	//*label_unchanged = false;
	if (min_neighbor != std::string(node) && nmembers != 1) {
		*label_unchanged = false;
	}
	//std::cerr << "node:" << node << " min:" << min_neighbor << " nmembers:" << nmembers << " completed:" << *label_unchanged << std::endl;
	int new_label_length = min_neighbor.length() + 1;
	char new_label[new_label_length];
	strcpy(new_label, min_neighbor.c_str());
	for (auto iter = neighbors.begin(); iter != neighbors.end(); iter++) {
		std::string neighbor = *iter;
	//	if (neighbor[0] != HTM_CC_LABEL_INDICATOR) {
			int neighbor_length = neighbor.length() + 1;
			char neighbor_array[neighbor_length];
			strcpy(neighbor_array, neighbor.c_str());
	//		if (min_neighbor != std::string(node)) {
			//send all your members to the new label
			//if (neighbor != min_neighbor) {
				kv -> add(new_label, new_label_length, neighbor_array, neighbor_length);
			//}
			//send the new label to your neigbors
	//		if (neighbor != std::string(node)) {
				kv -> add(neighbor_array , neighbor_length, new_label, new_label_length);
	//		}
			//repeat the neighbor for yourself
			//kv -> add(node, node_length, neighbor_array, neighbor_length);
	//	}
	}
	//min_neighbor = HTM_CC_LABEL_INDICATOR + min_neighbor;
	//int new_marked_label_length = min_neighbor.length() + 1;
	//char new_marked_label[new_marked_label_length];
	//strcpy(new_marked_label, min_neighbor.c_str());
	//send the new label to yourself
	//kv -> add(node, node_length, new_marked_label, new_marked_label_length);
}

void ShinglingClusterGen::find_label_and_min_neighbor(std::unordered_set<std::string>* neighbors, std::string* min_neighbor, 
							std::string* label) {
	*label = "";
	std::string current = *neighbors -> begin();
	if (current[0] == HTM_CC_LABEL_INDICATOR) {
		current = current.substr(1);
		*label = current;
	}
	*min_neighbor = current;
	for (auto iter = std::next(neighbors -> begin()); iter != neighbors -> end(); iter++) {
		current = *iter;
		if (current[0] == HTM_CC_LABEL_INDICATOR) {
			current = current.substr(1);
			*label = current;
		}
		if (current.compare(*min_neighbor) < 0) {
			*min_neighbor = current;
		}
	}
}

void ShinglingClusterGen::save_completed_components(char* node, int node_length, char* members, int nmembers, int* member_lengths,
							KeyValue* kv, void* ptr) {
	std::unordered_set<std::string> neighbors;
        int total_members_length = 0;
        for (int i = 0; i < nmembers; total_members_length += member_lengths[i], i++);
        std::vector<char> members_vector(members, members + total_members_length);
        MpiUtil::extract_strings(members_vector.begin(), members_vector.end(), &neighbors);
        std::string min_neighbor, current_label;
        find_label_and_min_neighbor(&neighbors, &min_neighbor, &current_label);
	if (min_neighbor == std::string(node)) {
		int marked_node_length = node_length + 1;
		char marked_node[marked_node_length];
		marked_node[0] = '@';
		memcpy(&marked_node[1], node, node_length);
		for (auto iter = neighbors.begin(); iter != neighbors.end(); iter++) {
			std::string neighbor = "/" + *iter;
			int neighbor_length = neighbor.length() + 1;
                        char neighbor_array[neighbor_length];
                        strcpy(neighbor_array, neighbor.c_str());
			kv -> add(marked_node, marked_node_length, neighbor_array, neighbor_length);			
		}
	}
/*	
	int label_index = -1, label_offset = -1;
	for (int i = 0, offset = 0; i < nmembers; offset += member_lengths[i], i++) {
		if (members[offset] == HTM_CC_LABEL_INDICATOR) {
			label_index = i;
			label_offset = offset;
			break;
		}
	}
	if (strcmp(node, members + label_offset + 1) == 0) {
		for (int i = 0, offset = 0; i < nmembers; offset += member_lengths[i], i++) {
			if (i != label_index) {
				kv -> add(node, node_length, members + offset, member_lengths[i]);
			}
		}
	}
*/
}

int ShinglingClusterGen::send_to_owner_proc(char* key, int keybytes) {
	return seq_owner_proc[std::string(key)];
}
/*
void ShinglingClusterGen::populate_seq_owners_proc(std::unordered_map<std::string, std::string>* local_seqs) {
	int local_seqs_count = local_seqs -> size();
	std::vector<std::string> local_seq_names;
	std::vector<std::string> global_seq_names;
	for (auto iter = local_seqs -> begin(); iter != local_seqs -> end(); iter++) {
		local_seq_names.push_back(iter -> first);
	}
	int nproc;
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	int local_counts[nproc];
	MPI_Allgather(&local_seqs_count, 1, MPI_INT, local_counts, 1, MPI_INT, MPI_COMM_WORLD);
	int proc_end_nums[nproc];
	int marginal_sum = 0;
	for (int i = 0; i < nproc; i++) {
		marginal_sum += local_counts[i];
		proc_end_nums[i] = marginal_sum;
	}
	MpiUtil::allgatherv_string_container(&local_seq_names, &global_seq_names);
	int proc_counter = 0;
	for (int i = 0; i < global_seq_names.size(); i++) {
		if (i == proc_end_nums[proc_counter]) {
			proc_counter++;
		}
		seq_owner_proc.insert(std::make_pair(global_seq_names[i], proc_counter));
	}
}
*/
void ShinglingClusterGen::establish_sequence_ownership(char* name, int name_length, char* sequence, int seq_length, void* local_seqs_ptr) {
	std::unordered_map<std::string, std::string>* local_seqs = (std::unordered_map<std::string, std::string>*) local_seqs_ptr;
	std::string seq_name = std::string(name);
	std::string seq = std::string(sequence);
	local_seqs -> insert(std::make_pair(seq_name, seq));
}

void ShinglingClusterGen::add_special_characters(char* cluster_label, int label_size, char* members, int num_members, int* member_lengths, 
							KeyValue* kv, void* ptr) {
	char marked_label[label_size + 1];
	marked_label[0] = '@';
	memcpy(marked_label + 1, cluster_label, label_size);
	for (int i = 0, offset = 0; i < num_members; offset += member_lengths[i], i++) {
		char marked_member[member_lengths[i] + 1];
		marked_member[0] = '/';
		memcpy(marked_member + 1, members + offset, member_lengths[i]);
		kv -> add(marked_label, label_size + 1, marked_member, member_lengths[i] + 1);
	}
}

void ShinglingClusterGen::remove_finalized_sequences(std::unordered_map<std::string, std::string>* local_sequences, 
							std::string finalized_seq_list_filename) {
	std::ifstream infile;
	//std::cerr << "going to read finalized seq from : " << finalized_seq_list_filename << std::endl;
	infile.open(&(finalized_seq_list_filename)[0], std::ios_base::app);
	if (infile.is_open())
	{
		std::string line;
		while ( infile.good() ) {
			getline (infile,line);
			int removed = local_sequences -> erase(line);
			//if (removed > 0) std::cerr << "removed: " << removed << " " << line << std::endl;
		}
		infile.close();
	} else {  
		std::cerr << "Unable to open file for finalized sequences!"; 
	}
}

std::unordered_map<std::string, int> ShinglingClusterGen::seq_owner_proc;
int ShinglingClusterGen::me = 0;
int ShinglingClusterGen::log_proc = 0 ;
int ShinglingClusterGen::log_scale = 0;
//int ShinglingClusterGen::latest_hashed = 0;
//std::unordered_map<std::string, std::vector<std::vector<unsigned int>>*> ShinglingClusterGen::persistent_hash_function_min_values;
//std::unordered_map<std::string, std::vector<std::vector<unsigned int>>*> ShinglingClusterGen::persistent_vertex_hashes;
//std::unordered_map<std::string, std::unordered_set<std::string>> ShinglingClusterGen::persistent_graph;
//unordered_map<std::string, std::vector<std::vector<unsigned int>>*>* ShinglingClusterGen::persistent_hash_function_min_values = 
//							new unordered_map<std::string, std::vector<std::vector<unsigned int>>*>;
//unordered_map<std::string, std::vector<std::vector<unsigned int>>*>* ShinglingClusterGen::persistent_vertex_hashes = 
//							new unordered_map<std::string, std::vector<std::vector<unsigned int>>*>;
