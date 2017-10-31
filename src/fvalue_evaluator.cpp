#include "fvalue_evaluator.h"

void ClusterEvaluator::compute_nmi(std::vector< std::set<std::string> > test_clusters, std::vector< std::set<std::string> > base_clusters, double* nmi_lfk, double* nmi_max, int print_match_stats_cluster_size_limit) {
	double H_X = 0, H_Y = 0, H_X_cond_Y = 0, H_Y_cond_X = 0;
	compute_entropies(test_clusters, base_clusters, &H_X, &H_X_cond_Y, print_match_stats_cluster_size_limit);
	compute_entropies(base_clusters, test_clusters, &H_Y, &H_Y_cond_X, -1);
	compute_nmi_from_entropies(H_X, H_Y, H_X_cond_Y, H_Y_cond_X, nmi_lfk, nmi_max);
}

//All processors should have the same collection of clusters (vector <unordered_set <string>>)
//Converts to MapReduce clusters and calls the other compute_nmi with MapReduce.
void ClusterEvaluator::compute_nmi(std::vector<std::unordered_set<std::string>>* test_clusters, 
					std::vector<std::unordered_set<std::string>>* base_clusters, double* nmi_lfk, double* nmi_max, 
					int min_cluster_size, int largest_clusters_count, int* largest_test_cluster_size,
					int* largest_base_cluster_size) {
	MapReduce* test_mr;
	MapReduce* base_mr;
	MrToolbox::set_mr_params(test_mr, 4196);
	MrToolbox::set_mr_params(base_mr, 4196);
	int me; MPI_Comm_rank(MPI_COMM_WORLD, &me);
	//if (me == 0) std::cerr << "converting test clusters to mr" << std::endl;
	convert_clusters_to_mr(test_mr, test_clusters, min_cluster_size, largest_clusters_count, largest_test_cluster_size);
	MPI_Barrier(MPI_COMM_WORLD);
	//if (me == 0) std::cerr << "converting base clusters to mr" << std::endl;
	convert_clusters_to_mr(base_mr, base_clusters, min_cluster_size, largest_clusters_count, largest_base_cluster_size);
	MPI_Barrier(MPI_COMM_WORLD);
	//if (me == 0) std::cerr << "computing the fvalue now" << std::endl;
	int root_processor = 0;
	compute_nmi(test_mr, base_mr, root_processor, nmi_lfk, nmi_max);
}

void ClusterEvaluator::convert_clusters_to_mr(MapReduce* mr, std::vector<std::unordered_set<std::string>>* clusters, int min_cluster_size, 
					int largest_clusters_count, int* largest_cluster_size) {
	//pick the largest n clusters
	int me; MPI_Comm_rank(MPI_COMM_WORLD, &me);
	//std::cerr << "largest clusters count " << largest_clusters_count << std::endl;
	if (largest_clusters_count > 0) {
		std::vector< std::unordered_set< std::string >> selected_clusters;
		find_largest_clusters(clusters, largest_clusters_count, &selected_clusters);
		*clusters = selected_clusters;
	}
	//std::cerr << me << " picking largest clusters over " << clusters -> size() << "\n" << std::endl;
	
	//generate_mr
	bool unlabeled = true;
	int size_limit = 1000;
	ClusterReadToolbox clusters_toolbox(clusters, unlabeled, &size_limit);
//	std::cerr << me << " " << clusters_toolbox.clusters -> size() << std::endl;
	MPI_Barrier(MPI_COMM_WORLD);
	//if (me == 1) std::cerr << me << " toolbox created" << std::endl;
	int nmap;
	MPI_Comm_size(MPI_COMM_WORLD, &nmap);
	MPI_Barrier(MPI_COMM_WORLD);
	//std::cerr << me << " calling distribute on " << nmap << " processes" << std::endl;
	MPI_Barrier(MPI_COMM_WORLD);
	mr -> map(32, distribute_clusters, &clusters_toolbox);
	int clusters_count = mr -> convert();
	//if (me == 0) std::cerr << "clusters distributed" << std::endl;
	
	//remove clusters of size < m away
	int local_largest_cluster_size = -1, global_largest_cluster_size = 0;
	ClusterSizeInspectionUtil util(min_cluster_size, &local_largest_cluster_size);
	mr -> reduce(remove_small_clusters, &util);
	mr -> convert();
	//if (me == 0) std::cerr << "small clusters removed" << std::endl;
	MPI_Allreduce(&local_largest_cluster_size, largest_cluster_size, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	//if (me == 0) std::cerr << "largest cluster size computed" << std::endl;
	//*largest_cluster_size = global_largest_cluster_size;
}

void ClusterEvaluator::compute_nmi(MapReduce* test_mr, MapReduce* base_mr, int root_processor, double* nmi_lfk, double* nmi_max, 
					bool rename, int min_cluster_size_limit, int* largest_cluster_size) {
	MapReduce* pair[2];
	//pair[0] = new MapReduce(MPI_COMM_WORLD);
	//pair[1] = new MapReduce(MPI_COMM_WORLD);
	//MrToolbox::set_mr_params(pair[0], 4196);
	//MrToolbox::set_mr_params(pair[1], 4196);
	pair[0] = test_mr -> copy();
	pair[1] = base_mr -> copy(); 
	int local_largest_cluster_size = -1;
	for (int i = 0; i < 2; i++) {
		int current_largest_cluster_size = -1;
		pair[i] -> reduce(rename_sequences, NULL);
		pair[i] -> convert();
		
		MrToolbox::remove_redundant(pair[i]);
		ClusterSizeInspectionUtil util(min_cluster_size_limit, &current_largest_cluster_size);
		pair[i] -> reduce(remove_small_clusters, &util);
		pair[i] -> convert();
		if (i == 0) local_largest_cluster_size = current_largest_cluster_size;
	}
	MPI_Allreduce(&local_largest_cluster_size, largest_cluster_size, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	compute_nmi(pair[0], pair[1], root_processor, nmi_lfk, nmi_max);
	delete pair[0];
	delete pair[1];
}

void ClusterEvaluator::compute_nmi(MapReduce* test_mr, MapReduce* base_mr, int root_processor, double* nmi_lfk, double* nmi_max) {
	//std::cerr << "computing nmi" << std::endl;
	NmiToolbox local_nmi_toolbox(seq_count);
	MapReduce* intersection_mr = new MapReduce(MPI_COMM_WORLD);
	MrToolbox::set_mr_params(intersection_mr, test_mr -> memsize);
	MapReduce* allclusters_mr = new MapReduce(MPI_COMM_WORLD);
	MrToolbox::set_mr_params(allclusters_mr, test_mr -> memsize);
	
	allclusters_mr -> open();
	local_nmi_toolbox.kv = allclusters_mr -> kv;
	int number_of_test_clusters = test_mr -> scan(compute_local_cluster_H_X, &local_nmi_toolbox); 
	allclusters_mr -> close();
	int total_memberships_in_test = 0;
	MPI_Reduce(&(local_nmi_toolbox.member_count_X), &total_memberships_in_test, 1, MPI_INT, MPI_SUM, root_processor, MPI_COMM_WORLD);

	//std::cerr << "local h_x computed " << local_nmi_toolbox.H_X << " " << number_of_test_clusters << std::endl;
	//test_mr -> kmv_stats(2);
	//std::cerr << "and here are the stats for base_mr: " << std::endl;
	//base_mr -> kmv_stats(2);
	//std::cerr << "and finally for intersection_mr: " << std::endl;
	//intersection_mr -> kv_stats(2);
	
	intersection_mr -> open();
	local_nmi_toolbox.kv = intersection_mr -> kv;
	test_mr -> scan(emit_cluster_members_as_key, &local_nmi_toolbox);
//	std::cerr << "hereis h_x : " << local_nmi_toolbox.H_X << std::endl;
	intersection_mr -> close();
	//intersection_mr -> kv_stats(2);
	//std::cerr << "test members intersected" << std::endl;
	
	local_nmi_toolbox.working_on_X = false;
	allclusters_mr -> open(1);
	local_nmi_toolbox.kv = allclusters_mr -> kv;
	int number_of_base_clusters = base_mr -> scan(compute_local_cluster_H_X, &local_nmi_toolbox);
	allclusters_mr -> close();
	int total_memberships_in_base = 0;
	MPI_Reduce(&(local_nmi_toolbox.member_count_Y), &total_memberships_in_base, 1, MPI_INT, MPI_SUM, root_processor, MPI_COMM_WORLD);

	//int num_clusters = cluster_pairs_mr -> gather(1);
	//num_clusters = cluster_pairs_mr -> broadcast(0);
	//cluster_pairs_mr -> convert();
	//int pairs_count = cluster_pairs_mr -> reduce(generate_all_pairs, NULL);
	int total_clusters = allclusters_mr -> convert();
	//std::cerr << "allpairscounted: " << total_clusters << " clusters." << std::endl;
	//std::cerr << "testclusterscount: " << number_of_test_clusters << " baseclusterscount: " << number_of_base_clusters << std::endl;
	
	intersection_mr -> open(1);
	local_nmi_toolbox.kv = intersection_mr -> kv;
	base_mr -> scan(emit_cluster_members_as_key, &local_nmi_toolbox);
	intersection_mr -> close();
	//std::cerr <<"base members intersected" << std::endl;

	intersection_mr -> collate(NULL);
	//intersection_mr -> kmv_stats(2); 
	//cluster_pairs_mr -> open(1);
	//intersection_mr -> reduce(couple_intersecting_clusters, cluster_pairs_mr -> kv); //generates <(label1, label2) , (size1, size2)>
	intersection_mr -> reduce(couple_intersecting_clusters, NULL); //generates <(label1, label2), (size1, size2)>
	//std::cerr << "ccoupled intersecting clusters" << std::endl;
	//intersection_mr -> kv_stats(2);
	//cluster_pairs_mr -> close();
	
	intersection_mr -> collate(NULL);
//	std::cerr << "intersecting clusters coupled and collated" << std::endl;
	allclusters_mr -> open(1);
	intersection_mr -> scan(send_intersection_sizes, allclusters_mr -> kv);
//	std::cerr << "sizes send" << std::endl;
	allclusters_mr -> close();
	
	delete intersection_mr;
	//if (log_scale > 2 && (log_proc == me || log_proc == -1)) std::cerr << "intersection_mr deleted" << std::endl;
	allclusters_mr -> collate(NULL);
	allclusters_mr -> scan(compute_cluster_entropies, &local_nmi_toolbox);
	//cluster_pairs_mr -> reduce(compute_pair_entropy, &seq_count);
	//cluster_pairs_mr -> collate(NULL);
	//cluster_pairs_mr -> scan(compute_single_cluster_entropy, &local_nmi_toolbox);
	delete allclusters_mr;
	int me;
	MPI_Comm_rank(MPI_COMM_WORLD, &me);
	double f1score;
	double f1score_based_on_shingling;
	MPI_Reduce(&(local_nmi_toolbox.H_Y_cond_X), &f1score, 1, MPI_DOUBLE, MPI_SUM, root_processor, MPI_COMM_WORLD);
	MPI_Reduce(&(local_nmi_toolbox.H_X_cond_Y), &f1score_based_on_shingling, 1, MPI_DOUBLE, MPI_SUM, root_processor, MPI_COMM_WORLD);
	//std::cerr << "me: " << me << f1score << std::endl;
	//f1score /= number_of_base_clusters;
	if (me == root_processor) {
	//	std::cerr << " total memberships: in test (shingling):" << total_memberships_in_test << " in base (pfam):" 
	//			<< total_memberships_in_base << std::endl;
	}
	f1score /= total_memberships_in_base;
	//f1score_based_on_shingling /= number_of_test_clusters;
	f1score_based_on_shingling /= total_memberships_in_test;
	double avg_f1 = (f1score + f1score_based_on_shingling) / 2;
//	aggregate_local_metrics(&local_nmi_toolbox, root_processor, me);
	if (me == root_processor) {
//		std::cerr << "before sendngto accum: " << local_nmi_toolbox.H_X << std::endl;
//		compute_nmi_from_entropies(local_nmi_toolbox.H_X, local_nmi_toolbox.H_Y, local_nmi_toolbox.H_X_cond_Y, 
//						local_nmi_toolbox.H_Y_cond_X, nmi_lfk, nmi_max
		*nmi_max = avg_f1;
		*nmi_lfk = f1score;
	}
}

void ClusterEvaluator::compute_entropies(std::vector< std::set<std::string> > X, 
						std::vector< std::set<std::string> > Y, 
						double* H_X, double* H_X_cond_Y, int print_match_stats_cluster_size_limit) {
	compute_H_X(X, H_X);
	compute_H_X_cond_Y(X, Y, H_X_cond_Y, print_match_stats_cluster_size_limit);
}

void ClusterEvaluator::compute_H_X(std::vector< std::set<std::string> > X, double* H_X) {
	*H_X = 0;
	int clusters_count = X.size();
	for (int cluster_index = 0; cluster_index < clusters_count; cluster_index++) {
		int cluster_size = X[cluster_index].size();
		*H_X += compute_h(cluster_size, seq_count);
		*H_X += compute_h(seq_count - cluster_size, seq_count);
	}
}

void ClusterEvaluator::compute_H_X_cond_Y(std::vector< std::set<std::string> > X, std::vector< std::set<std::string> > Y, double* H_X_cond_Y,
						int print_match_stats_cluster_size_limit) {
	int x_clusters_count = X.size();
	int y_clusters_count = Y.size();
	*H_X_cond_Y = 0;
	double H_Xi_cond_Y = 0, H_Xi_cond_Yj = 0, H_star_Xi_cond_Yj = 0;
	int X_clust_size = 0, Y_clust_size = 0;
	for (int x_index = 0; x_index < x_clusters_count; x_index++) {
		/*
		std::cerr << "Checking X[" << x_index <<"]: ";
		for ( auto iter = X[x_index].begin(); iter != X[x_index].end(); iter++) {
                        std::cerr << *iter << "\t";
                }
		std::cerr << std::endl;
		*/
		H_Xi_cond_Y = DBL_MAX;
		double best_H_Xi_cond_Yj = DBL_MAX;
		double best_a, best_b, best_c, best_d, best_Y_size;
		X_clust_size = X[x_index].size();
		std::set<std::string>::iterator x_begin_iter = X[x_index].begin();
		std::set<std::string>::iterator x_end_iter = X[x_index].end();
		for (int y_index = 0; y_index < y_clusters_count; y_index++) {
			/*
			std::cerr << "\twith Y[" << y_index <<"]: ";
                	for ( auto iter = Y[y_index].begin(); iter != Y[y_index].end(); iter++) {
                	        std::cerr << *iter << "\t";
                	}
                	std::cerr << std::endl;
			*/
			Y_clust_size = Y[y_index].size();
			std::set<std::string>::iterator y_begin_iter = Y[y_index].begin();
			std::set<std::string>::iterator y_end_iter = Y[y_index].end();
			std::vector<std::string> set_op_result(std::max(X_clust_size, Y_clust_size));
			std::vector<std::string>::iterator set_op_result_begin_iter = set_op_result.begin();
			//std::vector<std::string>::iterator result_end;	
			int d = std::distance(set_op_result_begin_iter, std::set_intersection(x_begin_iter, x_end_iter, 
												y_begin_iter, y_end_iter, 
												set_op_result_begin_iter)); 
			int c = std::distance(set_op_result_begin_iter, std::set_difference(x_begin_iter, x_end_iter, 
												y_begin_iter, y_end_iter, 
												set_op_result_begin_iter));
			int b = std::distance(set_op_result_begin_iter, std::set_difference(y_begin_iter, y_end_iter, 
												x_begin_iter, x_end_iter,
												set_op_result_begin_iter));
			int a = seq_count - b - c - d;
			/*
			std::cerr << "\t\ta=" << a << "\tb=" << b << "\tc=" << c << "\td=" << d << std::endl;
			*/
			H_star_Xi_cond_Yj = compute_H_star_Xi_cond_Yj(a, b, c, d, seq_count, &H_Xi_cond_Yj);
			//std::cerr << "paired_cluster 1 (index): " << x_index << " 2: " << y_index << " from the pair: " << X_clust_size << " " << Y_clust_size << " " << d << " " << c << " " << b << " " << a << " " << H_star_Xi_cond_Yj << std::endl;
			if (H_Xi_cond_Yj < best_H_Xi_cond_Yj) {
				best_H_Xi_cond_Yj = H_Xi_cond_Yj;
				best_a = a;
				best_b = b;
				best_c = c;
				best_d = d;
				best_Y_size = Y_clust_size;
			}
			H_Xi_cond_Y = H_star_Xi_cond_Yj < H_Xi_cond_Y ? H_star_Xi_cond_Yj : H_Xi_cond_Y;
			//if (X_clust_size == 50) {
			//	std::cerr << "candidate: " << Y_clust_size << " " << H_Xi_cond_Yj << " a=" << a << " b=" << b << " c=" << c << " d=" << d << std::endl;
			//}
			/*
			if (H_star_Xi_cond_Yj < H_Xi_cond_Y) std::cerr << "\t\tupdating H_Xi_cond_Y to " << H_Xi_cond_Y << std::endl;
			*/
		}
		/*
		std::cerr << "\tupdating H_X_cond_Y by " << H_Xi_cond_Y << std::endl;
		*/
		//if (print_match_stats_cluster_size_limit > 0 && X_clust_size > print_match_stats_cluster_size_limit) {
		//	std::cerr << "best match found for cluster of size: " << X_clust_size << " - match_size=" << best_Y_size <<
		//			" a=" << best_a << " b=" << best_b << " c=" << best_c << " d=" << best_d << std::endl;
		//}
		*H_X_cond_Y += H_Xi_cond_Y;
	}
}

double ClusterEvaluator::compute_h(int w, int n) {
	double result = w > 0 ? ( -w * log2((double) w / n) ) : 0;
	return result;
}

void ClusterEvaluator::compute_local_cluster_H_X(char* label, int label_length, char* members, int num_members, 
								int* member_lengths, void* nmi_toolbox_ptr) {
	NmiToolbox* nmi_toolbox = (NmiToolbox*) nmi_toolbox_ptr;
//	double local_H = compute_h(num_members, nmi_toolbox -> seq_count) +
//					compute_h(nmi_toolbox -> seq_count - num_members, nmi_toolbox -> seq_count);
	if (nmi_toolbox -> working_on_X) {
//		nmi_toolbox -> H_X += local_H;
		LabelAndSizeAsValue single_label(num_members, label);
		single_label.set_label_symbol(X_CLUSTER_LABEL);
		//nmi_toolbox -> kv -> add(NULL, 0, single_label.get_message(), single_label.get_message_length());
		nmi_toolbox -> kv -> add(single_label.get_message(), single_label.get_message_length(), NULL, 0);
		nmi_toolbox -> member_count_X += num_members;
	} else {
//		nmi_toolbox -> H_Y += local_H;
		LabelAndSizeAsValue single_label(num_members, label);
		//if (strcmp(label, "@0") == 0) std::cerr << " forcluster: " << label << " size: " << num_members << std::endl;
		single_label.set_label_symbol(Y_CLUSTER_LABEL);
		//nmi_toolbox -> kv -> add(NULL, 0, single_label.get_message(), single_label.get_message_length());
		nmi_toolbox -> kv -> add(single_label.get_message(), single_label.get_message_length(), NULL, 0);
		nmi_toolbox -> member_count_Y += num_members;
	}
}

double ClusterEvaluator::compute_H_star_Xi_cond_Yj(int a, int b, int c, int d, int seq_count, double* H_Xi_cond_Yj_ptr = NULL) {
	double H_star_Xi_cond_Yj = 0;
	double h_a = compute_h(a, seq_count);
	double h_b = compute_h(b, seq_count);
	double h_c = compute_h(c, seq_count);
	double h_d = compute_h(d, seq_count);
	double H_Xi_cond_Yj = 0;
	H_Xi_cond_Yj = h_a + h_b + h_c + h_d - compute_h( b + d, seq_count) - compute_h( a + c, seq_count);
	if (H_Xi_cond_Yj < DBL_EPSILON) {
		H_Xi_cond_Yj = 0;
	}
/*	if ( a == 94597 ) {
		std::cerr << "this is H_Xi_cond_Yj: " << H_Xi_cond_Yj << std::endl;
		std::cerr << "a= " << a << "h_a=" << h_a << '\n';
		std::cerr << "b= " << b << "h_b=" << h_b << '\n';
		std::cerr << "c= " << c << "h_c=" << h_c << '\n';
		std::cerr << "d= " << d << "h_d=" << h_d << '\n';
		std::cerr << "b+d gives " << compute_h( b + d, seq_count) << " and a+c gives " << compute_h( a + c, seq_count) << std::endl;
		if (h_a + h_d >= h_b + h_c) 
			std::cerr << "first condition" << std::endl;
		else {
			std::cerr << "otherwise: " << compute_h(c + d, seq_count) << " and " << compute_h(a + b, seq_count) << std::endl;
		}	
	}
*/	if (h_a + h_d >= h_b + h_c) {
		/*
		std::cerr << "\t\t\tcondition held; H_star_Xi_cond_Yj = H_Xi_cond_Yj = " << H_Xi_cond_Yj << std::endl;
		*/
		H_star_Xi_cond_Yj = H_Xi_cond_Yj;
	} else {
		H_star_Xi_cond_Yj = compute_h(c + d, seq_count) + compute_h(a + b, seq_count);
		/*
		std::cerr << "\t\t\tcondition not held; H_star_Xi_cond_Yj = H_Xi = " << H_star_Xi_cond_Yj << std::endl;
		*/
	}
	if(H_Xi_cond_Yj_ptr != NULL) {
		*H_Xi_cond_Yj_ptr = H_Xi_cond_Yj;
	}
/*	if (H_star_Xi_cond_Yj < -1) {
		std::cerr << "here it is negative: " << a << " " << b << " " << c << " " << d << " " << seq_count << " " << H_star_Xi_cond_Yj << std::endl;
	}
	std::cerr << "armprec: " << H_star_Xi_cond_Yj << std::endl;
*/	return H_star_Xi_cond_Yj;
}

void ClusterEvaluator::emit_cluster_members_as_key(char* label, int label_length, char* members, int num_members,
								int* member_lengths, void* nmi_toolbox_ptr) {
	NmiToolbox* nmi_toolbox = (NmiToolbox*) nmi_toolbox_ptr;
	//std::cerr << "extracting members" << std::endl;
	LabelAndSizeAsValue message(num_members, label);
	if (nmi_toolbox -> working_on_X) {
		//std::cerr << "Setting label to X_label " << X_CLUSTER_LABEL << std::endl;
		message.set_label_symbol(X_CLUSTER_LABEL);
	} else {
		//std::cerr << "Setting label to X_label" << Y_CLUSTER_LABEL << std::endl;
		message.set_label_symbol(Y_CLUSTER_LABEL);
	}
	char* label_w_cluster_size = message.get_message();
	int value_length = message.get_message_length();
	//std::cerr<< "for label: " << std::endl;
	for (int index = 0, offset = 0; index < num_members; offset += member_lengths[index], index++) {
		//std::cerr << "\temitting cluster member as key: " << &members[offset] << " " << *(int*)label_w_cluster_size << " " << label_w_cluster_size + sizeof(int) << " keysize: " << member_lengths[index] << " valuelength: " << value_length << std::endl;	
		nmi_toolbox -> kv -> add(&members[offset], member_lengths[index], label_w_cluster_size, value_length);
	}
}

void ClusterEvaluator::generate_all_pairs(char* null_key, int keybytes, char* all_labels, int clusters_count, 
								int* message_lengths, KeyValue* kv, void* ptr) {
	//int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	//if (rank == 0) { std::cerr << "generate_all_pairs is called " << clusters_count << std::endl;
	//for (int i = 0, offset1 = 0; i < clusters_count; offset1 += message_lengths[i], i++) {
	//	std::cerr << all_labels + offset1 + sizeof(int) << " ";
	//}
	//if (clusters_count > 0) std::cerr << std::endl;}
	//MPI_Barrier(MPI_COMM_WORLD);
	int me, nprocs;
	char default_first_set_symbol = all_labels[sizeof(int)];
	int default_first_set_count = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &me);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	for (int i = 0, offset1 = 0; i < clusters_count; offset1 += message_lengths[i], i++) {
		char current_label_symbol = all_labels[offset1 + sizeof(int)];
		if (default_first_set_count % nprocs == me && current_label_symbol == default_first_set_symbol) {
			LabelAndSizeAsValue label1(all_labels + offset1);
			for (int j = 1, offset2 = message_lengths[0]; j < clusters_count; offset2 += message_lengths[j], j++) {
				LabelAndSizeAsValue label2(all_labels + offset2);
				CoupledLabels coupled_labels(label1.get_label_name(), label2.get_label_name(), label1.get_size(), label2.get_size());
				if (!coupled_labels.are_from_same_side()) {
					kv -> add(coupled_labels.get_labels_message(), coupled_labels.get_labels_message_length(),
								coupled_labels.get_sizes_message(), coupled_labels.get_sizes_message_length());
					//std::cerr << "these clusters are coupled: " << label1.get_label_name() << " " << label2.get_label_name() << std::endl;
				}
			}
		}
		if (current_label_symbol == default_first_set_symbol) {
			default_first_set_count++;
		}
	}
}

void ClusterEvaluator::couple_intersecting_clusters(char* common_member, int member_length, char* labels, int num_labels, 
								int* label_lengths, KeyValue* kv, void* couples_mr_kv_ptr) {
	//KeyValue* couples_kv = (KeyValue*) couples_mr_kv_ptr;
	char *label1, *label2;
	int label1_members_count, label2_members_count;
	//std::cerr << "received common member: " << common_member << " shared between " << num_labels << " labels" << std::endl;
	for (int i = 0, offset1 = 0 ; i < (num_labels - 1); i++) {
		//std::cerr << "at offset1 " << offset1 << std::endl;
		//std::cerr << "\tlabel1: " << labels + offset1 + sizeof(int) << std::endl;
		LabelAndSizeAsValue value1(labels + offset1);
		label1 = value1.get_label_name();
		//std::cerr << "\t\tshould be equal to: " << label1 << std::endl;
		label1_members_count = value1.get_size();
		for (int j = i + 1, offset2 = offset1 + label_lengths[i]; j < num_labels; offset2 += label_lengths[j], j++) {
			LabelAndSizeAsValue value2(labels + offset2);
			label2 = value2.get_label_name();
			//std::cerr << "\t\tshould be equal to: " << label2 << std::endl;
			label2_members_count = value2.get_size();
			CoupledLabels coupled_labels(label1, label2, label1_members_count, label2_members_count);
			if (!coupled_labels.are_from_same_side()) {
				//couples_kv -> add(coupled_labels.get_labels_message(), coupled_labels.get_labels_message_length(),
				//		coupled_labels.get_sizes_message(), coupled_labels.get_sizes_message_length());
				kv -> add(coupled_labels.get_labels_message(), coupled_labels.get_labels_message_length(),
						coupled_labels.get_sizes_message(), coupled_labels.get_sizes_message_length());
				//std::cerr << "these clusters are coupled: " << coupled_labels.get_labels_message() << std::endl;
			}
			//std::cerr << "going up" << std::endl;
		}
		//std::cerr << "updating off1" << std::endl;
		offset1 += label_lengths[i];
		//std::cerr << "going for next i" << std::endl;
	}
	//std::cerr << "returning from couple intersecting clusters " << std::endl;
}

void ClusterEvaluator::compute_pair_entropy(char* paired_labels, int key_length, char* paired_cluster_sizes, int intersection_size_plus_one,
								int* valuebytes, KeyValue* kv, void* seq_count_ptr) {
	//std::cerr << " in computepairentropy" << paired_labels << std::endl; 
	//std::cerr << " " << paired_labels + strlen(paired_labels) << std::endl;
	int intersection_size = intersection_size_plus_one - 1;
	int seq_count = *(int*) seq_count_ptr;
	CoupledLabels coupled_labels(paired_labels, paired_cluster_sizes);
	char* label_x = coupled_labels.get_label(X_CLUSTER_LABEL);
	char* label_y = coupled_labels.get_label(Y_CLUSTER_LABEL);
	//std::cerr << "paired_cluster 1: " << label_x << " 2: " << label_y << std::endl;
	int x_cluster_size = coupled_labels.get_cluster_size(X_CLUSTER_LABEL);
	int y_cluster_size = coupled_labels.get_cluster_size(Y_CLUSTER_LABEL);
	int x_diff_y = x_cluster_size - intersection_size;
	int y_diff_x = y_cluster_size - intersection_size;
	int a = seq_count - intersection_size - x_diff_y - y_diff_x;
	double precision = (double) intersection_size / x_cluster_size;
	double recall = (double) intersection_size / y_cluster_size;
	double f1score = 2 * precision * recall / (precision + recall);
	//if (strcmp (label_x, "X1005_1") == 0) {
	//	std::cerr << "label_x: " << label_x << " " << x_cluster_size << " " << intersection_size << " " << y_cluster_size << " " << y_diff_x << " " << x_diff_y << " " << f1score << std::endl;
	//}
//	double H_star_Xi_cond_Yj = compute_H_star_Xi_cond_Yj(a, y_diff_x, x_diff_y, intersection_size, seq_count);
//      double H_star_Yi_cond_Xj = compute_H_star_Xi_cond_Yj(a, x_diff_y, y_diff_x, intersection_size, seq_count);
/*	if (y_cluster_size > 6000) std::cerr << "label_y: " << label_y << " " << y_cluster_size << " " << intersection_size << 
						" label_x:" << label_x << " " << x_cluster_size << " " << y_diff_x << " " << x_diff_y << " " 
						<< f1score << " " << recall << " " << precision << std::endl;
	if (x_cluster_size > 10000) std::cerr << "label_x: " << label_x << " " << x_cluster_size << " " << intersection_size <<  
						" label_y:" << label_y << " " << y_cluster_size << " " << y_diff_x << " " << x_diff_y << " " 
						<< f1score << " " << recall << " " << precision << std::endl;
*/	//if (x_cluster_size == largest_shingling_cluster_size) std::cerr << "label_x: " << label_x << " " << x_cluster_size << " " << intersection_size << " " << y_cluster_size << " " << y_diff_x << " " << x_diff_y << std::endl;
	//if (y_cluster_size == largest_pfam_cluster_size) std::cerr << "label_y: " << label_y << " " << y_cluster_size << " " << intersection_size << " " << x_cluster_size << " " << y_diff_x << " " << x_diff_y << std::endl;
	//std::cerr << "paired_cluster 1: " << label_x << " 2: " << label_y << " from the pair: " << x_cluster_size << " " << y_cluster_size << " " << intersection_size << " " << x_diff_y << " " << y_diff_x << " " << a << std::endl;
	//std::cerr << "paired_cluster 1: " << label_x << " 2: " << label_y << " from the pair: " << x_cluster_size << " " << y_cluster_size << " " << intersection_size << " " << x_diff_y << " " << y_diff_x << " " << a << " " << H_star_Xi_cond_Yj << " " << H_star_Yi_cond_Xj << std::endl;
/*	if (x_cluster_size < 6 && x_cluster_size > 2 && y_cluster_size > 2 && y_cluster_size < 8 && intersection_size > 1 && intersection_size < x_cluster_size && intersection_size < y_cluster_size) {
		std::cerr << " for " << label_x << " and " << label_y << ", values are " << intersection_size << " " << x_diff_y << " " << y_diff_x << " " << a << ". and: " << H_star_Xi_cond_Yj << " " << H_star_Yi_cond_Xj << std::endl;
	}
	if (strcmp(label_x, "X6156_1") == 0) {
		std::cerr << " herecheckthisonefor " << label_x << " and " << label_y << ", values are " << intersection_size << " " << x_diff_y << " " << y_diff_x << " " << a << ". and: " << H_star_Xi_cond_Yj << " " << H_star_Yi_cond_Xj << std::endl;
	}
*/	//std::cerr << "HstarXicondy: " << H_star_Xi_cond_Yj << " " << H_star_Yi_cond_Xj << std::endl;
	
//	kv -> add(label_x, strlen(label_x) + 1, (char*) &H_star_Xi_cond_Yj, sizeof(double));
//	kv -> add(label_y, strlen(label_y) + 1, (char*) &H_star_Yi_cond_Xj, sizeof(double)); 
	int combined_value_length = sizeof(double) + sizeof(int);
	char value_x[combined_value_length], value_y[combined_value_length];
	memcpy(value_x, &f1score, sizeof(double));
	memcpy(value_x + sizeof(double), &x_cluster_size, sizeof(int));
	memcpy(value_y, &f1score, sizeof(double));
	memcpy(value_y + sizeof(double), &y_cluster_size, sizeof(int));
	kv -> add(label_x, strlen(label_x) + 1, &value_x[0], combined_value_length);
	kv -> add(label_y, strlen(label_y) + 1, &value_y[0], combined_value_length);
	//kv -> add(label_y, strlen(label_y) + 1, (char*) &f1score, sizeof(double));
	//kv -> add(label_x, strlen(label_x) + 1, (char*) &f1score, sizeof(double));
	//std::cerr << "returning from compute_pair_entropy" << std::endl;
}

void ClusterEvaluator::compute_single_cluster_entropy(char* cluster_label, int label_length, char* pairwise_entropies, 
					int paired_clusters_count, int* size_of_double, void* nmi_toolbox_ptr) {
	NmiToolbox* nmi_toolbox = (NmiToolbox*) nmi_toolbox_ptr;
//	std::cerr << "computingsinglecluster: " << cluster_label << " " << paired_clusters_count << std::endl;
	char cluster_label_symbol = cluster_label[0];
/*	double H_Xi_cond_Y = DBL_MAX;
	double H_star_Xi_cond_Yj = 0;
	for (int j = 0, offset = 0; j < paired_clusters_count; offset += size_of_double[j], j++) {
		H_star_Xi_cond_Yj = *(double*) (pairwise_entropies + offset);
//		if (strcmp(cluster_label, "X6156_1") == 0) { std::cerr << "seewhatigot: " << cluster_label << " " << H_star_Xi_cond_Yj << std::endl;}
//		std::cerr << " received H_star for: " << cluster_label << ": " <<  H_star_Xi_cond_Yj << std::endl;
		H_Xi_cond_Y = H_star_Xi_cond_Yj < H_Xi_cond_Y ? H_star_Xi_cond_Yj : H_Xi_cond_Y;
	}
//	std::cerr << "computedsingleclusterentropy for label: " << cluster_label << " and symbol: " << cluster_label_symbol << "= " << H_Xi_cond_Y << std::endl;
//	std::cerr << "armselected: " << H_Xi_cond_Y << std::endl;
*/
	//double f1score = DBL_MIN;
	double f1score = 0;
	int best_cluster_size = 0;
	double current = 0;
	for (int j = 0, offset = 0; j < paired_clusters_count; offset += size_of_double[j], j++) {
		current = *(double*) (pairwise_entropies + offset);
		if (f1score < current) {
			f1score = current;
			best_cluster_size = *(int*)(pairwise_entropies + offset + sizeof(double));
		}
		//f1score = f1score < current ? current : f1score;
	}
	double unchanged_fscore = f1score;
	if (f1score > 0) f1score = f1score * best_cluster_size;
	if (cluster_label_symbol == X_CLUSTER_LABEL) {
		//if (best_cluster_size > 2000) {
//			std::cerr << "testcluster\t" << cluster_label << "\t" << best_cluster_size << "\t" << f1score << "\t" << unchanged_fscore << std::endl;
		//}
//		nmi_toolbox -> H_X_cond_Y += H_Xi_cond_Y;
		nmi_toolbox -> H_X_cond_Y += f1score;
	} else {
//		std::cerr << "basecluster\t" << cluster_label << "\t" << best_cluster_size << "\t" << f1score << "\t" << unchanged_fscore << std::endl;
//		nmi_toolbox -> H_Y_cond_X += H_Xi_cond_Y;
		nmi_toolbox -> H_Y_cond_X += f1score;
	}
	//if (strcmp (cluster_label, "X1005_1") == 0) std::cerr << cluster_label << " has max: " << f1score << std::endl;
	//kv -> add(cluster_label, label_length, (char*) &paired_clusters_count, sizeof(int));
	//kv -> add(&cluster_label_symbol, 1, (char*) &H_Xi_cond_Y, sizeof(double));
}

void ClusterEvaluator::compute_clusterings_conditional_entropy(char* clustering_symbol, int one, char* single_cluster_entropies, 
						int num_clusters, int* size_of_double, KeyValue* kv, void* root_ptr) {
	double H_X_cond_Y = 0;
	double H_Xi_cond_Y = 0;
	for (int i = 0, offset = 0; i < num_clusters; offset += size_of_double[i], i++) {
		H_Xi_cond_Y = *(double*)(single_cluster_entropies + offset);
		H_X_cond_Y += H_Xi_cond_Y;
	}
		
}

void ClusterEvaluator::aggregate_local_metrics(NmiToolbox* nmi_toolbox, int root, int me) {
	//std::cerr << "now aggregating the results" << std::endl;
	double local_nmi_metrics[4];
	local_nmi_metrics[0] = nmi_toolbox -> H_X;
	local_nmi_metrics[1] = nmi_toolbox -> H_Y;
	local_nmi_metrics[2] = nmi_toolbox -> H_X_cond_Y;
	local_nmi_metrics[3] = nmi_toolbox -> H_Y_cond_X;
	int nprocs;
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	double* rcvd_metrics = NULL;
	int numbers_to_rcv = nprocs * 4;
	if (me == root) {
		rcvd_metrics = new double[numbers_to_rcv];
	}
	MPI_Gather(local_nmi_metrics, 4, MPI_DOUBLE, rcvd_metrics, 4, MPI_DOUBLE, root, MPI_COMM_WORLD);
	if (me == root) {
/*		for (int i = 0; i < numbers_to_rcv; i++) {
			std::cerr << "rcvdval: " << i << ": " << rcvd_metrics[i] << std::endl;
		}
*/		nmi_toolbox -> H_X = 0;
		nmi_toolbox -> H_Y = 0;
		nmi_toolbox -> H_X_cond_Y = 0;
		nmi_toolbox -> H_Y_cond_X = 0;
		for (int i = 0; i < numbers_to_rcv; i += 4) {
//			std::cerr << " adding: " << rcvd_metrics[i] << " to " << nmi_toolbox -> H_X << std::endl;
			nmi_toolbox -> H_X += rcvd_metrics[i];
			nmi_toolbox -> H_Y += rcvd_metrics[i + 1];
			nmi_toolbox -> H_X_cond_Y += rcvd_metrics[i + 2];
			nmi_toolbox -> H_Y_cond_X += rcvd_metrics[i + 3];
		}
//		std::cerr << "in func " << nmi_toolbox -> H_X << std::endl;
		delete[] rcvd_metrics;
		rcvd_metrics = NULL;
	}
}

void ClusterEvaluator::compute_nmi_from_entropies(double H_X, double H_Y, double H_X_cond_Y, double H_Y_cond_X, 
							double* nmi_lfk, double* nmi_max) {
	//std::cerr << "summing up the metrics" << std::endl;
	double I_X_Y = 0.5 * (H_X - H_X_cond_Y + H_Y - H_Y_cond_X);
	//std::cerr << "H_X=" << H_X << "\tH_X_cond_Y=" << H_X_cond_Y << "\tH_Y=" << H_Y << "\tH_Y_cond_X=" << H_Y_cond_X << std::endl;
	*nmi_max = I_X_Y / fmax(H_X, H_Y);
	*nmi_lfk = 1 - 0.5 * ( H_X_cond_Y / H_X + H_Y_cond_X / H_Y );
//	std::cerr << "I_X_Y=" << I_X_Y << "\tnmi_lfk=" << *nmi_lfk << "\tnmi_max=" << *nmi_max << std::endl;
}

void ClusterEvaluator::rename_sequences(char* cluster_label, int label_length, char* seq_names, int num_members,
						int* seq_name_lengths, KeyValue* kv, void* ptr) {
	for (int index = 0, offset = 0; index < num_members; offset += seq_name_lengths[index], index++) {
		char* domained_seq_name = strtok(seq_names + offset + 1, "_");
		kv -> add(cluster_label, label_length, domained_seq_name, strlen(domained_seq_name) + 1);
	}
}

void ClusterEvaluator::remove_small_clusters(char* cluster_label, int label_length, char* member_names, int num_members,
						int* member_name_lengths, KeyValue* kv, void* util_ptr) {
	ClusterSizeInspectionUtil* util = (ClusterSizeInspectionUtil*) util_ptr;
	if (*(util -> largest_cluster_size) < num_members){
		*(util -> largest_cluster_size) = num_members;
	}
	int min_cluster_size_limit = util -> min_allowed_size;
	if (num_members >= min_cluster_size_limit) {
		for (int index = 0, offset = 0; index < num_members; offset += member_name_lengths[index], index++) {
			kv -> add(cluster_label, label_length, member_names + offset, member_name_lengths[index]);
		}
	}
}

void ClusterEvaluator::read_pfam_clusters(std::string filename, int* largest_cluster_size, int min_cluster_size_limit, 
						MapReduce* pfam_mr, int nmap, bool wholeseq) {
	//wholeseq is for overlapping clusters, when we disregard the domains; i.e. rename=true for shingling clusters
	std::vector< std::unordered_set<std::string> > pfam_clusters = read_pfam_clusters_on_all_proc(filename, largest_cluster_size, min_cluster_size_limit, wholeseq);
	//std::cout << "pfam clusters read:" << pfam_clusters.size() << std::endl;
	largest_pfam_cluster_size = *largest_cluster_size;
	bool unlabeled = true;
	int size_limit = 6000;
	ClusterReadToolbox clusters(&pfam_clusters, unlabeled, &size_limit);
	//std::cout << "base first cluster first member: " << *(pfam_clusters[0].begin()) << std::endl;
	pfam_mr -> map(nmap, distribute_clusters, &clusters);
	int pfam_clusters_count = pfam_mr -> convert();
	int me; MPI_Comm_rank(MPI_COMM_WORLD, &me);
	//if (me == 0) std::cerr << "pfam clusters count: " << pfam_clusters_count << std::endl;
}

std::vector< std::unordered_set<std::string> > ClusterEvaluator::read_pfam_clusters_on_all_proc(std::string filename, 
							int* largest_cluster_size, int min_cluster_size_limit, bool wholeseq) {
	*largest_cluster_size = 0;
	std::unordered_map< std::string, std::unordered_set<std::string> > named_clusters;
	std::unordered_map< std::string, int> domain_counter;
	std::ifstream pfam_file(filename);
	std::string line;
	while (std::getline(pfam_file, line)) {
		std::vector<std::string> tokens = split_string(line, '\t');
		std::string cluster_name = tokens[5];
		std::string seq_name = tokens[0];
		if (named_clusters.find(cluster_name) == named_clusters.end()) {
			std::unordered_set<std::string> empty_cluster;
			named_clusters[cluster_name] = empty_cluster;
		}
		if (domain_counter.find(seq_name) == domain_counter.end()) {
			domain_counter.insert(std::make_pair(seq_name, 0));
		}
		domain_counter[seq_name] = domain_counter[seq_name] + 1;
		if (wholeseq) {
			named_clusters[cluster_name].insert(seq_name);
		} else {
			named_clusters[cluster_name].insert("/" + seq_name + "_" + std::to_string(domain_counter[seq_name]));
			//named_clusters[cluster_name].insert("/" + seq_name);
		}
	}
	std::vector< std::unordered_set<std::string> > clusters;
	clusters.reserve(named_clusters.size());
	for (auto iter = named_clusters.begin(); iter != named_clusters.end(); iter++) {
		std::unordered_set<std::string> cluster = iter -> second;
		int cluster_size = cluster.size();
		if (cluster_size >= min_cluster_size_limit) {
			//std::cerr << "clname\t" << iter -> first << "\t" << clusters.size() << std::endl;
			*largest_cluster_size = cluster_size > *largest_cluster_size ? cluster_size : *largest_cluster_size;
			clusters.push_back(cluster);
		}
	}
	return clusters;
}

std::vector< std::string > ClusterEvaluator::split_string(std::string str, char delim) {
	std::vector< std::string > tokens;
	std::stringstream ss(str);
	std::string token;
	while (std::getline(ss, token, delim)) {
		tokens.push_back(token);
	}
	return tokens;
}

void ClusterEvaluator::distribute_clusters(int itask, KeyValue* kv, void* clusters_ptr) {
	int me, nprocs;
	MPI_Comm_rank(MPI_COMM_WORLD, &me);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	//std::cerr << me << " in distribute" << std::endl;
	MPI_Barrier(MPI_COMM_WORLD);
	ClusterReadToolbox* clusters = (ClusterReadToolbox*) clusters_ptr;
	//std::cerr << "distributing " << clusters -> clusters -> size() << " clusters" << std::endl;
	bool first_printed = false;
	int clusters_count = clusters -> clusters -> size();
	//std::cerr << "all clusters count: " << clusters_count << std::endl;
	int proc_share = ceil( (double) clusters_count / nprocs );
	int start_index = me * proc_share;
	int end_index = std::min((me + 1) * proc_share, clusters_count);
	int limit = *(clusters -> size_limit);
	//std::cerr << "proc" << me << " start: " << start_index << " end: " << end_index << " share: " << proc_share << " clusters_count: " << clusters_count << std::endl;
	for (int i = start_index; i < end_index; i++) {
		std::string label_str;
		if (clusters -> unlabeled == false) {
			//std::cerr << "clusters are labeled" << std::endl;
			label_str = clusters -> cluster_labels -> at(i);
		} else {
			//std::cerr << "clusters are unlabeled" << std::endl;
			label_str = '@' + std::to_string(i);
		}
		int label_length = label_str.length() + 1;
		char label[label_length];
		strcpy(label, label_str.c_str());
		bool write_to_file = false;
		//bool write_to_file = (label_str == "@204" || label_str == "@167" || label_str == "@165" || label_str == "@133" || 
		//	label_str == "@151" || label_str == "@152");
		std::string line = "";
		if (write_to_file) {
			line += "KMV pair: proc oo, nvalues oo, sizes oo oo, key " + label_str +", values ";
		}
//		if (clusters -> clusters -> at(i).size() > limit) {
//			std::cerr << "large cluster found: " << label << std::endl;
//		}
		//if (clusters -> unlabeled) {
		//	std::cerr << "cluster" << "\t" << i << "\t" << clusters -> clusters -> at(i).size() << std::endl;
		//}
		for (auto iter = clusters -> clusters -> at(i).begin(); iter != clusters -> clusters -> at(i).end(); iter++) {
			std::string member_str = *iter;
			if (write_to_file) {
				line += "/" + member_str + " ";
			}
			int member_length = member_str.length() + 1;
			char member[member_length];
			strcpy(member, member_str.c_str());
			//std::cerr << "adding kv now: " << label << " " << member << " by: " << me << std::endl;
			kv -> add(label, label_length, member, member_length);
			if (!first_printed && me == 0) {
			//	std::cerr << "sending: " << label << ", " << label_length << " and " << member << ", " << member_length << std::endl;
				first_printed = true;
			}
		}
		if (write_to_file) {
			std::string filename = "datasets/large70/" + label_str.substr(1);
			std::ofstream myfile (filename);
			myfile << line;
			myfile.close();
			line = "";
		}
	}
//	std::cerr << me << " distributed" << std::endl;
}

void ClusterEvaluator::read_shingling_clusters(std::string filename, int* largest_cluster_size, bool rename, int min_cluster_size, 
						MapReduce* shingling_mr, int nmap, int largest_clusters_count, std::string clustering_type) {
	std::vector<std::string> cluster_labels;
	//std::cerr << "reading clusters from: " << filename << std::endl;
	std::vector< std::unordered_set<std::string> > shingling_clusters = read_shingling_clusters_on_all_proc(filename, largest_cluster_size, rename, 
													min_cluster_size, &cluster_labels);
	//std::cerr << "clusters read: " << shingling_clusters.size() << std::endl;
	//std::cout << "largest clusters count: " << largest_clusters_count << std::endl;
	std::vector<int> selected_clusters_indices;
	if (largest_clusters_count > 0) {
		//std::cout << "finding largest clusters" << std::endl;
		find_largest_clusters(&shingling_clusters, &selected_clusters_indices, largest_clusters_count);
		std::vector< std::unordered_set<std::string> > temp_selected_clusters;
		std::vector< std::string > temp_selected_labels;
		for (int i = 0; i < selected_clusters_indices.size(); i++) {
			temp_selected_clusters.push_back(shingling_clusters[selected_clusters_indices[i]]);
			temp_selected_labels.push_back(cluster_labels[selected_clusters_indices[i]]);
		}
		//std::cout << "clisters found" << std::endl;
		shingling_clusters = temp_selected_clusters;
		cluster_labels = temp_selected_labels;
		int me = 0; MPI_Comm_rank(MPI_COMM_WORLD, &me);
//		if (me == 0) {
//			std::cerr << "selected cluster sizes: ";
//			for (int i = 0; i < shingling_clusters.size(); i++) {
//				std::cerr << shingling_clusters[i].size() << "\t" << cluster_labels[i] << "\t";
//			}
//			std::cerr << std::endl;
//		}
	}
	bool unlabeled;
	int size_limit = 1000;
	if (clustering_type == "as_test") {
		unlabeled = false;
		size_limit = 1000;
	}	
	else if (clustering_type == "as_base") {
		unlabeled = true;
		size_limit = 1000;
	}
	//std::cout << " and now... " << std::endl;
	//std::cout << "largest_cluster_size: " << *largest_cluster_size << std::endl;
	largest_shingling_cluster_size = *largest_cluster_size;
	//std::cout << "shingling first clusters first member: " << *(shingling_clusters[0].begin()) << std::endl;
	//std::cout << "going to distribute" << std::endl;
	ClusterReadToolbox clusters(&shingling_clusters, unlabeled, &size_limit, &cluster_labels);
	//std::cout << "calling map" << std::endl;
	shingling_mr -> map(nmap, distribute_clusters, &clusters);
	//std::cout << "map done" << std::endl;
	int shingling_clusters_count = shingling_mr -> convert();
	//std::cout << "converted" << std::endl;
	int me; MPI_Comm_rank(MPI_COMM_WORLD, &me);
//	if (me == 0)
//		std::cerr << "shingling clusters count: " << shingling_clusters_count << std::endl;
}

std::vector< std::unordered_set<std::string> > ClusterEvaluator::read_shingling_clusters_on_all_proc(std::string filename, 
								int* largest_cluster_size, bool rename, int min_cluster_limit,
								std::vector<std::string>* cluster_labels) {
	int me;
	MPI_Comm_rank(MPI_COMM_WORLD, &me);
	*largest_cluster_size = 0;
	//if (me == 0) std::cerr << "reading file: " << filename << std::endl;
	std::ifstream infile(filename);
	std::string line;
	std::vector< std::unordered_set<std::string> > clusters;
	std::unordered_set<std::string> labels_set;
	int index_counter = 0;
	while (std::getline(infile, line)) {
		std::unordered_set<std::string> cluster;
		std::vector< std::string > numbered_members = split_string(line, ' ');
		int members_size = numbered_members.size();
		std::string cluster_label = numbered_members[10].substr(0, numbered_members[10].length() - 1);
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
		//if (cluster_size == 5923) {
			//std::cerr << "this is the largest cluster: \n" << cluster_label << std::endl;
		//}
		if (cluster_size >= min_cluster_limit) {
			if (labels_set.find(cluster_label) == labels_set.end()) {
				labels_set.insert(cluster_label);
				*largest_cluster_size = cluster_size > *largest_cluster_size ? cluster_size : *largest_cluster_size;
				clusters.push_back(cluster);
				cluster_labels -> push_back(cluster_label);
				//std::cerr << "cluster label: " << numbered_members[10] << " in index " << index_counter << std::endl;
				index_counter++;
			} else {
//				std::cerr << "label already inserted: " << cluster_label << std::endl;
			}
		}
	}
	return clusters;
}

void ClusterEvaluator::find_largest_clusters(std::vector<std::unordered_set<std::string>>* shingling_clusters, 
						std::vector<int>* selected_clusters_indices, int count) {
	std::vector<int> selected_clusters_sizes;
	for (int i = 0; i < shingling_clusters -> size(); i++) {
		int current_size = shingling_clusters -> at(i).size();
		int j = 0;
		while (j < selected_clusters_indices -> size() && current_size < selected_clusters_sizes[j]) {
			j++;
		}
		if(j < count) {
			selected_clusters_indices -> insert(selected_clusters_indices -> begin() + j, i);
			selected_clusters_sizes.insert(selected_clusters_sizes.begin() + j, current_size);
			if (selected_clusters_indices -> size() > count) {
				selected_clusters_indices -> pop_back();
				selected_clusters_sizes.pop_back();
			}
		}
	}
}

void ClusterEvaluator::find_largest_clusters(std::vector<std::unordered_set<std::string>>* clusters, int count, 
						std::vector< std::unordered_set< std::string >>* selected_clusters) {
	std::vector<int> selected_clusters_sizes;
	std::vector<int> selected_clusters_indices;
	for (int i = 0; i < clusters -> size(); i++) {
		int current_size = clusters -> at(i).size();
		int j = 0;
		while (j < selected_clusters_indices.size() && current_size < selected_clusters_sizes[j]) {
			j++;
		}
		if(j < count) {
			selected_clusters_indices.insert(selected_clusters_indices.begin() + j, i);
			selected_clusters_sizes.insert(selected_clusters_sizes.begin() + j, current_size);
			if (selected_clusters_indices.size() > count) {
				selected_clusters_indices.pop_back();
				selected_clusters_sizes.pop_back();
			}
		}
	}
	for (int i = 0; i < selected_clusters_indices.size(); i++) {
		selected_clusters -> push_back(clusters -> at(selected_clusters_indices[i]) );
	}
}

void ClusterEvaluator::send_intersection_sizes(char* two_labels, int key_size, char* two_sizes_repeated, int intersection_size, 
							int* two_ints_repeated, void* allclusters_kv_ptr) {
	KeyValue* allclusters_kv = (KeyValue*) allclusters_kv_ptr;
	CoupledLabels coupled_labels(two_labels, two_sizes_repeated);
	char* label_x = coupled_labels.get_label(X_CLUSTER_LABEL);
	char* label_y = coupled_labels.get_label(Y_CLUSTER_LABEL);
	int x_cluster_size = coupled_labels.get_cluster_size(X_CLUSTER_LABEL);
	int y_cluster_size = coupled_labels.get_cluster_size(Y_CLUSTER_LABEL);
	LabelAndSizeAsValue x_key(x_cluster_size, label_x);
	LabelAndSizeAsValue y_key(y_cluster_size, label_y);
	MatchedCluster x_matched_y(intersection_size, y_cluster_size, label_y);
	MatchedCluster y_matched_x(intersection_size, x_cluster_size, label_x);
	//LabelAndSizeAsValue x_intersect_y(intersection_size, label_y);
	//LabelAndSizeAsValue y_intersect_x(intersection_size, label_x);
	allclusters_kv -> add(x_key.get_message(), x_key.get_message_length(), 
					x_matched_y.get_message(), x_matched_y.get_message_length());
	allclusters_kv -> add(y_key.get_message(), y_key.get_message_length(),
					y_matched_x.get_message(), y_matched_x.get_message_length());
}

void ClusterEvaluator::compute_cluster_entropies(char* cluster_label_and_size, int key_length, char* matched_clusters_and_sizes, 
							int num_matched_clusters, int* match_lengths, void* local_nmi_toolbox_ptr) {
	NmiToolbox* nmi_toolbox = (NmiToolbox*) local_nmi_toolbox_ptr;
	LabelAndSizeAsValue label_and_size(cluster_label_and_size);
	char* key_label = label_and_size.get_label_name();
	int key_cluster_size = label_and_size.get_size();
	double best_fvalue = 0;
	char* best_matched_cluster_label = "@none";
	int best_matched_cluster_size = 0;
	int best_intersection = 0;
	for (int i = 0, offset = 0; i < num_matched_clusters; offset += match_lengths[i], i++) {
		MatchedCluster matched_cluster(matched_clusters_and_sizes + offset);
		char* matched_cluster_label = matched_cluster.get_label_name();
		int matched_cluster_size = matched_cluster.get_cluster_size();
		int intersection = matched_cluster.get_intersection_size();
		double precision = (double) intersection / key_cluster_size;
		double recall = (double) intersection / matched_cluster_size;
		double current_fvalue = 2 * precision * recall / (precision + recall);
		//std::cerr << "svet\t" << key_label << "\t" << matched_cluster_label << "\t" << key_cluster_size << "\t" << matched_cluster_size << "\t" << intersection << std::endl;
		//if (strcmp(key_label, "Y151") == 0) {
		//	std::cerr << "nowfory151: " << matched_cluster_label << " " << intersection << " " << matched_cluster_size << " " << current_fvalue << std::endl;
		//}
		if (current_fvalue > best_fvalue) {
			best_fvalue = current_fvalue;
			best_matched_cluster_label = matched_cluster_label;
			best_matched_cluster_size = matched_cluster_size;
			best_intersection = intersection;
		}
	}
	double weighted_fvalue = key_cluster_size * best_fvalue;
	//if (strcmp(key_label, "Y151") == 0) {
//		std::cerr << key_label << "\t" << key_cluster_size << "\t" << best_intersection << "\t" << best_matched_cluster_label << "\t" 
//				<< best_matched_cluster_size << "\t" << best_fvalue << "\t" << weighted_fvalue << "\n";
	//}
	if (key_label[0] == X_CLUSTER_LABEL) {
		nmi_toolbox -> H_X_cond_Y += weighted_fvalue;
	} else if (key_label[0] == Y_CLUSTER_LABEL) {
		nmi_toolbox -> H_Y_cond_X += weighted_fvalue;
	}
}

void ClusterEvaluator::read_pclust_clusters(std::string filename, int* largest_cluster_size, int min_cluster_limit, 
						MapReduce* pclust_mr, int nmap, bool rename) {
	std::vector< std::unordered_set<std::string> > pclust_clusters;
	read_pclust_clusters_on_all_procs(filename, min_cluster_limit, largest_cluster_size, &pclust_clusters, rename);

	bool unlabeled = true;
	int size_limit = 6000;
	ClusterReadToolbox clusters(&pclust_clusters, unlabeled, &size_limit);

	pclust_mr -> map(nmap, distribute_clusters, &clusters);
	int pclust_clusters_count = pclust_mr -> convert();
	int me; MPI_Comm_rank(MPI_COMM_WORLD, &me);
//	if (me == 0)
//		std::cerr << "pclust clusters count: " << pclust_clusters_count << std::endl;
}

void ClusterEvaluator::read_pclust_clusters_on_all_procs(std::string filename, int min_cluster_limit, int* largest_cluster_size,
								std::vector< std::unordered_set<std::string> >* clusters, bool rename) {
	std::ifstream infile(filename);
	std::string line;
	std::unordered_set<std::string> cluster;
	int largest = 0;
	while(std::getline(infile, line)) {
		if (line.find("Cluster:") != std::string::npos) {
			if (cluster.size() > 0) {
				if (cluster.size() >= min_cluster_limit) {
					clusters -> push_back(cluster);
				}
				if (cluster.size() > largest) {
					largest = cluster.size();
				}
				cluster.clear();
			}
		} else {
			std::string member = line.substr(line.find(" ") + 1);
			if (member[0] == '>') {
				member = member.substr(1);
			}
			if (rename) {
				member = member.substr(0, member.find("_"));
			}
			if (rename) {
				if (cluster.find(member) == cluster.end()) {
					cluster.insert(member);
				}
			} else {
				cluster.insert("/" + member);
			}
		}
	}
	if (cluster.size() > 0) {
		if (cluster.size() >= min_cluster_limit) {
			clusters -> push_back(cluster);
		}
		if (cluster.size() > largest) {
			largest = cluster.size();
		}
	}
	*largest_cluster_size = largest;
}

int ClusterEvaluator::largest_shingling_cluster_size = -1;
int ClusterEvaluator::largest_pfam_cluster_size = -1;
