#include "graph_formater.h"

void GraphFormater::generate_metis_from_mapreduce_graph(MapReduce* mr, std::string metis_filename, std::string com_detection_method) {
	int num_vertices = 0, num_edges = 0;
	//std::cout << "gathering kmv" << std::endl;
	MrToolbox::gather_kmv(mr, 1, true, &num_vertices, &num_edges);
	//std::cout << "kmv gather done" << std::endl;
	std::unordered_map<std::string, int> name_to_num;
	std::vector<std::string> num_to_name(num_vertices, "");
	std::vector<std::set<int>> edge_lists(num_vertices, std::set<int>());
	VertexName2Id graph_container(num_vertices, num_edges, &name_to_num, &num_to_name, &edge_lists);
	int me; MPI_Comm_rank(MPI_COMM_WORLD, &me);
	//if (me == 0) std::cout << "calling the scan method" << std::endl;
	mr -> scan(convert_graph_to_ids, &graph_container);
	//if (me == 0) std::cout << "conversion done!" << std::endl;
	if (me == 0) {
		std::cout << "here: " << num_vertices << "\t" << graph_container.name_to_num -> size() << "\t" << graph_container.num_to_name -> size() << std::endl;
/*		// printing for debug purposes
		std::set<std::string> dict_names;
		std::stringstream dict_out;
		for (auto iter = graph_container.name_to_num -> begin(); iter != graph_container.name_to_num -> end(); iter++) {
			dict_out << iter -> first << "\n";
			dict_names.insert(iter -> first);
		}
		std::ofstream d_file("/home/armen.abnousi/shingle/temp_dict_output");
		d_file << dict_out.str();
		d_file.close();
		std::set<std::string> mr_names;
		mr -> scan(extract_key_names, &mr_names);
		mr -> print("/home/armen.abnousi/shingle/temp_mr_output_compare_to_dict", 0, -1, 1, 5, 0);
		std::vector<std::string> v(num_vertices);                      // 0  0  0  0  0  0  0  0  0  0
		std::vector<std::string>::iterator it;
		it = std::set_difference(mr_names.begin(), mr_names.end(), dict_names.begin(), dict_names.end(), v.begin());
		v.resize(it - v.begin());
		std::cout << "mr names size: " << mr_names.size() << " dict_names size: " << dict_names.size() << "difference: " << v.size() << std::endl;
		for (auto iter = v.begin(); iter != v.end(); iter++) {
			std::cout << "v member: " << *iter << std::endl;
		}
*/		//print done
		assert(num_vertices == graph_container.name_to_num -> size());
		assert(num_vertices == graph_container.edge_lists -> size());
		write_graph_to_metis_file(metis_filename, &graph_container, com_detection_method);
	}
}

void GraphFormater::extract_key_names(char* key, int key_len, char* vals, int nvals, int* val_lens, void* set_ptr) {
	std::set<std::string>* names = (std::set<std::string>*) set_ptr;
	names -> insert((std::string) key);
}

void GraphFormater::generate_mapreduce_formated_file_from_grappolo_communities(std::string dict_filename, std::string com_filename,
									std::string output_filename) {
	int me; MPI_Comm_rank(MPI_COMM_WORLD, &me);
	if (me == 0) {
		std::unordered_map<int, std::string> id2name;
		std::vector<std::unordered_set<std::string>> communities;
		read_node_dictionary(dict_filename, &id2name);
        	std::cerr << "graph read" << std::endl;
        	read_communities(com_filename, &id2name, &communities);
        	std::cerr << "communities read" << std::endl;
        	write_communities_mrmpi_format(output_filename, &communities);
	}
}

//All keys and members should have a prefix label (e.g. '@' or '/' when entering this function!
void GraphFormater::convert_graph_to_ids(char* key, int keybytes, char* members, int nmembers, int* member_lens, void* vertex_dict_ptr) {
	VertexName2Id* vertex_name2id = (VertexName2Id*) vertex_dict_ptr;
	int key_id = vertex_name2id -> generate_id(std::string(key + 1));
	for (int i = 0, offset = 0; i < nmembers; offset += member_lens[i], i++) {
		int member_id = vertex_name2id -> generate_id(std::string(members + offset + 1));
		vertex_name2id -> edge_lists -> at(key_id).insert(member_id);
	}
}

void GraphFormater::write_graph_to_metis_file(std::string filename, VertexName2Id* graph_container, std::string com_detection_method) {
	std::stringstream output;
	std::stringstream dictionary_output;
	std::stringstream edges;
	bool edge_file_requested = false;
	if (com_detection_method == "grappolo") {
		output << std::to_string(graph_container -> num_vertices) << " " << std::to_string(graph_container -> num_edges) 
			<< " " << 0 << "\n";
		edge_file_requested = false;
	} else if (com_detection_method == "usc_louvain") {
		output << std::to_string(graph_container -> num_vertices) << " " << std::to_string(graph_container -> num_edges) << std::endl;
		edge_file_requested = true;
	}
	for (int i = 0; i < graph_container -> num_vertices; i++) {
		std::string current_node = std::to_string(i + 1);
		dictionary_output << "#" << current_node << " " << graph_container -> num_to_name -> at(i) << "\n";
		for (auto member = graph_container -> edge_lists -> at(i).begin(); 
			member != graph_container -> edge_lists -> at(i).end(); member++) {
			std::string neighbor = std::to_string(*member + 1);
			output << neighbor << " ";
			if (edge_file_requested) {
				edges << current_node << " " << neighbor << "\n";
			}
		}
		output << "\n";
	}

	std::string dictionary_filename = filename + "_dictionary";
	std::ofstream dict_file(dictionary_filename);
	dict_file << dictionary_output.str();
	dict_file.close();

	std::ofstream metis_file(filename);
	metis_file << output.str();
	metis_file.close();

	if (edge_file_requested) {
		std::string edge_filename = filename + "_edges";
		std::ofstream edge_file(edge_filename);
		edge_file << edges.str();
		edge_file.close();
	}
}

void GraphFormater::read_node_dictionary(std::string filename, std::unordered_map<int, std::string>* dict, 
						std::unordered_map<std::string, int>* name2num) {
	//std::cerr << "trying to open file" << std::endl;
	std::ifstream infile(filename);
	std::string line;
	//std::cerr << "starting the loop" << std::endl;
	while (std::getline(infile, line)) {
		//std::cerr << "line read: " << line << std::endl;
		int id;
		char name[100];
		sscanf(&line[0], "#%d %s", &id, name);
		//std::cerr << " read: " << id << " " << name << std::endl;
		dict -> insert(std::make_pair(id, std::string(name)));
		if (name2num) {
			name2num -> insert(std::make_pair(std::string(name), id));
		}
		//std::cerr << "added to dict" << std::endl;
	}
}

void GraphFormater::read_communities(std::string filename, std::unordered_map<int, std::string>* dict, 
					std::vector<std::unordered_set<std::string>>* communities) {
	int line_number = 1;
	std::ifstream infile(filename);
	std::string line;
	std::unordered_map<int, std::unordered_set<std::string>> com_map;
	while (std::getline(infile, line)) {
		int community;
		sscanf(&line[0], "%d", &community);
		if (com_map.find(community) == com_map.end()) {
			//std::unordered_set<std::string> empty_set;
			com_map.insert(std::make_pair(community, std::unordered_set<std::string>()));
		}
		com_map[community].insert("/" + dict -> at(line_number));
		line_number++;
	}
	//std::cerr << "communities read: " << com_map.size() << std::endl;
	for (auto iter = com_map.begin(); iter != com_map.end(); iter++) {
		communities -> push_back(iter -> second);
	}
}

void GraphFormater::write_communities_mrmpi_format(std::string filename, std::vector<std::unordered_set<std::string>>* communities) {
	std::stringstream output;
	int num_communities = communities -> size();
	for (int i = 0; i < num_communities; i++) {
		output << "KMV pair: proc 0, nvalues " << communities -> at(i).size()  << ", sizes 9 4679, key @" << std::to_string(i) << ", values";
		for (auto iter = communities -> at(i).begin(); iter != communities -> at(i).end(); iter++) {
			output << " " << *iter;
		}
		output << "\n";
	}
	
	std::ofstream outfile(filename);
	outfile << output.str();
	outfile.close();
}

void GraphFormater::generate_metis_from_serial_graph(std::unordered_map<std::string, std::vector<std::string>>* named_graph, int num_vertices,
							int num_edges, std::string metis_filename, std::string com_detection_method) {
	std::unordered_map<std::string, int> name_to_num;
	std::vector<std::string> num_to_name(num_vertices, "");
	//std::vector<std::string> num_to_name;
	std::vector<std::set<int>> edge_lists(num_vertices, std::set<int>());
	/*for (auto iter = named_graph -> begin(); iter != named_graph -> end(); iter++) {
		std::cout << "node " << iter -> first << std::endl;
		for (auto it = iter -> second.begin(); it != iter -> second.end(); it++) {
			std::cout << "\t" << *it << std::endl;
		}
	}*/
	//read_node_dictionary(dict_filename, &num_to_name, &name_to_num);
        VertexName2Id graph_container(num_vertices, num_edges, &name_to_num, &num_to_name, &edge_lists);
	std::cout << "going_to_convert" << std::endl;
	convert_graph_to_ids_serial(named_graph, &graph_container);
	std::cout << "converted" << std::endl;
	graph_container.num_vertices = num_vertices;
	graph_container.num_edges = num_edges;
	std::cout << "going to write" << std::endl;
	write_graph_to_metis_file(metis_filename, &graph_container, com_detection_method);
}

void GraphFormater::convert_graph_to_ids_serial(std::unordered_map<std::string, std::vector<std::string>>* named_graph, 
						VertexName2Id* vertex_name2id) {
	for (auto iter = named_graph -> begin(); iter != named_graph -> end(); iter++) {
		int key_id = vertex_name2id -> generate_id(iter -> first);
		for (auto neighbor = iter -> second.begin(); neighbor != iter -> second.end(); neighbor++) {
			int neighbor_id = vertex_name2id -> generate_id(*neighbor);
			vertex_name2id -> edge_lists -> at(key_id).insert(neighbor_id);
		}
	}
}

void GraphFormater::read_mapreduce_edges(std::string filename, std::unordered_map<std::string, std::vector<std::string>>* graph,
						int* num_vertices, int* num_edges, bool convert_names_to_ids, std::string dict_filename) {
	FastaToolbox fasta_util;
	if (convert_names_to_ids) {
		fasta_util.read_conversion_table(dict_filename);
	}
	std::ifstream infile(filename);
	std::string line;
	int edge_count = 0;
	while(std::getline(infile, line)) {
		edge_count++;
		std::vector<std::string> tokens = split_string(line, ' ');
		std::string source = tokens[8].substr(0, tokens[8].length() - 1);
		std::string dest = tokens[10];
		if (convert_names_to_ids) {
			int source_dom_pos = source.find("_");
			std::string source_seq_id = source.substr(0, source_dom_pos);
			std::string source_dom_id = source.substr(source_dom_pos);
			int dest_dom_pos = dest.find("_");
			std::string dest_seq_id = dest.substr(0, dest_dom_pos);
			std::string dest_dom_id = dest.substr(dest_dom_pos);
			std::string source_id = std::to_string(fasta_util.name_to_num[source_seq_id]) + source_dom_id ;
			std::string dest_id = std::to_string(fasta_util.name_to_num[dest_seq_id]) + dest_dom_id;
			source = source_id;
			dest = dest_id;
		}
		if (graph -> find(source) == graph -> end()) {
			std::vector<std::string> edges;
			graph -> insert(std::make_pair(source, edges));
		}
		graph -> at(source).push_back(dest);
	}
	assert(edge_count % 2 == 0);
	edge_count /= 2;
	*num_vertices = graph -> size();
	*num_edges = edge_count;
}

std::vector< std::string > GraphFormater::split_string(std::string str, char delim) {
        std::vector< std::string > tokens;
        std::stringstream ss(str);
        std::string token;
        while (std::getline(ss, token, delim)) {
                tokens.push_back(token);
        }
        return tokens;
}
