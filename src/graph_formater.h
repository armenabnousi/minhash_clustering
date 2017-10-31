#pragma once
#include "mpi.h"
#include "mapreduce.h"
#include "keyvalue.h"
#include "mr_util.h"
#include "assert.h"
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <set>
#include <sstream>
#include <fstream>
#include <algorithm>
#include "fasta_util.h"
using namespace MAPREDUCE_NS;

typedef struct VertexName2Id {
	int num_vertices;
	int num_edges;
	std::unordered_map<std::string, int>* name_to_num;
	std::vector<std::string>* num_to_name; //index_to_name
        std::vector<std::set<int>>* edge_lists;
	VertexName2Id(int num_vertices, int num_edges, std::unordered_map<std::string, int>* name_to_num, 
				std::vector<std::string>* num_to_name, std::vector<std::set<int>>* edge_lists) : 
				edge_lists(edge_lists), name_to_num(name_to_num), num_to_name(num_to_name),
				num_vertices(num_vertices), num_edges(num_edges) { };
	int generate_id(std::string name) {
		int id;
		//std::cout << "generating id" << std::endl;
		if (name_to_num -> find(name) == name_to_num -> end()) {
                	id = name_to_num -> size();
                	name_to_num -> insert(std::make_pair(name, id));
			num_to_name -> at(id) = name;
		} else {
			id = name_to_num -> at(name);
		}
		//std::cout << "id generated" << std::endl;
		return id;
        }
} VertexName2Id;

class GraphFormater {
public:
	GraphFormater() {}
	void generate_metis_from_mapreduce_graph(MapReduce* mr, std::string filename, std::string com_detection_method = "grappolo");
	void read_node_dictionary(std::string filename, std::unordered_map<int, std::string>* dict, 
					std::unordered_map<std::string, int>* name2num = NULL);
	void read_communities(std::string filename, std::unordered_map<int, std::string>* dict, 
				std::vector<std::unordered_set<std::string>>* communities);
	void write_communities_mrmpi_format(std::string filename, std::vector<std::unordered_set<std::string>>* communities);
	void generate_mapreduce_formated_file_from_grappolo_communities(std::string dict_filename, std::string com_file, 
										std::string output_filename);
	void generate_metis_from_serial_graph(std::unordered_map<std::string, std::vector<std::string>>* named_graph, int num_vertices,
										int num_edges, std::string metis_filename, 
										std::string com_detection_method = "grappolo"); 
	void read_mapreduce_edges(std::string filename, std::unordered_map<std::string, std::vector<std::string>>* graph, int* num_vertices,
										int* num_edges, bool convert_names_to_ids, 
										std::string dict_filename = NULL);
private:
	std::vector< std::string > split_string(std::string str, char delim);
	void write_graph_to_metis_file(std::string filename, VertexName2Id* graph_container, std::string com_detection_method);
	void static convert_graph_to_ids(char* key, int keybytes, char* values, int nvalues, int* valuebytes, void* ptr);
	void convert_graph_to_ids_serial(std::unordered_map<std::string, std::vector<std::string>>* named_graph,
						VertexName2Id* vertex_name2id);
	void static extract_key_names(char* key, int key_len, char* vals, int nvals, int* val_lens, void* set_ptr);
};
