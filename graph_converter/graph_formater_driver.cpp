#include "../src/graph_formater.h"

void parse_args_mapreduce_edges_to_metis(char** argv, std::string* mrmpi_filename, std::string* metis_filename, std::string* dict_filename);
void parse_args(char** argv, std::string* operation_type);
void parse_args_grappolo_to_mrmpi(char** argv, std::string* dict_filename, std::string* com_filename, std::string* output_filename);
void grappolo_to_mrmpi_formated_file(char** argv);
void mapreduce_edges_to_metis(char** argv);

int main(int argc, char** argv) {
	std::string operation_type;
	parse_args(argv, &operation_type);
	if (operation_type == "grappolo_to_mrmpi_formatted_file") {
		grappolo_to_mrmpi_formated_file(argv);
	} else if (operation_type == "mapreduce_edges_to_metis") {
		mapreduce_edges_to_metis(argv);
	}
}

void grappolo_to_mrmpi_formated_file(char** argv) {
	//std::cout << "entering the program" << std::endl;
	std::string dict_filename = "", com_filename = "", output_filename = "";
	parse_args_grappolo_to_mrmpi(argv, &dict_filename, &com_filename, &output_filename);
	GraphFormater grapher;
	std::unordered_map<int, std::string> id2name;
	std::vector<std::unordered_set<std::string>> communities;
	//std::cout << "going to read dictionary " << std::endl;
	grapher.read_node_dictionary(dict_filename, &id2name);
	//std::cout << "graph read" << std::endl;
	grapher.read_communities(com_filename, &id2name, &communities);
	//std::cout << "communities read" << std::endl;
	grapher.write_communities_mrmpi_format(output_filename, &communities);
}

void mapreduce_edges_to_metis(char** argv) {
	std::string mrmpi_filename, metis_filename, dict_filename;
	parse_args_mapreduce_edges_to_metis(argv, &mrmpi_filename, &metis_filename, &dict_filename);
	GraphFormater grapher;
	std::unordered_map<std::string, std::vector<std::string>> graph;
	int num_vertices = 0, num_edges = 0;
	bool convert_names_to_ids = true;
	grapher.read_mapreduce_edges(mrmpi_filename, &graph, &num_vertices, &num_edges, convert_names_to_ids, dict_filename);
	std::cout << "graph read" << std::endl;
	grapher.generate_metis_from_serial_graph(&graph, num_vertices, num_edges, metis_filename, "usc_louvain");
}

void parse_args(char** argv, std::string* operation_type) {
	*operation_type = argv[1];
}

void parse_args_grappolo_to_mrmpi(char** argv, std::string* dict_filename, std::string* com_filename, std::string* output_filename) {
	//std::cout << "here: " << argv[2] << std::endl;
	//std::cout << "here: " << argv[3] << std::endl;
	*dict_filename = argv[2];
	*com_filename = argv[3];
	*output_filename = argv[4];
}

void parse_args_mapreduce_edges_to_metis(char** argv, std::string* mrmpi_filename, std::string* metis_filename, std::string* dict_filename) {
	*mrmpi_filename = argv[2];
	*metis_filename = argv[3];
	*dict_filename = argv[4];
}

