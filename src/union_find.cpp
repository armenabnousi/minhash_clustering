#include "union_find.h"
#include <iostream>
using std::string;
using std::cout;
using std::endl;

UnionFind::UnionFind(){
	node_count = 0;
}

UnionFind::~UnionFind(){};

void UnionFind::run_connected_component(std::vector<std::set<std::string>>* clusters, 
				std::unordered_map<std::string, std::vector<std::string>>* connected_components) {
	for (int i = 0; i < clusters -> size(); i++) {
		std::set<std::string> cluster = clusters -> at(i);
		for (auto node1 = cluster.begin(); node1 != std::prev(cluster.end()); node1++) {
			std::string root1 = find(*node1);
			/*
			if (node1 -> length() < 3 || root1.length() < 3) {
				std::cout << "error in rcvd" << *node1 << " and " << root1 << std::endl;
				exit(0);
			}
			*/
			for (auto node2 = std::next(node1); node2 != cluster.end(); node2++) {
				//since each edge is present twice in the list, this 'if' avoids extra work
				if (node1 -> compare(*node2) < 0) {
					std::string root2 = find(*node2);
					/*
					if (node2 -> length() < 3 || root2.length() < 3) {
                                		std::cout << "error in rcvd" << *node2 << " and " << root2 << std::endl;
						exit(0);
                        		}
					*/
					if (root1 != root2) {
						union_trees(root1, root2);
					}
				}
			}
		}
	}

	std::unordered_map<std::string, std::string> connected_component_labels;
	std::unordered_map<std::string, std::vector<std::string>> temp_connected_components;
	for (auto iter = name2id.begin(); iter != name2id.end(); iter++) {
		std::string seq = iter -> first;
		std::string cluster_ufo_label = find(seq);
		if (temp_connected_components.find(cluster_ufo_label) == temp_connected_components.end()) {
			std::vector<std::string> cluster;
			/*
			if (cluster_ufo_label.length() < 3) {
				std::cout << "found an error " << cluster_ufo_label << " for seq " << seq << std::endl;
				exit(0);
			}
			*/
			temp_connected_components.insert(std::make_pair(cluster_ufo_label, cluster));
			connected_component_labels.insert(std::make_pair(cluster_ufo_label, seq));
		}
		(temp_connected_components.at(cluster_ufo_label)).push_back(seq);
		if (seq.compare(connected_component_labels[cluster_ufo_label]) < 0 ) {
			connected_component_labels[cluster_ufo_label] = seq;
		}
	}

	for (auto iter = temp_connected_components.begin(); iter != temp_connected_components.end(); iter++) {
		/*
		if (iter -> first.length() < 3) {
			std::cout <<"found label is short: " << iter -> first << std::endl;
		}
		*/
		std::string new_label = connected_component_labels.at(iter -> first);
		/*
		if (new_label.length() < 3) {
			std::cout << "here is an error " << iter -> first << " and " << new_label << std::endl;
			exit(0);
		}
		*/
		connected_components -> insert(std::make_pair(new_label, iter -> second));
	}
}

string UnionFind::find(string str){
	size_t str_id = get_id(str);
	while ( str_id != parent_id[str_id]){
		parent_id[ str_id ] = parent_id[ parent_id[str_id] ];
		str_id = parent_id [ str_id ];
	}
	string root_name = id2name[str_id];
	return root_name;
}

void UnionFind::union_trees(string s1, string s2){
	size_t root1_id = get_id(s1);
	size_t root2_id = get_id(s2);
	if ( subtree_size[root1_id] < subtree_size[root2_id] ){
		parent_id[root1_id] = root2_id;
	} else {
		parent_id[root2_id] = root1_id;
	}
}

size_t UnionFind::get_id(string str){
	size_t str_id = 0;
	auto str_hash = name2id.find(str);
	if (str_hash == name2id.end()){
		str_id = node_count;
		name2id[str] = str_id;
		id2name[str_id] = str;
		parent_id[str_id] = str_id;
		subtree_size[str_id] = 1;
		node_count++;
	} else {
		str_id = name2id[str];
	}
	return str_id;
}

unordered_map<string, size_t> UnionFind::get_name2id(){
	return name2id;
}
