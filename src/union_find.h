#pragma once
#include <unordered_map>
#include <vector>
#include <set>
#include <string>
using std::unordered_map;
using std::string;

class UnionFind{
	unordered_map<size_t, size_t> parent_id;
	unordered_map<size_t, size_t> subtree_size;
	unordered_map<string, size_t> name2id;
	unordered_map<size_t, string> id2name;
	int node_count;

public:
	UnionFind();
	~UnionFind();
	void run_connected_component(std::vector<std::set<std::string>>* clusters,
					std::unordered_map<std::string, std::vector<std::string>>* connected_components);
	string find(string str);
	void union_trees(string s1, string s2);
	unordered_map<string, size_t> get_name2id();
private:
	size_t get_id(string str);
};
