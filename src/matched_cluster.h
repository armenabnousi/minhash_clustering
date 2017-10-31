#pragma once
#include <cstddef>
#include <string.h>
#include <iostream>

#define X_CLUSTER_LABEL 'X'
#define Y_CLUSTER_LABEL 'Y'

class MatchedCluster {
	char* message = NULL;
	int message_length;
	char* label_name = NULL;
	char label_symbol;
	int size;
	int intersection;
public:
	MatchedCluster();
	MatchedCluster(char* value);
	MatchedCluster(int intersection, int size, char* label_name);
	~MatchedCluster();
	char* get_message();
	int get_message_length();
	void set_label_symbol(char symbol);
	int get_cluster_size();
	int get_intersection_size();
	char* get_label_name();
};
