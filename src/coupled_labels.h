#pragma once
#include <cstddef>	//for NULL
#include <string.h>
#include "label_as_value.h" //for label symbol macros

class CoupledLabels {
	char* label_x = NULL;
	char* label_y = NULL;
	int label_x_cluster_size;
	int label_y_cluster_size;
	char* labels_message = NULL;
	int labels_message_length;
	char* sizes_message = NULL;
	int sizes_message_length;
	bool from_same_side;
public:
	CoupledLabels();
	CoupledLabels(char* label1, char* label2, int label1_members_count, int label2_members_count);
	CoupledLabels(char* coupled_labels, char* coupled_cluster_sizes);
	~CoupledLabels();
	char* get_labels_message();
	int get_labels_message_length();
	char* get_sizes_message();
	int get_sizes_message_length();
	char* get_label(char label_symbol);
	int get_cluster_size(char label_symbol);
	bool are_from_same_side();
};
