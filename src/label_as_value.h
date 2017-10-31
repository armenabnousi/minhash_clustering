#pragma once
#include <cstddef>
#include <string.h>
#include <iostream>

#define X_CLUSTER_LABEL 'X'
#define Y_CLUSTER_LABEL 'Y'

class LabelAndSizeAsValue {
	char* message = NULL;
	int message_length;
	char* label_name = NULL;
	char label_symbol;
	int size;
public:
	LabelAndSizeAsValue();
	LabelAndSizeAsValue(char* value);
	LabelAndSizeAsValue(int size, char* label_name);
	~LabelAndSizeAsValue();
	char* get_message();
	int get_message_length();
	void set_label_symbol(char symbol);
	int get_size();
	char* get_label_name();
};
