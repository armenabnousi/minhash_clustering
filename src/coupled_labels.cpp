#include "coupled_labels.h"

CoupledLabels::CoupledLabels() : labels_message(NULL), labels_message_length(0), sizes_message(NULL), sizes_message_length(0) { }

CoupledLabels::CoupledLabels(char* label1, char* label2, int label1_members_count, int label2_members_count) : 
						label_x(label1), label_y(label2), label_x_cluster_size(label1_members_count),
						label_y_cluster_size(label2_members_count) {
	if (label_x[0] == label_y[0]) {
		from_same_side = true;
	} else {
		from_same_side = false;
		if (label_x[0] != X_CLUSTER_LABEL) {
			label_x = label2;
			label_y = label1;
			label_x_cluster_size = label2_members_count;
			label_y_cluster_size = label1_members_count;
		}
	}
	
	//generating the labels message
	int label_x_length = strlen(label_x) + 1;
	labels_message_length = label_x_length + strlen(label_y) + 1;
	labels_message = new char[labels_message_length];
	strcpy(labels_message, label_x);
	strcpy(&labels_message[label_x_length], label_y);

	//generating the sizes message
	int size_of_int = sizeof(int);
	sizes_message_length = size_of_int * 2;
	sizes_message = new char[sizes_message_length];
	memcpy(sizes_message, (char*) &label_x_cluster_size, size_of_int);
	memcpy(sizes_message + size_of_int, (char*) &label_y_cluster_size, size_of_int);
}

CoupledLabels::CoupledLabels(char* coupled_labels, char* coupled_cluster_sizes) : from_same_side(true), 
										labels_message(NULL), sizes_message(NULL) {
	label_x = coupled_labels;
	label_y = coupled_labels + strlen(label_x) + 1;
	label_x_cluster_size = *(int*) coupled_cluster_sizes;
	label_y_cluster_size = *(int*) (coupled_cluster_sizes + sizeof(int));
	if (label_x[0] != label_y[0]) from_same_side = false;
}

CoupledLabels::~CoupledLabels() {
	if (labels_message) {
		delete[] labels_message;
		delete[] sizes_message;
	}
}

char* CoupledLabels::get_labels_message() {
	return labels_message;
}

int CoupledLabels::get_labels_message_length() {
	return labels_message_length;
}

char* CoupledLabels::get_sizes_message() {
	return sizes_message;
}

int CoupledLabels::get_sizes_message_length() {
	return sizes_message_length;
}

char* CoupledLabels::get_label(char label_symbol) {
	char* result;
	if (label_symbol == X_CLUSTER_LABEL) {
		result = label_x;
	} else {
		result = label_y;
	}
	return result;
}

int CoupledLabels::get_cluster_size(char label_symbol) {
	int result;
	if (label_symbol == X_CLUSTER_LABEL) {
		result = label_x_cluster_size;
	} else {
		result = label_y_cluster_size;
	}
	return result;
}

bool CoupledLabels::are_from_same_side() {
	return from_same_side;
}
