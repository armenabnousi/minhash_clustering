#include "matched_cluster.h"

MatchedCluster::MatchedCluster() : message(NULL), message_length(0) { }

MatchedCluster::MatchedCluster(char* value) {
	intersection = *(int*) value;
	size = *(int*) (value + sizeof(int));
	label_name = value + (2 * sizeof(int));
	label_symbol = label_name[0];
	message_length = 2 * sizeof(int) + strlen(label_name) + 1;
	message = new char[message_length];
	memcpy(message, value, message_length);
}

MatchedCluster::MatchedCluster(int intersection, int size, char* label_name) : size(size), label_name(label_name) { 
	message_length = 2 * sizeof(int) + strlen(label_name) + 1;
	//message_length = strlen(label_name) + 1;
	message = new char[message_length];
	memcpy(message, (char*) &intersection, sizeof(int));
	memcpy(message + sizeof(int), (char*) &size, sizeof(int));
	strcpy(&message[2 * sizeof(int)], &label_name[0]);
	//strcpy(&message[0], label_name);
	label_symbol = label_name[0];
}

MatchedCluster::~MatchedCluster() {
	if (message) {
		delete[] message;
	}
}

char* MatchedCluster::get_message() {
	return message;
}

int MatchedCluster::get_message_length() {
	return message_length;
}

void MatchedCluster::set_label_symbol(char symbol) {
//	std::cout << "setting label to " << symbol << std::endl;
	message[sizeof(int)] = symbol;
	label_symbol = symbol;
//	std::cout << "new message is:" << message + sizeof(int) << std::endl;
}

int MatchedCluster::get_cluster_size() {
	return size;
}

char* MatchedCluster::get_label_name() {
	return label_name;
}

int MatchedCluster::get_intersection_size() {
	return intersection;
}
