#include "label_as_value.h"

LabelAndSizeAsValue::LabelAndSizeAsValue() : message(NULL), message_length(0) { }

LabelAndSizeAsValue::LabelAndSizeAsValue(char* value) {
	size = *(int*) value;
	label_name = value + sizeof(int);
	label_symbol = label_name[0];
	message_length = sizeof(int) + strlen(label_name) + 1;
	message = new char[message_length];
	memcpy(message, value, message_length);
}

LabelAndSizeAsValue::LabelAndSizeAsValue(int size, char* label_name) : size(size), label_name(label_name) { 
	message_length = sizeof(int) + strlen(label_name) + 1;
	//message_length = strlen(label_name) + 1;
	message = new char[message_length];
	memcpy(message, (char*) &size, sizeof(int));
	strcpy(&message[sizeof(int)], &label_name[0]);
	//strcpy(&message[0], label_name);
	label_symbol = label_name[0];
}

LabelAndSizeAsValue::~LabelAndSizeAsValue() {
	if (message) {
		delete[] message;
	}
}

char* LabelAndSizeAsValue::get_message() {
	return message;
}

int LabelAndSizeAsValue::get_message_length() {
	return message_length;
}

void LabelAndSizeAsValue::set_label_symbol(char symbol) {
//	std::cout << "setting label to " << symbol << std::endl;
	message[sizeof(int)] = symbol;
	label_symbol = symbol;
//	std::cout << "new message is:" << message + sizeof(int) << std::endl;
}

int LabelAndSizeAsValue::get_size() {
	return size;
}

char* LabelAndSizeAsValue::get_label_name() {
	return label_name;
}
