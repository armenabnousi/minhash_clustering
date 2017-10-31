#include "message_value.h"
//#include <iostream>

MessageValue::MessageValue() : message_array(NULL), sequence_name(NULL) {}

MessageValue::MessageValue(char type, char* name) : message_type(type) {
	sequence_name = name;
	length = 1 + strlen(sequence_name) + 1;
	if (message_array != NULL) {
		delete[] message_array;
		message_array = NULL;
	}
	message_array = new char[length];
	message_array[0] = message_type;
	strcpy(&message_array[1], sequence_name);
}

MessageValue::MessageValue(char* message) {
	if (message[0] == CURRENT_LABEL || message[0] == INCOMING_LABEL || message[0] == NEIGHBOR_NODE) {
		message_type = message[0];
		sequence_name = &message[1];
	} else {
		message_type = NEIGHBOR_NODE;
		sequence_name = message;
	}
	length = 1 + strlen(sequence_name) + 1;
	if (message_array != NULL) {
		delete[] message_array;
		message_array = NULL;
	}
	message_array = new char[length];
	message_array[0] = message_type;
	strcpy(&message_array[1], sequence_name);
}

MessageValue::~MessageValue() {
	if (message_array) {
		//std::cout << "deleting array " << std::endl;
		//std::cout << message_array << std::endl;
		delete[] message_array;
		message_array = NULL;
		//std::cout << "it's now deleted" << std::endl;
	}
	//std::cout << "exiting destructor" << std::endl;
}

MessageValue::MessageValue(const MessageValue& rhs) {
	this -> length = rhs.length;
	this -> sequence_name = rhs.sequence_name;
	this -> message_type = rhs.message_type;
	if (this -> message_array != NULL) {
		delete[] (this -> message_array);
		this -> message_array = NULL;
	}
	this -> message_array = new char[length];
	memcpy(this -> message_array, rhs.message_array, length);
}

char* MessageValue::get_message_array() {
	return message_array;
}

void MessageValue::set_name(char* name) {
	strcpy(sequence_name, name);
}

char* MessageValue::get_name() {
	return sequence_name;
}

void MessageValue::set_type(char type) {
	message_type = type;
	if (message_array != NULL) {
		message_array[0] = type;
	}
}
char MessageValue::get_type() {
	return message_type;
}

MessageValue& MessageValue::operator=( const MessageValue& rhs ) {
	if (this -> message_array != NULL) {
		delete[] (this -> message_array);
		this -> message_array = NULL;
	}
	this -> length = rhs.length;
	this -> sequence_name = rhs.sequence_name;
	this -> message_type = rhs.message_type;
	this -> message_array = new char[length];
	memcpy(this -> message_array, rhs.message_array, length);
}

bool MessageValue::operator==(MessageValue& rhs) {
	return ( !strcmp( this -> sequence_name, rhs.get_name() ) && 
		this -> message_type == rhs.get_type() );
}

void MessageValue::assign(MessageValue* rhs) {
	this -> message_type = rhs -> get_type();
	this -> sequence_name = rhs -> get_name();
	if (this -> message_array != NULL) {
		delete[] message_array;
		message_array = NULL;
	}
	length = rhs -> get_length();
	message_array = new char[length];
	strcpy( this -> message_array, rhs -> get_message_array());
}

int MessageValue::get_length() {
	return length;
}
