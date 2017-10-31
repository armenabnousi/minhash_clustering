#pragma once
#include <string.h>
#define CURRENT_LABEL '@'
#define INCOMING_LABEL '^'
#define NEIGHBOR_NODE '/'

class MessageValue{
private:
	char message_type;
	char* sequence_name = NULL;
	char* message_array = NULL;
	int length;
public:
	MessageValue();
	MessageValue(char type, char* name);
	MessageValue(char* message);
	~MessageValue();
	MessageValue(const MessageValue& rhs);
	char* get_message_array();
	void set_name(char*);
	char* get_name();
	void set_type(char type);
	char get_type();
	bool operator==(MessageValue& rhs);
	MessageValue& operator=( const MessageValue& rhs );
	void assign(MessageValue* rhs);
	int get_length();
	
};
