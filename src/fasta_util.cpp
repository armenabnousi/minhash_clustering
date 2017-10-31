#include "fasta_util.h"
#ifndef FASTA_HEADER_LENGTH
	#define FASTA_HEADER_LENGTH 200
#endif
#ifndef LOG_SCALE
	#define LOG_SCALE 0
#endif
#ifndef LOG_RANK
	#define LOG_RANK -1
#endif

std::list<Protein> FastaToolbox::read_fasta(char* text, long long int text_size, int proc) {
	std::list<Protein> proteins;
	bool serial = false;
	char* text_begin = text;
	while (text[0] == '\n') {
		text++;
		text_size--;
		if(LOG_SCALE > 1 && (LOG_RANK < 0 || LOG_RANK == proc)) std::cout << proc << ": forwarding text. next is:" << (int)text[0] << std::endl;
	}
	//std::cout << "text is: " << text << std::endl;
	if (text_size == -1) {
		serial = true;	
		read_fasta_serial(text, &proteins);
		//std::cout << "read proteins: " << proteins.size() << " ex: " << (*(proteins.begin())).name << " " << (*(proteins.begin())).sequence << std::endl;
		//std::cout << "last: " << proteins.size() << " ex: " << (*(proteins.rbegin())).name << " " << (*(proteins.rbegin())).sequence << std::endl;
		//text_size = strlen(text);
	} else {
		if(LOG_SCALE > 1 && (LOG_RANK < 0 || LOG_RANK == proc)) {
			std::cout << "\tTASK " << proc << ": COUNTING SEQUENCE KMERS..." << std::endl;
		}
	
		//FIRST READING: COUNTING KMERS
		//reading the sequences start from here:
		if (text[0] == '>' || text[0] == '\0'){	//this is only for task=0, because strtok automatically moves
					//the start pointer but since the program doesn't know, we move the
					//text pointer later in the loop based on the strlen of the first token
					//which will be off by one.
			if(LOG_SCALE > 1 && (LOG_RANK < 0 || LOG_RANK == proc)) {
	         		std::cout << "\tTASK " << proc << ": COUNTING SEQUENCE KMERS..." << std::endl;
	        	}
			text++;
			text_size--;
		}
		if(LOG_SCALE > 1 && (LOG_RANK < 0 || LOG_RANK == proc)) {
			std::cout << "\tTASK " << proc << ": COUNTING SEQUENCE KMERS..." << std::endl;
	        }
		char* fasta_pack = strtok(text, ">");
		if (LOG_SCALE > 1 && (LOG_RANK < 0 || LOG_RANK == proc)) std::cout<<"checkthis:"<<(int) text[0]<<"\tand\t" << strlen(fasta_pack) <<std::endl;
		int name_length = FASTA_HEADER_LENGTH;
		int next_offset = 0;
		int moved_offset = 0;
		while (text_size > 0 && fasta_pack != NULL) {
			moved_offset = next_offset;
			if (LOG_SCALE > 1 && (LOG_RANK < 0 || LOG_RANK == proc)) std::cout << "first" << std::endl;
			next_offset = strlen(fasta_pack) + 1;
			if (LOG_SCALE > 1 && (LOG_RANK < 0 || LOG_RANK == proc)) std::cout << next_offset << " " << text_size << std::endl;
			text += next_offset;
			text_size -= next_offset;
			if (LOG_SCALE > 1 && (LOG_RANK < 0 || LOG_RANK == proc)) std::cout << "here done with this!" << std::endl;
			if (LOG_SCALE > 1 && (LOG_RANK < 0 || LOG_RANK == proc)) std::cout << "text[0] = " << text[0] << std::endl;
			//std::cout << "processing: " << fasta_pack << std::endl;
			char sequence_name[name_length + 1], sequence[strlen(fasta_pack)];
			parse_fasta_sequence(fasta_pack, sequence_name, sequence, name_length, proc);
			if (LOG_SCALE > 1 && (LOG_RANK < 0 || LOG_RANK == proc)) std::cout << "<proc " << proc << " seq read: " << sequence_name << "\t" << std::endl;
			Protein protein((std::string) sequence_name, (std::string) sequence);
			proteins.push_back(protein);
			fasta_pack[-1] = '>';
			if (LOG_SCALE > 1 && (LOG_RANK < 0 || LOG_RANK == proc)) std::cout << "second" << std::endl;
			fasta_pack = strtok(text, ">");
			if (LOG_SCALE > 1 && (LOG_RANK < 0 || LOG_RANK == proc)) std::cout << "third" << std::endl;
		}
		text = text_begin;
		if (LOG_SCALE > 1 && (LOG_RANK < 0 || LOG_RANK == proc)) std::cout << "fourth" << std::endl;
		//if(serial) delete[] text;
		if (LOG_SCALE > 3 && (LOG_RANK < 0 || LOG_RANK == proc)) std::cout << "fifth" << std::endl;
	}
	return proteins;
}

void FastaToolbox::parse_fasta_sequence(char* fasta_pack, char* seq_name, char* sequence, int name_length, int proc = 0){
	char* temp_name = strtok(fasta_pack, "\n");
	strncpy(seq_name, temp_name, name_length);
	seq_name[name_length - 1] = '\0';
	temp_name[strlen(temp_name)] = '\n';
	if (LOG_SCALE > 3 && LOG_RANK == proc) std::cout<<"name read: "<< seq_name << std::endl;
	char* seq_line = strtok(NULL, "\n");
	if (LOG_SCALE > 3 && LOG_RANK == proc) std::cout << "check1 " << seq_line << std::endl;
	int offset = 0;
	while(seq_line != NULL && seq_line[0] != '>'){
		strcpy(sequence + offset, seq_line);
		offset += strlen(seq_line);
		seq_line[strlen(seq_line)] = '\n';
		if (LOG_SCALE > 3 && LOG_RANK == proc) std::cout << "seqsofar: " << sequence << std::endl;
		seq_line = strtok(NULL, "\n");
	}
	sequence[offset] = '\0';
}

void FastaToolbox::convert_name_to_num(std::list<Protein> *proteins, bool domained) {
	if (!domained) {
		int seq_id = 1;
		for (auto protein = proteins -> begin(); protein != proteins -> end(); protein++) {
			if ( name_to_num.find(protein -> name) == name_to_num.end() ) {
				protein -> name_id = seq_id;
				name_to_num[protein -> name] = seq_id;
				seq_id++;
			} else {
				protein -> name_id = name_to_num[protein->name];
			}
		}
	} else if (domained) {
		int seq_id = 1;
		for (auto protein = proteins -> begin(); protein != proteins -> end(); protein++) {
			int dom_id;
			std::string name = parse_protein_name(protein -> name, &dom_id);
			if (name_to_num.find(name) == name_to_num.end()) {
				protein -> name_id = seq_id;
				protein -> domain_id = dom_id;
				name_to_num[name] = seq_id;
				seq_id++;
			} else {
				protein -> name_id = name_to_num[name];
				protein -> domain_id = dom_id;
			}
		}
	}
}

std::string FastaToolbox::parse_protein_name(std::string name, int* domain_id) {
	std::string seq_name;
	int sep_index = name.find_last_of('_');
	seq_name = name.substr(0, sep_index);
	std::string domain_name = name.substr(sep_index + 1);
	*domain_id = std::stoi(domain_name);
	return seq_name;
}

void FastaToolbox::write_conversion_table(std::string file_name) {
	std::string output = "";
	for (auto protein = name_to_num.begin(); protein != name_to_num.end(); protein++) {
		output += protein -> first + '\t' + std::to_string((long long int) protein -> second) + '\n';
	}
	std::ofstream my_file;
	my_file.open(file_name);
	my_file << output;
	my_file.close();
}

void FastaToolbox::read_conversion_table(std::string file_name) {
	std::ifstream my_file(file_name);
	std::string line;
	if (my_file.is_open()) {
		while( getline(my_file, line) ) {
			int name_id;
			std::string name;
			int name_end = line.find('\t');
			name = line.substr(0, name_end);
			int name_id_end = line.find(name_end + 1, '\t');
			if (name_id_end != std::string::npos) {
				name_id = stoi( line.substr(name_end + 1, name_id_end) );
			} else {
				name_id = stoi( line.substr(name_end + 1) );
			}
			name_to_num[name] = name_id;
		}
		my_file.close();
	}
}

void FastaToolbox::convert_pfam_file_to_numbered(std::string input_filename, std::string output_filename) {
	std::ifstream my_file(input_filename);
	std::string line;
	std::string output = "";
	while( getline(my_file, line) ) {
		if (line[0] == '#' || line[0] == '\n') {
			continue;
		}
		int name_end = line.find('\t');
		std::string name = line.substr(0, name_end);
		int name_id;
		if (name_to_num.find(name) == name_to_num.end()) {
			name_id = name_to_num.size() + 1;
			name_to_num[name] = name_id;
		} else {
			name_id = name_to_num[name];
		}
		output += std::to_string( (long long int) name_id) + line.substr(name_end) + "\n";
	}
	std::ofstream outfile;
        outfile.open(output_filename);
        outfile << output;
        outfile.close();
}

char* FastaToolbox::read_file_raw(char* filename) {
	std::ifstream myfile;
	myfile.open(filename);
	myfile.seekg(0, std::ios::end);
	int length = myfile.tellg();
	myfile.seekg(0, std::ios::beg);
	char* text = new char[length];
	myfile.read(text,length);
	myfile.close();
	return text;
}

void FastaToolbox::write_fasta(std::list<Protein> proteins, std::string file_name, bool write_numbered, bool write_domained) {
	std::string output = "";
	for (auto protein = proteins.begin(); protein != proteins.end(); protein++) {
		if (write_numbered) {
			output += ">" + std::to_string((long long int) protein -> name_id);
			if (write_domained) {
				output += "_" + std::to_string((long long int) protein -> domain_id);
			}
		} else {
			output += ">" + protein -> name;
		}
		output += "\n" + protein -> sequence + "\n";
	}
	std::ofstream my_file;
	my_file.open(file_name);
	my_file << output;
	my_file.close();
}

void FastaToolbox::read_fasta_serial(char* filename, std::list<Protein>* proteins) {
	std::ifstream my_file(filename);
	std::string line = "";
	std::string sequence = "";
	std::string name = "";
	std::stringstream seq_stream;
	while( getline(my_file, line) ) {
		if (line[0] == '>') {
			sequence = seq_stream.str();
			if (sequence != "") {
				Protein protein(name, sequence);
				proteins -> push_back(protein);
			}
			name = line.substr(1);
			seq_stream.str("");
			seq_stream.clear();
			sequence = "";
		} else {
			seq_stream << line;
		}
	}
	sequence = seq_stream.str();
	Protein protein(name, sequence);
	proteins -> push_back(protein);
}
