#pragma once
#include <string>
#include <string.h>
#include <list>
#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>

class Protein {
public:
	std::string name;
	int name_id;
	int domain_id;
	std::string sequence;
	Protein(std::string prot_name, std::string prot_seq) : name(prot_name), sequence(prot_seq) { };
};

class FastaToolbox {
	
public:
	std::unordered_map<std::string, int> name_to_num;
	std::list<Protein> read_fasta(char* text, long long int text_size = -1, int proc = 0);
	void convert_name_to_num(std::list<Protein> *proteins, bool domained);
	std::list<Protein> convert_num_to_name(std::list<Protein> numbered_proteins, bool domained);
	void write_fasta(std::list<Protein> proteins, std::string file_name, bool write_numbered, bool write_domained);
	void write_conversion_table(std::string file_name);
	void read_conversion_table(std::string file_name);
	void convert_pfam_file_to_numbered(std::string input_filename, std::string output_filename);
	std::string parse_protein_name(std::string name, int* domain_id);
private:
	void parse_fasta_sequence(char* fasta_pack, char* seq_name, char* sequence, int name_length, int proc); 
	char* read_file_raw(char* file_name);
	void read_fasta_serial(char* filename, std::list<Protein>* proteins);
};
