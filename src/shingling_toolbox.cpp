#include "shingling_toolbox.h"

ShinglingToolbox::ShinglingToolbox(int prime_num, int num_hash, int shingle_size, int kmer_len, char* random_filename, int max_hash, 
		bool use_persistent_hashing) :
		Prime_num(prime_num), num_hash(num_hash), shingle_size(shingle_size), kmer_len(kmer_len), 
		max_hash(max_hash), use_persistent_hashing(use_persistent_hashing) {
	hash_randoms = new int[max_hash * 2];
	set_hash_randoms(hash_randoms, max_hash, random_filename);
}

ShinglingToolbox::~ShinglingToolbox() {
	if (hash_randoms) delete[] hash_randoms;
}

void ShinglingToolbox::set_hash_randoms(int* rand_array, int max_hash, char* filename) {
	std::string line;
	std::ifstream imyfile (filename);
	int required = max_hash * 2;
	if (imyfile.is_open())
 	{
		//cout<<"ifile open!" << std::endl;
		int counter = 0;
    		while ( getline (imyfile,line) && counter < required)
    		{
			rand_array[counter] = stoi(line);
			counter++;
      			//cout << line << '\t' << stoi(line) << std::endl;
    		}
    		imyfile.close();
  	}
	//std::cerr << "chkmarands: " << rand_array[0] << " " << rand_array[1] << std::endl;
}

