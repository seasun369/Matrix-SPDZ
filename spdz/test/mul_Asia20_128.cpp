#include "emp-tool/emp-tool.h"
#include "emp-ot/emp-ot.h"
#include "../network/network.h"
#include "../online/20Asiacrypt_128.h"
#include "emp-tool/utils/ThreadPool.h"

using namespace emp;
using namespace std;

int main(int argc, char **argv) {
  	int port, party;
	parse_party_and_port(argv, &party, &port);

	const static int nP = 2;
	NetIOMP<nP> io(party, port);
	NetIOMP<nP> io2(party, port+2*(nP+1)*(nP+1)+1);
	NetIOMP<nP> *ios[2] = {&io, &io2};
	

  	uint64_t *x, *y, *mac_x, *mac_y, *output, *output_mac;
  	int num_mul=16384;
  	x = new uint64_t[num_mul];
  	mac_x = new uint64_t[num_mul];
  	y = new uint64_t[num_mul];
  	mac_y = new uint64_t[num_mul];
 	output = new uint64_t[num_mul];
  	output_mac = new uint64_t[num_mul];

	ThreadPool pool(4);	
  	std::ifstream fin1,fin2,fin3,fin4;
	fin1.open("/home/jackie/spdz/pre_data/sextuple_Asiacrypt20/128/input_x.txt",std::ios::in);
	fin2.open("/home/jackie/spdz/pre_data/sextuple_Asiacrypt20/128/input_y.txt",std::ios::in);
	fin3.open("/home/jackie/spdz/pre_data/sextuple_Asiacrypt20/128/mac_x.txt",std::ios::in);
	fin4.open("/home/jackie/spdz/pre_data/sextuple_Asiacrypt20/128/mac_y.txt",std::ios::in);

  	if(!(fin1.is_open() && fin2.is_open() && fin3.is_open() && fin4.is_open()))
	{
		std::cerr<<"cannot open the file";
	}
	
		char line1[1024000] = { 0 };
		std::vector<triple> Triple;
		while (fin1.getline(line1, sizeof(line1)))
		{
			triple t;
			std::stringstream word(line1);
			word >> t.party;
			uint64_t num;
			while (word >> num)
				t.value.push_back(num);
			Triple.push_back(t);
		}

		for (int i = 0; i < nP; i++) {
			if (party == stoi(Triple[i].party)) {
				for (int j = 0; j < num_mul; j++) x[j] = Triple[i].value[j];
			}
		}

		Triple.clear();
		memset(line1, 0, sizeof(line1));

		while (fin2.getline(line1, sizeof(line1)))
		{
			triple t;
			std::stringstream word(line1);
			word >> t.party;
			uint64_t num;
			while (word >> num)
				t.value.push_back(num);
			Triple.push_back(t);
		}

		for (int i = 0; i < nP; i++) {
			if (party == stoi(Triple[i].party)) {
				for (int j = 0; j < num_mul; j++) y[j] = Triple[i].value[j];
			}
		}

		Triple.clear();
		memset(line1, 0, sizeof(line1));

		while (fin3.getline(line1, sizeof(line1)))
		{
			triple t;
			std::stringstream word(line1);
			word >> t.party;
			uint64_t num;
			while (word >> num)
				t.value.push_back(num);
			Triple.push_back(t);
		}

		for (int i = 0; i < nP; i++) {
			if (party == stoi(Triple[i].party)) {
				for (int j = 0; j < num_mul; j++) mac_x[j] = Triple[i].value[j];
			}
		}

		Triple.clear();
		memset(line1, 0, sizeof(line1));

		while (fin4.getline(line1, sizeof(line1)))
		{
			triple t;
			std::stringstream word(line1);
			word >> t.party;
			uint64_t num;
			while (word >> num)
				t.value.push_back(num);
			Triple.push_back(t);
		}

		for (int i = 0; i < nP; i++) {
			if (party == stoi(Triple[i].party)) {
				for (int j = 0; j < num_mul; j++) mac_y[j] = Triple[i].value[j];
			}
		}

		Triple.clear();
		memset(line1, 0, sizeof(line1));

  	std::cout << std::endl
            << "------------ BASE MUL ------------" << std::endl
            << std::endl;
  	;

	double sum_time = 0;
	for(int i=0;i<1;i++){
		Asiacrypt20_SPDZ<nP>* mpc = new Asiacrypt20_SPDZ<nP>(ios, &pool, party);
  		//cout <<"Setup:\t"<<party<<"\n";
  		mpc->get_triple();
  
  		//start counting
  		auto start = clock_start();
  		mpc->Online_mul(x,y,mac_x,mac_y,output,output_mac);
  		mpc->mac_check(output,output_mac);
  		sum_time += time_from(start) * 1000;

  		delete mpc;
	}

	std::cout << "BASE MUL: " << sum_time/1
                << " ns" << std::endl;
  	
  	return 0;
}