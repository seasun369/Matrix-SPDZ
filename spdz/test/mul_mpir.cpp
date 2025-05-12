//one element mul test
#include "emp-tool/emp-tool.h"
#include "emp-ot/emp-ot.h"
#include "../network/network.h"
#include "../online/spdz_mpir.h"
#include "emp-tool/utils/ThreadPool.h"
#include "../mpir/mpir.h"

using namespace emp;
using namespace std;

int main(int argc, char **argv) {
  	int port, party;
	parse_party_and_port(argv, &party, &port);

	const static int nP = 2;
	NetIOMP<nP> io(party, port);
	NetIOMP<nP> io2(party, port+2*(nP+1)*(nP+1)+1);
	NetIOMP<nP> *ios[2] = {&io, &io2};
	std::cout << "To party : sent bytes = "
              << io.count() <<"B"<< std::endl;
	ThreadPool pool(4);	

  	mpz_t x,  y,  mac_x,  mac_y, output, output_mac;
	mpz_init(output);
	mpz_init(output_mac);

  	std::ifstream fin1,fin2,fin3,fin4;
	fin1.open("/home/jackie/spdz/pre_data/predata_spdz/input_x.txt",std::ios::in);
	fin2.open("/home/jackie/spdz/pre_data/predata_spdz/input_y.txt",std::ios::in);
	fin3.open("/home/jackie/spdz/pre_data/predata_spdz/mac_x.txt",std::ios::in);
	fin4.open("/home/jackie/spdz/pre_data/predata_spdz/mac_y.txt", std::ios::in);

  	if(!(fin1.is_open() && fin2.is_open() && fin3.is_open() && fin4.is_open()))
	{
		std::cerr<<"cannot open the file";
	}
	char line1[1024] = { 0 };
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
			mpz_init_set_ui(x, Triple[i].value[0]);
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
			mpz_init_set_ui(y, Triple[i].value[0]);
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
			mpz_init_set_ui(mac_x, Triple[i].value[0]);
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
			mpz_init_set_ui(mac_y, Triple[i].value[0]);
		}
	}

	Triple.clear();
	memset(line1, 0, sizeof(line1));

  	std::cout << std::endl
            << "------------ BASE MUL ------------" << std::endl
            << std::endl;
  	;
  	double sum_time = 0;
	
	for (int i = 0; i < 1; i++) {
		SPDZ<nP>* mpc = new SPDZ<nP>(ios, &pool, party);
		//cout <<"Setup:\t"<<party<<"\n";
		mpc->get_triple();
		//start counting
		auto start = clock_start();
		mpc->Online_mul(x, y, mac_x, mac_y, output, output_mac);
		mpc->mac_check(output, output_mac);
		sum_time += time_from(start) * 1000;

		delete mpc;
	}
	

	std::cout << "BASE MUL: " << sum_time / 1
		<< " ns" << std::endl;
		std::cout << "To party : sent bytes = "
              << io.count() <<"B"<< std::endl;
  	return 0;
}
