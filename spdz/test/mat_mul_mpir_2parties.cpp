// matrix mul test
#include "emp-tool/emp-tool.h"
#include "emp-ot/emp-ot.h"
#include "../network/network.h"
#include "../online/m_spdz_mpir.h"
#include "../mpir/mpir.h"

using namespace emp;
using namespace std;

int main(int argc, char **argv)
{
	int port, party;
	parse_party_and_port(argv, &party, &port);

	const static int nP = 2;
	NetIOMP<nP> io(party, port);
	NetIOMP<nP> io2(party, port + 2 * (nP + 1) * (nP + 1) + 1);
	NetIOMP<nP> *ios[2] = {&io, &io2};
	ThreadPool pool(4);

	int mat_sz = 128;
	int sz = mat_sz * mat_sz;
	mpz_t *x = new mpz_t[sz];
	mpz_t *y = new mpz_t[sz];
	mpz_t *output = new mpz_t[sz];
	mpz_t *mac_x = new mpz_t[mat_sz];
	mpz_t *mac_y = new mpz_t[mat_sz];
	mpz_t *output_mac = new mpz_t[mat_sz];

	std::ifstream fin1, fin2, fin3, fin4;
	fin1.open("/home/jackie/spdz/pre_data/predata_mspdz/128/input_x.txt", std::ios::in);
	fin2.open("/home/jackie/spdz/pre_data/predata_mspdz/128/input_y.txt", std::ios::in);
	fin3.open("/home/jackie/spdz/pre_data/predata_mspdz/128/mac_x.txt", std::ios::in);
	fin4.open("/home/jackie/spdz/pre_data/predata_mspdz/128/mac_y.txt", std::ios::in);

	if (!(fin1.is_open() && fin2.is_open() && fin3.is_open() && fin4.is_open()))
	{
		std::cerr << "cannot open the file";
	}

	//char line1[4096000] = {0};
	const size_t arraySize = 4096000*4*4;
    char* line1 = new(std::nothrow) char[arraySize]; 

    if (line1 == nullptr) {
        std::cerr << "Memory allocation failed" << std::endl;
        //return 1;
    }

	std::vector<triple> Triple;
	while (fin1.getline(line1, arraySize))
	{
		triple t;
		std::stringstream word(line1);
		word >> t.party;
		uint64_t num;
		while (word >> num)
			t.value.push_back(num);
		Triple.push_back(t);
	}

	for (int i = 0; i < nP; i++)
	{
		if (party == stoi(Triple[i].party))
		{
			for (int j = 0; j < sz; j++)
			{
				mpz_init_set_ui(x[j], Triple[i].value[j]);
			}
		}
	}

	Triple.clear();
	delete[] line1;
	fin1.close();
    char* line2 = new(std::nothrow) char[arraySize];

    if (line2 == nullptr) {
        std::cerr << "Memory allocation failed" << std::endl;
        //return 1;
    }

	while (fin2.getline(line2, arraySize))
	{
		triple t;
		std::stringstream word(line2);
		word >> t.party;
		uint64_t num;
		while (word >> num)
			t.value.push_back(num);
		Triple.push_back(t);
	}

	for (int i = 0; i < nP; i++)
	{
		if (party == stoi(Triple[i].party))
		{
			for (int j = 0; j < sz; j++)
			{
				mpz_init_set_ui(y[j], Triple[i].value[j]);
			}
		}
	}

	Triple.clear();
	delete[] line2;
	fin2.close();
    char* line3 = new(std::nothrow) char[arraySize];

	while (fin3.getline(line3, arraySize))
	{
		triple t;
		std::stringstream word(line3);
		word >> t.party;
		uint64_t num;
		while (word >> num)
			t.value.push_back(num);
		Triple.push_back(t);
	}

	for (int i = 0; i < nP; i++)
	{
		if (party == stoi(Triple[i].party))
		{
			for (int j = 0; j < mat_sz; j++)
			{
				mpz_init_set_ui(mac_x[j], Triple[i].value[j]);
			}
		}
	}

	Triple.clear();
	delete[] line3;
	fin3.close();
    char* line4 = new(std::nothrow) char[arraySize];

	while (fin4.getline(line4, arraySize))
	{
		triple t;
		std::stringstream word(line4);
		word >> t.party;
		uint64_t num;
		while (word >> num)
			t.value.push_back(num);
		Triple.push_back(t);
	}

	for (int i = 0; i < nP; i++)
	{
		if (party == stoi(Triple[i].party))
		{
			for (int j = 0; j < mat_sz; j++)
			{
				mpz_init_set_ui(mac_y[j], Triple[i].value[j]);
			}
		}
	}

	Triple.clear();
	delete[] line4;
	fin4.close();

	//delete[] line1;

	std::cout << std::endl
			  << "------------ MATRIX MUL ------------" << std::endl
			  << std::endl;
	;

	double sum_time = 0;
	for (int i = 0; i < 20; i++)
	{
		MSPDZ<nP> *mpc = new MSPDZ<nP>(ios, &pool, party);
		mpc->get_triple();

		auto start = clock_start();
		mpc->Online_mul(x, y, mac_x, mac_y, output, output_mac);
		mpc->mac_check(output, output_mac);

		sum_time += time_from(start) * 1000;
		delete mpc;
	}

	std::cout << "MATRIX MUL: " << sum_time / 20
			  << " ns" << std::endl;

	for (int i = 0; i < sz; i++)
	{
		mpz_clear(x[i]);
		mpz_clear(y[i]);
		mpz_clear(output[i]);
	}

	for (int i = 0; i < mat_sz; i++)
	{
		mpz_clear(mac_x[i]);
		mpz_clear(mac_y[i]);
		mpz_clear(output_mac[i]);
	}

	delete[] x;
	delete[] y;
	delete[] output;
	delete[] mac_x;
	delete[] mac_y;
	delete[] output_mac;

	return 0;
}
