// matrix mul test
#include "emp-tool/emp-tool.h"
#include "emp-ot/emp-ot.h"
#include "../network/network.h"
#include "../online/DNN.h"
#include "../mpir/mpir.h"
#include <algorithm>
#include <omp.h>

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

	int m = 6000;
	int n = 784;
	int l = 1;
	int inputnodes = 784;
	int hiddennodes = 128;
	int outputnodes = 10;
	int epoch = 10;
	int batch_size = 128;
	int batch_size_bit = 7;
	int learningrate_DNN_bit = 3;
	int key_length = std::max({inputnodes, batch_size, m, n, l});
	int max_m_n = std::max({m, n});
	mpz_t *x = new mpz_t[m * n];
	mpz_t *y = new mpz_t[m * l];
	mpz_t *w1 = new mpz_t[inputnodes * hiddennodes];
	mpz_t *w2 = new mpz_t[hiddennodes * outputnodes];
	mpz_t *b1 = new mpz_t[hiddennodes];
	mpz_t *b2 = new mpz_t[outputnodes];
	mpz_t *output = new mpz_t[m * l];
	mpz_t *mac_x = new mpz_t[m];
	mpz_t *mac_y = new mpz_t[m];
	mpz_t *output_mac = new mpz_t[m];
	mpz_t *mac_w1 = new mpz_t[inputnodes];
	mpz_t *mac_w2 = new mpz_t[hiddennodes];
	mpz_t *mac_b1 = new mpz_t[1];
	mpz_t *mac_b2 = new mpz_t[1];

	/*
	for (int i = 0; i < n * l; i++)
	{
		mpz_init(w[i]);
		mpz_set_ui(w[i], 0);
	}

	for (int i = 0; i < n * l; i++)
	{
		mpz_init(mac_w[i]);
		mpz_set_ui(mac_w[i], 0);
	}
	*/

	std::ifstream fin1, fin2, fin3, fin4, fin5, fin6, fin7, fin8;
	fin1.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/input_x.txt", std::ios::in);
	fin2.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/input_y.txt", std::ios::in);
	fin3.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/mac_x.txt", std::ios::in);
	fin4.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/mac_y.txt", std::ios::in);
	fin5.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/w1.txt", std::ios::in);
	fin6.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/w2.txt", std::ios::in);
	fin7.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/mac_w1.txt", std::ios::in);
	fin8.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/mac_w2.txt", std::ios::in);

	if (!(fin1.is_open() && fin2.is_open() && fin3.is_open() && fin4.is_open() && fin5.is_open() && fin6.is_open() && fin7.is_open() && fin8.is_open()))
	{
		std::cerr << "cannot open the file";
	}

	// char line1[4096000] = {0};
	const size_t arraySize = 4096000 * 4 * 4;
	char *line1 = new (std::nothrow) char[arraySize];

	if (line1 == nullptr)
	{
		std::cerr << "Memory allocation failed" << std::endl;
		// return 1;
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
			for (int j = 0; j < m * n; j++)
			{
				mpz_init_set_ui(x[j], Triple[i].value[j]);
			}
		}
	}

	Triple.clear();
	delete[] line1;
	fin1.close();
	char *line2 = new (std::nothrow) char[arraySize];

	if (line2 == nullptr)
	{
		std::cerr << "Memory allocation failed" << std::endl;
		// return 1;
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
			for (int j = 0; j < m * l; j++)
			{
				mpz_init_set_ui(y[j], Triple[i].value[j]);
			}
		}
	}

	Triple.clear();
	delete[] line2;
	fin2.close();
	char *line3 = new (std::nothrow) char[arraySize];

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
			for (int j = 0; j < m; j++)
			{
				mpz_init_set_ui(mac_x[j], Triple[i].value[j]);
			}
		}
	}

	Triple.clear();
	delete[] line3;
	fin3.close();
	char *line4 = new (std::nothrow) char[arraySize];

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
			for (int j = 0; j < m; j++)
			{
				mpz_init_set_ui(mac_y[j], Triple[i].value[j]);
			}
		}
	}

	Triple.clear();
	delete[] line4;
	fin4.close();
	char *line5 = new (std::nothrow) char[arraySize];

	while (fin5.getline(line5, arraySize))
	{
		triple t;
		std::stringstream word(line5);
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
			for (int j = 0; j < inputnodes * hiddennodes; j++)
			{
				mpz_init_set_ui(w1[j], Triple[i].value[j]);
			}
		}
	}

	Triple.clear();
	delete[] line5;
	fin5.close();
	char *line6 = new (std::nothrow) char[arraySize];

	while (fin6.getline(line6, arraySize))
	{
		triple t;
		std::stringstream word(line6);
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
			for (int j = 0; j < hiddennodes * outputnodes; j++)
			{
				mpz_init_set_ui(w2[j], Triple[i].value[j]);
			}
		}
	}

	Triple.clear();
	delete[] line6;
	fin6.close();
	char *line7 = new (std::nothrow) char[arraySize];

	while (fin7.getline(line7, arraySize))
	{
		triple t;
		std::stringstream word(line7);
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
			for (int j = 0; j < inputnodes; j++)
			{
				mpz_init_set_ui(mac_w1[j], Triple[i].value[j]);
			}
		}
	}

	Triple.clear();
	delete[] line7;
	fin7.close();
	char *line8 = new (std::nothrow) char[arraySize];

	while (fin8.getline(line8, arraySize))
	{
		triple t;
		std::stringstream word(line8);
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
			for (int j = 0; j < hiddennodes; j++)
			{
				mpz_init_set_ui(mac_w2[j], Triple[i].value[j]);
			}
		}
	}

	Triple.clear();
	delete[] line8;
	fin8.close();

	// delete[] line1;

	std::cout << std::endl
			  << "------------ MATRIX MUL ------------" << std::endl
			  << std::endl;
	;

	double sum_time = 0;
	for (int i = 0; i < 1; i++)
	{
		MSPDZ<nP> *mpc = new MSPDZ<nP>(ios, &pool, party, m, n, l, inputnodes, hiddennodes, outputnodes, batch_size, key_length, max_m_n); // revise
		std::cout << "To party : sent bytes = "
				  << io.count() << "B" << std::endl;
		mpc->m = m;
		mpc->n = n;
		mpc->l = l;
		mpc->key_length = key_length;
		mpc->max_m_n = max_m_n;
		mpc->inputnodes = inputnodes;
		mpc->hiddennodes = hiddennodes;
		mpc->outputnodes = outputnodes;
		mpc->epoch = epoch;
		mpc->batch_size = batch_size; // batch size of DNN
		mpc->batch_size_bit = batch_size_bit;
		mpc->learningrate_DNN_bit = learningrate_DNN_bit;

		// mpc->get_triple(m, n, l, key_length);
		// mpc->mac_check(x, mac_x, m, n);
		// mpc->mac_check(w, mac_w, n, l);
		auto start = clock_start();

		// mpc->Online_nonsqure_mul(x, w, mac_x, mac_w, output, output_mac, m, n, l, key_length);

		// mpc->mac_check(output, output_mac, m, l);
		for (int i = 0; i < hiddennodes; i++)
		{
			mpz_init_set_ui(b1[i], 0);
		}

		for (int i = 0; i < outputnodes; i++)
		{
			mpz_init_set_ui(b2[i], 0);
		}

		mpz_init_set_ui(mac_b1[0], 0);
		mpz_init_set_ui(mac_b2[0], 0);

		mpc->Gradient_Descent_DNN(x, w1, b1, w2, b2, y, mac_x, mac_w1, mac_b1, mac_w2, mac_b2, mac_y, m, key_length);

		sum_time += time_from(start) * 1000;
		delete mpc;
	}

	std::cout << "MATRIX MUL: " << sum_time / 1
			  << " ns" << std::endl;
	std::cout << "To party : sent bytes = "
			  << io.count() << "B" << std::endl;

	for (int i = 0; i < m * n; i++)
	{
		mpz_clear(x[i]);
	}

	for (int i = 0; i < m * l; i++)
	{
		mpz_clear(y[i]);
		mpz_clear(output[i]);
	}

	for (int i = 0; i < inputnodes * hiddennodes; i++)
	{
		mpz_clear(w1[i]);
	}

	for (int i = 0; i < hiddennodes * outputnodes; i++)
	{
		mpz_clear(w2[i]);
	}

	for (int i = 0; i < m; i++)
	{
		mpz_clear(mac_x[i]);
		mpz_clear(output_mac[i]);
		mpz_clear(mac_y[i]);
	}

	for (int i = 0; i < inputnodes; i++)
	{
		mpz_clear(mac_w1[i]);
	}

	for (int i = 0; i < hiddennodes; i++)
	{
		mpz_clear(mac_w2[i]);
		mpz_clear(b1[i]);
	}

	for (int i = 0; i < outputnodes; i++)
	{
		mpz_clear(b2[i]);
	}

	delete[] x;
	delete[] y;
	delete[] output;
	delete[] mac_x;
	delete[] mac_y;
	delete[] output_mac;
	delete[] w1;
	delete[] w2;
	delete[] b1;
	delete[] b2;
	delete[] mac_w1;
	delete[] mac_w2;
	delete[] mac_b1;
	delete[] mac_b2;

	return 0;
}