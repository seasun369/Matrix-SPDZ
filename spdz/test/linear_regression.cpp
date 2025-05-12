// matrix mul test
#include "emp-tool/emp-tool.h"
#include "emp-ot/emp-ot.h"
#include "../network/network.h"
#include "../online/linear_regression.h"
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

	int m = 1913;
	int n = 14;
	int l = 1;
	int key_length = std::max({m, n, l});
	int max_m_n = std::max({m, n});
	mpz_t *x = new mpz_t[m * n];
	mpz_t *y = new mpz_t[m * l];
	mpz_t *w = new mpz_t[n * l];
	mpz_t *output = new mpz_t[m * l];
	mpz_t *mac_x = new mpz_t[m];
	mpz_t *mac_y = new mpz_t[m];
	mpz_t *output_mac = new mpz_t[m];
	mpz_t *mac_w = new mpz_t[n];

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

	std::ifstream fin1, fin2, fin3, fin4;
	fin1.open("/home/jackie/spdz/pre_data/predata_mspdz/nonsqure/input_x.txt", std::ios::in);
	fin2.open("/home/jackie/spdz/pre_data/predata_mspdz/nonsqure/input_y.txt", std::ios::in);
	fin3.open("/home/jackie/spdz/pre_data/predata_mspdz/nonsqure/mac_x.txt", std::ios::in);
	fin4.open("/home/jackie/spdz/pre_data/predata_mspdz/nonsqure/mac_y.txt", std::ios::in);

	if (!(fin1.is_open() && fin2.is_open() && fin3.is_open() && fin4.is_open()))
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

	// delete[] line1;

	std::cout << std::endl
			  << "------------ MATRIX MUL ------------" << std::endl
			  << std::endl;
	;

	double sum_time = 0;
	for (int i = 0; i < 3; i++)
	{
		MSPDZ<nP> *mpc = new MSPDZ<nP>(ios, &pool, party, m, n, l, key_length, max_m_n);
		mpc->m = m;
		mpc->n = n;
		mpc->l = l;
		mpc->key_length = key_length;
		mpc->max_m_n = max_m_n;

		// mpc->get_triple(m, n, l, key_length);
		// mpc->mac_check(x, mac_x, m, n);
		// mpc->mac_check(w, mac_w, n, l);
		auto start = clock_start();

		// mpc->Online_nonsqure_mul(x, w, mac_x, mac_w, output, output_mac, m, n, l, key_length);

		// mpc->mac_check(output, output_mac, m, l);
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

		mpc->Gradient_Descent_linear_regression(x, w, y, mac_x, mac_w, mac_y, m, n, l, key_length);

		sum_time += time_from(start) * 1000;
		delete mpc;
	}

	std::cout << "MATRIX MUL: " << sum_time / 3
			  << " ns" << std::endl;

	for (int i = 0; i < m * n; i++)
	{
		mpz_clear(x[i]);
	}

	for (int i = 0; i < m * l; i++)
	{
		mpz_clear(y[i]);
		mpz_clear(output[i]);
	}

	for (int i = 0; i < n * l; i++)
	{
		mpz_clear(w[i]);
	}

	for (int i = 0; i < m; i++)
	{
		mpz_clear(mac_x[i]);
		mpz_clear(output_mac[i]);
		mpz_clear(mac_y[i]);
	}

	for (int i = 0; i < n; i++)
	{
		mpz_clear(mac_w[i]);
	}

	delete[] x;
	delete[] y;
	delete[] output;
	delete[] mac_x;
	delete[] mac_y;
	delete[] output_mac;
	delete[] w;
	delete[] mac_w;

	return 0;
}
