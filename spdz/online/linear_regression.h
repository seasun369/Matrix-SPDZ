#ifndef MSPDZ_H__
#define MSPDZ_H__
#define NUM_THREADS 8

#include "../network/network.h"
#include <emp-tool/emp-tool.h>
#include "utility.h"
#include "../network/helper.h"
#include <omp.h>
#include "../mpir/mpir.h"
#include <algorithm>
#include <vector>
using namespace emp;

// Matrix SPDZ
template <int nP>
class MSPDZ
{
public:
	mpz_t *key;

	mpz_t *a, *a_t, *b, *c, *r, *r_t, *s, *s_bit, *t;
	mpz_t *mac_a, *mac_a_t, *mac_b, *mac_c, *mac_r, *mac_r_t, *mac_s_spdz, *mac_s_bit_spdz, *mac_t, *mac_t_spdz;
	vector<vector<int>> new_matrix = {
		{1, 2, 3},
		{4, 5, 6},
		{7, 8, 9}};

	NetIOMP<nP> *io;
	int party, total_pre, ssp;

	int m;
	int n;
	int l;
	int key_length;
	int max_m_n;
	int precision = 15;
	int l_bit = 61;
	ThreadPool *pool;

	PRP prp;
	MSPDZ(NetIOMP<nP> *io[2], ThreadPool *pool, int party, int m, int n, int l, int key_length, int max_m_n, bool *_delta = nullptr, int ssp = 40)
	{
		this->party = party;
		this->io = io[0];
		this->ssp = ssp;
		this->pool = pool;

		// A more flexible strategy is needed for setting the array size
		a = new mpz_t[m * n];
		a_t = new mpz_t[n * m];
		b = new mpz_t[m * n];
		c = new mpz_t[m * l];
		r = new mpz_t[m * l];
		r_t = new mpz_t[m * l];
		s = new mpz_t[m * n];
		s_bit = new mpz_t[m * n * l_bit];
		t = new mpz_t[m * n];

		key = new mpz_t[key_length];
		mac_a = new mpz_t[max_m_n];
		mac_a_t = new mpz_t[max_m_n];
		mac_b = new mpz_t[max_m_n];
		mac_c = new mpz_t[max_m_n];
		mac_r = new mpz_t[max_m_n];
		mac_r_t = new mpz_t[max_m_n];
		mac_s_spdz = new mpz_t[m * n];
		mac_s_bit_spdz = new mpz_t[m * n * l_bit];
		mac_t = new mpz_t[max_m_n];
		mac_t_spdz = new mpz_t[m * n];

		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < m * n; i++)
		{
			mpz_init(a[i]);
			mpz_init(a_t[i]);
			mpz_init(b[i]);
			mpz_init(s[i]);
			mpz_init(t[i]);
			mpz_init(mac_t_spdz[i]);
			mpz_init(mac_s_spdz[i]);
		}
		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < m * l; i++)
		{
			mpz_init(c[i]);
			mpz_init(r[i]);
			mpz_init(r_t[i]);
		}
		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < m * n * l_bit; i++)
		{
			mpz_init(s_bit[i]);
			mpz_init(mac_s_bit_spdz[i]);
		}
		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < key_length; i++)
		{
			mpz_init(key[i]);
		}
		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < max_m_n; i++)
		{
			mpz_init(mac_a[i]);
			mpz_init(mac_a_t[i]);
			mpz_init(mac_b[i]);
			mpz_init(mac_c[i]);
			mpz_init(mac_r[i]);
			mpz_init(mac_r_t[i]);
			mpz_init(mac_t[i]);
		}
	}
	~MSPDZ()
	{
		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < m * n; i++)
		{
			mpz_clear(a[i]);
			mpz_clear(a_t[i]);
			mpz_clear(b[i]);
			mpz_clear(s[i]);
			mpz_clear(t[i]);
			mpz_clear(mac_t_spdz[i]);
			mpz_clear(mac_s_spdz[i]);
		}
		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < m * l; i++)
		{
			mpz_clear(c[i]);
			mpz_clear(r[i]);
			mpz_clear(r_t[i]);
		}
		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < m * n * l_bit; i++)
		{
			mpz_clear(s_bit[i]);
			mpz_clear(mac_s_bit_spdz[i]);
		}
		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < key_length; i++)
		{
			mpz_clear(key[i]);
		}
		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < max_m_n; i++)
		{
			mpz_clear(mac_a[i]);
			mpz_clear(mac_a_t[i]);
			mpz_clear(mac_b[i]);
			mpz_clear(mac_c[i]);
			mpz_clear(mac_r[i]);
			mpz_clear(mac_r_t[i]);
			mpz_clear(mac_t[i]);
		}

		delete[] a;
		delete[] a_t;
		delete[] b;
		delete[] c;
		delete[] r;
		delete[] r_t;
		delete[] s;
		delete[] s_bit;
		delete[] t;
		delete[] key;
		delete[] mac_a;
		delete[] mac_a_t;
		delete[] mac_b;
		delete[] mac_c;
		delete[] mac_r;
		delete[] mac_r_t;
		delete[] mac_s_spdz;
		delete[] mac_s_bit_spdz;
		delete[] mac_t;
		delete[] mac_t_spdz;
	}
	PRG prg;

	// it should be implemented by offline.
	void get_transform_tuple(int l1, int l3, int key_length)
	{
		std::ifstream fin1;
		std::ifstream fin_1, fin_2;
		if ((l1 == 14) && (l3 == 14))
		{
			fin1.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/transform_predata/1/t.txt", std::ios::in);
			fin_1.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/transform_predata/1/mac_t.txt", std::ios::in);
			fin_2.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/transform_predata/1/mac_t_spdz.txt", std::ios::in);
		}
		else if ((l1 == 14) && (l3 == 1))
		{
			fin1.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/transform_predata/2/t.txt", std::ios::in);
			fin_1.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/transform_predata/2/mac_t.txt", std::ios::in);
			fin_2.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/transform_predata/2/mac_t_spdz.txt", std::ios::in);
		}
		else if ((l1 == 1913) && (l3 == 1))
		{
			fin1.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/transform_predata/4/t.txt", std::ios::in);
			fin_1.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/transform_predata/4/mac_t.txt", std::ios::in);
			fin_2.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/transform_predata/4/mac_t_spdz.txt", std::ios::in);
		}

		if (!(fin1.is_open() && fin_1.is_open() && fin_2.is_open()))
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
				for (int j = 0; j < l1 * l3; j++)
				{
					mpz_set_ui(t[j], Triple[i].value[j]);
				}
			}
		}

		Triple.clear();
		delete[] line1;
		fin1.close();
		char *line_1 = new (std::nothrow) char[arraySize];

		while (fin_1.getline(line_1, arraySize))
		{
			triple t;
			std::stringstream word(line_1);
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
				for (int j = 0; j < l1; j++)
					mpz_set_ui(mac_t[j], Triple[i].value[j]);
			}
		}

		Triple.clear();
		delete[] line_1;
		fin_1.close();
		char *line_2 = new (std::nothrow) char[arraySize];

		while (fin_2.getline(line_2, arraySize))
		{
			triple t;
			std::stringstream word(line_2);
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
				for (int j = 0; j < l1 * l3; j++)
					mpz_set_ui(mac_t_spdz[j], Triple[i].value[j]);
			}
		}

		Triple.clear();
		delete[] line_2;
		fin_2.close();
	}

	void get_truncate_tuple(int l1, int l3, int key_length)
	{
		std::ifstream fin1, fin2;
		std::ifstream fin_1, fin_2;
		if ((l1 == 14) && (l3 == 14))
		{
			fin1.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/truncate_predata/1/s.txt", std::ios::in);
			fin2.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/truncate_predata/1/s_bit.txt", std::ios::in);
			fin_1.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/truncate_predata/1/mac_s_spdz.txt", std::ios::in);
			fin_2.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/truncate_predata/1/mac_s_bit_spdz.txt", std::ios::in);
		}
		else if ((l1 == 14) && (l3 == 1))
		{
			fin1.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/truncate_predata/2/s.txt", std::ios::in);
			fin2.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/truncate_predata/2/s_bit.txt", std::ios::in);
			fin_1.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/truncate_predata/2/mac_s_spdz.txt", std::ios::in);
			fin_2.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/truncate_predata/2/mac_s_bit_spdz.txt", std::ios::in);
		}
		else if ((l1 == 1913) && (l3 == 1))
		{
			fin1.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/truncate_predata/3/s.txt", std::ios::in);
			fin2.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/truncate_predata/3/s_bit.txt", std::ios::in);
			fin_1.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/truncate_predata/3/mac_s_spdz.txt", std::ios::in);
			fin_2.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/truncate_predata/3/mac_s_bit_spdz.txt", std::ios::in);
		}

		if (!(fin1.is_open() && fin2.is_open() && fin_1.is_open() && fin_2.is_open()))
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
				for (int j = 0; j < l1 * l3; j++)
				{
					mpz_set_ui(s[j], Triple[i].value[j]);
				}
			}
		}

		Triple.clear();
		delete[] line1;
		fin1.close();
		char *line2 = new (std::nothrow) char[arraySize];

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
				for (int j = 0; j < l1 * l3 * l_bit; j++)
					mpz_set_ui(s_bit[j], Triple[i].value[j]);
			}
		}

		Triple.clear();
		delete[] line2;
		fin2.close();
		char *line_1 = new (std::nothrow) char[arraySize];

		while (fin_1.getline(line_1, arraySize))
		{
			triple t;
			std::stringstream word(line_1);
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
				for (int j = 0; j < l1 * l3; j++)
					mpz_set_ui(mac_s_spdz[j], Triple[i].value[j]);
			}
		}

		Triple.clear();
		delete[] line_1;
		fin_1.close();
		char *line_2 = new (std::nothrow) char[arraySize];

		while (fin_2.getline(line_2, arraySize))
		{
			triple t;
			std::stringstream word(line_2);
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
				for (int j = 0; j < l1 * l3 * l_bit; j++)
					mpz_set_ui(mac_s_bit_spdz[j], Triple[i].value[j]);
			}
		}

		Triple.clear();
		delete[] line_2;
		fin_2.close();
	}

	void get_triple(int l1, int l2, int l3, int key_length)
	{
		std::ifstream fin1, fin2, fin3, fin4, fin5, fin6;
		std::ifstream fin_1, fin_2, fin_3, fin_4, fin_5, fin_6;
		if ((l1 == 14) && (l2 == 1913) && (l3 == 14))
		{
			fin1.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/triple_predata/1/a.txt", std::ios::in);
			fin2.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/triple_predata/1/a_t.txt", std::ios::in);
			fin3.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/triple_predata/1/b.txt", std::ios::in);
			fin4.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/triple_predata/1/c.txt", std::ios::in);
			fin5.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/triple_predata/1/r.txt", std::ios::in);
			fin6.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/triple_predata/1/r_t.txt", std::ios::in);
			fin_1.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/triple_predata/1/mac_a.txt", std::ios::in);
			fin_2.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/triple_predata/1/mac_a_t.txt", std::ios::in);
			fin_3.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/triple_predata/1/mac_b.txt", std::ios::in);
			fin_4.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/triple_predata/1/mac_c.txt", std::ios::in);
			fin_5.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/triple_predata/1/mac_r.txt", std::ios::in);
			fin_6.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/triple_predata/1/mac_r_t.txt", std::ios::in);
		}
		else if ((l1 == 14) && (l2 == 1913) && (l3 == 1))
		{
			fin1.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/triple_predata/2/a.txt", std::ios::in);
			fin2.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/triple_predata/2/a_t.txt", std::ios::in);
			fin3.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/triple_predata/2/b.txt", std::ios::in);
			fin4.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/triple_predata/2/c.txt", std::ios::in);
			fin5.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/triple_predata/2/r.txt", std::ios::in);
			fin6.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/triple_predata/2/r_t.txt", std::ios::in);
			fin_1.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/triple_predata/2/mac_a.txt", std::ios::in);
			fin_2.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/triple_predata/2/mac_a_t.txt", std::ios::in);
			fin_3.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/triple_predata/2/mac_b.txt", std::ios::in);
			fin_4.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/triple_predata/2/mac_c.txt", std::ios::in);
			fin_5.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/triple_predata/2/mac_r.txt", std::ios::in);
			fin_6.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/triple_predata/2/mac_r_t.txt", std::ios::in);
		}
		else if ((l1 == 14) && (l2 == 14) && (l3 == 1))
		{
			fin1.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/triple_predata/3/a.txt", std::ios::in);
			fin2.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/triple_predata/3/a_t.txt", std::ios::in);
			fin3.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/triple_predata/3/b.txt", std::ios::in);
			fin4.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/triple_predata/3/c.txt", std::ios::in);
			fin5.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/triple_predata/3/r.txt", std::ios::in);
			fin6.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/triple_predata/3/r_t.txt", std::ios::in);
			fin_1.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/triple_predata/3/mac_a.txt", std::ios::in);
			fin_2.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/triple_predata/3/mac_a_t.txt", std::ios::in);
			fin_3.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/triple_predata/3/mac_b.txt", std::ios::in);
			fin_4.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/triple_predata/3/mac_c.txt", std::ios::in);
			fin_5.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/triple_predata/3/mac_r.txt", std::ios::in);
			fin_6.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/triple_predata/3/mac_r_t.txt", std::ios::in);
		}
		else if ((l1 == 1913) && (l2 == 14) && (l3 == 1))
		{
			fin1.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/triple_predata/4/a.txt", std::ios::in);
			fin2.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/triple_predata/4/a_t.txt", std::ios::in);
			fin3.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/triple_predata/4/b.txt", std::ios::in);
			fin4.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/triple_predata/4/c.txt", std::ios::in);
			fin5.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/triple_predata/4/r.txt", std::ios::in);
			fin6.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/triple_predata/4/r_t.txt", std::ios::in);
			fin_1.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/triple_predata/4/mac_a.txt", std::ios::in);
			fin_2.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/triple_predata/4/mac_a_t.txt", std::ios::in);
			fin_3.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/triple_predata/4/mac_b.txt", std::ios::in);
			fin_4.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/triple_predata/4/mac_c.txt", std::ios::in);
			fin_5.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/triple_predata/4/mac_r.txt", std::ios::in);
			fin_6.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/triple_predata/4/mac_r_t.txt", std::ios::in);
		}

		if (!(fin1.is_open() && fin2.is_open() && fin3.is_open() && fin4.is_open() && fin5.is_open() && fin6.is_open() && fin_1.is_open() && fin_2.is_open() && fin_3.is_open() && fin_4.is_open() && fin_5.is_open() && fin_6.is_open()))
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
				for (int j = 0; j < l1 * l2; j++)
				{
					mpz_set_ui(a[j], Triple[i].value[j]);
				}
			}
		}

		Triple.clear();
		delete[] line1;
		fin1.close();
		char *line2 = new (std::nothrow) char[arraySize];

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
				for (int j = 0; j < l2 * l1; j++)
					mpz_set_ui(a_t[j], Triple[i].value[j]);
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
				for (int j = 0; j < l2 * l3; j++)
					mpz_set_ui(b[j], Triple[i].value[j]);
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
				for (int j = 0; j < l1 * l3; j++)
					mpz_set_ui(c[j], Triple[i].value[j]);
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
				for (int j = 0; j < l1 * l3; j++)
					mpz_set_ui(r[j], Triple[i].value[j]);
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
				for (int j = 0; j < l3 * l1; j++)
					mpz_set_ui(r_t[j], Triple[i].value[j]);
			}
		}

		Triple.clear();
		delete[] line6;
		fin6.close();

		char *line_1 = new (std::nothrow) char[arraySize];

		while (fin_1.getline(line_1, arraySize))
		{
			triple t;
			std::stringstream word(line_1);
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
				for (int j = 0; j < l1; j++)
					mpz_set_ui(mac_a[j], Triple[i].value[j]);
			}
		}

		Triple.clear();
		delete[] line_1;
		fin_1.close();
		char *line_2 = new (std::nothrow) char[arraySize];

		while (fin_2.getline(line_2, arraySize))
		{
			triple t;
			std::stringstream word(line_2);
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
				for (int j = 0; j < l2; j++)
					mpz_set_ui(mac_a_t[j], Triple[i].value[j]);
			}
		}

		Triple.clear();
		delete[] line_2;
		fin_2.close();
		char *line_3 = new (std::nothrow) char[arraySize];

		while (fin_3.getline(line_3, arraySize))
		{
			triple t;
			std::stringstream word(line_3);
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
				for (int j = 0; j < l2; j++)
					mpz_set_ui(mac_b[j], Triple[i].value[j]);
			}
		}

		Triple.clear();
		delete[] line_3;
		fin_3.close();
		char *line_4 = new (std::nothrow) char[arraySize];

		while (fin_4.getline(line1, arraySize))
		{
			triple t;
			std::stringstream word(line_4);
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
				for (int j = 0; j < l1; j++)
					mpz_set_ui(mac_c[j], Triple[i].value[j]);
			}
		}

		Triple.clear();
		delete[] line_4;
		fin_4.close();
		char *line_5 = new (std::nothrow) char[arraySize];

		while (fin_5.getline(line_5, arraySize))
		{
			triple t;
			std::stringstream word(line_5);
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
				for (int j = 0; j < l1; j++)
					mpz_set_ui(mac_r[j], Triple[i].value[j]);
			}
		}

		Triple.clear();
		delete[] line_5;
		fin_5.close();
		char *line_6 = new (std::nothrow) char[arraySize];

		while (fin_6.getline(line_6, arraySize))
		{
			triple t;
			std::stringstream word(line_6);
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
				for (int j = 0; j < l3; j++)
					mpz_set_ui(mac_r_t[j], Triple[i].value[j]);
			}
		}

		Triple.clear();
		delete[] line_6;
		fin_6.close();
	}

	void get_key(int key_length)
	{
		std::ifstream fin1;
		fin1.open("/home/jackie/spdz/pre_data/predata_mspdz/linear_regression/triple_predata/1/key.txt", std::ios::in);
		if (!(fin1.is_open()))
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
				for (int j = 0; j < key_length; j++)
				{
					mpz_set_ui(key[j], Triple[i].value[j]);
				}
			}
		}

		Triple.clear();
		delete[] line1;
		fin1.close();
	}

	void Online_nonsqure_mul(mpz_t *x, mpz_t *y, mpz_t *mac_x, mpz_t *mac_y, mpz_t *output, mpz_t *output_mac, int l1, int l2, int l3, int key_length)
	{
		get_triple(l1, l2, l3, key_length);

		mpz_t ppp;
		mpz_init(ppp);
		mpz_set_ui(ppp, 2305843009213693951UL);
		mpz_t *d = new mpz_t[l1 * l2];
		mpz_t *e = new mpz_t[l2 * l3];
		mpz_t *f = new mpz_t[l3 * l1];
		uint64_t *d_int = new uint64_t[l1 * l2], *e_int = new uint64_t[l2 * l3], *f_int = new uint64_t[l3 * l1];

		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < l1 * l2; i++)
		{
			mpz_init(d[i]);
			mpz_sub(d[i], x[i], a[i]);
			mpz_mod(d[i], d[i], ppp);
			d_int[i] = mpz_get_ui(d[i]);
		}

		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < l2 * l3; i++)
		{
			mpz_init(e[i]);
			mpz_sub(e[i], y[i], b[i]);
			mpz_mod(e[i], e[i], ppp);
			e_int[i] = mpz_get_ui(e[i]);
		}

		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < l3 * l1; i++)
		{
			mpz_init(f[i]);
			mpz_set_ui(f[i], 0);
		}

		if (party != 1)
		{
			io->send_data(1, d_int, l1 * l2 * sizeof(uint64_t));
			io->flush(1);
			io->recv_data(1, d_int, l1 * l2 * sizeof(uint64_t));
			omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
			for (int i = 0; i < l1 * l2; i++)
			{
				mpz_set_ui(d[i], d_int[i]);
			}
		}
		else
		{
			uint64_t *tmp[nP + 1];
			for (int i = 1; i <= nP; ++i)
			{
				tmp[i] = new uint64_t[l1 * l2];
			}
			vector<future<void>> res;
			int party2 = 2;
			res.push_back(pool->enqueue([this, tmp, party2, l1, l2]()
										{ io->recv_data(party2, tmp[party2], l1 * l2 * sizeof(uint64_t)); }));
			joinNclean(res);

			mpz_t *tmp_i = new mpz_t[l1 * l2];
			omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
			for (int i = 0; i < l1 * l2; ++i)
			{
				mpz_init(tmp_i[i]);
				mpz_set_ui(tmp_i[i], tmp[2][i]);
				mpz_add(d[i], d[i], tmp_i[i]);
				mpz_mod(d[i], d[i], ppp);
				d_int[i] = mpz_get_ui(d[i]);
				mpz_clear(tmp_i[i]);
			}
			delete[] tmp_i;

			res.push_back(pool->enqueue([this, d_int, party2, l1, l2]()
										{
					io->send_data(party2, d_int, l1 * l2 * sizeof(uint64_t));
					io->flush(party2); }));
			joinNclean(res);

			for (int i = 1; i <= nP; ++i)
			{
				delete[] tmp[i];
			}
		}

		if (party != 1)
		{
			io->send_data(1, e_int, l2 * l3 * sizeof(uint64_t));
			io->flush(1);
			io->recv_data(1, e_int, l2 * l3 * sizeof(uint64_t));
			omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
			for (int i = 0; i < l2 * l3; i++)
			{
				mpz_set_ui(e[i], e_int[i]);
			}
		}
		else
		{
			uint64_t *tmp[nP + 1];
			for (int i = 1; i <= nP; ++i)
			{
				tmp[i] = new uint64_t[l2 * l3];
			}
			vector<future<void>> res;
			int party2 = 2;
			res.push_back(pool->enqueue([this, tmp, party2, l2, l3]()
										{ io->recv_data(party2, tmp[party2], l2 * l3 * sizeof(uint64_t)); }));
			joinNclean(res);

			mpz_t *tmp_i = new mpz_t[l2 * l3];
			omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
			for (int i = 0; i < l2 * l3; ++i)
			{
				mpz_init(tmp_i[i]);
				mpz_set_ui(tmp_i[i], tmp[2][i]);
				mpz_add(e[i], e[i], tmp_i[i]);
				mpz_mod(e[i], e[i], ppp);
				e_int[i] = mpz_get_ui(e[i]);
				mpz_clear(tmp_i[i]);
			}
			delete[] tmp_i;

			res.push_back(pool->enqueue([this, e_int, party2, l2, l3]()
										{
					io->send_data(party2, e_int, l2 * l3 * sizeof(uint64_t));
					io->flush(party2); }));
			joinNclean(res);
			for (int i = 1; i <= nP; ++i)
			{
				delete[] tmp[i];
			}
		}

		mpz_t *multc1;
		multc1 = new mpz_t[l3 * l1];

		// compute F = E^T * A^T - R^T
		for (int i = 0; i < l3; i++)
		{
			for (int j = 0; j < l1; j++)
			{
				int index1 = i * l1 + j;
				mpz_init(multc1[index1]);
				for (int k = 0; k < l2; k++)
				{
					mpz_mul(multc1[index1], e[k * l3 + i], a_t[k * l1 + j]);
					mpz_add(f[index1], f[index1], multc1[index1]);
				}
				mpz_sub(f[index1], f[index1], r_t[index1]);
				mpz_mod(f[index1], f[index1], ppp);
				f_int[index1] = mpz_get_ui(f[index1]);
			}
		}

		if (party != 1)
		{
			io->send_data(1, f_int, l3 * l1 * sizeof(uint64_t));
			io->flush(1);
			io->recv_data(1, f_int, l3 * l1 * sizeof(uint64_t));
			omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
			for (int i = 0; i < l3 * l1; i++)
			{
				mpz_set_ui(f[i], f_int[i]);
			}
		}
		else
		{
			uint64_t *tmp[nP + 1];
			for (int i = 1; i <= nP; ++i)
			{
				tmp[i] = new uint64_t[l3 * l1];
			}
			vector<future<void>> res;
			int party2 = 2;
			res.push_back(pool->enqueue([this, tmp, party2, l1, l3]()
										{ io->recv_data(party2, tmp[party2], l3 * l1 * sizeof(uint64_t)); }));
			joinNclean(res);
			mpz_t *tmp_i = new mpz_t[l3 * l1];
			omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
			for (int i = 0; i < l3 * l1; ++i)
			{
				mpz_init(tmp_i[i]);
				mpz_set_ui(tmp_i[i], tmp[2][i]);
				mpz_add(f[i], f[i], tmp_i[i]);
				mpz_mod(f[i], f[i], ppp);
				f_int[i] = mpz_get_ui(f[i]);
				mpz_clear(tmp_i[i]);
			}
			delete[] tmp_i;

			res.push_back(pool->enqueue([this, f_int, party2, l1, l3]()
										{
					io->send_data(party2, f_int, l3 * l1 * sizeof(uint64_t));
					io->flush(party2); }));
			joinNclean(res);

			for (int i = 1; i <= nP; ++i)
			{
				delete[] tmp[i];
			}
		}

		mpz_t *f_t = new mpz_t[l1 * l3];
		mpz_t *db = new mpz_t[l1 * l3];
		mpz_t *de = new mpz_t[l1 * l3];
		mpz_t *mac_db = new mpz_t[l1];
		mpz_t *mac_de = new mpz_t[l1];
		mpz_t *mac_f_t = new mpz_t[l1];
		mpz_t *multc_1 = new mpz_t[l1];
		mpz_t *multc_2 = new mpz_t[l1 * l3];
		mpz_t *multc_3 = new mpz_t[l1 * l3];

		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < l1; i++)
		{
			mpz_init(multc_1[i]);
		}

		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < l1 * l3; i++)
		{
			mpz_init(multc_2[i]);
			mpz_init(multc_3[i]);
		}

		if (party == 1)
		{
			// compute db[sz],de[sz],mac_db[mat_sz],f_t[sz],mac_de[mat_sz],mac_f_t[mat_sz]
			omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
			for (int i = 0; i < l1; i++)
			{
				mpz_init(mac_db[i]);
				mpz_set_ui(mac_db[i], 0);
				mpz_init(mac_de[i]);
				mpz_set_ui(mac_de[i], 0);
				mpz_init(mac_f_t[i]);
				mpz_set_ui(mac_f_t[i], 0);
				for (int j = 0; j < l3; j++)
				{
					int index2 = i * l3 + j;
					mpz_init(f_t[index2]);
					mpz_set(f_t[index2], f[j * l1 + i]);
					mpz_init(db[index2]);
					mpz_set_ui(db[index2], 0);
					mpz_init(de[index2]);
					mpz_set_ui(de[index2], 0);
					for (int k = 0; k < l2; k++)
					{
						mpz_mul(multc_2[index2], d[i * l2 + k], b[k * l3 + j]);
						mpz_add(db[index2], db[index2], multc_2[index2]);
						mpz_mul(multc_3[index2], d[i * l2 + k], e[k * l3 + j]);
						mpz_add(de[index2], de[index2], multc_3[index2]);
					}
					mpz_mod(db[index2], db[index2], ppp);
					mpz_mod(de[index2], de[index2], ppp);
					mpz_mul(multc_1[i], de[i * l3 + j], key[j]);
					mpz_add(mac_de[i], mac_de[i], multc_1[i]);
					mpz_mul(multc_1[i], f_t[i * l3 + j], key[j]);
					mpz_add(mac_f_t[i], mac_f_t[i], multc_1[i]);
				}
				for (int j = 0; j < l2; j++)
				{
					mpz_mul(multc_1[i], d[i * l2 + j], mac_b[j]);
					mpz_add(mac_db[i], mac_db[i], multc_1[i]);
				}
				mpz_mod(mac_db[i], mac_db[i], ppp);
				mpz_mod(mac_de[i], mac_de[i], ppp);
				mpz_mod(mac_f_t[i], mac_f_t[i], ppp);
			}
		}
		else
		{
			// compute db[sz],de[sz],mac_db[mat_sz],f_t[sz],mac_de[mat_sz],mac_f_t[mat_sz]
			omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
			for (int i = 0; i < l1; i++)
			{
				mpz_init(mac_db[i]);
				mpz_set_ui(mac_db[i], 0);
				mpz_init(mac_de[i]);
				mpz_set_ui(mac_de[i], 0);
				mpz_init(mac_f_t[i]);
				mpz_set_ui(mac_f_t[i], 0);
				for (int j = 0; j < l3; j++)
				{
					int index2 = i * l3 + j;
					mpz_init(f_t[index2]);
					mpz_set(f_t[index2], f[j * l1 + i]);
					mpz_init(db[index2]);
					mpz_set_ui(db[index2], 0);
					mpz_init(de[index2]);
					mpz_set_ui(de[index2], 0);
					for (int k = 0; k < l2; k++)
					{
						mpz_mul(multc_2[index2], d[i * l2 + k], b[k * l3 + j]);
						mpz_add(db[index2], db[index2], multc_2[index2]);
						mpz_mul(multc_3[index2], d[i * l2 + k], e[k * l3 + j]);
						mpz_add(de[index2], de[index2], multc_3[index2]);
					}
					mpz_mod(db[index2], db[index2], ppp);
					mpz_mod(de[index2], de[index2], ppp);
					mpz_mul(multc_1[i], de[i * l3 + j], key[j]);
					mpz_add(mac_de[i], mac_de[i], multc_1[i]);
					mpz_mul(multc_1[i], f_t[i * l3 + j], key[j]);
					mpz_add(mac_f_t[i], mac_f_t[i], multc_1[i]);
				}
				for (int j = 0; j < l2; j++)
				{
					mpz_mul(multc_1[i], d[i * l2 + j], mac_b[j]);
					mpz_add(mac_db[i], mac_db[i], multc_1[i]);
				}
				mpz_mod(mac_db[i], mac_db[i], ppp);
				mpz_mod(mac_de[i], mac_de[i], ppp);
				mpz_mod(mac_f_t[i], mac_f_t[i], ppp);
			}
		}

		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < l1; i++)
		{
			mpz_clear(multc_1[i]);
		}

		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < l1 * l3; ++i)
		{
			mpz_clear(multc_2[i]);
			mpz_clear(multc_3[i]);
			if (party == 1)
			{
				mpz_add(output[i], c[i], db[i]);
				mpz_add(output[i], output[i], r[i]);
				mpz_add(output[i], output[i], de[i]);
				mpz_add(output[i], output[i], f_t[i]);
				mpz_mod(output[i], output[i], ppp);
			}
			else
			{
				mpz_add(output[i], c[i], db[i]);
				mpz_add(output[i], output[i], r[i]);
				mpz_mod(output[i], output[i], ppp);
			}
		}

		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < l1; i++)
		{
			mpz_add(output_mac[i], mac_c[i], mac_db[i]);
			mpz_add(output_mac[i], output_mac[i], mac_r[i]);
			mpz_add(output_mac[i], output_mac[i], mac_de[i]);
			mpz_add(output_mac[i], output_mac[i], mac_f_t[i]);
			mpz_mod(output_mac[i], output_mac[i], ppp);
		}

		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < l1 * l2; i++)
		{
			mpz_clear(d[i]);
		}

		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < l2 * l3; i++)
		{
			mpz_clear(e[i]);
		}

		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < l1 * l3; i++)
		{
			mpz_clear(f[i]);
			mpz_clear(f_t[i]);
			mpz_clear(db[i]);
			mpz_clear(de[i]);
		}

		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < l1; i++)
		{
			mpz_clear(mac_db[i]);
			mpz_clear(mac_de[i]);
			mpz_clear(mac_f_t[i]);
		}

		mpz_clear(ppp);

		delete[] multc_1;
		delete[] multc_2;
		delete[] multc_3;
		delete[] d;
		delete[] d_int;
		delete[] e;
		delete[] e_int;
		delete[] f;
		delete[] f_int;
		delete[] f_t;
		delete[] db;
		delete[] de;
		delete[] mac_db;
		delete[] mac_de;
		delete[] mac_f_t;
	}
	
	void mac_check(mpz_t *x, mpz_t *mac_x, int l1, int l3)
	{
		mpz_t ppp;
		mpz_init(ppp);
		mpz_set_ui(ppp, 2305843009213693951UL);
		// gmp_printf("Value of x[0]: %Zd\n", x[0]);

		uint64_t *x_int = new uint64_t[l1 * l3];
		mpz_t *x1 = new mpz_t[l1 * l3];

		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < l1 * l3; i++)
		{
			x_int[i] = mpz_get_ui(x[i]);
			mpz_init(x1[i]);
		}

		if (party != 1)
		{
			io->send_data(1, x_int, l1 * l3 * sizeof(uint64_t));
			io->flush(1);
			io->recv_data(1, x_int, l1 * l3 * sizeof(uint64_t));
			omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
			for (int i = 0; i < l1 * l3; i++)
			{
				mpz_set_ui(x1[i], x_int[i]);
			}
		}
		else
		{
			uint64_t *tmp[nP + 1];
			for (int i = 1; i <= nP; ++i)
			{
				tmp[i] = new uint64_t[l1 * l3];
			}
			vector<future<void>> res;
			int party2 = 2;
			res.push_back(pool->enqueue([this, tmp, party2, l1, l3]()
										{ io->recv_data(party2, tmp[party2], l1 * l3 * sizeof(uint64_t)); }));
			joinNclean(res);

			mpz_t *tmp_i = new mpz_t[l1 * l3];
			omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
			for (int i = 0; i < l1 * l3; ++i)
			{
				mpz_init(tmp_i[i]);
				mpz_set_ui(tmp_i[i], tmp[2][i]);
				mpz_add(x1[i], x[i], tmp_i[i]);
				mpz_mod(x1[i], x1[i], ppp);
				x_int[i] = mpz_get_ui(x1[i]);
				mpz_clear(tmp_i[i]);
			}
			delete[] tmp_i;

			res.push_back(pool->enqueue([this, x_int, party2, l1, l3]()
										{
					io->send_data(party2, x_int, l1 * l3 * sizeof(uint64_t));
					io->flush(party2); }));
			joinNclean(res);

			for (int i = 1; i <= nP; ++i)
			{
				delete[] tmp[i];
			}
		}
		// gmp_printf("Value of x1[0]: %Zd\n", x1[0]);

		mpz_t *sigma_x = new mpz_t[l1];
		mpz_t *kk = new mpz_t[l1];

		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < l1; i++)
		{
			mpz_init(kk[i]);
		}
		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < l1; i++)
		{
			mpz_init(sigma_x[i]);
			mpz_set_ui(sigma_x[i], 0);
			for (int j = 0; j < l3; j++)
			{
				mpz_mul(kk[i], x1[i * l3 + j], key[j]);
				mpz_add(sigma_x[i], sigma_x[i], kk[i]);
			}
			mpz_sub(sigma_x[i], mac_x[i], sigma_x[i]);
			mpz_mod(sigma_x[i], sigma_x[i], ppp);
		}

		uint64_t *sigma_x_int = new uint64_t[m];
		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < l1; i++)
		{
			sigma_x_int[i] = mpz_get_ui(sigma_x[i]);
		}
		if (party != 1)
		{
			io->send_data(1, sigma_x_int, l1 * sizeof(uint64_t));
			io->flush(1);
		}
		else
		{
			uint64_t *tmp[nP + 1];
			for (int i = 1; i <= nP; ++i)
			{
				tmp[i] = new uint64_t[m];
			}
			vector<future<void>> res;
			int party2 = 2;
			res.push_back(pool->enqueue([this, tmp, party2, l1]()
										{ io->recv_data(party2, tmp[party2], l1 * sizeof(uint64_t)); }));
			joinNclean(res);

			uint64_t sum_check = 0;
			mpz_t *tmp_i = new mpz_t[l1];
			omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
			for (int i = 0; i < l1; ++i)
			{
				for (int j = 2; j <= nP; ++j)
				{
					mpz_init(tmp_i[i]);
					mpz_set_ui(tmp_i[i], tmp[j][i]);
					mpz_add(sigma_x[i], sigma_x[i], tmp_i[i]);
				}
				mpz_mod(sigma_x[i], sigma_x[i], ppp);
				if (mpz_sgn(sigma_x[i]) != 0)
				{
					std::cout << "check error" << endl;
					std::cout << i << endl;
					sum_check = sum_check + 1;
				}
			}
			if (sum_check == 0)
				std::cout << "Correct!" << endl;
			mpz_clear(ppp);
			omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
			for (int i = 0; i < l1; i++)
			{
				mpz_clear(tmp_i[i]);
				mpz_clear(kk[i]);
				mpz_clear(sigma_x[i]);
			}
			delete[] tmp_i;
			for (int i = 1; i <= nP; ++i)
			{
				delete[] tmp[i];
			}
		}
		delete[] x_int;
		delete[] kk;
		delete[] sigma_x;
	}

	void Gradient_Descent_linear_regression(mpz_t *x, mpz_t *w, mpz_t *y, mpz_t *mac_x, mpz_t *mac_w, mpz_t *mac_y, int m, int n, int l, int key_length)
	{
		mpz_t ppp;
		mpz_init(ppp);
		mpz_set_ui(ppp, 2305843009213693951UL);
		int accuracy = 15;

		get_key(key_length);

		mpz_t *middle = new mpz_t[n * l];
		mpz_t *x_t = new mpz_t[n * m];
		mpz_t *x_tx = new mpz_t[n * n];
		mpz_t *x_ty = new mpz_t[n * l];
		mpz_t *middle_mac = new mpz_t[n];
		mpz_t *mac_x_t = new mpz_t[n];
		mpz_t *mac_x_tx = new mpz_t[n];
		mpz_t *mac_x_ty = new mpz_t[n];
		mpz_t *derivative = new mpz_t[n * l];
		mpz_t *derivative_mac = new mpz_t[n];

		// init
		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < n * l; i++)
		{
			mpz_init(middle[i]);
		}

		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < n * l; i++)
		{
			mpz_init(derivative[i]);
		}

		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < n; i++)
		{
			mpz_init(middle_mac[i]);
		}

		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < n; i++)
		{
			mpz_init(derivative_mac[i]);
		}

		mpz_t *d = new mpz_t[m * n];
		uint64_t *d_int = new uint64_t[m * n];
		get_triple(m, n, l, key_length);
		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < m * n; i++)
		{
			mpz_init(d[i]);
			mpz_sub(d[i], x[i], a[i]);
			mpz_mod(d[i], d[i], ppp);
			d_int[i] = mpz_get_ui(d[i]);
		}

		if (party != 1)
		{
			io->send_data(1, d_int, m * n * sizeof(uint64_t));
			io->flush(1);
			io->recv_data(1, d_int, m * n * sizeof(uint64_t));
			omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
			for (int i = 0; i < m * n; i++)
			{
				mpz_set_ui(d[i], d_int[i]);
			}
		}
		else
		{
			uint64_t *tmp[nP + 1];
			for (int i = 1; i <= nP; ++i)
			{
				tmp[i] = new uint64_t[m * n];
			}
			vector<future<void>> res;
			int party2 = 2;
			res.push_back(pool->enqueue([this, tmp, party2, m, n]()
										{ io->recv_data(party2, tmp[party2], m * n * sizeof(uint64_t)); }));
			joinNclean(res);
			mpz_t *tmp_i = new mpz_t[m * n];
			omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
			for (int i = 0; i < m * n; ++i)
			{
				mpz_init(tmp_i[i]);
				mpz_set_ui(tmp_i[i], tmp[2][i]);
				mpz_add(d[i], d[i], tmp_i[i]);
				mpz_mod(d[i], d[i], ppp);
				d_int[i] = mpz_get_ui(d[i]);
			}
			res.push_back(pool->enqueue([this, d_int, party2, m, n]()
										{
					io->send_data(party2, d_int, m * n * sizeof(uint64_t));
					io->flush(party2); }));
			joinNclean(res);
			for (int i = 1; i <= nP; ++i)
			{
				delete[] tmp[i];
			}
		}

		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < n * m; i++)
		{
			mpz_init(x_t[i]);
		}

		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < n * n; i++)
		{
			mpz_init(x_tx[i]);
			mpz_set_ui(x_tx[i], 0);
		}

		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < n * l; i++)
		{
			mpz_init(x_ty[i]);
			mpz_set_ui(x_ty[i], 0);
		}

		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < n; i++)
		{
			mpz_init(mac_x_t[i]);
			mpz_init(mac_x_tx[i]);
			mpz_init(mac_x_ty[i]);
			mpz_set_ui(mac_x_tx[i], 0);
			mpz_set_ui(mac_x_ty[i], 0);
		}

		// compute X^T and mac_x_t
		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < n; i++)
		{
			mpz_set(mac_x_t[i], mac_a_t[i]);
			for (int j = 0; j < m; j++)
			{
				mpz_set(x_t[i * m + j], x[j * n + i]);
				mpz_mul(derivative_mac[i], d[j * n + i], key[j]);
				mpz_add(mac_x_t[i], mac_x_t[i], derivative_mac[i]);
			}
			mpz_mod(mac_x_t[i], mac_x_t[i], ppp);
		}
		// compute X^T*X and mac_x_tx
		Online_nonsqure_mul(x_t, x, mac_x_t, mac_x, x_tx, mac_x_tx, n, m, n, key_length);
		mpz_t *x_tx_mac_spdz = new mpz_t[n * n];
		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for(int i = 0;i<n*n;i++)
		{
			mpz_init(x_tx_mac_spdz[i]);
		}
		transform_vectormac_to_spdzmac(x_tx, mac_x_tx, x_tx_mac_spdz, n, m, n, key_length);
		truncate(x_tx, mac_x_tx, x_tx_mac_spdz, accuracy, n, m, n, key_length);
		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for(int i = 0;i<n*n;i++)
		{
			mpz_clear(x_tx_mac_spdz[i]);
		}
		delete []x_tx_mac_spdz;
		// gmp_printf("Value of x_tx: %Zd\n", x_tx[0]);

		// compute X^T*y and mac_x_ty
		Online_nonsqure_mul(x_t, y, mac_x_t, mac_y, x_ty, mac_x_ty, n, m, l, key_length);
		mpz_t *x_ty_mac_spdz = new mpz_t[n * l];
		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for(int i = 0;i<n*l;i++)
		{
			mpz_init(x_ty_mac_spdz[i]);
		}
		transform_vectormac_to_spdzmac(x_ty, mac_x_ty, x_ty_mac_spdz, n, m, l, key_length);
		truncate(x_ty, mac_x_ty, x_ty_mac_spdz, accuracy, n, m, l, key_length);
		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for(int i = 0;i<n*l;i++)
		{
			mpz_clear(x_ty_mac_spdz[i]);
		}
		delete []x_ty_mac_spdz;
		// gmp_printf("Value of x_ty: %Zd\n", x_ty[0]);		

		// train 50 times
		for (int i = 0; i < 50; i++)
		{
			// compute (X^T*X)w
			//std::cout << "To party : sent bytes = "
              //<< io->count() <<"B"<< std::endl;
			Online_nonsqure_mul(x_tx, w, mac_x_tx, mac_w, middle, middle_mac, n, n, l, key_length);
			mpz_t *middle_mac_spdz = new mpz_t[n * l];
			omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
			for(int i = 0;i<n*l;i++)
			{
				mpz_init(middle_mac_spdz[i]);
			}
			transform_vectormac_to_spdzmac(middle, middle_mac, middle_mac_spdz, n, n, l, key_length);
			truncate(middle, middle_mac, middle_mac_spdz, accuracy, n, n, l, key_length);

			// compute (X^T*X)w-X^T*y
			omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
			for (int j = 0; j < n * l; j++)
			{
				mpz_sub(middle[j], middle[j], x_ty[j]);
				mpz_mod(middle[j], middle[j], ppp);
			}
			omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
			for (int j = 0; j < n; j++)
			{
				mpz_sub(middle_mac[j], middle_mac[j], mac_x_ty[j]);
				mpz_mod(middle_mac[j], middle_mac[j], ppp);
			}
			// mac_check(middle, middle_mac, m, l);
			// gmp_printf("Value of middle[0]: %Zd\n", middle[0]);

			// learning rate = 2392*2^{-13} = 0.292
			transform_vectormac_to_spdzmac(middle, middle_mac, middle_mac_spdz, n, n, l, key_length);
			truncate(middle, middle_mac, middle_mac_spdz, 13, n, n, l, key_length);
			omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
			for(int i = 0;i<n*l;i++)
			{
				mpz_clear(middle_mac_spdz[i]);
			}
			delete []middle_mac_spdz;

			for (int j = 0; j < n * l; j++)
			{
				mpz_mod(middle[j], middle[j], ppp);
				mpz_sub(w[j], w[j], middle[j]);
				mpz_mod(w[j], w[j], ppp);
			}
			// gmp_printf("Value of derivative[0]: %Zd\n", derivative[0]);
			for (int j = 0; j < n; j++)
			{
				// mpz_mul(derivative_mac[j],derivative_mac[j],inverse_accuracy);
				mpz_mod(middle_mac[j], middle_mac[j], ppp);
				mpz_sub(mac_w[j], mac_w[j], middle_mac[j]);
				mpz_mod(mac_w[j], mac_w[j], ppp);
			}
		}
		mac_check(w, mac_w, n, l);

		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < l * n; i++)
		{
			d_int[i] = mpz_get_ui(w[i]);
		}

		if (party != 1)
		{
			io->send_data(1, d_int, l * n * sizeof(uint64_t));
			io->flush(1);
			io->recv_data(1, d_int, l * n * sizeof(uint64_t));
			omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
			for (int i = 0; i < l * n; i++)
			{
				mpz_set_ui(w[i], d_int[i]);
			}
		}
		else
		{
			uint64_t *tmp[nP + 1];
			for (int i = 1; i <= nP; ++i)
			{
				tmp[i] = new uint64_t[l * n];
			}
			vector<future<void>> res;
			int party2 = 2;
			res.push_back(pool->enqueue([this, tmp, party2, l, n]()
										{ io->recv_data(party2, tmp[party2], l * n * sizeof(uint64_t)); }));
			joinNclean(res);
			mpz_t *tmp_i = new mpz_t[l * n];
			omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
			for (int i = 0; i < l * n; ++i)
			{
				mpz_init(tmp_i[i]);
				mpz_set_ui(tmp_i[i], tmp[2][i]);
				mpz_add(w[i], w[i], tmp_i[i]);
				mpz_mod(w[i], w[i], ppp);
				d_int[i] = mpz_get_ui(w[i]);
			}
			res.push_back(pool->enqueue([this, d_int, party2, l, n]()
										{
					io->send_data(party2, d_int, l * n * sizeof(uint64_t));
					io->flush(party2); }));
			joinNclean(res);
			for (int i = 1; i <= nP; ++i)
			{
				delete[] tmp[i];
			}
		}
		for (int i = 0; i < n; i++)
		{
			cout << "Value of w[";
			cout << i;
			gmp_printf("]: %Zd\n", w[i]);
		}
	}

	void transform_vectormac_to_spdzmac(mpz_t *x, mpz_t *x_mac, mpz_t *x_mac_spdz, int l1, int l2, int l3, int key_length)
	{
		mpz_t ppp;
		mpz_init(ppp);
		mpz_set_ui(ppp, 2305843009213693951UL);
		get_transform_tuple(l1, l3, key_length);

		mpz_t *transform = new mpz_t[l1 * l3];
		uint64_t *transform_int = new uint64_t[l1 * l3];
		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < l1 * l3; i++)
		{
			mpz_init(transform[i]);
			mpz_sub(transform[i], x[i], t[i]);
			transform_int[i] = mpz_get_ui(transform[i]);
		}

		// open transform
		if (party != 1)
		{
			io->send_data(1, transform_int, l1 * l3 * sizeof(uint64_t));
			io->flush(1);
			io->recv_data(1, transform_int, l1 * l3 * sizeof(uint64_t));
			omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
			for (int i = 0; i < l1 * l3; i++)
			{
				mpz_set_ui(transform[i], transform_int[i]);
			}
		}
		else
		{
			uint64_t *tmp[nP + 1];
			for (int i = 1; i <= nP; ++i)
			{
				tmp[i] = new uint64_t[l1 * l3];
			}
			vector<future<void>> res;
			int party2 = 2;
			res.push_back(pool->enqueue([this, tmp, party2, l1, l3]()
										{ io->recv_data(party2, tmp[party2], l1 * l3 * sizeof(uint64_t)); }));
			joinNclean(res);

			mpz_t *tmp_i = new mpz_t[l1 * l3];
			omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
			for (int i = 0; i < l1 * l3; ++i)
			{
				mpz_init(tmp_i[i]);
				mpz_set_ui(tmp_i[i], tmp[2][i]);
				mpz_add(transform[i], transform[i], tmp_i[i]);
				mpz_mod(transform[i], transform[i], ppp);
				transform_int[i] = mpz_get_ui(transform[i]);
				mpz_clear(tmp_i[i]);
			}
			delete[] tmp_i;

			res.push_back(pool->enqueue([this, transform_int, party2, l1, l3]()
										{
					io->send_data(party2, transform_int, l1 * l3 * sizeof(uint64_t));
					io->flush(party2); }));
			joinNclean(res);
			for (int i = 1; i <= nP; ++i)
			{
				delete[] tmp[i];
			}
		}

		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < l1; i++)
		{
			for (int j = 0; j < l3; j++)
			{
				mpz_mul(x_mac_spdz[i * l3 + j], transform[i * l3 + j], key[j]);
				mpz_add(x_mac_spdz[i * l3 + j], x_mac_spdz[i * l3 + j], mac_t_spdz[i * l3 + j]);
				mpz_mod(x_mac_spdz[i * l3 + j], x_mac_spdz[i * l3 + j], ppp);
			}
		}
		// gmp_printf("Value of x_mac_spdz[0]: %Zd\n", x_mac_spdz[0]);
	}

	void truncate(mpz_t *x, mpz_t *x_mac, mpz_t *x_mac_spdz, int accuracy, int l1, int l2, int l3, int key_length)
	{
		mpz_t ppp;
		mpz_t *k_const = new mpz_t[l1 * l3];
		mpz_t *s_msb = new mpz_t[l1 * l3];
		mpz_t *s1 = new mpz_t[l1 * l3];
		mpz_t *s_msb_mac_spdz = new mpz_t[l1 * l3];
		mpz_t *s1_mac_spdz = new mpz_t[l1 * l3];
		mpz_init(ppp);
		mpz_set_ui(ppp, 2305843009213693951UL);
		// mac_check(x, x_mac, l1, l3);

		get_truncate_tuple(l1, l3, key_length);

		// compute s_msb and s1
		// gmp_printf("Value of s[0]: %Zd\n", s[0]);
		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < l1 * l3; i++)
		{
			mpz_init(k_const[i]);
			mpz_init(s_msb[i]);
			mpz_init(s1[i]);
			mpz_set(s_msb[i], s_bit[l1 * l3 * (l_bit - 1) + i]);
			mpz_set_ui(s1[i], 0);
			for (int j = accuracy; j < l_bit; j++)
			{
				mpz_ui_pow_ui(k_const[i], 2, j - accuracy);
				mpz_mul(k_const[i], k_const[i], s_bit[j * l1 * l3 + i]);
				mpz_add(s1[i], s1[i], k_const[i]);
				mpz_mod(s1[i], s1[i], ppp);
			}
			for (int j = l_bit - accuracy; j < l_bit; j++)
			{
				mpz_ui_pow_ui(k_const[i], 2, j);
				mpz_mul(k_const[i], k_const[i], s_msb[i]);
				mpz_add(s1[i], s1[i], k_const[i]);
				mpz_mod(s1[i], s1[i], ppp);
			}
		}
		// gmp_printf("Value of s_msb[0]: %Zd\n", s_msb[0]);
		// gmp_printf("Value of s1[0]: %Zd\n", s1[0]);

		// compute s_msb_mac_spdz and s1_mac_spdz
		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < l1 * l3; i++)
		{
			mpz_init(s_msb_mac_spdz[i]);
			mpz_init(s1_mac_spdz[i]);
			mpz_set_ui(s1_mac_spdz[i], 0);
			mpz_set(s_msb_mac_spdz[i], mac_s_bit_spdz[l1 * l3 * (l_bit - 1) + i]);
			for (int j = accuracy; j < l_bit; j++)
			{
				mpz_ui_pow_ui(k_const[i], 2, j - accuracy);
				mpz_mul(k_const[i], k_const[i], mac_s_bit_spdz[j * l1 * l3 + i]);
				mpz_add(s1_mac_spdz[i], s1_mac_spdz[i], k_const[i]);
				mpz_mod(s1_mac_spdz[i], s1_mac_spdz[i], ppp);
			}
			for (int j = l_bit - accuracy; j < l_bit; j++)
			{
				mpz_ui_pow_ui(k_const[i], 2, j);
				mpz_mul(k_const[i], k_const[i], s_msb_mac_spdz[i]);
				mpz_add(s1_mac_spdz[i], s1_mac_spdz[i], k_const[i]);
				mpz_mod(s1_mac_spdz[i], s1_mac_spdz[i], ppp);
			}
		}
		// mac_check(s_msb, s_msb_mac, l1, l3);
		// mac_check(s1, s1_mac, l1, l3);

		// compute d=x+2^{l-2}
		mpz_t *d = new mpz_t[l1 * l3];
		if (party == 1)
		{
			omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
			for (int i = 0; i < l1 * l3; i++)
			{
				mpz_init(d[i]);
				mpz_ui_pow_ui(k_const[i], 2, l_bit - 2);
				mpz_add(d[i], x[i], k_const[i]);
				mpz_mod(d[i], d[i], ppp);
			}
		}
		else
		{
			omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
			for (int i = 0; i < l1 * l3; i++)
			{
				mpz_init(d[i]);
				mpz_set(d[i], x[i]);
			}
		}

		// compute d_mac
		mpz_t *d_mac_spdz = new mpz_t[l1 * l3];
		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < l1 * l3; i++)
		{
			mpz_init(d_mac_spdz[i]);
			mpz_ui_pow_ui(k_const[i], 2, l_bit - 2);
			mpz_mul(k_const[i], k_const[i], key[i]);
			mpz_add(d_mac_spdz[i], x_mac_spdz[i], k_const[i]);
			mpz_mod(d_mac_spdz[i], d_mac_spdz[i], ppp);
		}
		// gmp_printf("Value of x[0]: %Zd\n", x[0]);
		// gmp_printf("Value of d[0]: %Zd\n", d[0]);
		// mac_check(d, d_mac, l1, l3);

		// compute e=d+s
		mpz_t *e = new mpz_t[l1 * l3];
		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < l1 * l3; i++)
		{
			mpz_init(e[i]);
			mpz_add(e[i], d[i], s[i]);
			mpz_mod(e[i], e[i], ppp);
		}
		// gmp_printf("Value of e[0]: %Zd\n", e[0]);

		// compute e_mac_spdz
		mpz_t *e_mac_spdz = new mpz_t[l1 * l3];
		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < l1 * l3; i++)
		{
			mpz_init(e_mac_spdz[i]);
			mpz_add(e_mac_spdz[i], d_mac_spdz[i], mac_s_spdz[i]);
			mpz_mod(e_mac_spdz[i], e_mac_spdz[i], ppp);
		}
		// mac_check(e, e_mac, l1, l3);

		uint64_t *e_int = new uint64_t[l1 * l3];

		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < l1 * l3; i++)
		{
			e_int[i] = mpz_get_ui(e[i]);
		}
		// cout << "Value of e[0]:" << e_int[0] << endl;

		// open e
		if (party != 1)
		{
			io->send_data(1, e_int, l1 * l3 * sizeof(uint64_t));
			io->flush(1);
			io->recv_data(1, e_int, l1 * l3 * sizeof(uint64_t));
			omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
			for (int i = 0; i < l1 * l3; i++)
			{
				mpz_set_ui(e[i], e_int[i]);
			}
		}
		else
		{
			uint64_t *tmp[nP + 1];
			for (int i = 1; i <= nP; ++i)
			{
				tmp[i] = new uint64_t[l1 * l3];
			}
			vector<future<void>> res;
			int party2 = 2;
			res.push_back(pool->enqueue([this, tmp, party2, l1, l3]()
										{ io->recv_data(party2, tmp[party2], l1 * l3 * sizeof(uint64_t)); }));
			joinNclean(res);

			mpz_t *tmp_i = new mpz_t[l1 * l3];
			omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
			for (int i = 0; i < l1 * l3; ++i)
			{
				mpz_init(tmp_i[i]);
				mpz_set_ui(tmp_i[i], tmp[2][i]);
				mpz_add(e[i], e[i], tmp_i[i]);
				mpz_mod(e[i], e[i], ppp);
				e_int[i] = mpz_get_ui(e[i]);
				mpz_clear(tmp_i[i]);
			}
			delete[] tmp_i;

			res.push_back(pool->enqueue([this, e_int, party2, l1, l3]()
										{
					io->send_data(party2, e_int, l1 * l3 * sizeof(uint64_t));
					io->flush(party2); }));
			joinNclean(res);
			for (int i = 1; i <= nP; ++i)
			{
				delete[] tmp[i];
			}
		}
		// gmp_printf("Value of e[0]: %Zd\n", e[0]);

		// compute trunc_e and e_msb
		mpz_t *trunc_e = new mpz_t[l1 * l3];
		mpz_t *e_msb = new mpz_t[l1 * l3];
		mpz_t compare_const;
		mpz_init(compare_const);
		mpz_sub_ui(compare_const, ppp, 1);
		mpz_fdiv_q_ui(compare_const, compare_const, 2);

		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < l1 * l3; i++)
		{
			mpz_init(trunc_e[i]);
			mpz_init(e_msb[i]);
			if (mpz_cmp_si(e[i], 0) >= 0 && mpz_cmp(e[i], compare_const) <= 0)
			{
				mpz_fdiv_q_2exp(trunc_e[i], e[i], accuracy);
			}
			else
			{
				mpz_sub(trunc_e[i], ppp, e[i]);
				mpz_fdiv_q_2exp(trunc_e[i], trunc_e[i], accuracy);
				mpz_sub(trunc_e[i], ppp, trunc_e[i]);
			}
			mpz_fdiv_q_2exp(e_msb[i], e[i], l_bit - 1);
		}
		// gmp_printf("Value of trunc_e[0]: %Zd\n", trunc_e[0]);
		// gmp_printf("Value of e_msb[0]: %Zd\n", e_msb[0]);

		// compute f=(1-s_msb)e_msb
		mpz_t *f = new mpz_t[l1 * l3];
		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < l1 * l3; i++)
		{
			mpz_init(f[i]);
			mpz_mul(f[i], e_msb[i], s_msb[i]);
			mpz_sub(f[i], ppp, f[i]);
			if (party == 1)
			{
				mpz_add(f[i], e_msb[i], f[i]);
			}
			mpz_mod(f[i], f[i], ppp);
		}

		// compute f_mac
		mpz_t *f_mac_spdz = new mpz_t[l1 * l3];
		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < l1; i++)
		{
			for (int j = 0; j < l3; j++)
			{
				mpz_init(f_mac_spdz[i * l3 + j]);
				mpz_sub(s_msb_mac_spdz[i * l3 + j], key[j], s_msb_mac_spdz[i * l3 + j]);
				mpz_mod(s_msb_mac_spdz[i * l3 + j], s_msb_mac_spdz[i * l3 + j], ppp);
				mpz_mul(f_mac_spdz[i * l3 + j], s_msb_mac_spdz[i * l3 + j], e_msb[i * l3 + j]);
			}
		}

		// compute g=trunc_e-s1+f*(2^{l-d}-1) and g_mac
		mpz_t *g = new mpz_t[l1 * l3];
		mpz_t *g_mac_spdz = new mpz_t[l1 * l3];
		mpz_ui_pow_ui(k_const[0], 2, l_bit - accuracy);
		mpz_t x_const;
		mpz_t m_const;
		mpz_init(x_const);
		mpz_init(m_const);
		mpz_set_si(x_const, -1);
		mpz_add(k_const[0], k_const[0], x_const);		 // k_const[0]=2^{l_bit-accuracy}-1
		mpz_ui_pow_ui(m_const, 2, l_bit - 2 - accuracy); // m_const=2^{l_bit-accuracy-2}

		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < l1 * l3; i++)
		{
			mpz_init(g[i]);
			mpz_mul(g[i], f[i], k_const[0]);
			mpz_sub(x[i], g[i], s1[i]);
			if (party == 1)
			{
				mpz_add(x[i], x[i], trunc_e[i]);
				mpz_sub(x[i], x[i], m_const);
			}
			mpz_mod(x[i], x[i], ppp);
		}

		mpz_t *count1 = new mpz_t[l1 * l3];
		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < l1; i++)
		{
			for (int j = 0; j < l3; j++)
			{
				mpz_init(g_mac_spdz[i * l3 + j]);
				mpz_mul(g_mac_spdz[i * l3 + j], f_mac_spdz[i * l3 + j], k_const[0]);
				mpz_sub(g_mac_spdz[i * l3 + j], g_mac_spdz[i * l3 + j], s1_mac_spdz[i * l3 + j]);
				mpz_mod(g_mac_spdz[i * l3 + j], g_mac_spdz[i * l3 + j], ppp);
				mpz_init(count1[i * l3 + j]);
				mpz_sub(count1[i * l3 + j], trunc_e[i * l3 + j], m_const);
				mpz_mul(count1[i * l3 + j], count1[i * l3 + j], key[j]);
				mpz_add(x_mac_spdz[i * l3 + j], g_mac_spdz[i * l3 + j], count1[i * l3 + j]);
				mpz_mod(x_mac_spdz[i * l3 + j], x_mac_spdz[i * l3 + j], ppp);
			}
		}
		// mac_check(x, x_mac, l1, l3);

		// transform spdz mac to vector mac
		for (int i = 0; i < l1; i++)
		{
			mpz_set_ui(x_mac[i], 0);
			for (int j = 0; j < l3; j++)
			{
				mpz_add(x_mac[i], x_mac[i], x_mac_spdz[i * l3 + j]);
			}
			mpz_mod(x_mac[i], x_mac[i], ppp);
		}

		mpz_clear(m_const);

		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < l1 * l3; i++)
		{
			mpz_clear(count1[i]);
			mpz_clear(k_const[i]);
			mpz_clear(s_msb[i]);
			mpz_clear(s1[i]);
			mpz_clear(d[i]);
			mpz_clear(e[i]);
			mpz_clear(trunc_e[i]);
			mpz_clear(e_msb[i]);
			mpz_clear(f[i]);
			mpz_clear(g[i]);
			mpz_clear(f_mac_spdz[i]);
			mpz_clear(d_mac_spdz[i]);
			mpz_clear(e_mac_spdz[i]);
			mpz_clear(g_mac_spdz[i]);
			mpz_clear(s1_mac_spdz[i]);
			mpz_clear(s_msb_mac_spdz[i]);
		}

		delete[] k_const;
		delete[] s_msb;
		delete[] s1;
		delete[] d;
		delete[] e;
		delete[] trunc_e;
		delete[] e_msb;
		delete[] f;
		delete[] g;
		delete[] count1;
		delete[] d_mac_spdz;
		delete[] e_mac_spdz;
		delete[] f_mac_spdz;
		delete[] s1_mac_spdz;
		delete[] s_msb_mac_spdz;
	}
};
#endif // MSPDZ_H__