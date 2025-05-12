#ifndef MSPDZ_H__
#define MSPDZ_H__
#define NUM_THREADS 8

#include "../network/network.h"
#include <emp-tool/emp-tool.h>
#include "utility.h"
#include "../network/helper.h"
#include <omp.h>
#include "../mpir/mpir.h"
using namespace emp;

// Matrix SPDZ
template <int nP>
class MSPDZ
{
public:
	mpz_t *key;

	mpz_t *a, *a_t, *b, *c, *r, *r_t;
	mpz_t *mac_a, *mac_a_t, *mac_b, *mac_c, *mac_r, *mac_r_t;

	NetIOMP<nP> *io;
	int party, total_pre, ssp;

	int mat_sz = 128;
	int sz = mat_sz * mat_sz;
	ThreadPool *pool;

	PRP prp;
	MSPDZ(NetIOMP<nP> *io[2], ThreadPool *pool, int party, bool *_delta = nullptr, int ssp = 40)
	{
		this->party = party;
		this->io = io[0];
		this->ssp = ssp;
		this->pool = pool;

		a = new mpz_t[sz];
		a_t = new mpz_t[sz];
		b = new mpz_t[sz];
		c = new mpz_t[sz];
		r = new mpz_t[sz];
		r_t = new mpz_t[sz];
		key = new mpz_t[mat_sz];

		mac_a = new mpz_t[mat_sz];
		mac_a_t = new mpz_t[mat_sz];
		mac_b = new mpz_t[mat_sz];
		mac_c = new mpz_t[mat_sz];
		mac_r = new mpz_t[mat_sz];
		mac_r_t = new mpz_t[mat_sz];

		for (int i = 0; i < sz; i++)
		{
			mpz_init(a[i]);
			mpz_init(a_t[i]);
			mpz_init(b[i]);
			mpz_init(c[i]);
			mpz_init(r[i]);
			mpz_init(r_t[i]);
		}
		for (int i = 0; i < mat_sz; i++)
		{
			mpz_init(key[i]);
			mpz_init(mac_a[i]);
			mpz_init(mac_a_t[i]);
			mpz_init(mac_b[i]);
			mpz_init(mac_c[i]);
			mpz_init(mac_r[i]);
			mpz_init(mac_r_t[i]);
		}
	}
	~MSPDZ()
	{
		for (int i = 0; i < sz; i++)
		{
			mpz_clear(a[i]);
			mpz_clear(a_t[i]);
			mpz_clear(b[i]);
			mpz_clear(c[i]);
			mpz_clear(r[i]);
			mpz_clear(r_t[i]);
		}
		for (int i = 0; i < mat_sz; i++)
		{
			mpz_clear(key[i]);
			mpz_clear(mac_a[i]);
			mpz_clear(mac_a_t[i]);
			mpz_clear(mac_b[i]);
			mpz_clear(mac_c[i]);
			mpz_clear(mac_r[i]);
			mpz_clear(mac_r_t[i]);
		}
		delete[] a;
		delete[] a_t;
		delete[] b;
		delete[] c;
		delete[] r;
		delete[] r_t;
		delete[] key;
		delete[] mac_a;
		delete[] mac_a_t;
		delete[] mac_b;
		delete[] mac_c;
		delete[] mac_r;
		delete[] mac_r_t;
	}
	PRG prg;

	// it should be implemented by offline.
	void get_triple()
	{
		std::ifstream fin1, fin2, fin3, fin4, fin5, fin6, fin7;
		std::ifstream fin_1, fin_2, fin_3, fin_4, fin_5, fin_6;
		fin1.open("/home/jackie/spdz/pre_data/predata_mspdz/128/a.txt", std::ios::in);
		fin2.open("/home/jackie/spdz/pre_data/predata_mspdz/128/a_t.txt", std::ios::in);
		fin3.open("/home/jackie/spdz/pre_data/predata_mspdz/128/b.txt", std::ios::in);
		fin4.open("/home/jackie/spdz/pre_data/predata_mspdz/128/c.txt", std::ios::in);
		fin5.open("/home/jackie/spdz/pre_data/predata_mspdz/128/r.txt", std::ios::in);
		fin6.open("/home/jackie/spdz/pre_data/predata_mspdz/128/r_t.txt", std::ios::in);
		fin7.open("/home/jackie/spdz/pre_data/predata_mspdz/128/key.txt", std::ios::in);
		fin_1.open("/home/jackie/spdz/pre_data/predata_mspdz/128/mac_a.txt", std::ios::in);
		fin_2.open("/home/jackie/spdz/pre_data/predata_mspdz/128/mac_a_t.txt", std::ios::in);
		fin_3.open("/home/jackie/spdz/pre_data/predata_mspdz/128/mac_b.txt", std::ios::in);
		fin_4.open("/home/jackie/spdz/pre_data/predata_mspdz/128/mac_c.txt", std::ios::in);
		fin_5.open("/home/jackie/spdz/pre_data/predata_mspdz/128/mac_r.txt", std::ios::in);
		fin_6.open("/home/jackie/spdz/pre_data/predata_mspdz/128/mac_r_t.txt", std::ios::in);

		if (!(fin1.is_open() && fin2.is_open() && fin3.is_open() && fin4.is_open() && fin5.is_open() && fin6.is_open() && fin7.is_open() && fin_1.is_open() && fin_2.is_open() && fin_3.is_open() && fin_4.is_open() && fin_5.is_open() && fin_6.is_open()))
		{
			std::cerr << "cannot open the file";
		}

		char line1[512000] = {0};
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

		for (int i = 0; i < nP; i++)
		{
			if (party == stoi(Triple[i].party))
			{
				for (int j = 0; j < sz; j++)
				{
					mpz_set_ui(a[j], Triple[i].value[j]);
				}
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

		for (int i = 0; i < nP; i++)
		{
			if (party == stoi(Triple[i].party))
			{
				for (int j = 0; j < sz; j++)
					mpz_set_ui(a_t[j], Triple[i].value[j]);
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

		for (int i = 0; i < nP; i++)
		{
			if (party == stoi(Triple[i].party))
			{
				for (int j = 0; j < sz; j++)
					mpz_set_ui(b[j], Triple[i].value[j]);
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

		for (int i = 0; i < nP; i++)
		{
			if (party == stoi(Triple[i].party))
			{
				for (int j = 0; j < sz; j++)
					mpz_set_ui(c[j], Triple[i].value[j]);
			}
		}

		Triple.clear();
		memset(line1, 0, sizeof(line1));

		while (fin5.getline(line1, sizeof(line1)))
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
					mpz_set_ui(r[j], Triple[i].value[j]);
			}
		}

		Triple.clear();
		memset(line1, 0, sizeof(line1));

		while (fin6.getline(line1, sizeof(line1)))
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
					mpz_set_ui(r_t[j], Triple[i].value[j]);
			}
		}

		Triple.clear();
		memset(line1, 0, sizeof(line1));

		while (fin7.getline(line1, sizeof(line1)))
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
				for (int j = 0; j < mat_sz; j++)
					mpz_set_ui(key[j], Triple[i].value[j]);
			}
		}

		Triple.clear();
		memset(line1, 0, sizeof(line1));

		while (fin_1.getline(line1, sizeof(line1)))
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
				for (int j = 0; j < mat_sz; j++)
					mpz_set_ui(mac_a[j], Triple[i].value[j]);
			}
		}

		Triple.clear();
		memset(line1, 0, sizeof(line1));

		while (fin_2.getline(line1, sizeof(line1)))
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
				for (int j = 0; j < mat_sz; j++)
					mpz_set_ui(mac_a_t[j], Triple[i].value[j]);
			}
		}

		Triple.clear();
		memset(line1, 0, sizeof(line1));

		while (fin_3.getline(line1, sizeof(line1)))
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
				for (int j = 0; j < mat_sz; j++)
					mpz_set_ui(mac_b[j], Triple[i].value[j]);
			}
		}

		Triple.clear();
		memset(line1, 0, sizeof(line1));

		while (fin_4.getline(line1, sizeof(line1)))
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
				for (int j = 0; j < mat_sz; j++)
					mpz_set_ui(mac_c[j], Triple[i].value[j]);
			}
		}

		Triple.clear();
		memset(line1, 0, sizeof(line1));

		while (fin_5.getline(line1, sizeof(line1)))
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
				for (int j = 0; j < mat_sz; j++)
					mpz_set_ui(mac_r[j], Triple[i].value[j]);
			}
		}

		Triple.clear();
		memset(line1, 0, sizeof(line1));

		while (fin_6.getline(line1, sizeof(line1)))
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
				for (int j = 0; j < mat_sz; j++)
					mpz_set_ui(mac_r_t[j], Triple[i].value[j]);
			}
		}

		Triple.clear();
		memset(line1, 0, sizeof(line1));
	}

	void Online_mul(mpz_t *x, mpz_t *y, mpz_t *mac_x, mpz_t *mac_y, mpz_t *output, mpz_t *output_mac)
	{
		mpz_t ppp;
		mpz_init(ppp);
		mpz_set_ui(ppp, 2305843009213693951UL);

		mpz_t *d = new mpz_t[sz];
		mpz_t *e = new mpz_t[sz];
		mpz_t *f = new mpz_t[sz];

		omp_set_num_threads(NUM_THREADS);
		#pragma omp parallel for
		for (int i = 0; i < sz; i++)
		{
			mpz_init(d[i]);
			mpz_init(e[i]);
			mpz_init(f[i]);
			mpz_set_ui(f[i], 0);
		}

		omp_set_num_threads(NUM_THREADS);
		#pragma omp parallel for
		for (int i = 0; i < sz; ++i)
		{
			mpz_sub(d[i], x[i], a[i]);
			mpz_mod(d[i], d[i], ppp);
			mpz_sub(e[i], y[i], b[i]);
			mpz_mod(e[i], e[i], ppp);
		}

		uint64_t *d_int;
		d_int = new uint64_t[sz];
		omp_set_num_threads(NUM_THREADS);
		#pragma omp parallel for
		for (int i = 0; i < sz; i++)
		{
			d_int[i]=mpz_get_ui(d[i]);
		}
		if (party != 1)
		{
			io->send_data(1, d_int, sz * sizeof(uint64_t));
			io->flush(1);
			io->recv_data(1, d_int, sz * sizeof(uint64_t));
			omp_set_num_threads(NUM_THREADS);
			#pragma omp parallel for
			for(int i=0;i<sz;i++){
				mpz_set_ui(d[i],d_int[i]);
			}
		}
		else
		{
			uint64_t *tmp[nP + 1];
			for (int i = 1; i <= nP; ++i)
			{
				tmp[i] = new uint64_t[sz];
			}
			vector<future<void>> res;
			for (int i = 2; i <= nP; ++i)
			{
				int party2 = i;
				res.push_back(pool->enqueue([this, tmp, party2]()
											{ io->recv_data(party2, tmp[party2], sz * sizeof(uint64_t)); }));
			}
			joinNclean(res);
			mpz_t tmp_ji;
			mpz_init(tmp_ji);
			for (int i = 0; i < sz; ++i)
			{
				for (int j = 2; j <= nP; ++j)
				{
					mpz_set_ui(tmp_ji,tmp[j][i]);
					mpz_add(d[i],d[i],tmp_ji);
				}
				mpz_mod(d[i], d[i], ppp);
			}
			mpz_clear(tmp_ji);
			omp_set_num_threads(NUM_THREADS);
			#pragma omp parallel for
			for (int i = 0; i < sz; i++)
			{
				d_int[i]=mpz_get_ui(d[i]);
			}
			for (int i = 2; i <= nP; ++i)
			{
				int party2 = i;
				res.push_back(pool->enqueue([this, d_int, party2]()
											{
					io->send_data(party2, d_int, sz * sizeof(uint64_t));
					io->flush(party2); }));
			}
			joinNclean(res);
			for (int i = 1; i <= nP; ++i)
			{
				delete[] tmp[i];
			}
		}
		
		uint64_t *e_int;
		e_int = new uint64_t[sz];
		omp_set_num_threads(NUM_THREADS);
		#pragma omp parallel for
		for (int i = 0; i < sz; i++)
		{
			e_int[i]=mpz_get_ui(e[i]);
		}
		if (party != 1)
		{
			io->send_data(1, e_int, sz * sizeof(uint64_t));
			io->flush(1);
			io->recv_data(1, e_int, sz * sizeof(uint64_t));
			omp_set_num_threads(NUM_THREADS);
			#pragma omp parallel for
			for(int i=0;i<sz;i++){
				mpz_set_ui(e[i],e_int[i]);
			}
		}
		else
		{
			uint64_t *tmp[nP + 1];
			for (int i = 1; i <= nP; ++i)
			{
				tmp[i] = new uint64_t[sz];
			}
			vector<future<void>> res1;
			for (int i = 2; i <= nP; ++i)
			{
				int party2 = i;
				res1.push_back(pool->enqueue([this, tmp, party2]()
											 { io->recv_data(party2, tmp[party2], sz * sizeof(uint64_t)); }));
			}
			joinNclean(res1);
			mpz_t *tmp_ji;
			tmp_ji=new mpz_t[sz];
			omp_set_num_threads(NUM_THREADS);
			#pragma omp parallel for
			for(int i=0;i<sz;i++){
				mpz_init(tmp_ji[i]);
			}
			omp_set_num_threads(NUM_THREADS);
			#pragma omp parallel for
			for (int i = 0; i < sz; ++i)
			{
				for (int j = 2; j <= nP; ++j)
				{
					mpz_set_ui(tmp_ji[i],tmp[j][i]);
					mpz_add(e[i],e[i],tmp_ji[i]);
				}
				mpz_mod(e[i], e[i], ppp);
			}
			omp_set_num_threads(NUM_THREADS);
			#pragma omp parallel for
			for(int i=0;i<sz;i++){
				mpz_clear(tmp_ji[i]);
			}
			omp_set_num_threads(NUM_THREADS);
			#pragma omp parallel for
			for (int i = 0; i < sz; i++)
			{
				e_int[i]=mpz_get_ui(e[i]);
			}
			for (int i = 2; i <= nP; ++i)
			{
				int party2 = i;
				res1.push_back(pool->enqueue([this, e_int, party2]()
											 {
					io->send_data(party2, e_int, sz * sizeof(uint64_t));
					io->flush(party2); }));
			}
			joinNclean(res1);
			for (int i = 1; i <= nP; ++i)
			{
				delete[] tmp[i];
			}
		}
		
		mpz_t *multc1;
		multc1=new mpz_t[sz];
		omp_set_num_threads(NUM_THREADS);
		#pragma omp parallel for
		for(int i=0;i<sz;i++){
			mpz_init(multc1[i]);
		}
		omp_set_num_threads(NUM_THREADS);
		#pragma omp parallel for
		for (int i = 0; i < mat_sz; i++)
		{
			for (int j = 0; j < mat_sz; j++)
			{
				int index1 = i * mat_sz + j;
				for (int k = 0; k < mat_sz; k++)
				{
					mpz_mul(multc1[index1], e[k * mat_sz + i], a_t[k * mat_sz + j]);
					mpz_add(f[index1], f[index1], multc1[index1]);
				}
				mpz_sub(f[index1], f[index1], r_t[index1]);
				mpz_mod(f[index1], f[index1], ppp);
			}
		}
						
		uint64_t *f_int;
		f_int = new uint64_t[sz];
		omp_set_num_threads(NUM_THREADS);
		#pragma omp parallel for
		for (int i = 0; i < sz; i++)
		{
			f_int[i]=mpz_get_ui(f[i]);
		}
		if (party != 1)
		{
			io->send_data(1, f_int, sz * sizeof(uint64_t));
			io->flush(1);
			io->recv_data(1, f_int, sz * sizeof(uint64_t));
			omp_set_num_threads(NUM_THREADS);
			#pragma omp parallel for
			for(int i=0;i<sz;i++){
				mpz_set_ui(f[i],f_int[i]);
			}
		}
		else
		{
			uint64_t *tmp[nP + 1];
			for (int i = 1; i <= nP; ++i)
			{
				tmp[i] = new uint64_t[sz];
			}
			vector<future<void>> res;
			for (int i = 2; i <= nP; ++i)
			{
				int party2 = i;
				res.push_back(pool->enqueue([this, tmp, party2]()
											{ io->recv_data(party2, tmp[party2], sz * sizeof(uint64_t)); }));
			}
			joinNclean(res);
			mpz_t *tmp_ji;
			tmp_ji=new mpz_t[sz];
			omp_set_num_threads(NUM_THREADS);
			#pragma omp parallel for
			for(int i=0;i<sz;i++){
				mpz_init(tmp_ji[i]);
			}
			omp_set_num_threads(NUM_THREADS);
			#pragma omp parallel for
			for (int i = 0; i < sz; ++i)
			{
				for (int j = 2; j <= nP; ++j)
				{
					mpz_set_ui(tmp_ji[i],tmp[j][i]);
					mpz_add(f[i],f[i],tmp_ji[i]);
				}
				mpz_mod(f[i], f[i], ppp);
			}
			omp_set_num_threads(NUM_THREADS);
			#pragma omp parallel for
			for (int i = 0; i < sz; i++)
			{
				mpz_clear(tmp_ji[i]);
				f_int[i]=mpz_get_ui(f[i]);
			}
			for (int i = 2; i <= nP; ++i)
			{
				int party2 = i;
				res.push_back(pool->enqueue([this, f_int, party2]()
											{
					io->send_data(party2, f_int, sz * sizeof(uint64_t));
					io->flush(party2); }));
			}
			joinNclean(res);
			for (int i = 1; i <= nP; ++i)
			{
				delete[] tmp[i];
			}
		}

		mpz_t *f_t = new mpz_t[sz];
		mpz_t *db = new mpz_t[sz];
		mpz_t *de = new mpz_t[sz];
		mpz_t *mac_db = new mpz_t[mat_sz];
		mpz_t *mac_de = new mpz_t[mat_sz];
		mpz_t *mac_f_t = new mpz_t[mat_sz];
		mpz_t *multc_1;
		multc_1 = new mpz_t[sz];
		mpz_t *multc_2;
		multc_2 = new mpz_t[sz];
		mpz_t *multc_3;
		multc_3 = new mpz_t[sz];
		omp_set_num_threads(NUM_THREADS);
		#pragma omp parallel for
		for(int i=0;i<mat_sz;i++){
			mpz_init(multc_1[i]);
		}
		omp_set_num_threads(NUM_THREADS);
		#pragma omp parallel for
		for(int i=0;i<sz;i++){
			mpz_init(multc_1[i]);
			mpz_init(multc_2[i]);
			mpz_init(multc_3[i]);
		}

		if (party == 1)
		{
			omp_set_num_threads(NUM_THREADS);
			#pragma omp parallel for
			for (int i = 0; i < mat_sz; i++)
			{
				mpz_init(mac_db[i]);
				mpz_set_ui(mac_db[i], 0);
				for (int j = 0; j < mat_sz; j++)
				{
					int index2 = i * mat_sz + j;
					mpz_init(f_t[index2]);
					mpz_set(f_t[index2], f[j * mat_sz + i]);
					mpz_init(db[index2]);
					mpz_set_ui(db[index2], 0);
					mpz_init(de[index2]);
					mpz_set_ui(de[index2], 0);
					mpz_mul(multc_1[i], d[i * mat_sz + j], mac_b[j]);
					mpz_add(mac_db[i], mac_db[i], multc_1[i]);
					for (int k = 0; k < mat_sz; k++)
					{
						mpz_mul(multc_2[index2], d[i * mat_sz + k], b[k * mat_sz + j]);
						mpz_add(db[index2], db[index2], multc_2[index2]);
						mpz_mul(multc_3[index2], d[i * mat_sz + k], e[k * mat_sz + j]);
						mpz_add(de[index2], de[index2], multc_3[index2]);
					}
					mpz_mod(db[index2], db[index2], ppp);
					mpz_mod(de[index2], de[index2], ppp);
				}
				mpz_mod(mac_db[i], mac_db[i], ppp);
			}
			omp_set_num_threads(NUM_THREADS);
			#pragma omp parallel for
			for (int i = 0; i < mat_sz; i++)
			{
				mpz_init(mac_de[i]);
				mpz_set_ui(mac_de[i], 0);
				mpz_init(mac_f_t[i]);
				mpz_set_ui(mac_f_t[i], 0);
				for (int j = 0; j < mat_sz; j++)
				{
					mpz_mul(multc_1[i], de[i * mat_sz + j], key[j]);
					mpz_add(mac_de[i], mac_de[i], multc_1[i]);
					mpz_mul(multc_2[i], f_t[i * mat_sz + j], key[j]);
					mpz_add(mac_f_t[i], mac_f_t[i], multc_2[i]);
				}
				mpz_mod(mac_de[i], mac_de[i], ppp);
				mpz_mod(mac_f_t[i], mac_f_t[i], ppp);
			}
		}
		else
		{
			omp_set_num_threads(NUM_THREADS);
			#pragma omp parallel for
			for (int i = 0; i < mat_sz; i++)
			{
				mpz_init(mac_db[i]);
				mpz_set_ui(mac_db[i], 0);
				for (int j = 0; j < mat_sz; j++)
				{
					int index2 = i * mat_sz + j;
					mpz_init(f_t[index2]);
					mpz_set(f_t[index2], f[j * mat_sz + i]);
					mpz_init(db[index2]);
					mpz_set_ui(db[index2], 0);
					mpz_init(de[index2]);
					mpz_set_ui(de[index2], 0);
					mpz_mul(multc_1[i], d[i * mat_sz + j], mac_b[j]);
					mpz_add(mac_db[i], mac_db[i], multc_1[i]);
					for (int k = 0; k < mat_sz; k++)
					{
						mpz_mul(multc_2[index2], d[i * mat_sz + k], b[k * mat_sz + j]);
						mpz_add(db[index2], db[index2], multc_2[index2]);
						mpz_mul(multc_3[index2], d[i * mat_sz + k], e[k * mat_sz + j]);
						mpz_add(de[index2], de[index2], multc_3[index2]);
					}
					mpz_mod(db[index2], db[index2], ppp);
					mpz_mod(de[index2], de[index2], ppp);
				}
				mpz_mod(mac_db[i], mac_db[i], ppp);
			}
			omp_set_num_threads(NUM_THREADS);
			#pragma omp parallel for
			for (int i = 0; i < mat_sz; i++)
			{
				mpz_init(mac_de[i]);
				mpz_set_ui(mac_de[i], 0);
				mpz_init(mac_f_t[i]);
				mpz_set_ui(mac_f_t[i], 0);
				for (int j = 0; j < mat_sz; j++)
				{
					mpz_mul(multc_1[i], de[i * mat_sz + j], key[j]);
					mpz_add(mac_de[i], mac_de[i], multc_1[i]);
					mpz_mul(multc_2[i], f_t[i * mat_sz + j], key[j]);
					mpz_add(mac_f_t[i], mac_f_t[i], multc_2[i]);
				}
				mpz_mod(mac_de[i], mac_de[i], ppp);
				mpz_mod(mac_f_t[i], mac_f_t[i], ppp);
			}
		}
		omp_set_num_threads(NUM_THREADS);
		#pragma omp parallel for
		for(int i=0;i<mat_sz;i++){
			mpz_clear(multc_1[i]);
		}
		omp_set_num_threads(NUM_THREADS);
		#pragma omp parallel for
		for(int i=0;i<sz;i++){
			mpz_clear(multc_2[i]);
			mpz_clear(multc_3[i]);
		}

		omp_set_num_threads(NUM_THREADS);
		#pragma omp parallel for
		for (int i = 0; i < sz; ++i)
		{
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
		for (int i = 0; i < mat_sz; i++)
		{
			mpz_add(output_mac[i], mac_c[i], mac_db[i]);
			mpz_add(output_mac[i], output_mac[i], mac_r[i]);
			mpz_add(output_mac[i], output_mac[i], mac_de[i]);
			mpz_add(output_mac[i], output_mac[i], mac_f_t[i]);
			mpz_mod(output_mac[i], output_mac[i], ppp);
		}

		omp_set_num_threads(NUM_THREADS);
		#pragma omp parallel for
		for (int i = 0; i < sz; i++)
		{
			mpz_clear(d[i]);
			mpz_clear(e[i]);
			mpz_clear(f[i]);
			mpz_clear(f_t[i]);
			mpz_clear(db[i]);
			mpz_clear(de[i]);
			mpz_clear(multc1[i]);
		}
		omp_set_num_threads(NUM_THREADS);
		#pragma omp parallel for
		for (int i = 0; i < mat_sz; i++)
		{
			mpz_clear(mac_db[i]);
			mpz_clear(mac_de[i]);
			mpz_clear(mac_f_t[i]);
		}
		mpz_clear(ppp);
		delete[] multc1;
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

	void mac_check(mpz_t *x, mpz_t *mac_x)
	{
		mpz_t ppp;
		mpz_init(ppp);
		mpz_set_ui(ppp, 2305843009213693951UL);

		uint64_t *x_int;
		x_int = new uint64_t[sz];
		omp_set_num_threads(NUM_THREADS);
		#pragma omp parallel for
		for (int i = 0; i < sz; i++)
		{
			x_int[i] = mpz_get_ui(x[i]);
		}
		if (party != 1)
		{
			io->send_data(1, x_int, sz * sizeof(uint64_t));
			io->flush(1);
			io->recv_data(1, x_int, sz * sizeof(uint64_t));
			omp_set_num_threads(NUM_THREADS);
			#pragma omp parallel for
			for(int i=0;i<sz;i++){
				mpz_set_ui(x[i],x_int[i]);
			}
		}
		else
		{
			uint64_t *tmp[nP + 1];
			for (int i = 1; i <= nP; ++i)
			{
				tmp[i] = new uint64_t[sz];
			}
			vector<future<void>> res;
			for (int i = 2; i <= nP; ++i)
			{
				int party2 = i;
				res.push_back(pool->enqueue([this, tmp, party2]()
											{ io->recv_data(party2, tmp[party2], sz * sizeof(uint64_t)); }));
			}
			joinNclean(res);
			mpz_t *tmp_ji;
			tmp_ji=new mpz_t[sz];
			omp_set_num_threads(NUM_THREADS);
			#pragma omp parallel for
			for(int i=0;i<sz;i++){
				mpz_init(tmp_ji[i]);
			}
			omp_set_num_threads(NUM_THREADS);
			#pragma omp parallel for
			for (int i = 0; i < sz; ++i)
			{
				for (int j = 2; j <= nP; ++j)
				{
					mpz_set_ui(tmp_ji[i],tmp[j][i]);
					mpz_add(x[i],x[i],tmp_ji[i]);
				}
				mpz_mod(x[i], x[i], ppp);
			}
			omp_set_num_threads(NUM_THREADS);
			#pragma omp parallel for
			for (int i = 0; i < sz; i++)
			{
				mpz_clear(tmp_ji[i]);
				x_int[i]=mpz_get_ui(x[i]);
			}
			for (int i = 2; i <= nP; ++i)
			{
				int party2 = i;
				res.push_back(pool->enqueue([this, x_int, party2]()
											{
					io->send_data(party2, x_int, sz * sizeof(uint64_t));
					io->flush(party2); }));
			}
			joinNclean(res);
			for (int i = 1; i <= nP; ++i)
			{
				delete[] tmp[i];
			}
		}

		mpz_t *sigma_x = new mpz_t[mat_sz];
		mpz_t *kk;
		kk =new mpz_t[mat_sz];
		omp_set_num_threads(NUM_THREADS);
		#pragma omp parallel for
		for(int i=0;i<mat_sz;i++){
			mpz_init(kk[i]);
		}
		omp_set_num_threads(NUM_THREADS);
		#pragma omp parallel for
		for (int i = 0; i < mat_sz; i++)
		{
			mpz_init(sigma_x[i]);
			mpz_set_ui(sigma_x[i], 0);
			for (int j = 0; j < mat_sz; j++)
			{
				mpz_mul(kk[i], x[i * mat_sz + j], key[j]);
				mpz_add(sigma_x[i], sigma_x[i], kk[i]);
			}
			mpz_sub(sigma_x[i], mac_x[i], sigma_x[i]);
			mpz_mod(sigma_x[i], sigma_x[i], ppp);
		}

		uint64_t *sigma_x_int;
		sigma_x_int = new uint64_t[mat_sz];
		omp_set_num_threads(NUM_THREADS);
		#pragma omp parallel for
		for(int i = 0; i < mat_sz; i++){
			sigma_x_int[i]=mpz_get_ui(sigma_x[i]);
		}
		if (party != 1)
		{
			io->send_data(1, sigma_x_int, mat_sz * sizeof(uint64_t));
			io->flush(1);
		}
		else
		{
			uint64_t *tmp[nP + 1];
			for (int i = 1; i <= nP; ++i)
			{
				tmp[i] = new uint64_t[mat_sz];
			}
			vector<future<void>> res;
			for (int i = 2; i <= nP; ++i)
			{
				int party2 = i;
				res.push_back(pool->enqueue([this, tmp, party2]()
											{ io->recv_data(party2, tmp[party2], mat_sz * sizeof(uint64_t)); }));
			}
			joinNclean(res);
			mpz_t *tmp_ji;
			tmp_ji =new mpz_t[mat_sz];
			omp_set_num_threads(NUM_THREADS);
			#pragma omp parallel for
			for(int i=0;i<mat_sz;i++){
				mpz_init(tmp_ji[i]);
			}
			uint64_t sum_check = 0;
			omp_set_num_threads(NUM_THREADS);
			#pragma omp parallel for
			for (int i = 0; i < mat_sz; ++i)
			{
				for (int j = 2; j <= nP; ++j)
				{
					mpz_set_ui(tmp_ji[i],tmp[j][i]);
					mpz_add(sigma_x[i],sigma_x[i],tmp_ji[i]);
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
			for(int i=0;i<mat_sz;i++){
				mpz_clear(tmp_ji[i]);
				mpz_clear(kk[i]);
				mpz_clear(sigma_x[i]);
			}
			delete[] sigma_x;
			for (int i = 1; i <= nP; ++i)
			{
				delete[] tmp[i];
			}
		}
		delete[] x_int;
		delete[] kk;
	}
};
#endif // MSPDZ_H__