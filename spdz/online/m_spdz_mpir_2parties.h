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
					mpz_set_ui(a[j], Triple[i].value[j]);
				}
			}
		}

		Triple.clear();
		delete[] line1;
		fin1.close();
    	char* line2 = new(std::nothrow) char[arraySize];

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
					mpz_set_ui(a_t[j], Triple[i].value[j]);
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
				for (int j = 0; j < sz; j++)
					mpz_set_ui(b[j], Triple[i].value[j]);
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
				for (int j = 0; j < sz; j++)
					mpz_set_ui(c[j], Triple[i].value[j]);
			}
		}

		Triple.clear();
		delete[] line4;
		fin4.close();
    	char* line5 = new(std::nothrow) char[arraySize];

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
				for (int j = 0; j < sz; j++)
					mpz_set_ui(r[j], Triple[i].value[j]);
			}
		}

		Triple.clear();
		delete[] line5;
		fin5.close();
    	char* line6 = new(std::nothrow) char[arraySize];

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
				for (int j = 0; j < sz; j++)
					mpz_set_ui(r_t[j], Triple[i].value[j]);
			}
		}

		Triple.clear();
		delete[] line6;
		fin6.close();
    	char* line7 = new(std::nothrow) char[arraySize];

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
				for (int j = 0; j < mat_sz; j++)
					mpz_set_ui(key[j], Triple[i].value[j]);
			}
		}

		Triple.clear();
		delete[] line7;
		fin7.close();
    	char* line_1 = new(std::nothrow) char[arraySize];

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
				for (int j = 0; j < mat_sz; j++)
					mpz_set_ui(mac_a[j], Triple[i].value[j]);
			}
		}

		Triple.clear();
		delete[] line_1;
		fin_1.close();
    	char* line_2 = new(std::nothrow) char[arraySize];

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
				for (int j = 0; j < mat_sz; j++)
					mpz_set_ui(mac_a_t[j], Triple[i].value[j]);
			}
		}

		Triple.clear();
		delete[] line_2;
		fin_2.close();
    	char* line_3 = new(std::nothrow) char[arraySize];

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
				for (int j = 0; j < mat_sz; j++)
					mpz_set_ui(mac_b[j], Triple[i].value[j]);
			}
		}

		Triple.clear();
		delete[] line_3;
		fin_3.close();
    	char* line_4 = new(std::nothrow) char[arraySize];

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
				for (int j = 0; j < mat_sz; j++)
					mpz_set_ui(mac_c[j], Triple[i].value[j]);
			}
		}

		Triple.clear();
		delete[] line_4;
		fin_4.close();
    	char* line_5 = new(std::nothrow) char[arraySize];

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
				for (int j = 0; j < mat_sz; j++)
					mpz_set_ui(mac_r[j], Triple[i].value[j]);
			}
		}

		Triple.clear();
		delete[] line_5;
		fin_5.close();
    	char* line_6 = new(std::nothrow) char[arraySize];

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
				for (int j = 0; j < mat_sz; j++)
					mpz_set_ui(mac_r_t[j], Triple[i].value[j]);
			}
		}

		Triple.clear();
		delete[] line_6;
		fin_6.close();

		//delete[] line1; 
	}

	void Online_mul(mpz_t *x, mpz_t *y, mpz_t *mac_x, mpz_t *mac_y, mpz_t *output, mpz_t *output_mac)
	{
		mpz_t ppp;
		mpz_init(ppp);
		mpz_set_ui(ppp, 2305843009213693951UL);

		mpz_t *d = new mpz_t[sz];
		mpz_t *e = new mpz_t[sz];
		mpz_t *f = new mpz_t[sz];
		uint64_t *d_int = new uint64_t[sz], *e_int = new uint64_t[sz], *f_int = new uint64_t[sz];

		omp_set_num_threads(NUM_THREADS);
		#pragma omp parallel for
		for (int i = 0; i < sz; i++)
		{
			mpz_init(d[i]);
			mpz_sub(d[i], x[i], a[i]);
			mpz_mod(d[i], d[i], ppp);
			d_int[i]=mpz_get_ui(d[i]);
			mpz_init(e[i]);
			mpz_sub(e[i], y[i], b[i]);
			mpz_mod(e[i], e[i], ppp);
			e_int[i]=mpz_get_ui(e[i]);
			mpz_init(f[i]);
			mpz_set_ui(f[i], 0);
		}
		
		if (party != 1)
		{
			io->send_data(1, d_int, sz * sizeof(uint64_t));
			io->flush(1);
			io->recv_data(1, d_int, sz * sizeof(uint64_t));
			omp_set_num_threads(NUM_THREADS);
			#pragma omp parallel for
			for(int i = 0; i < sz; i++){
				mpz_set_ui(d[i], d_int[i]);
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
			int party2 = 2;
			res.push_back(pool->enqueue([this, tmp, party2]()
										{ io->recv_data(party2, tmp[party2], sz * sizeof(uint64_t)); }));
			joinNclean(res);

			mpz_t *tmp_i = new mpz_t[sz];
			omp_set_num_threads(NUM_THREADS);
			#pragma omp parallel for
			for (int i = 0; i < sz; ++i)
			{
				mpz_init(tmp_i[i]);
				mpz_set_ui(tmp_i[i],tmp[2][i]);
				mpz_add(d[i],d[i],tmp_i[i]);
				mpz_mod(d[i], d[i], ppp);
				d_int[i]=mpz_get_ui(d[i]);
				//mpz_clear(tmp_i[i]);
			}
			//delete[] tmp_i;

			res.push_back(pool->enqueue([this, d_int, party2]()
											{
					io->send_data(party2, d_int, sz * sizeof(uint64_t));
					io->flush(party2); }));
			joinNclean(res);

			for (int i = 1; i <= nP; ++i)
			{
				delete[] tmp[i];
			}
		}
		
		if (party != 1)
		{
			io->send_data(1, e_int, sz * sizeof(uint64_t));
			io->flush(1);
			io->recv_data(1, e_int, sz * sizeof(uint64_t));
			omp_set_num_threads(NUM_THREADS);
			#pragma omp parallel for
			for(int i = 0; i < sz; i++){
				mpz_set_ui(e[i], e_int[i]);
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
			int party2 = 2;
			res.push_back(pool->enqueue([this, tmp, party2]()
										{ io->recv_data(party2, tmp[party2], sz * sizeof(uint64_t)); }));
			joinNclean(res);

			//mpz_t *tmp_i = new mpz_t[sz];
			omp_set_num_threads(NUM_THREADS);
			#pragma omp parallel for
			for (int i = 0; i < sz; ++i)
			{
				//mpz_init(tmp_i[i]);
				mpz_set_ui(tmp_i[i],tmp[2][i]);
				mpz_add(e[i],e[i],tmp_i[i]);
				mpz_mod(e[i], e[i], ppp);
				e_int[i]=mpz_get_ui(e[i]);
				//mpz_clear(tmp_i[i]);
			}
			//delete[] tmp_i;

			res.push_back(pool->enqueue([this, e_int, party2]()
											{
					io->send_data(party2, e_int, sz * sizeof(uint64_t));
					io->flush(party2); }));
			joinNclean(res);
			
			for (int i = 1; i <= nP; ++i)
			{
				delete[] tmp[i];
			}
		}
		
		mpz_t *multc1;
		multc1=new mpz_t[sz];

		//compute f[sz]
		omp_set_num_threads(NUM_THREADS);
		#pragma omp parallel for
		for (int i = 0; i < mat_sz; i++)
		{
			for (int j = 0; j < mat_sz; j++)
			{
				int index1 = i * mat_sz + j;
				mpz_init(multc1[index1]);
				for (int k = 0; k < mat_sz; k++)
				{
					mpz_mul(multc1[index1], e[k * mat_sz + i], a_t[k * mat_sz + j]);
					mpz_add(f[index1], f[index1], multc1[index1]);
				}
				mpz_sub(f[index1], f[index1], r_t[index1]);
				mpz_mod(f[index1], f[index1], ppp);
				f_int[index1]=mpz_get_ui(f[index1]);
			}
		}

		if (party != 1)
		{
			io->send_data(1, f_int, sz * sizeof(uint64_t));
			io->flush(1);
			io->recv_data(1, f_int, sz * sizeof(uint64_t));
			omp_set_num_threads(NUM_THREADS);
			#pragma omp parallel for
			for(int i = 0; i < sz; i++){
				mpz_set_ui(f[i], f_int[i]);
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
			int party2 = 2;
			res.push_back(pool->enqueue([this, tmp, party2]()
										{ io->recv_data(party2, tmp[party2], sz * sizeof(uint64_t)); }));
			joinNclean(res);

			//mpz_t *tmp_i = new mpz_t[sz];
			omp_set_num_threads(NUM_THREADS);
			#pragma omp parallel for
			for (int i = 0; i < sz; ++i)
			{
				//mpz_init(tmp_i[i]);
				mpz_set_ui(tmp_i[i],tmp[2][i]);
				mpz_add(f[i],f[i],tmp_i[i]);
				mpz_mod(f[i], f[i], ppp);
				f_int[i]=mpz_get_ui(f[i]);
				mpz_clear(tmp_i[i]);
			}
			delete[] tmp_i;

			res.push_back(pool->enqueue([this, f_int, party2]()
											{
					io->send_data(party2, f_int, sz * sizeof(uint64_t));
					io->flush(party2); }));
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
		mpz_t *multc_1 = new mpz_t[mat_sz];
		mpz_t *multc_2 = new mpz_t[sz];
		mpz_t *multc_3 = new mpz_t[sz];

		omp_set_num_threads(NUM_THREADS);
		#pragma omp parallel for
		for(int i = 0; i < mat_sz; i++){
			mpz_init(multc_1[i]);
		}

		omp_set_num_threads(NUM_THREADS);
		#pragma omp parallel for
		for(int i=0;i<sz;i++){
			mpz_init(multc_2[i]);
			mpz_init(multc_3[i]);
		}

		if (party == 1)
		{
			//compute db[sz],de[sz],mac_db[mat_sz],f_t[sz],mac_de[mat_sz],mac_f_t[mat_sz]
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
					mpz_init(mac_de[i]);
					mpz_set_ui(mac_de[i], 0);
					mpz_init(mac_f_t[i]);
					mpz_set_ui(mac_f_t[i], 0);
					for (int k = 0; k < mat_sz; k++)
					{
						mpz_mul(multc_2[index2], d[i * mat_sz + k], b[k * mat_sz + j]);
						mpz_add(db[index2], db[index2], multc_2[index2]);
						mpz_mul(multc_3[index2], d[i * mat_sz + k], e[k * mat_sz + j]);
						mpz_add(de[index2], de[index2], multc_3[index2]);
					}
					mpz_mod(db[index2], db[index2], ppp);
					mpz_mod(de[index2], de[index2], ppp);
					mpz_mul(multc_1[i], de[i * mat_sz + j], key[j]);
					mpz_add(mac_de[i], mac_de[i], multc_1[i]);
					mpz_mul(multc_1[i], f_t[i * mat_sz + j], key[j]);
					mpz_add(mac_f_t[i], mac_f_t[i], multc_1[i]);
				}
				mpz_mod(mac_db[i], mac_db[i], ppp);
				mpz_mod(mac_de[i], mac_de[i], ppp);
				mpz_mod(mac_f_t[i], mac_f_t[i], ppp);
			}
		}
		else
		{
			//compute db[sz],de[sz],mac_db[mat_sz],f_t[sz],mac_de[mat_sz],mac_f_t[mat_sz]
			omp_set_num_threads(NUM_THREADS);
			#pragma omp parallel for
			for (int i = 0; i < mat_sz; i++)
			{
				mpz_init(mac_db[i]);
				mpz_set_ui(mac_db[i], 0);
				mpz_init(mac_de[i]);
				mpz_set_ui(mac_de[i], 0);
				mpz_init(mac_f_t[i]);
				mpz_set_ui(mac_f_t[i], 0);
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
					mpz_mul(multc_1[i], de[i * mat_sz + j], key[j]);
					mpz_add(mac_de[i], mac_de[i], multc_1[i]);
					mpz_mul(multc_1[i], f_t[i * mat_sz + j], key[j]);
					mpz_add(mac_f_t[i], mac_f_t[i], multc_1[i]);
				}
				mpz_mod(mac_db[i], mac_db[i], ppp);
				mpz_mod(mac_de[i], mac_de[i], ppp);
				mpz_mod(mac_f_t[i], mac_f_t[i], ppp);
			}
		}

		omp_set_num_threads(NUM_THREADS);
		#pragma omp parallel for
		for(int i = 0; i < mat_sz; i++){
			mpz_clear(multc_1[i]);
		}

		omp_set_num_threads(NUM_THREADS);
		#pragma omp parallel for
		for (int i = 0; i < sz; ++i)
		{
			if (party == 1)
			{
				mpz_clear(multc_2[i]);
				mpz_clear(multc_3[i]);
				mpz_add(output[i], c[i], db[i]);
				mpz_add(output[i], output[i], r[i]);
				mpz_add(output[i], output[i], de[i]);
				mpz_add(output[i], output[i], f_t[i]);
				mpz_mod(output[i], output[i], ppp);
			}
			else
			{
				mpz_clear(multc_2[i]);
				mpz_clear(multc_3[i]);
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
		delete[] multc2;
		delete[] multc3;
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

		uint64_t *x_int= new uint64_t[sz];

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
			for(int i = 0; i < sz; i++){
				mpz_set_ui(x[i], x_int[i]);
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
			int party2 = 2;
			res.push_back(pool->enqueue([this, tmp, party2]()
										{ io->recv_data(party2, tmp[party2], sz * sizeof(uint64_t)); }));
			joinNclean(res);

			mpz_t *tmp_i = new mpz_t[sz];
			omp_set_num_threads(NUM_THREADS);
			#pragma omp parallel for
			for (int i = 0; i < sz; ++i)
			{
				mpz_init(tmp_i[i]);
				mpz_set_ui(tmp_i[i],tmp[2][i]);
				mpz_add(x[i],x[i],tmp_i[i]);
				mpz_mod(x[i], x[i], ppp);
				x_int[i]=mpz_get_ui(x[i]);
				mpz_clear(tmp_i[i]);
			}
			//delete[] tmp_i;

			res.push_back(pool->enqueue([this, x_int, party2]()
											{
					io->send_data(party2, x_int, sz * sizeof(uint64_t));
					io->flush(party2); }));
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
		for(int i = 0; i < mat_sz; i++){
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

		uint64_t *sigma_x_int= new uint64_t[mat_sz];
		omp_set_num_threads(NUM_THREADS);
		#pragma omp parallel for
		for(int i = 0; i < mat_sz; i++){
			sigma_x_int[i] = mpz_get_ui(sigma_x[i]);
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
			int party2 = 2;
			res.push_back(pool->enqueue([this, tmp, party2]()
										{ io->recv_data(party2, tmp[party2], mat_sz * sizeof(uint64_t)); }));
			joinNclean(res);

			uint64_t sum_check = 0;
			omp_set_num_threads(NUM_THREADS);
			#pragma omp parallel for
			for (int i = 0; i < mat_sz; ++i)
			{
				for (int j = 2; j <= nP; ++j)
				{
					mpz_set_ui(tmp_i[i],tmp[j][i]);
					mpz_add(sigma_x[i],sigma_x[i],tmp_i[i]);
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
			for(int i = 0; i < mat_sz; i++){
				mpz_clear(tmp_i[i]);
				mpz_clear(kk[i]);
				mpz_clear(sigma_x[i]);
			}
			for (int i = 1; i <= nP; ++i)
			{
				delete[] tmp[i];
			}
		}
		delete[] x_int;
		delete[] kk;
		delete[] sigma_x;
		delete[] tmp_i;
	}
};
#endif // MSPDZ_H__