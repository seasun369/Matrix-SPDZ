#ifndef SPDZ_H__
#define SPDZ_H__
#define NUM_THREADS 8

#include "../network/network.h"
#include "utility.h"
#include <emp-tool/emp-tool.h>
#include "../network/helper.h"
#include "../mpir/mpir.h"
using namespace emp;

// base SPDZ protocol with scalar elements
template<int nP>
class SPDZ { public:
	const block MASK = makeBlock(0x0ULL, 0xFFFFFULL);
	//const uint64_t PR = 2305843009213693951;
	mpz_t key;
    mpz_t a,b,c;
	mpz_t mac_a,mac_b,mac_c;

	NetIOMP<nP> * io;
	int party, total_pre, ssp;
	ThreadPool * pool;
	int sz=1;

	PRP prp;
	SPDZ(NetIOMP<nP> * io[2], ThreadPool * pool, int party, bool * _delta = nullptr, int ssp = 40) {
		this->party = party;
		this->io = io[0];
		this->ssp = ssp;
		this->pool = pool;

		mpz_init(key);
		mpz_init(a);
		mpz_init(b);
		mpz_init(c);
		mpz_init(mac_a);
		mpz_init(mac_b);
		mpz_init(mac_c);
	}
	~SPDZ () {
		mpz_clear(a);
		mpz_clear(b);
		mpz_clear(c);
		mpz_clear(mac_a);
		mpz_clear(mac_b);
		mpz_clear(mac_c);
		mpz_clear(key);
	}
	PRG prg;

    // it should be implemented by offline.
    void get_triple(){
		std::ifstream fin1,fin2,fin3,fin4;
		std::ifstream fin_1,fin_2,fin_3;
		fin1.open("/home/jackie/spdz/pre_data/predata_spdz/a.txt",std::ios::in);
		fin2.open("/home/jackie/spdz/pre_data/predata_spdz/b.txt",std::ios::in);
		fin3.open("/home/jackie/spdz/pre_data/predata_spdz/c.txt",std::ios::in);
		fin4.open("/home/jackie/spdz/pre_data/predata_spdz/key.txt", std::ios::in);

		fin_1.open("/home/jackie/spdz/pre_data/predata_spdz/mac_a.txt", std::ios::in);
		fin_2.open("/home/jackie/spdz/pre_data/predata_spdz/mac_b.txt", std::ios::in);
		fin_3.open("/home/jackie/spdz/pre_data/predata_spdz/mac_c.txt", std::ios::in);

		if(!(fin1.is_open() && fin2.is_open() && fin3.is_open() && fin4.is_open() && fin_1.is_open() && fin_2.is_open() && fin_3.is_open()))
		{
		    std::cerr<<"cannot open the file";
		}

		char line[1024] = { 0 };
		std::vector<triple> Triple;
		while (fin1.getline(line, sizeof(line)))
		{
			triple t;
			std::stringstream word(line);
			word >> t.party;
			uint64_t num;
			while (word >> num)
				t.value.push_back(num);
			Triple.push_back(t);
		}

		for (int i = 0; i < nP; i++) {
			if (party == stoi(Triple[i].party)) {
				mpz_set_ui(a, Triple[i].value[0]);
			}
		}

		Triple.clear();
		memset(line, 0, sizeof(line));

		while (fin2.getline(line, sizeof(line)))
		{
			triple t;
			std::stringstream word(line);
			word >> t.party;
			uint64_t num;
			while (word >> num)
				t.value.push_back(num);
			Triple.push_back(t);
		}

		for (int i = 0; i < nP; i++) {
			if (party == stoi(Triple[i].party)) {
				mpz_set_ui(b, Triple[i].value[0]);
			}
		}

		Triple.clear();
		memset(line, 0, sizeof(line));

		while (fin3.getline(line, sizeof(line)))
		{
			triple t;
			std::stringstream word(line);
			word >> t.party;
			uint64_t num;
			while (word >> num)
				t.value.push_back(num);
			Triple.push_back(t);
		}

		for (int i = 0; i < nP; i++) {
			if (party == stoi(Triple[i].party)) {
				mpz_set_ui(c, Triple[i].value[0]);
			}
		}

		Triple.clear();
		memset(line, 0, sizeof(line));

		while (fin4.getline(line, sizeof(line)))
		{
			triple t;
			std::stringstream word(line);
			word >> t.party;
			uint64_t num;
			while (word >> num)
				t.value.push_back(num);
			Triple.push_back(t);
		}

		for (int i = 0; i < nP; i++) {
			if (party == stoi(Triple[i].party)) {
				mpz_set_ui(key, Triple[i].value[0]);
			}
		}

		Triple.clear();
		memset(line, 0, sizeof(line));

		while (fin_1.getline(line, sizeof(line)))
		{
			triple t;
			std::stringstream word(line);
			word >> t.party;
			uint64_t num;
			while (word >> num)
				t.value.push_back(num);
			Triple.push_back(t);
		}

		for (int i = 0; i < nP; i++) {
			if (party == stoi(Triple[i].party)) {
				mpz_set_ui(mac_a, Triple[i].value[0]);
			}
		}

		Triple.clear();
		memset(line, 0, sizeof(line));

		while (fin_2.getline(line, sizeof(line)))
		{
			triple t;
			std::stringstream word(line);
			word >> t.party;
			uint64_t num;
			while (word >> num)
				t.value.push_back(num);
			Triple.push_back(t);
		}

		for (int i = 0; i < nP; i++) {
			if (party == stoi(Triple[i].party)) {
				mpz_set_ui(mac_b, Triple[i].value[0]);
			}
		}

		Triple.clear();
		memset(line, 0, sizeof(line));

		while (fin_3.getline(line, sizeof(line)))
		{
			triple t;
			std::stringstream word(line);
			word >> t.party;
			uint64_t num;
			while (word >> num)
				t.value.push_back(num);
			Triple.push_back(t);
		}

		for (int i = 0; i < nP; i++) {
			if (party == stoi(Triple[i].party)) {
				mpz_set_ui(mac_c, Triple[i].value[0]);
			}
		}

		Triple.clear();
		memset(line, 0, sizeof(line));
	}

	void Online_mul (mpz_t x, mpz_t y, mpz_t mac_x, mpz_t mac_y, mpz_t output, mpz_t output_mac) {
		mpz_t ppp;
		mpz_init(ppp);
		mpz_set_ui(ppp, 2305843009213693951UL);
		
		mpz_t d;
		mpz_t e;

		mpz_init(d);
		mpz_init(e);
		mpz_sub(d, x, a);
		mpz_mod(d, d, ppp);
		mpz_sub(e, y, b);
		mpz_mod(e, e, ppp);

		uint64_t *d_int;
		d_int = new uint64_t[sz];
		for (int i = 0; i < sz; i++)
		{
			d_int[i]=mpz_get_ui(d);
		}

		if(party != 1) {
			io->send_data(1, d_int, sz*sizeof(uint64_t));
			io->flush(1);
			io->recv_data(1, d_int, sz*sizeof(uint64_t));
			mpz_set_ui(d,d_int[0]);
		} else {
			uint64_t * tmp[nP+1];
			for (int i = 1; i <= nP; ++i)
			{
				tmp[i] = new uint64_t[sz];
			}
			vector<future<void>> res;
			for(int i = 2; i <= nP; ++i) {
				int party2 = i;
				res.push_back(pool->enqueue([this, tmp, party2]() {
					io->recv_data(party2, tmp[party2], sz*sizeof(uint64_t));
				}));
			}
			joinNclean(res);
			mpz_t tmp_ji;
			mpz_init(tmp_ji);
			for (int j = 2; j <= nP; ++j)
			{
				mpz_set_ui(tmp_ji,tmp[j][0]);
				mpz_add(d,d,tmp_ji);
			}
			mpz_mod(d, d, ppp);
			mpz_clear(tmp_ji);
			d_int[0]=mpz_get_ui(d);
			for(int i = 2; i <= nP; ++i) {
				int party2 = i;
				res.push_back(pool->enqueue([this, d_int, party2]() {
					io->send_data(party2, d_int, sz*sizeof(uint64_t));
					io->flush(party2);
				}));
			}
			joinNclean(res);
			for (int i = 1; i <= nP; ++i)
			{
				delete[] tmp[i];
			}
		}
		delete[] d_int;
		
		uint64_t *e_int;
		e_int = new uint64_t[sz];
		for (int i = 0; i < sz; i++)
		{
			e_int[i]=mpz_get_ui(e);
		}

		if(party != 1) {
			io->send_data(1, e_int, sz*sizeof(uint64_t));
			io->flush(1);
			io->recv_data(1, e_int, sz*sizeof(uint64_t));
			mpz_set_ui(e,e_int[0]);
		} else {
			uint64_t * tmp[nP+1];
			for (int i = 1; i <= nP; ++i)
			{
				tmp[i] = new uint64_t[sz];
			}
			vector<future<void>> res;
			for(int i = 2; i <= nP; ++i) {
				int party2 = i;
				res.push_back(pool->enqueue([this, tmp, party2]() {
					io->recv_data(party2, tmp[party2], sz*sizeof(uint64_t));
				}));
			}
			joinNclean(res);
			mpz_t tmp_ji;
			mpz_init(tmp_ji);
			for (int j = 2; j <= nP; ++j)
			{
				mpz_set_ui(tmp_ji,tmp[j][0]);
				mpz_add(e,e,tmp_ji);
			}
			mpz_mod(e, e, ppp);
			mpz_clear(tmp_ji);
			e_int[0]=mpz_get_ui(e);
			for(int i = 2; i <= nP; ++i) {
				int party2 = i;
				res.push_back(pool->enqueue([this, e_int, party2]() {
					io->send_data(party2, e_int, sz*sizeof(uint64_t));
					io->flush(party2);
				}));
			}
			joinNclean(res);
			for (int i = 1; i <= nP; ++i)
			{
				delete[] tmp[i];
			}
		}
		delete[] e_int;
		
		mpz_t db;
		mpz_t ae;
		mpz_t de;
		mpz_t mac_db;
		mpz_t mac_ae;
		mpz_t mac_de;

		mpz_init(db);
		mpz_init(ae);
		mpz_init(de);
		mpz_init(mac_db);
		mpz_init(mac_ae);
		mpz_init(mac_de);

		if(party==1){
			mpz_mul(de,d,e);
			mpz_mul(db,d,b);
			mpz_mul(ae,a,e);
			mpz_mul(mac_de,de,key);
			mpz_mul(mac_db,d,mac_b);
			mpz_mul(mac_ae,mac_a,e);
			mpz_add(output,de,db);
			mpz_add(output,output,ae);
			mpz_add(output,output,c);
			mpz_mod(output,output,ppp);
			mpz_add(output_mac,mac_de,mac_db);
			mpz_add(output_mac,output_mac,mac_ae);
			mpz_add(output_mac,output_mac,mac_c);
			mpz_mod(output_mac,output_mac,ppp);
		}else{
			mpz_mul(de,d,e);
			mpz_mul(db,d,b);
			mpz_mul(ae,a,e);
			mpz_mul(mac_de,de,key);
			mpz_mul(mac_db,d,mac_b);
			mpz_mul(mac_ae,mac_a,e);
			mpz_add(output,db,ae);
			mpz_add(output,output,c);
			mpz_mod(output,output,ppp);
			mpz_add(output_mac,mac_de,mac_db);
			mpz_add(output_mac,output_mac,mac_ae);
			mpz_add(output_mac,output_mac,mac_c);
			mpz_mod(output_mac,output_mac,ppp);
		}
		
		mpz_clear(db);
		mpz_clear(ae);
		mpz_clear(de);
		mpz_clear(mac_db);
		mpz_clear(mac_ae);
		mpz_clear(mac_de);
		mpz_clear(d);
		mpz_clear(e);
	}	
	
    void mac_check(mpz_t x, mpz_t mac_x){
		mpz_t ppp;
		mpz_init(ppp);
		mpz_set_ui(ppp, 2305843009213693951UL);

		uint64_t *x_int;
		x_int = new uint64_t[sz];
		for (int i = 0; i < sz; i++)
		{
			x_int[i]=mpz_get_ui(x);
		}
		if (party != 1)
		{
			io->send_data(1, x_int, sz * sizeof(uint64_t));
			io->flush(1);
			io->recv_data(1, x_int, sz * sizeof(uint64_t));
			mpz_set_ui(x,x_int[0]);
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
			for (int j = 2; j <= nP; ++j)
			{
				mpz_set_ui(tmp_ji,tmp[j][0]);
				mpz_add(x,x,tmp_ji);
			}
			mpz_mod(x, x, ppp);
			mpz_clear(tmp_ji);
			x_int[0]=mpz_get_ui(x);
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
		delete[] x_int;

        mpz_t sigma_x;
		mpz_init(sigma_x);
		mpz_mul(sigma_x, x,key);
		mpz_sub(sigma_x,mac_x,sigma_x);
		mpz_mod(sigma_x, sigma_x, ppp);

		uint64_t *sigma_x_int;
		sigma_x_int = new uint64_t[sz];
		for (int i = 0; i < sz; i++)
		{
			sigma_x_int[i]=mpz_get_ui(sigma_x);
		}
		
        if(party != 1) {
			io->send_data(1, sigma_x_int, sz*sizeof(uint64_t));
			io->flush(1);
		} else {
			uint64_t * tmp[nP+1];
			for (int i = 1; i <= nP; ++i)
			{
				tmp[i] = new uint64_t[sz];
			}
			vector<future<void>> res;
			for(int i = 2; i <= nP; ++i) {
				int party2 = i;
				res.push_back(pool->enqueue([this, tmp, party2]() {
					io->recv_data(party2, tmp[party2], sz*sizeof(uint64_t));
				}));
			}
			joinNclean(res);
			mpz_t tmp_ji;
			mpz_init(tmp_ji);
			uint64_t sum_check = 0;
			for (int j = 2; j <= nP; ++j){
				mpz_set_ui(tmp_ji,tmp[j][0]);
				mpz_add(sigma_x,sigma_x,tmp_ji);
			}
			mpz_mod(sigma_x, sigma_x, ppp);
			if (mpz_sgn(sigma_x) != 0) {
				std::cout << "check error" << endl;
				sum_check = sum_check + 1;
			}
			if (sum_check == 0) std::cout << "Correct!" << endl;
			mpz_clear(ppp);
			mpz_clear(tmp_ji);
			mpz_clear(sigma_x);
			for (int i = 1; i <= nP; ++i) delete[] tmp[i];
		}
    }
};
#endif// SPDZ_H__