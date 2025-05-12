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
    mpz_t *a,*b,*c;
	mpz_t *mac_a,*mac_b,*mac_c;

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

		//mpz_t key;
        a = new mpz_t[sz];
        b = new mpz_t[sz];
        c = new mpz_t[sz];
		mac_a = new mpz_t[sz];
        mac_b = new mpz_t[sz];
        mac_c = new mpz_t[sz];

		mpz_init(key);
		for (int i = 0; i < sz; i++)
		{
			mpz_init(a[i]);
			mpz_init(b[i]);
			mpz_init(c[i]);
			mpz_init(mac_a[i]);
			mpz_init(mac_b[i]);
			mpz_init(mac_c[i]);
		}
	}
	~SPDZ () {
		for (int i = 0; i < sz; i++)
		{
			mpz_clear(a[i]);
			mpz_clear(b[i]);
			mpz_clear(c[i]);
			mpz_clear(mac_a[i]);
			mpz_clear(mac_b[i]);
			mpz_clear(mac_c[i]);
		}
		mpz_clear(key);
		delete[] a;
		delete[] b;
		delete[] c;
		delete[] mac_a;
		delete[] mac_b;
		delete[] mac_c;
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

		char line[10240] = { 0 };
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
				for (int j = 0; j < sz; j++) mpz_set_ui(a[j], Triple[i].value[j]);
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
				for (int j = 0; j < sz; j++) mpz_set_ui(b[j], Triple[i].value[j]);
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
				for (int j = 0; j < sz; j++) mpz_set_ui(c[j], Triple[i].value[j]);
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
				for (int j = 0; j < sz; j++) mpz_set_ui(mac_a[j], Triple[i].value[j]);
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
				for (int j = 0; j < sz; j++) mpz_set_ui(mac_b[j], Triple[i].value[j]);
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
				for (int j = 0; j < sz; j++) mpz_set_ui(mac_c[j], Triple[i].value[j]);
			}
		}

		Triple.clear();
		memset(line, 0, sizeof(line));

	}

	void Online_mul (mpz_t* x, mpz_t* y, mpz_t* mac_x, mpz_t* mac_y, mpz_t* output, mpz_t* output_mac) {
		mpz_t ppp;
		mpz_init(ppp);
		mpz_set_ui(ppp, 2305843009213693951UL);
		
		mpz_t *d = new mpz_t[sz];
		mpz_t *e = new mpz_t[sz];

		for (int i = 0; i < sz; i++)
		{
			mpz_init(d[i]);
			mpz_init(e[i]);
			mpz_sub(d[i], x[i], a[i]);
			mpz_mod(d[i], d[i], ppp);
			mpz_sub(e[i], y[i], b[i]);
			mpz_mod(e[i], e[i], ppp);
		}


		uint64_t *d_int;
		d_int = new uint64_t[sz];
		for (int i = 0; i < sz; i++)
		{
			d_int[i]=mpz_get_ui(d[i]);
		}
		if(party != 1) {
			io->send_data(1, d_int, sz*sizeof(uint64_t));
			io->flush(1);
			io->recv_data(1, d_int, sz*sizeof(uint64_t));
			for(int i=0;i<sz;i++){
				mpz_set_ui(d[i],d_int[i]);
			}
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
			mpz_t *tmp_ji;
			tmp_ji=new mpz_t[sz];
			for(int i=0;i<sz;i++){
				mpz_init(tmp_ji[i]);
			}
			for (int i = 0; i < sz; ++i)
			{
				for (int j = 2; j <= nP; ++j)
				{
					mpz_set_ui(tmp_ji[i],tmp[j][i]);
					mpz_add(d[i],d[i],tmp_ji[i]);
				}
				mpz_mod(d[i], d[i], ppp);
			}
			for (int i = 0; i < sz; i++)
			{
				mpz_clear(tmp_ji[i]);
				d_int[i]=mpz_get_ui(d[i]);
			}
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

		uint64_t *e_int;
		e_int = new uint64_t[sz];
		for (int i = 0; i < sz; i++)
		{
			e_int[i]=mpz_get_ui(e[i]);
		}
		if(party != 1) {
			io->send_data(1, e_int, sz*sizeof(uint64_t));
			io->flush(1);
			io->recv_data(1, e_int, sz*sizeof(uint64_t));
			for(int i=0;i<sz;i++){
				mpz_set_ui(e[i],e_int[i]);
			}
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
			mpz_t *tmp_ji;
			tmp_ji=new mpz_t[sz];
			for(int i=0;i<sz;i++){
				mpz_init(tmp_ji[i]);
			}
			for (int i = 0; i < sz; ++i)
			{
				for (int j = 2; j <= nP; ++j)
				{
					mpz_set_ui(tmp_ji[i],tmp[j][i]);
					mpz_add(e[i],e[i],tmp_ji[i]);
				}
				mpz_mod(e[i], e[i], ppp);
			}
			for (int i = 0; i < sz; i++)
			{
				mpz_clear(tmp_ji[i]);
				e_int[i]=mpz_get_ui(e[i]);
			}
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
			mpz_mul(de,d[0],e[0]);
			mpz_mul(db,d[0],b[0]);
			mpz_mul(ae,a[0],e[0]);
			mpz_mul(mac_de,de,key);
			mpz_mul(mac_db,d[0],mac_b[0]);
			mpz_mul(mac_ae,mac_a[0],e[0]);
			mpz_add(output[0],de,db);
			mpz_add(output[0],output[0],ae);
			mpz_add(output[0],output[0],c[0]);
			mpz_mod(output[0],output[0],ppp);
			mpz_add(output_mac[0],mac_de,mac_db);
			mpz_add(output_mac[0],output_mac[0],mac_ae);
			mpz_add(output_mac[0],output_mac[0],mac_c[0]);
			mpz_mod(output_mac[0],output_mac[0],ppp);
		}else{
			mpz_mul(de,d[0],e[0]);
			mpz_mul(db,d[0],b[0]);
			mpz_mul(ae,a[0],e[0]);
			mpz_mul(mac_de,de,key);
			mpz_mul(mac_db,d[0],mac_b[0]);
			mpz_mul(mac_ae,mac_a[0],e[0]);
			mpz_add(output[0],db,ae);
			mpz_add(output[0],output[0],c[0]);
			mpz_mod(output[0],output[0],ppp);
			mpz_add(output_mac[0],mac_de,mac_db);
			mpz_add(output_mac[0],output_mac[0],mac_ae);
			mpz_add(output_mac[0],output_mac[0],mac_c[0]);
			mpz_mod(output_mac[0],output_mac[0],ppp);
		}
		
		mpz_clear(db);
		mpz_clear(ae);
		mpz_clear(de);
		mpz_clear(mac_db);
		mpz_clear(mac_ae);
		mpz_clear(mac_de);
		mpz_clear(d[0]);
		mpz_clear(e[0]);
		delete[] d;
		delete[] e;

	}	
	
    void mac_check(mpz_t* x, mpz_t* mac_x){
		mpz_t ppp;
		mpz_init(ppp);
		mpz_set_ui(ppp, 2305843009213693951UL);

        uint64_t *x_int;
		x_int = new uint64_t[sz];
		for (int i = 0; i < sz; i++)
		{
			x_int[i] = mpz_get_ui(x[i]);
		}
		if (party != 1)
		{
			io->send_data(1, x_int, sz * sizeof(uint64_t));
			io->flush(1);
			io->recv_data(1, x_int, sz * sizeof(uint64_t));
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
			for(int i=0;i<sz;i++){
				mpz_init(tmp_ji[i]);
			}
			for (int i = 0; i < sz; ++i)
			{
				for (int j = 2; j <= nP; ++j)
				{
					mpz_set_ui(tmp_ji[i],tmp[j][i]);
					mpz_add(x[i],x[i],tmp_ji[i]);
				}
				mpz_mod(x[i], x[i], ppp);
			}
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

        mpz_t* sigma_x = new mpz_t[sz];
		for(int i=0;i<sz;i++){
			mpz_init(sigma_x[i]);
			mpz_mul(sigma_x[i], x[i],key);
			mpz_sub(sigma_x[i],mac_x[i],sigma_x[i]);
			mpz_mod(sigma_x[i], sigma_x[i], ppp);
		}

		uint64_t *sigma_x_int;
		sigma_x_int = new uint64_t[sz];
		for(int i = 0; i < sz; i++){
			sigma_x_int[i]=mpz_get_ui(sigma_x[i]);
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
			//for (int i = 0; i < sz; ++i) tmp[1][i] = sigma_x[i];
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
			for (int i = 0; i < sz; ++i) {
				for (int j = 2; j <= nP; ++j){
					mpz_set_ui(tmp_ji,tmp[j][i]);
					mpz_add(sigma_x[i],sigma_x[i],tmp_ji);
				}
				mpz_mod(sigma_x[i], sigma_x[i], ppp);
				if (mpz_sgn(sigma_x[i]) != 0) {
					std::cout << "check error" << endl;
					sum_check = sum_check + 1;
				}
			}
			if (sum_check == 0) std::cout << "Correct!" << endl;
			mpz_clear(ppp);
			mpz_clear(tmp_ji);
			for(int i=0;i<sz;i++){
				mpz_clear(sigma_x[i]);
			}
			delete[] sigma_x;
			for (int i = 1; i <= nP; ++i) delete[] tmp[i];
			delete[] x_int;
		}
    }
};
#endif// SPDZ_H__