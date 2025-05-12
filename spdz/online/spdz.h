#ifndef SPDZ_H__
#define SPDZ_H__

#include "../network/network.h"
#include "utility.h"
#include <emp-tool/emp-tool.h>
#include "../network/helper.h"
using namespace emp;

// base SPDZ protocol with scalar elements
template<int nP>
class SPDZ { public:
	const block MASK = makeBlock(0x0ULL, 0xFFFFFULL);
	const uint64_t PR = 2305843009213693951;
	uint64_t key;
    uint64_t *a,*b,*c;
	uint64_t *mac_a,*mac_b,*mac_c;

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

		uint64_t key;
        a = new uint64_t[sz];
        b = new uint64_t[sz];
        c = new uint64_t[sz];
		mac_a = new uint64_t[sz];
        mac_b = new uint64_t[sz];
        mac_c = new uint64_t[sz];
	
	}
	~SPDZ () {
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
				for (int j = 0; j < sz; j++) a[j] = Triple[i].value[j];
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
				for (int j = 0; j < sz; j++) b[j] = Triple[i].value[j];
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
				for (int j = 0; j < sz; j++) c[j] = Triple[i].value[j];
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
				key = Triple[i].value[0];
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
				for (int j = 0; j < sz; j++) mac_a[j] = Triple[i].value[j];
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
				for (int j = 0; j < sz; j++) mac_b[j] = Triple[i].value[j];
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
				for (int j = 0; j < sz; j++) mac_c[j] = Triple[i].value[j];
			}
		}

		Triple.clear();
		memset(line, 0, sizeof(line));

	}

	void Online_mul (uint64_t* x, uint64_t* y, uint64_t* mac_x, uint64_t* mac_y, uint64_t* output, uint64_t* output_mac) {
		uint64_t* d = new uint64_t[sz];
		uint64_t* e = new uint64_t[sz];

		for(int i=0;i<sz;i++){
			d[i] = mod_sub1(x[i] , a[i]);
        	e[i] = mod_sub1(y[i] , b[i]);
		}
        
		uint64_t* D;
		D=new uint64_t[sz];
		if(party != 1) {
			io->send_data(1, d, sz*sizeof(uint64_t));
			io->flush(1);
			io->recv_data(1, d, sz*sizeof(uint64_t));
		} else {
			uint64_t * tmp[nP+1];
			for (int i = 1; i <= nP; ++i) tmp[i] = new uint64_t[sz];
			for (int i = 0; i < sz; ++i) tmp[1][i] = d[i];
			vector<future<void>> res;
			for(int i = 2; i <= nP; ++i) {
				int party2 = i;
				res.push_back(pool->enqueue([this, tmp, party2]() {
					io->recv_data(party2, tmp[party2], sz*sizeof(uint64_t));
				}));
			}
			joinNclean(res);
			for (int i = 0; i < sz; ++i) {
				for (int j = 2; j <= nP; ++j) {
					tmp[1][i] = mod_add(tmp[1][i], tmp[j][i]);
				}
			}
            d = tmp[1];
			for (int i = 0; i < sz; i++) {
				D[i] = d[i];
			}
			for(int i = 2; i <= nP; ++i) {
				int party2 = i;
				res.push_back(pool->enqueue([this, d, party2]() {
					io->send_data(party2, d, sz*sizeof(uint64_t));
					io->flush(party2);
				}));
			}
			joinNclean(res);
			for (int i = 1; i <= nP; ++i) delete[] tmp[i];
		}

		uint64_t* E;
		E = new uint64_t[sz];
        if(party != 1) {
			io->send_data(1, e, sz*sizeof(uint64_t));
			io->flush(1);
			io->recv_data(1, e, sz*sizeof(uint64_t));
		} else {
			uint64_t * tmp[nP+1];
			for (int i = 1; i <= nP; ++i) tmp[i] = new uint64_t[sz];
			for (int i = 0; i < sz; ++i) tmp[1][i] = e[i];
			vector<future<void>> res1;
			for(int i = 2; i <= nP; ++i) {
				int party2 = i;
				res1.push_back(pool->enqueue([this, tmp, party2]() {
					io->recv_data(party2, tmp[party2], sz*sizeof(uint64_t));
				}));
			}
			joinNclean(res1);
			for (int i = 0; i < sz; ++i) {
				for (int j = 2; j <= nP; ++j)
					tmp[1][i] = mod_add(tmp[1][i], tmp[j][i]);
			}
			e = tmp[1];
			for (int i = 0; i < sz; i++) {
				E[i] = e[i];
			}
			for(int i = 2; i <= nP; ++i) {
				int party2 = i;
				res1.push_back(pool->enqueue([this, e, party2]() {
					io->send_data(party2, e, sz*sizeof(uint64_t));
					io->flush(party2);
				}));
			}
			joinNclean(res1);
			for (int i = 1; i <= nP; ++i) delete[] tmp[i];
		}

		uint64_t db;
		uint64_t ae;
		uint64_t de;
		uint64_t mac_db;
		uint64_t mac_ae;
		uint64_t mac_de;

		for(int i=0;i<sz;i++){
			if(party==1){
				de=mod_mul(D[0],E[0]);
				db=mod_mul(D[0],b[0]);
				ae=mod_mul(a[0],E[0]);
				mac_de=mod_mul(de,key);
				mac_db=mod_mul(D[0],mac_b[0]);
				mac_ae=mod_mul(mac_a[0],E[0]);
				output[0]=mod_add(mod_add(mod_add(de,db),ae),c[0]);
				output_mac[0]=mod_add(mod_add(mod_add(mac_de,mac_db),mac_ae),mac_c[0]);
			}else{
				de=mod_mul(d[0],e[0]);
				db=mod_mul(d[0],b[0]);
				ae=mod_mul(a[0],e[0]);
				mac_de=mod_mul(de,key);
				mac_db=mod_mul(d[0],mac_b[0]);
				mac_ae=mod_mul(mac_a[0],e[0]);
				output[0]=mod_add(mod_add(db,ae),c[0]);
				output_mac[0]=mod_add(mod_add(mod_add(mac_de,mac_db),mac_ae),mac_c[0]);
			}
		}
				
	}	
	
    void mac_check(uint64_t* x, uint64_t* mac_x){
		uint64_t* X;
		X = new uint64_t[sz];
        if(party != 1) {
			io->send_data(1, x, sz*sizeof(uint64_t));
			io->flush(1);
			io->recv_data(1, x, sz*sizeof(uint64_t));
		} else {
			uint64_t * tmp[nP+1];
			for (int i = 1; i <= nP; ++i) tmp[i] = new uint64_t[sz];
			for (int i = 0; i < sz; i++) tmp[1][i] = x[i];
			vector<future<void>> res;
			for(int i = 2; i <= nP; ++i) {
				int party2 = i;
				res.push_back(pool->enqueue([this, tmp, party2]() {
					io->recv_data(party2, tmp[party2], sz*sizeof(uint64_t));
				}));
			}
			joinNclean(res);
			for (int i = 0; i < sz; ++i)
				for (int j = 2; j <= nP; ++j)
					tmp[1][i] = mod_add(tmp[1][i], tmp[j][i]);
			x = tmp[1];
			for (int i = 0; i < sz; i++) {
				X[i] = x[i];
			}
			for(int i = 2; i <= nP; ++i) {
				int party2 = i;
				res.push_back(pool->enqueue([this, x, party2]() {
					io->send_data(party2, x, sz*sizeof(uint64_t));
					io->flush(party2);
				}));
			}
			joinNclean(res);
			for (int i = 1; i <= nP; ++i) delete[] tmp[i];
		}

        uint64_t* sigma_x = new uint64_t[sz];
		if(party==1){
			for (int i = 0; i < sz; i++) sigma_x[i]=mod_sub1(mac_x[i],mod_mul(X[i],key));
		}else{
			for (int i = 0; i < sz; i++) sigma_x[i]=mod_sub1(mac_x[i],mod_mul(x[i],key));
		}

        if(party != 1) {
			io->send_data(1, sigma_x, sz*sizeof(uint64_t));
			io->flush(1);
		} else {
			uint64_t * tmp[nP+1];
			for (int i = 1; i <= nP; ++i) tmp[i] = new uint64_t[sz];
			for (int i = 0; i < sz; ++i) tmp[1][i] = sigma_x[i];
			vector<future<void>> res;
			for(int i = 2; i <= nP; ++i) {
				int party2 = i;
				res.push_back(pool->enqueue([this, tmp, party2]() {
					io->recv_data(party2, tmp[party2], sz*sizeof(uint64_t));
				}));
			}
			joinNclean(res);
			uint64_t sum_check = 0;
			for (int i = 0; i < sz; ++i) {
				for (int j = 2; j <= nP; ++j)
					tmp[1][i] = mod_add(tmp[1][i], tmp[j][i]);
				if (tmp[1][i] != 0) {
					std::cout << "check error" << endl;
					sum_check = sum_check + 1;
				}
			}
			if (sum_check == 0) std::cout << "Correct!" << endl;
			for (int i = 1; i <= nP; ++i) delete[] tmp[i];
		}
    }

	uint64_t mod_add(uint64_t a, uint64_t b) {
		if (a >= PR - b) {
			return a - (PR - b);
		}
		else {
			return a + b;
		}
	}

	uint64_t mod_sub1(uint64_t a, uint64_t b) {
		if (a < b) {
			return PR - (b - a);
		}
		else {
			return a - b;
		}
	}

	uint64_t mod_mul(uint64_t a, uint64_t b) {
		uint64_t result = 0;
		while (b > 0) {
			if (b & 1) {
				result = mod_add(result, a);
			}
			a = mod_add(a, a);
			b >>= 1;
		}
		return result;
	}
};
#endif// SPDZ_H__