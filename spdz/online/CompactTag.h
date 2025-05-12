#ifndef CompactTag_SPDZ_H__
#define CompactTag_SPDZ_H__

#include "../network/network.h"
#include <emp-tool/emp-tool.h>
#include "utility.h"
#include "../network/helper.h"
using namespace emp;

template<int nP>
class CompactTag_SPDZ {
public:
	//const static int SSP = 5;
	//const uint64_t MASK = makeuint64_t(0x0ULL, 0xFFFFFULL);
	//Offline<nP>* fpre = nullptr; //not write
	//uint64_t* mac[nP+1];
	uint64_t key;
	//uint64_t* value[nP+1];

	//triple
	uint64_t* a, * b, * c, * r, * d, * e, * f;
	uint64_t* mac_a, * mac_b, * mac_c, * mac_r;

	//uint64_t * labels;
	//BristolFormat * cf;
	NetIOMP<nP>* io;
	//int num_ands = 0, num_in;
	int party, total_pre, ssp;
	int mat_sz = 2;
	int sz = mat_sz * mat_sz;
	ThreadPool* pool;
	//uint64_t Delta;

	PRP prp;
	CompactTag_SPDZ(NetIOMP<nP>* io[2], ThreadPool* pool, int party, bool* _delta = nullptr, int ssp = 40) {
		this->party = party;
		this->io = io[0];
		this->ssp = ssp;
		this->pool = pool;

		//num_in = cf->n1+cf->n2;
		//total_pre = num_in + num_ands + 3*ssp;
		//fpre = new FpreMP<nP>(io, pool, party, _delta, ssp);

		a = new uint64_t[sz];
		b = new uint64_t[sz];
		c = new uint64_t[sz];
		r = new uint64_t[sz];
		d = new uint64_t[sz];
		e = new uint64_t[sz];
		f = new uint64_t[sz];
		mac_a = new uint64_t[sz];
		mac_b = new uint64_t[sz];
		mac_c = new uint64_t[sz];
		mac_r = new uint64_t[sz];

		//key is a scalar
		uint64_t key;

		uint64_t test[3] = { 7, 13, 17 };

		if (party == 1) key = test[0];
		if (party == 2) key = test[1];
		if (party == 3) key = test[2];

	}
	~CompactTag_SPDZ() {
		delete[] a;
		delete[] b;
		delete[] c;
		delete[] r;
		delete[] mac_a;
		delete[] mac_b;
		delete[] mac_c;
		delete[] mac_r;
	}

	PRG prg;

	// it should be implemented by offline.
	void get_triple() {
		std::ifstream fin1, fin2, fin3, fin4;
		std::ifstream fin_1, fin_2, fin_3, fin_4;
		fin1.open("/home/jackie/spdz/pre_data/predata_CompactTag/a.txt", std::ios::in);
		fin2.open("/home/jackie/spdz/pre_data/predata_CompactTag/b.txt", std::ios::in);
		fin3.open("/home/jackie/spdz/pre_data/predata_CompactTag/c.txt", std::ios::in);
		fin4.open("/home/jackie/spdz/pre_data/predata_CompactTag/c.txt", std::ios::in);

		fin_1.open("/home/jackie/spdz/pre_data/predata_CompactTag/mac_a.txt", std::ios::in);
		fin_2.open("/home/jackie/spdz/pre_data/predata_CompactTag/mac_b.txt", std::ios::in);
		fin_3.open("/home/jackie/spdz/pre_data/predata_CompactTag/mac_c.txt", std::ios::in);
		fin_4.open("/home/jackie/spdz/pre_data/predata_CompactTag/mac_r.txt", std::ios::in);

		if (!(fin1.is_open() && fin2.is_open() && fin3.is_open() && fin4.is_open() && fin_1.is_open() && fin_2.is_open() && fin_3.is_open() && fin_4.is_open()))
		{
			std::cerr << "cannot open the file";
		}

		char line[1024] = { 0 };
		std::vector<triple> Triple;
		//从文件中提取“行”
		while (fin1.getline(line, sizeof(line)))
		{
			//定义局部变量
			triple t;
			//从“行”中提取“单词”
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
			//定义局部变量
			triple t;
			//从“行”中提取“单词”
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
			//定义局部变量
			triple t;
			//从“行”中提取“单词”
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
			//定义局部变量
			triple t;
			//从“行”中提取“单词”
			std::stringstream word(line);
			word >> t.party;
			uint64_t num;
			while (word >> num)
				t.value.push_back(num);
			Triple.push_back(t);
		}

		for (int i = 0; i < nP; i++) {
			if (party == stoi(Triple[i].party)) {
				for (int j = 0; j < sz; j++) r[j] = Triple[i].value[j];
			}
		}

		Triple.clear();
		memset(line, 0, sizeof(line));

		while (fin_1.getline(line, sizeof(line)))
		{
			//定义局部变量
			triple t;
			//从“行”中提取“单词”
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
			//定义局部变量
			triple t;
			//从“行”中提取“单词”
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
			//定义局部变量
			triple t;
			//从“行”中提取“单词”
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

		while (fin_4.getline(line, sizeof(line)))
		{
			//定义局部变量
			triple t;
			//从“行”中提取“单词”
			std::stringstream word(line);
			word >> t.party;
			uint64_t num;
			while (word >> num)
				t.value.push_back(num);
			Triple.push_back(t);
		}
		for (int i = 0; i < nP; i++) {
			if (party == stoi(Triple[i].party)) {
				for (int j = 0; j < sz; j++) mac_r[j] = Triple[i].value[j];
			}
		}

		Triple.clear();
		memset(line, 0, sizeof(line));
	}

	void Online_mul(uint64_t* x, uint64_t* y, uint64_t* mac_x, uint64_t* mac_y, uint64_t* output, uint64_t* output_mac) {
		uint64_t* d = new uint64_t[sz];
		uint64_t* e = new uint64_t[sz];

		for (int i = 0; i < sz; ++i) {
			d[i] = x[i] - a[i];
			e[i] = y[i] - b[i];
		}

		//open d
		if (party != 1) {
			io->send_data(1, d, sz);
			io->flush(1);
			io->recv_data(1, d, sz);
		}
		else {
			uint64_t* tmp[nP + 1];
			for (int i = 1; i <= nP; ++i) tmp[i] = new uint64_t[sz];
			vector<future<void>> res;
			for (int i = 2; i <= nP; ++i) {
				int party2 = i;
				res.push_back(pool->enqueue([this, tmp, party2]() {
					io->recv_data(party2, tmp[party2], sz);
					}));
			}
			joinNclean(res);
			for (int i = 0; i < sz; ++i)
				for (int j = 2; j <= nP; ++j)
					tmp[1][i] = add_mod(tmp[1][i], tmp[j][i]);
			d = tmp[1];
			for (int i = 2; i <= nP; ++i) {
				int party2 = i;
				res.push_back(pool->enqueue([this, d, party2]() {
					io->send_data(party2, d, sz);
					io->flush(party2);
					}));
			}
			joinNclean(res);
			for (int i = 1; i <= nP; ++i) delete[] tmp[i];
		}

		//open e
		if (party != 1) {
			io->send_data(1, e, sz);
			io->flush(1);
			io->recv_data(1, e, sz);
		}
		else {
			uint64_t* tmp[nP + 1];
			for (int i = 1; i <= nP; ++i) tmp[i] = new uint64_t[sz];
			vector<future<void>> res;
			for (int i = 2; i <= nP; ++i) {
				int party2 = i;
				res.push_back(pool->enqueue([this, tmp, party2]() {
					io->recv_data(party2, tmp[party2], sz);
					}));
			}
			joinNclean(res);
			for (int i = 0; i < sz; ++i)
				for (int j = 2; j <= nP; ++j)
					tmp[1][i] = add_mod(tmp[1][i], tmp[j][i]);
			e = tmp[1];
			for (int i = 2; i <= nP; ++i) {
				int party2 = i;
				res.push_back(pool->enqueue([this, e, party2]() {
					io->send_data(party2, e, sz);
					io->flush(party2);
					}));
			}
			joinNclean(res);
			for (int i = 1; i <= nP; ++i) delete[] tmp[i];
		}

		uint64_t* f = new uint64_t[sz];
		uint64_t* db = new uint64_t[sz];
		uint64_t* ae = new uint64_t[sz];
		uint64_t* de = new uint64_t[sz];

		db = mat_mult_mod(d, b, mat_sz);
		ae = mat_mult_mod(a, e, mat_sz);
		de = mat_mult_mod(d, e, mat_sz);

		for (int i = 0; i < sz; ++i) {
			output[i] = mod(c[i] + db[i] + ae[i] + de[i]);
			f[i] = mod(output[i]-r[i]);
		}

		//open f
		if (party != 1) {
			io->send_data(1, f, sz);
			io->flush(1);
			io->recv_data(1, f, sz);
		}
		else {
			uint64_t* tmp[nP + 1];
			for (int i = 1; i <= nP; ++i) tmp[i] = new uint64_t[sz];
			vector<future<void>> res;
			for (int i = 2; i <= nP; ++i) {
				int party2 = i;
				res.push_back(pool->enqueue([this, tmp, party2]() {
					io->recv_data(party2, tmp[party2], sz);
					}));
			}
			joinNclean(res);
			for (int i = 0; i < sz; ++i)
				for (int j = 2; j <= nP; ++j)
					tmp[1][i] = add_mod(tmp[1][i], tmp[j][i]);
			f = tmp[1];
			for (int i = 2; i <= nP; ++i) {
				int party2 = i;
				res.push_back(pool->enqueue([this, f, party2]() {
					io->send_data(party2, f, sz);
					io->flush(party2);
					}));
			}
			joinNclean(res);
			for (int i = 1; i <= nP; ++i) delete[] tmp[i];
		}

		//generate mac_z
		uint64_t* mac_f = new uint64_t[sz];
		uint64_t* key1 = new uint64_t[sz];

		for (int i = 0; i < mat_sz; ++i) {
			for (int j = 0; j < mat_sz; ++j) {
				if (i != j) {
					key1[i * mat_sz + j] = 0;
				}
				else {
					key1[i * mat_sz + j] = key;
				}
			}
		}
		mac_f = mat_mult_mod(f, key1, mat_sz);

		for (int i = 0; i < sz; ++i) {
			output_mac[i] = mod(mac_r[i] + mac_f[i]);
		}
	}

	void mac_check(uint64_t* z, uint64_t* mac_z) {
		//generate random vector v
		uint64_t* v = new uint64_t[mat_sz];
		std::ifstream fin5;
		fin5.open("/home/jackie/spdz/pre_data/predata_CompactTag/v.txt", std::ios::in);
		if (!(fin5.is_open()))
		{
			std::cerr << "cannot open the file";
		}

		char line[1024] = { 0 };
		std::vector<triple> Triple;
		//从文件中提取“行”
		while (fin5.getline(line, sizeof(line)))
		{
			//定义局部变量
			triple t;
			//从“行”中提取“单词”
			std::stringstream word(line);
			word >> t.party;
			uint64_t num;
			while (word >> num)
				t.value.push_back(num);
			Triple.push_back(t);
		}

		for (int i = 0; i < nP; i++) {
			if (party == stoi(Triple[i].party)) {
				for (int j = 0; j < mat_sz; j++) v[j] = Triple[i].value[j];
			}
		}

		//open vector v
		if (party != 1) {
			io->send_data(1, v, mat_sz);
			io->flush(1);
			io->recv_data(1, v, mat_sz);
		}
		else {
			uint64_t* tmp[nP + 1];
			for (int i = 1; i <= nP; ++i) tmp[i] = new uint64_t[mat_sz];
			vector<future<void>> res;
			for (int i = 2; i <= nP; ++i) {
				int party2 = i;
				res.push_back(pool->enqueue([this, tmp, party2]() {
					io->recv_data(party2, tmp[party2], mat_sz);
					}));
			}
			joinNclean(res);
			for (int i = 0; i < mat_sz; ++i)
				for (int j = 2; j <= nP; ++j)
					tmp[1][i] = add_mod(tmp[1][i], tmp[j][i]);
			v = tmp[1];
			for (int i = 2; i <= nP; ++i) {
				int party2 = i;
				res.push_back(pool->enqueue([this, v, party2]() {
					io->send_data(party2, v, mat_sz);
					io->flush(party2);
				}));
			}
			joinNclean(res);
			for (int i = 1; i <= nP; ++i) delete[] tmp[i];
		}

		Triple.clear();
		memset(line, 0, sizeof(line));

		//locally compute 
		uint64_t* fv = new uint64_t[mat_sz];
		uint64_t* de = new uint64_t[sz];	
		uint64_t* mac_cv = new uint64_t[mat_sz];
		uint64_t* mac_aev = new uint64_t[mat_sz];
		uint64_t* mac_dbv = new uint64_t[mat_sz];
		uint64_t* mac_dev = new uint64_t[mat_sz];
		uint64_t* mac_rv = new uint64_t[mat_sz];
		uint64_t* mac_zv = new uint64_t[mat_sz];
		uint64_t* mac_fv = new uint64_t[mat_sz];
		uint64_t* key1 = new uint64_t[sz];

		for (int i = 0; i < mat_sz; ++i) {
			for (int j = 0; j < mat_sz; ++j) {
				if (i != j) {
					key1[i * mat_sz + j] = 0;
				}
				else {
					key1[i * mat_sz + j] = key;
				}
			}
		}
		fv = mat_vec_mult(f, v, mat_sz);
		de= mat_mult_mod(d, e, mat_sz);
		mac_cv = mat_vec_mult(mac_c, v, mat_sz);
		mac_rv = mat_vec_mult(mac_r, v, mat_sz);
		mac_aev = mat_vec_mult(mat_mult_mod(mac_a, e, mat_sz), v, mat_sz);
		mac_dbv = mat_vec_mult(mat_mult_mod(d, mac_b, mat_sz), v, mat_sz);
		mac_dev = mat_vec_mult(mat_mult_mod(de, key1, mat_sz), v, mat_sz);	

		for(int i = 0; i < mat_sz; ++i){
            mac_zv[i] = mod(mac_cv[i]+mac_aev[i]+mac_dbv[i]+mac_dev[i]);
			mac_fv[i] = mod(mac_cv[i]+mac_aev[i]+mac_dbv[i]+mac_dev[i]-mac_rv[i]);
        }

		//inner product
		uint64_t* g = new uint64_t[mat_sz];
		std::ifstream fin6;
		fin6.open("/home/jackie/spdz/pre_data/predata_CompactTag/g.txt", std::ios::in);
		if (!(fin6.is_open()))
		{
			std::cerr << "cannot open the file";
		}
		
		
		//从文件中提取“行”
		while (fin6.getline(line, sizeof(line)))
		{
			//定义局部变量
			triple t;
			//从“行”中提取“单词”
			std::stringstream word(line);
			word >> t.party;
			uint64_t num;
			while (word >> num)
				t.value.push_back(num);
			Triple.push_back(t);
		}

		for (int i = 0; i < nP; i++) {
			if (party == stoi(Triple[i].party)) {
				for (int j = 0; j < mat_sz; j++) g[j] = Triple[i].value[j];
			}
		}
		
		if (party != 1) {
			io->send_data(1, g, mat_sz);
			io->flush(1);
			io->recv_data(1, g, mat_sz);
		}
		else {
			uint64_t* tmp[nP + 1];
			for (int i = 1; i <= nP; ++i) tmp[i] = new uint64_t[mat_sz];
			vector<future<void>> res;
			for (int i = 2; i <= nP; ++i) {
				int party2 = i;
				res.push_back(pool->enqueue([this, tmp, party2]() {
					io->recv_data(party2, tmp[party2], mat_sz);
					}));
			}
			joinNclean(res);
			for (int i = 0; i < mat_sz; ++i)
				for (int j = 2; j <= nP; ++j)
					tmp[1][i] = add_mod(tmp[1][i], tmp[j][i]);
			g = tmp[1];
			for (int i = 2; i <= nP; ++i) {
				int party2 = i;
				res.push_back(pool->enqueue([this, g, party2]() {
					io->send_data(party2, g, mat_sz);
					io->flush(party2);
					}));
			}
			joinNclean(res);
			for (int i = 1; i <= nP; ++i) delete[] tmp[i];
		}

		Triple.clear();
		memset(line, 0, sizeof(line));

		uint64_t inn_prdt_fv;
		uint64_t mac_inn_prdt_fv;
		
		inn_prdt_fv = vector_inn_prdt_sum_red(fv, g, mat_sz);
		mac_inn_prdt_fv = vector_inn_prdt_sum_red(mac_fv, g, mat_sz);

		//check
		uint64_t* sigma_fv = new uint64_t[1];
		sigma_fv[0] = mac_inn_prdt_fv - mult_mod(inn_prdt_fv, key);

		if (party != 1) {
			io->send_data(1, sigma_fv, 1);
			io->flush(1);
		}
		else {
			uint64_t* tmp[nP + 1];
			for (int i = 1; i <= nP; ++i) tmp[i] = new uint64_t[1];
			vector<future<void>> res;
			for (int i = 2; i <= nP; ++i) {
				int party2 = i;
				res.push_back(pool->enqueue([this, tmp, party2]() {
					io->recv_data(party2, tmp[party2], 1);
					}));
			}
			joinNclean(res);
			for (int j = 2; j <= nP; ++j)
					tmp[1][0] = add_mod(tmp[1][0], tmp[j][0]);
			if (tmp[1][0] != 0) cout << "check error" << endl;
			for (int i = 1; i <= nP; ++i) delete[] tmp[i];
		}
	}
};
#endif// CompactTag_SPDZ_H__