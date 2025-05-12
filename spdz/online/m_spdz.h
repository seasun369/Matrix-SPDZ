#ifndef MSPDZ_H__
#define MSPDZ_H__
#define NUM_THREADS 12

#include "../network/network.h"
#include <emp-tool/emp-tool.h>
#include "utility.h"
#include "../network/helper.h"
#include <omp.h>
using namespace emp;

// Matrix SPDZ
template<int nP>
class MSPDZ { 
public:
	uint64_t* key;

    uint64_t *a,*a_t,*b,*c,*r,*r_t;
	uint64_t *mac_a,*mac_a_t,*mac_b,*mac_c,*mac_r,*mac_r_t;

	NetIOMP<nP> * io;
	int party, total_pre, ssp;
	const uint64_t PR = 2305843009213693951ULL;
	int mat_sz = 128;
	int sz = mat_sz * mat_sz;
	ThreadPool * pool;
	
		
	PRP prp;
	MSPDZ(NetIOMP<nP> * io[2], ThreadPool * pool, int party, bool * _delta = nullptr, int ssp = 40) {
		this->party = party;
		this->io = io[0];
		this->ssp = ssp;
		this->pool = pool;

        a = new uint64_t[sz];
		a_t = new uint64_t[sz];
        b = new uint64_t[sz];
        c = new uint64_t[sz];
        r = new uint64_t[sz];
        r_t = new uint64_t[sz];
		key = new uint64_t[mat_sz];

		mac_a = new uint64_t[mat_sz];
		mac_a_t = new uint64_t[mat_sz];
        mac_b = new uint64_t[mat_sz];
        mac_c = new uint64_t[mat_sz];
        mac_r = new uint64_t[mat_sz];
        mac_r_t = new uint64_t[mat_sz]; 

	}
	~MSPDZ() {
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
    void get_triple(){
		std::ifstream fin1,fin2,fin3,fin4,fin5,fin6,fin7;
		std::ifstream fin_1,fin_2,fin_3,fin_4,fin_5,fin_6;
		fin1.open("/home/jackie/spdz/pre_data/predata_mspdz/128/a.txt",std::ios::in);
		fin2.open("/home/jackie/spdz/pre_data/predata_mspdz/128/a_t.txt",std::ios::in);
		fin3.open("/home/jackie/spdz/pre_data/predata_mspdz/128/b.txt",std::ios::in);
		fin4.open("/home/jackie/spdz/pre_data/predata_mspdz/128/c.txt",std::ios::in);
		fin5.open("/home/jackie/spdz/pre_data/predata_mspdz/128/r.txt",std::ios::in);
		fin6.open("/home/jackie/spdz/pre_data/predata_mspdz/128/r_t.txt",std::ios::in);
		fin7.open("/home/jackie/spdz/pre_data/predata_mspdz/128/key.txt",std::ios::in);
		fin_1.open("/home/jackie/spdz/pre_data/predata_mspdz/128/mac_a.txt",std::ios::in);
		fin_2.open("/home/jackie/spdz/pre_data/predata_mspdz/128/mac_a_t.txt",std::ios::in);
		fin_3.open("/home/jackie/spdz/pre_data/predata_mspdz/128/mac_b.txt",std::ios::in);
		fin_4.open("/home/jackie/spdz/pre_data/predata_mspdz/128/mac_c.txt",std::ios::in);
		fin_5.open("/home/jackie/spdz/pre_data/predata_mspdz/128/mac_r.txt",std::ios::in);
		fin_6.open("/home/jackie/spdz/pre_data/predata_mspdz/128/mac_r_t.txt",std::ios::in);

		if(!(fin1.is_open() && fin2.is_open() && fin3.is_open()&& fin4.is_open() && fin5.is_open()&& fin6.is_open() && fin7.is_open()&& fin_1.is_open() && fin_2.is_open()&& fin_3.is_open() && fin_4.is_open()&& fin_5.is_open() && fin_6.is_open()))
		{
		    std::cerr<<"cannot open the file";
		}

		char line1[512000]={0};
		std::vector<triple> Triple;
		while(fin1.getline(line1,sizeof(line1)))
		{
		    triple t;
		    std::stringstream word(line1);
		    word>>t.party;
		    uint64_t num;
		    while(word>>num)
		        t.value.push_back(num);
		    Triple.push_back(t);
		}

		for(int i=0; i< nP; i++){
			if(party == stoi(Triple[i].party)){
				for(int j=0; j< sz; j++) a[j] = Triple[i].value[j];
			}
		}

		Triple.clear();
		memset(line1, 0, sizeof(line1));

		while(fin2.getline(line1,sizeof(line1)))
		{
		    triple t;
		    std::stringstream word(line1);
		    word>>t.party;
		    uint64_t num;
		    while(word>>num)
		        t.value.push_back(num);
		    Triple.push_back(t);
		}

		for(int i=0; i< nP; i++){
			if(party == stoi(Triple[i].party)){
				for(int j=0; j< sz; j++) a_t[j] = Triple[i].value[j];
			}
		}

		Triple.clear();
		memset(line1, 0, sizeof(line1));

		while(fin3.getline(line1,sizeof(line1)))
		{
		    triple t;
		    std::stringstream word(line1);
		    word>>t.party;
		    uint64_t num;
		    while(word>>num)
		        t.value.push_back(num);
		    Triple.push_back(t);
		}

		for(int i=0; i< nP; i++){
			if(party == stoi(Triple[i].party)){
				for(int j=0; j< sz; j++) b[j] = Triple[i].value[j];
			}
		}

		Triple.clear();
		memset(line1, 0, sizeof(line1));

		while(fin4.getline(line1,sizeof(line1)))
		{
		    triple t;
		    std::stringstream word(line1);
		    word>>t.party;
		    uint64_t num;
		    while(word>>num)
		        t.value.push_back(num);
		    Triple.push_back(t);
		}

		for(int i=0; i< nP; i++){
			if(party == stoi(Triple[i].party)){
				for(int j=0; j< sz; j++) c[j] = Triple[i].value[j];
			}
		}

		Triple.clear();
		memset(line1, 0, sizeof(line1));

		while(fin5.getline(line1,sizeof(line1)))
		{
		    triple t;
		    std::stringstream word(line1);
		    word>>t.party;
		    uint64_t num;
		    while(word>>num)
		        t.value.push_back(num);
		    Triple.push_back(t);
		}

		for(int i=0; i< nP; i++){
			if(party == stoi(Triple[i].party)){
				for(int j=0; j< sz; j++) r[j] = Triple[i].value[j];
			}
		}

		Triple.clear();
		memset(line1, 0, sizeof(line1));

		while(fin6.getline(line1,sizeof(line1)))
		{
		    triple t;
		    std::stringstream word(line1);
		    word>>t.party;
		    uint64_t num;
		    while(word>>num)
		        t.value.push_back(num);
		    Triple.push_back(t);
		}

		for(int i=0; i< nP; i++){
			if(party == stoi(Triple[i].party)){
				for(int j=0; j< sz; j++) r_t[j] = Triple[i].value[j];
			}
		}

		Triple.clear();
		memset(line1, 0, sizeof(line1));

		while(fin7.getline(line1,sizeof(line1)))
		{
		    triple t;
		    std::stringstream word(line1);
		    word>>t.party;
		    uint64_t num;
		    while(word>>num)
		        t.value.push_back(num);
		    Triple.push_back(t);
		}

		for(int i=0; i< nP; i++){
			if(party == stoi(Triple[i].party)){
				for(int j=0; j< mat_sz; j++) key[j] = Triple[i].value[j];
			}
		}

		Triple.clear();
		memset(line1, 0, sizeof(line1));

		while(fin_1.getline(line1,sizeof(line1)))
		{
		    triple t;
		    std::stringstream word(line1);
		    word>>t.party;
		    uint64_t num;
		    while(word>>num)
		        t.value.push_back(num);
		    Triple.push_back(t);
		}

		for(int i=0; i< nP; i++){
			if(party == stoi(Triple[i].party)){
				for(int j=0; j< mat_sz; j++) mac_a[j] = Triple[i].value[j];
			}
		}

		Triple.clear();
		memset(line1, 0, sizeof(line1));

		while(fin_2.getline(line1,sizeof(line1)))
		{
		    triple t;
		    std::stringstream word(line1);
		    word>>t.party;
		    uint64_t num;
		    while(word>>num)
		        t.value.push_back(num);
		    Triple.push_back(t);
		}

		for(int i=0; i< nP; i++){
			if(party == stoi(Triple[i].party)){
				for(int j=0; j< mat_sz; j++) mac_a_t[j] = Triple[i].value[j];
			}
		}

		Triple.clear();
		memset(line1, 0, sizeof(line1));

		while(fin_3.getline(line1,sizeof(line1)))
		{
		    triple t;
		    std::stringstream word(line1);
		    word>>t.party;
		    uint64_t num;
		    while(word>>num)
		        t.value.push_back(num);
		    Triple.push_back(t);
		}

		for(int i=0; i< nP; i++){
			if(party == stoi(Triple[i].party)){
				for(int j=0; j< mat_sz; j++) mac_b[j] = Triple[i].value[j];
			}
		}

		Triple.clear();
		memset(line1, 0, sizeof(line1));

		while(fin_4.getline(line1,sizeof(line1)))
		{
		    triple t;
		    std::stringstream word(line1);
		    word>>t.party;
		    uint64_t num;
		    while(word>>num)
		        t.value.push_back(num);
		    Triple.push_back(t);
		}

		for(int i=0; i< nP; i++){
			if(party == stoi(Triple[i].party)){
				for(int j=0; j< mat_sz; j++) mac_c[j] = Triple[i].value[j];
			}
		}

		Triple.clear();
		memset(line1, 0, sizeof(line1));

		while(fin_5.getline(line1,sizeof(line1)))
		{
		    triple t;
		    std::stringstream word(line1);
		    word>>t.party;
		    uint64_t num;
		    while(word>>num)
		        t.value.push_back(num);
		    Triple.push_back(t);
		}

		for(int i=0; i< nP; i++){
			if(party == stoi(Triple[i].party)){
				for(int j=0; j< mat_sz; j++) mac_r[j] = Triple[i].value[j];
			}
		}

		Triple.clear();
		memset(line1, 0, sizeof(line1));

		while(fin_6.getline(line1,sizeof(line1)))
		{
		    triple t;
		    std::stringstream word(line1);
		    word>>t.party;
		    uint64_t num;
		    while(word>>num)
		        t.value.push_back(num);
		    Triple.push_back(t);
		}

		for(int i=0; i< nP; i++){
			if(party == stoi(Triple[i].party)){
				for(int j=0; j< mat_sz; j++) mac_r_t[j] = Triple[i].value[j];
			}
		}

		Triple.clear();
		memset(line1, 0, sizeof(line1));
	}

	void Online_mul(uint64_t * x, uint64_t * y, uint64_t * mac_x, uint64_t * mac_y, uint64_t * output, uint64_t *output_mac) {
		uint64_t *d = new uint64_t[sz];
        uint64_t *e = new uint64_t[sz];
        uint64_t *f = new uint64_t[sz];

		omp_set_num_threads(NUM_THREADS);
		#pragma omp parallel for
        for(int i=0; i<sz; ++i) {
            d[i]=mod_sub1(x[i],a[i]);
            e[i]=mod_sub1(y[i],b[i]);
        }

		uint64_t* D;
		D = new uint64_t[sz];
		if(party != 1) {
			io->send_data(1, d, sz * sizeof(uint64_t));
			io->flush(1);
			io->recv_data(1, d, sz * sizeof(uint64_t));
		} else {
			uint64_t * tmp[nP+1];
			for(int i = 1; i <= nP; ++i) tmp[i] = new uint64_t[sz];
			omp_set_num_threads(NUM_THREADS);
			#pragma omp parallel for
			for(int i = 0 ; i < sz; i++) tmp[1][i] = d[i];
			vector<future<void>> res;
			for(int i = 2; i <= nP; ++i) {
				int party2 = i;
				res.push_back(pool->enqueue([this, tmp, party2]() {
					io->recv_data(party2, tmp[party2], sz * sizeof(uint64_t));
				}));
			}
			joinNclean(res);
			omp_set_num_threads(NUM_THREADS);
			#pragma omp parallel for
			for(int i = 0; i < sz; ++i){
				for(int j = 2; j <= nP; ++j){
					tmp[1][i] = mod_add(tmp[1][i] , tmp[j][i]);
				}
			}
            d = tmp[1];
			omp_set_num_threads(NUM_THREADS);
			#pragma omp parallel for
			for(int i =0;i < sz; i++){
				D[i] = d[i];
			}
			for(int i = 2; i <= nP; ++i) {
				int party2 = i;
				res.push_back(pool->enqueue([this, d, party2]() {
					io->send_data(party2, d, sz * sizeof(uint64_t));
					io->flush(party2);
				}));
			}
			joinNclean(res);
			for(int i = 1; i <= nP; ++i) delete[] tmp[i];
		}		

		uint64_t* E;
		E = new uint64_t[sz];
        if(party != 1) {
			io->send_data(1, e, sz * sizeof(uint64_t));
			io->flush(1);
			io->recv_data(1, e, sz * sizeof(uint64_t));
		} else {
			uint64_t * tmp[nP+1];
			for(int i = 1; i <= nP; ++i) tmp[i] = new uint64_t[sz];
			omp_set_num_threads(NUM_THREADS);
			#pragma omp parallel for
			for(int i = 0; i < sz; i++) tmp[1][i] = e[i];
			vector<future<void>> res1;
			for(int i = 2; i <= nP; ++i) {
				int party2 = i;
				res1.push_back(pool->enqueue([this, tmp, party2]() {
					io->recv_data(party2, tmp[party2], sz * sizeof(uint64_t));
				}));
			}
			joinNclean(res1);
			omp_set_num_threads(NUM_THREADS);
			#pragma omp parallel for
			for(int i = 0; i < sz; ++i)
				for(int j = 2; j <= nP; ++j)
					tmp[1][i] = mod_add(tmp[1][i] , tmp[j][i]);
            e = tmp[1];
			omp_set_num_threads(NUM_THREADS);
			#pragma omp parallel for
			for(int i = 0; i < sz; i++){
				E[i] = e[i];
			}
			for(int i = 2; i <= nP; ++i) {
				int party2 = i;
				res1.push_back(pool->enqueue([this, e, party2]() {
					io->send_data(party2, e, sz * sizeof(uint64_t));
					io->flush(party2);
				}));
			}
			joinNclean(res1);
			for(int i = 1; i <= nP; ++i) delete[] tmp[i];
		}

		if(party==1){
			f = mult_mod(trans_mat(E, mat_sz), a_t, mat_sz);
		}else{
			f = mult_mod(trans_mat(e, mat_sz), a_t, mat_sz);
		}

		omp_set_num_threads(NUM_THREADS);
		#pragma omp parallel for
        for(int i=0; i < sz ; ++i){
            f[i]=mod_sub1(f[i],r_t[i]);
        }

		uint64_t* F;
		F = new uint64_t[sz];
        if(party != 1) {
			io->send_data(1, f, sz * sizeof(uint64_t));
			io->flush(1);
			io->recv_data(1, f, sz * sizeof(uint64_t));
		} else {
			uint64_t * tmp[nP+1];
			for(int i = 1; i <= nP; ++i) tmp[i] = new uint64_t[sz];
			omp_set_num_threads(NUM_THREADS);
			#pragma omp parallel for
			for(int i=0;i<sz;i++) tmp[1][i]=f[i];
			vector<future<void>> res;
			for(int i = 2; i <= nP; ++i) {
				int party2 = i;
				res.push_back(pool->enqueue([this, tmp, party2]() {
					io->recv_data(party2, tmp[party2], sz * sizeof(uint64_t));
				}));
			}
			joinNclean(res);
			omp_set_num_threads(NUM_THREADS);
			#pragma omp parallel for
			for(int i = 0; i < sz; ++i)
				#pragma omp parallel for
				for(int j = 2; j <= nP; ++j)
					tmp[1][i] = mod_add(tmp[1][i] , tmp[j][i]);
            f=tmp[1];
			omp_set_num_threads(NUM_THREADS);
			#pragma omp parallel for
			for(int i=0;i<sz;i++){
				F[i]=f[i];
			}
			for(int i = 2; i <= nP; ++i) {
				int party2 = i;
				res.push_back(pool->enqueue([this, f, party2]() {
					io->send_data(party2, f, sz * sizeof(uint64_t));
					io->flush(party2);
				}));
			}
			joinNclean(res);
			for(int i = 1; i <= nP; ++i) delete[] tmp[i];
		}

        uint64_t *f_t = new uint64_t[sz];
        uint64_t *db = new uint64_t[sz];
        uint64_t *de = new uint64_t[sz];
        uint64_t *mac_db = new uint64_t[mat_sz];
        uint64_t *mac_de = new uint64_t[mat_sz];
		uint64_t *mac_f_t = new uint64_t[mat_sz];
		
		if(party==1){
			f_t = trans_mat(F, mat_sz);
        	db = mult_mod(D, b, mat_sz);
        	de = mult_mod(D, E, mat_sz);
        	mac_db = mat_vec(D, mac_b, mat_sz);
        	mac_de = mat_vec(de, key, mat_sz);
			mac_f_t = mat_vec(f_t, key, mat_sz);
		}else{
			f_t = trans_mat(f, mat_sz);
        	db = mult_mod(d, b, mat_sz);
        	de = mult_mod(d, e, mat_sz);
        	mac_db = mat_vec(d, mac_b, mat_sz);
        	mac_de = mat_vec(de, key, mat_sz);
			mac_f_t = mat_vec(f_t, key, mat_sz);
		}

		omp_set_num_threads(NUM_THREADS);
		#pragma omp parallel for
        for(int i = 0; i < sz; ++i){
			if(party==1){
				output[i] = mod_add(mod_add(mod_add(c[i],db[i]),mod_add(r[i],de[i])),f_t[i]);
			}else{
				output[i] = mod_add(mod_add(c[i],db[i]),r[i]);
			}
        }
		omp_set_num_threads(NUM_THREADS);
		#pragma omp parallel for
		for(int i=0;i<mat_sz;i++){
			output_mac[i] = mod_add(mod_add(mod_add(mac_c[i],mac_db[i]),mod_add(mac_r[i],mac_de[i])),mac_f_t[i]);
		}
	}
	
    void mac_check(uint64_t * x, uint64_t * mac_x){
		uint64_t* X;
		X=new uint64_t[sz];

        if(party != 1) {
			io->send_data(1, x, sz * sizeof(uint64_t));
			io->flush(1);
			io->recv_data(1, x, sz * sizeof(uint64_t));
		} else {
			uint64_t * tmp[nP+1];
			for(int i = 1; i <= nP; ++i) tmp[i] = new uint64_t[sz];
			omp_set_num_threads(NUM_THREADS);
			#pragma omp parallel for
			for(int i=0;i<sz;i++) tmp[1][i]=x[i];
			vector<future<void>> res;
			for(int i = 2; i <= nP; ++i) {
				int party2 = i;
				res.push_back(pool->enqueue([this, tmp, party2]() {
					io->recv_data(party2, tmp[party2], sz * sizeof(uint64_t));
				}));
			}
			joinNclean(res);
			omp_set_num_threads(NUM_THREADS);
			#pragma omp parallel for
			for(int i = 0; i < sz; ++i)
				for(int j = 2; j <= nP; ++j)
					tmp[1][i] = mod_add(tmp[1][i] , tmp[j][i]);
            x=tmp[1];
			omp_set_num_threads(NUM_THREADS);
			#pragma omp parallel for
			for(int i=0;i<sz;i++){
				X[i]=x[i];
			}
			for(int i = 2; i <= nP; ++i) {
				int party2 = i;
				res.push_back(pool->enqueue([this, x, party2]() {
					io->send_data(party2, x, sz * sizeof(uint64_t));
					io->flush(party2);
				}));
			}
			joinNclean(res);
			for(int i = 1; i <= nP; ++i) delete[] tmp[i];
		}

        uint64_t * sigma_x = new uint64_t[mat_sz];
		if(party==1){
			sigma_x = mat_vec(X,key,mat_sz);
		}else{
			sigma_x = mat_vec(x,key,mat_sz);
		}
		omp_set_num_threads(NUM_THREADS);
		#pragma omp parallel for
        for(int i = 0; i < mat_sz; ++i) {
            sigma_x[i] = mod_sub1(mac_x[i] , sigma_x[i]);
        }

        if(party != 1) {
			io->send_data(1, sigma_x, mat_sz * sizeof(uint64_t));
			io->flush(1);
		} else {
			uint64_t * tmp[nP+1];
			for(int i = 1; i <= nP; ++i) tmp[i] = new uint64_t[mat_sz];
			omp_set_num_threads(NUM_THREADS);
			#pragma omp parallel for
			for(int i=0;i<mat_sz;i++) tmp[1][i]=sigma_x[i];
			vector<future<void>> res;
			for(int i = 2; i <= nP; ++i) {
				int party2 = i;
				res.push_back(pool->enqueue([this, tmp, party2]() {
					io->recv_data(party2, tmp[party2], mat_sz * sizeof(uint64_t));
				}));
			}
			joinNclean(res);
			uint64_t sum_check=0;
			omp_set_num_threads(NUM_THREADS);
			#pragma omp parallel for
			for(int i = 0; i < mat_sz; ++i){
				for(int j = 2; j <= nP; ++j)
					tmp[1][i] = mod_add(tmp[1][i] , tmp[j][i]);
                if(tmp[1][i] != 0){
					std::cout <<"check error"<<endl;
					std::cout<<i<<endl;
					sum_check=sum_check+1;
				} 
            }
			if(sum_check==0) std::cout<<"Correct!"<<endl;
			for(int i = 1; i <= nP; ++i) delete[] tmp[i];
		}
    }

	uint64_t* mult_mod(uint64_t* a,uint64_t* b, int mat_sz){
		uint64_t* res = new uint64_t[sz];
		omp_set_num_threads(NUM_THREADS);
		#pragma omp parallel for
		for(int i = 0;i<sz;i++) res[i]=0;
		omp_set_num_threads(NUM_THREADS);
		#pragma omp parallel for
		for(int i=0;i<mat_sz;i++){
			#pragma omp parallel for
			for(int j=0;j<mat_sz;j++){
				for(int k=0;k<mat_sz;k++){
					res[i*mat_sz+k]=mod_add(res[i*mat_sz+k],mod_mul(a[i*mat_sz+j],b[j*mat_sz+k]));
				}
			}
		}
		return res;
	}

	uint64_t* trans_mat(uint64_t* A, int mat_sz){
  		uint64_t *A_T = new uint64_t[sz];
		omp_set_num_threads(NUM_THREADS);
		#pragma omp parallel for
  		for(int i=0; i < mat_sz; ++i){
			A_T[i*mat_sz+i]=A[i*mat_sz+i];
    		for(int j =i+1; j < mat_sz; ++j){
      			A_T[i*mat_sz+j]=A[j*mat_sz+i];
      			A_T[j*mat_sz+i]=A[i*mat_sz+j];
    		}
  		}
  		return A_T;
	}

	uint64_t* mat_vec(uint64_t* a, uint64_t* b, int mat_sz){
  		uint64_t * res = new uint64_t[mat_sz];
		omp_set_num_threads(NUM_THREADS);
		#pragma omp parallel for
		for(int i = 0 ; i < mat_sz;i++) res[i]=0;
		omp_set_num_threads(NUM_THREADS);
		#pragma omp parallel for
  		for (int i = 0; i < mat_sz; ++i) {
            for (int j = 0; j < mat_sz; ++j) {  
              res[i] = mod_add(res[i] , mod_mul(a[i * mat_sz + j] , b[j]));
            }
        } 
  		return res;
	}

	uint64_t mod_add(uint64_t a, uint64_t b) {
    	if (a >= PR - b) {
        	return a - (PR - b);
    	} else {
        	return a + b;
    	}
	}

	uint64_t mod_sub1(uint64_t a, uint64_t b) {
    	if (a < b) { 
        	return PR - (b - a);
    	} else {
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
#endif// MSPDZ_H__
