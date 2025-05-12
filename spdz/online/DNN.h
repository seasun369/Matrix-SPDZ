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
	mpz_t *a_bit, *b_bit, *c_bit, *a_bit_mac, *b_bit_mac, *c_bit_mac;

	NetIOMP<nP> *io;
	int party, total_pre, ssp;

	int m;
	int n;
	int l;
	int key_length;
	int max_m_n;
	int accuracy = 15; // bit of accuracy
	int l_bit = 61;
	int inputnodes;
	int hiddennodes;
	int outputnodes;
	int epoch;
	int batch_size; // batch size of DNN
	int batch_size_bit;
	int learningrate_DNN_bit; // learning rate of DNN equals to 0.125
	ThreadPool *pool;

	PRP prp;
	MSPDZ(NetIOMP<nP> *io[2], ThreadPool *pool, int party, int m, int n, int l, int inputnodes, int hiddennodes, int outputnodes, int batch_size, int key_length, int max_m_n, bool *_delta = nullptr, int ssp = 40)
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

		a_bit = new mpz_t[m*n];
		b_bit = new mpz_t[m*n];
		c_bit = new mpz_t[m*n];
		a_bit_mac = new mpz_t[m*n];
		b_bit_mac = new mpz_t[m*n];
		c_bit_mac = new mpz_t[m*n];

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
			mpz_init(a_bit[i]);
			mpz_init(b_bit[i]);
			mpz_init(c_bit[i]);
			mpz_init(a_bit_mac[i]);
			mpz_init(b_bit_mac[i]);
			mpz_init(c_bit_mac[i]);
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
			mpz_clear(a_bit[i]);
			mpz_clear(b_bit[i]);
			mpz_clear(c_bit[i]);
			mpz_clear(a_bit_mac[i]);
			mpz_clear(b_bit_mac[i]);
			mpz_clear(c_bit_mac[i]);
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
	void get_bitwise_triple(int l1, int l2, int equal){

		std::filesystem::path baseDir = (equal == 0) ?  "/mnt/extra_space/ztx/spdz-copy/pre_data/predata_mspdz/bitwise/0" : "/mnt/extra_space/ztx/spdz-copy/pre_data/predata_mspdz/bitwise/1";
		std::vector<std::string> filenames = (equal == 0) ?
    		std::vector<std::string>{"a_bitwise_mul.txt", "b_bitwise_mul.txt", "c_bitwise_mul.txt", "mac_a_bitwise_mul_spdz.txt","mac_b_bitwise_mul_spdz.txt","mac_c_bitwise_mul_spdz.txt"} :
    		std::vector<std::string>{"a_bitwise_mul.txt", "square_a_bitwise_mul.txt", "mac_a_bitwise_mul_spdz.txt","mac_square_a_bitwise_mul_spdz.txt"};
		
    	// 创建一个 std::vector 来存储 std::ifstream 对象
    	std::vector<std::ifstream> fileStreams;

    	// 预留空间可以提高效率（可选）
   	 	fileStreams.reserve(filenames.size());
		//cout << filenames.size() << endl;

    	//std::cout << "尝试打开文件..." << std::endl;
		for (const std::string& filename_ : filenames) {
			try {
				// 使用 emplace_back 直接在 vector 的末尾构造 ifstream 对象
				std::filesystem::path filename = baseDir / filename_;
				fileStreams.emplace_back(filename);
	
				// 检查刚刚添加的文件流是否成功打开
				if (!fileStreams.back()) { // 使用 operator bool() 检查状态
					 // 注意：如果构造函数本身抛出异常（某些实现可能这样做），这里可能不会执行
					 // 但通常，如果文件打不开，流对象会被创建，但处于错误状态
					std::cerr << "警告: 文件 '" << filename << "' 打开失败或流状态错误,但对象已添加到vector中。" << std::endl;
					// 你可以选择在这里移除这个无效的流，或者保留它并在后续处理中跳过
					// 例如：fileStreams.pop_back(); // 如果想移除失败的流
				} else {
					//std::cout << "文件 '" << filename << "' 成功添加到vector并打开。" << std::endl;
				}
			} catch (const std::ifstream::failure& e) {
				// 捕获可能的构造或打开异常 (如果设置了 exception mask)
				 std::cerr << "打开文件 '" << filename_ << "' 时发生异常: " << e.what() << std::endl;
			} catch (const std::exception& e) {
				 // 捕获其他可能的异常 (例如内存分配失败)
				 std::cerr << "处理文件 '" << filename_ << "' 时发生一般异常: " << e.what() << std::endl;
			}
		}

		std::vector<std::vector<triple>> allTriples;

		for (size_t i = 0; i < fileStreams.size(); ++i) {
			// 再次检查流是否处于良好状态，特别是如果之前没有移除失败的流
			if (fileStreams[i]) {
				std::string line;
    			std::vector<triple> Triples;

				while (std::getline(fileStreams[i], line)) {
					triple currentData;
					std::stringstream lineStream(line);
			
					// --- Parsing Logic ---
					// Assuming the first element is the party identifier (as an integer)
					if (!(lineStream >> currentData.party)) {
						std::cerr << "Warning: Failed to parse party ID from line: \"" << line << "\". Skipping line.\n";
						continue; // Skip this line if party ID is missing/invalid
					}
			
					// Read all subsequent numbers into the values vector
					uint64_t num;
					while (lineStream >> num) {
						currentData.value.push_back(num);
					}
			
					// Optional: Check if any non-numeric data was encountered after the party ID
					if (lineStream.fail() && !lineStream.eof()) {
						std::cerr << "Warning: Invalid numeric data encountered in line: \"" << line << "\". Data read so far kept.\n";
					}
			
					// Store the parsed data (using move semantics for efficiency)
					Triples.push_back(std::move(currentData));
				}
				allTriples.push_back(std::move(Triples));
				// 不需要手动 close()，当 fileStreams[i] 离开作用域或 vector 被销毁时，
				// ifstream 的析构函数会自动关闭文件（RAII）。
			} else {
				 std::cout << "\n--- 跳过文件: " << filenames[i] << " (之前打开失败或状态错误) ---" << std::endl;
			}
		}
		if(equal==0){
			for (int j = 0; j < l1*l2; j++){
				mpz_init_set_ui(a_bit[j], allTriples[0][party-1].value[j]);
				mpz_init_set_ui(b_bit[j], allTriples[1][party-1].value[j]);
				mpz_init_set_ui(c_bit[j], allTriples[2][party-1].value[j]);
				mpz_init_set_ui(a_bit_mac[j], allTriples[3][party-1].value[j]);
				mpz_init_set_ui(b_bit_mac[j], allTriples[4][party-1].value[j]);
				mpz_init_set_ui(c_bit_mac[j], allTriples[5][party-1].value[j]);
			}
		} else{
			for (int j = 0; j < l1*l2; j++){
				mpz_init_set_ui(a_bit[j], allTriples[0][party-1].value[j]);
				mpz_init_set(b_bit[j], a_bit[j]);
				mpz_init_set_ui(c_bit[j], allTriples[1][party-1].value[j]);
				mpz_init_set_ui(a_bit_mac[j], allTriples[2][party-1].value[j]);
				mpz_init_set(b_bit_mac[j], a_bit_mac[j]);
				mpz_init_set_ui(c_bit_mac[j], allTriples[3][party-1].value[j]);
			}
		}
		
	}
	
	void get_transform_tuple(int l1, int l3)
	{
		std::ifstream fin1;
		std::ifstream fin_1, fin_2;
		if ((l1 == 14) && (l3 == 14))
		{
			fin1.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/transform_predata/1/t.txt", std::ios::in);
			fin_1.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/transform_predata/1/mac_t.txt", std::ios::in);
			fin_2.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/transform_predata/1/mac_t_spdz.txt", std::ios::in);
		}
		else if ((l1 == 14) && (l3 == 1))
		{
			fin1.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/transform_predata/2/t.txt", std::ios::in);
			fin_1.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/transform_predata/2/mac_t.txt", std::ios::in);
			fin_2.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/transform_predata/2/mac_t_spdz.txt", std::ios::in);
		}
		else if ((l1 == 1913) && (l3 == 1))
		{
			fin1.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/transform_predata/4/t.txt", std::ios::in);
			fin_1.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/transform_predata/4/mac_t.txt", std::ios::in);
			fin_2.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/transform_predata/4/mac_t_spdz.txt", std::ios::in);
		}
		else if((l1 == 128) && (l3 == 10))
		{
			fin1.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/transform_predata/5/t.txt", std::ios::in);
			fin_1.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/transform_predata/5/mac_t.txt", std::ios::in);
			fin_2.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/transform_predata/5/mac_t_spdz.txt", std::ios::in);
		}
		else if((l1 == 128) && (l3 == 128))
		{
			fin1.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/transform_predata/6/t.txt", std::ios::in);
			fin_1.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/transform_predata/6/mac_t.txt", std::ios::in);
			fin_2.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/transform_predata/6/mac_t_spdz.txt", std::ios::in);
		}
		else if((l1 == 784) && (l3 == 128))
		{
			fin1.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/transform_predata/7/t.txt", std::ios::in);
			fin_1.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/transform_predata/7/mac_t.txt", std::ios::in);
			fin_2.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/transform_predata/7/mac_t_spdz.txt", std::ios::in);
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

	void get_transpose(int l1, int l2)
	{
		std::ifstream fin1, fin2;
		std::ifstream fin_1, fin_2;
		if ((l1 == 14) && (l2 == 1913))
		{
			fin1.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/transpose_predata/1/a.txt", std::ios::in);
			fin2.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/transpose_predata/1/a_t.txt", std::ios::in);
			fin_1.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/transpose_predata/1/mac_a.txt", std::ios::in);
			fin_2.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/transpose_predata/1/mac_a_t.txt", std::ios::in);
		}
		else if ((l1 == 14) && (l2 == 14))
		{
			fin1.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/transpose_predata/3/a.txt", std::ios::in);
			fin2.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/transpose_predata/3/a_t.txt", std::ios::in);
			fin_1.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/transpose_predata/3/mac_a.txt", std::ios::in);
			fin_2.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/transpose_predata/3/mac_a_t.txt", std::ios::in);
		}
		else if ((l1 == 1913) && (l2 == 14))
		{
			fin1.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/transpose_predata/4/a.txt", std::ios::in);
			fin2.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/transpose_predata/4/a_t.txt", std::ios::in);
			fin_1.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/transpose_predata/4/mac_a.txt", std::ios::in);
			fin_2.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/transpose_predata/4/mac_a_t.txt", std::ios::in);
		}
		else if ((l1 == 128) && (l2 == 10))
		{
			fin1.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/transpose_predata/5/a.txt", std::ios::in);
			fin2.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/transpose_predata/5/a_t.txt", std::ios::in);
			fin_1.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/transpose_predata/5/mac_a.txt", std::ios::in);
			fin_2.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/transpose_predata/5/mac_a_t.txt", std::ios::in);
		}
		else if ((l1 == 128) && (l2 == 128))
		{
			fin1.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/transpose_predata/6/a.txt", std::ios::in);
			fin2.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/transpose_predata/6/a_t.txt", std::ios::in);
			fin_1.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/transpose_predata/6/mac_a.txt", std::ios::in);
			fin_2.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/transpose_predata/6/mac_a_t.txt", std::ios::in);
		}
		else if ((l1 == 128) && (l2 == 784))
		{
			fin1.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/transpose_predata/7/a.txt", std::ios::in);
			fin2.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/transpose_predata/7/a_t.txt", std::ios::in);
			fin_1.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/transpose_predata/7/mac_a.txt", std::ios::in);
			fin_2.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/transpose_predata/7/mac_a_t.txt", std::ios::in);
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
				for (int j = 0; j < l1 * l2; j++)
					mpz_set_ui(a_t[j], Triple[i].value[j]);
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
	}

	void get_truncate_tuple(int l1, int l3)
	{
		std::ifstream fin1, fin2;
		std::ifstream fin_1, fin_2;
		if ((l1 == 14) && (l3 == 14))
		{
			fin1.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/truncate_predata/1/s.txt", std::ios::in);
			fin2.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/truncate_predata/1/s_bit.txt", std::ios::in);
			fin_1.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/truncate_predata/1/mac_s_spdz.txt", std::ios::in);
			fin_2.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/truncate_predata/1/mac_s_bit_spdz.txt", std::ios::in);
		}
		else if ((l1 == 14) && (l3 == 1))
		{
			fin1.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/truncate_predata/2/s.txt", std::ios::in);
			fin2.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/truncate_predata/2/s_bit.txt", std::ios::in);
			fin_1.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/truncate_predata/2/mac_s_spdz.txt", std::ios::in);
			fin_2.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/truncate_predata/2/mac_s_bit_spdz.txt", std::ios::in);
		}
		else if ((l1 == 1913) && (l3 == 1))
		{
			fin1.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/truncate_predata/3/s.txt", std::ios::in);
			fin2.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/truncate_predata/3/s_bit.txt", std::ios::in);
			fin_1.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/truncate_predata/3/mac_s_spdz.txt", std::ios::in);
			fin_2.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/truncate_predata/3/mac_s_bit_spdz.txt", std::ios::in);
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
			fin1.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/triple_predata/1/a.txt", std::ios::in);
			fin2.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/triple_predata/1/a_t.txt", std::ios::in);
			fin3.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/triple_predata/1/b.txt", std::ios::in);
			fin4.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/triple_predata/1/c.txt", std::ios::in);
			fin5.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/triple_predata/1/r.txt", std::ios::in);
			fin6.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/triple_predata/1/r_t.txt", std::ios::in);
			fin_1.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/triple_predata/1/mac_a.txt", std::ios::in);
			fin_2.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/triple_predata/1/mac_a_t.txt", std::ios::in);
			fin_3.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/triple_predata/1/mac_b.txt", std::ios::in);
			fin_4.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/triple_predata/1/mac_c.txt", std::ios::in);
			fin_5.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/triple_predata/1/mac_r.txt", std::ios::in);
			fin_6.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/triple_predata/1/mac_r_t.txt", std::ios::in);
		}
		else if ((l1 == 14) && (l2 == 1913) && (l3 == 1))
		{
			fin1.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/triple_predata/2/a.txt", std::ios::in);
			fin2.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/triple_predata/2/a_t.txt", std::ios::in);
			fin3.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/triple_predata/2/b.txt", std::ios::in);
			fin4.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/triple_predata/2/c.txt", std::ios::in);
			fin5.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/triple_predata/2/r.txt", std::ios::in);
			fin6.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/triple_predata/2/r_t.txt", std::ios::in);
			fin_1.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/triple_predata/2/mac_a.txt", std::ios::in);
			fin_2.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/triple_predata/2/mac_a_t.txt", std::ios::in);
			fin_3.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/triple_predata/2/mac_b.txt", std::ios::in);
			fin_4.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/triple_predata/2/mac_c.txt", std::ios::in);
			fin_5.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/triple_predata/2/mac_r.txt", std::ios::in);
			fin_6.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/triple_predata/2/mac_r_t.txt", std::ios::in);
		}
		else if ((l1 == 14) && (l2 == 14) && (l3 == 1))
		{
			fin1.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/triple_predata/3/a.txt", std::ios::in);
			fin2.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/triple_predata/3/a_t.txt", std::ios::in);
			fin3.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/triple_predata/3/b.txt", std::ios::in);
			fin4.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/triple_predata/3/c.txt", std::ios::in);
			fin5.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/triple_predata/3/r.txt", std::ios::in);
			fin6.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/triple_predata/3/r_t.txt", std::ios::in);
			fin_1.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/triple_predata/3/mac_a.txt", std::ios::in);
			fin_2.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/triple_predata/3/mac_a_t.txt", std::ios::in);
			fin_3.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/triple_predata/3/mac_b.txt", std::ios::in);
			fin_4.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/triple_predata/3/mac_c.txt", std::ios::in);
			fin_5.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/triple_predata/3/mac_r.txt", std::ios::in);
			fin_6.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/triple_predata/3/mac_r_t.txt", std::ios::in);
		}
		else if ((l1 == 1913) && (l2 == 14) && (l3 == 1))
		{
			fin1.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/triple_predata/4/a.txt", std::ios::in);
			fin2.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/triple_predata/4/a_t.txt", std::ios::in);
			fin3.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/triple_predata/4/b.txt", std::ios::in);
			fin4.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/triple_predata/4/c.txt", std::ios::in);
			fin5.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/triple_predata/4/r.txt", std::ios::in);
			fin6.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/triple_predata/4/r_t.txt", std::ios::in);
			fin_1.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/triple_predata/4/mac_a.txt", std::ios::in);
			fin_2.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/triple_predata/4/mac_a_t.txt", std::ios::in);
			fin_3.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/triple_predata/4/mac_b.txt", std::ios::in);
			fin_4.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/triple_predata/4/mac_c.txt", std::ios::in);
			fin_5.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/triple_predata/4/mac_r.txt", std::ios::in);
			fin_6.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/triple_predata/4/mac_r_t.txt", std::ios::in);
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
		fin1.open("/home/jackie/spdz/pre_data/predata_mspdz/DNN/triple_predata/1/key.txt", std::ios::in);
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

	void Online_bitwise_mul(mpz_t *x, mpz_t *y, mpz_t *output, mpz_t *output_mac, int l1, int l2, int equal)
	{
		get_bitwise_triple(l1,l2,equal);
		gmp_printf("]: %Zd\n", key[0]);

		mpz_t ppp;
		mpz_init(ppp);
		mpz_set_ui(ppp, 2305843009213693951UL);

		mpz_t *d = new mpz_t[l1 * l2];
		mpz_t *e = new mpz_t[l1 * l2];
		//mpz_t *f = new mpz_t[l3 * l1];
		uint64_t *d_int = new uint64_t[l1 * l2], *e_int = new uint64_t[l1 * l2];//*f_int = new uint64_t[l3 * l1];

		omp_set_num_threads(NUM_THREADS);
		#pragma omp parallel for
		for (int i = 0; i < l1 * l2; i++)
		{
			mpz_init(d[i]);
			mpz_sub(d[i], x[i], a_bit[i]);
			mpz_mod(d[i], d[i], ppp);
			d_int[i] = mpz_get_ui(d[i]);
		}

		omp_set_num_threads(NUM_THREADS);
		#pragma omp parallel for
		for (int i = 0; i < l1 * l2; i++)
		{
			mpz_init(e[i]);
			mpz_sub(e[i], y[i], b_bit[i]);
			mpz_mod(e[i], e[i], ppp);
			e_int[i] = mpz_get_ui(e[i]);
		}
		
		// open d
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
		
		// open e
		if (party != 1)
		{
			io->send_data(1, e_int, l1 * l2 * sizeof(uint64_t));
			io->flush(1);
			io->recv_data(1, e_int, l1 * l2 * sizeof(uint64_t));
			omp_set_num_threads(NUM_THREADS);
			#pragma omp parallel for
			for (int i = 0; i < l1 * l2; i++)
			{
				mpz_set_ui(e[i], e_int[i]);
			}
		}
		else
		{
			uint64_t *tmp[nP + 1];
			for (int i = 1; i <= nP; ++i)
			{
				tmp[i] = new uint64_t[l2 * l1];
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
				mpz_add(e[i], e[i], tmp_i[i]);
				mpz_mod(e[i], e[i], ppp);
				e_int[i] = mpz_get_ui(e[i]);
				mpz_clear(tmp_i[i]);
			}
			delete[] tmp_i;

			res.push_back(pool->enqueue([this, e_int, party2, l1, l2]()
										{
					io->send_data(party2, e_int, l1 * l2 * sizeof(uint64_t));
					io->flush(party2); }));
			joinNclean(res);
			for (int i = 1; i <= nP; ++i)
			{
				delete[] tmp[i];
			}
		}

		gmp_printf("a_bit : %Zd\n", a_bit[0]);
		gmp_printf("b_bit : %Zd\n", b_bit[0]);
		gmp_printf("d : %Zd\n", d[0]);
		gmp_printf("e : %Zd\n", e[0]);
		
		mpz_t *de = new mpz_t[l1 * l2];
		mpz_t *ae = new mpz_t[l1 * l2];
		mpz_t *db = new mpz_t[l1 * l2];

		mpz_t *de_spdz_mac = new mpz_t[l1 * l2];
		mpz_t *ae_spdz_mac = new mpz_t[l1 * l2];
		mpz_t *db_spdz_mac = new mpz_t[l1 * l2];

		omp_set_num_threads(NUM_THREADS);
		#pragma omp parallel for
		for (int i = 0; i < l1; i++){
			for( int j = 0; j < l2; j++){
				int idx = i * l2 + j;
				mpz_init(de[idx]);
				mpz_set_ui(de[idx], 0);
				mpz_init(ae[idx]);
				mpz_set_ui(ae[idx], 0);
				mpz_init(db[idx]);
				mpz_set_ui(db[idx], 0);
				mpz_init(de_spdz_mac[idx]);
				mpz_set_ui(de_spdz_mac[idx], 0);
				mpz_init(ae_spdz_mac[idx]);
				mpz_set_ui(ae_spdz_mac[idx], 0);
				mpz_init(db_spdz_mac[idx]);
				mpz_set_ui(db_spdz_mac[idx], 0);

				mpz_mul(de[idx], d[idx], e[idx]);
				mpz_mul(ae[idx], a_bit[idx], e[idx]);
				mpz_mul(db[idx], d[idx], b_bit[idx]);
				mpz_mul(de_spdz_mac[idx], de[idx], key[j]);
				mpz_mul(ae_spdz_mac[idx], a_bit_mac[idx], e[idx]);
				mpz_mul(db_spdz_mac[idx], d[idx], b_bit_mac[idx]);

				mpz_mod(de[idx], de[idx], ppp);
				mpz_mod(ae[idx], ae[idx], ppp);
				mpz_mod(db[idx], db[idx], ppp);
				mpz_mod(de_spdz_mac[idx], de_spdz_mac[idx], ppp);
				mpz_mod(ae_spdz_mac[idx], ae_spdz_mac[idx], ppp);
				mpz_mod(db_spdz_mac[idx], db_spdz_mac[idx], ppp);
			}
		}

		omp_set_num_threads(NUM_THREADS);
		#pragma omp parallel for
		for (int i = 0; i < l1 * l2; ++i)
		{
			if (party == 1)
			{
				mpz_add(output[i], c_bit[i], db[i]);
				mpz_add(output[i], output[i], de[i]);
				mpz_add(output[i], output[i], ae[i]);
				mpz_mod(output[i], output[i], ppp);
			}
			else
			{
				mpz_add(output[i], c_bit[i], db[i]);
				mpz_add(output[i], output[i], ae[i]);
				mpz_mod(output[i], output[i], ppp);
			}
		}

		omp_set_num_threads(NUM_THREADS);
		#pragma omp parallel for
		for (int i = 0; i < l1 * l2; i++)
		{
			mpz_add(output_mac[i], c_bit_mac[i], db_spdz_mac[i]);
			mpz_add(output_mac[i], output_mac[i], de_spdz_mac[i]);
			mpz_add(output_mac[i], output_mac[i], ae_spdz_mac[i]);
			mpz_mod(output_mac[i], output_mac[i], ppp);
		}

		omp_set_num_threads(NUM_THREADS);
		#pragma omp parallel for
		for (int i = 0; i < l1 * l2; i++)
		{
			mpz_clear(d[i]);
			mpz_clear(e[i]);
			mpz_clear(de[i]);
			mpz_clear(db[i]);
			mpz_clear(ae[i]);
			mpz_clear(de_spdz_mac[i]);
			mpz_clear(db_spdz_mac[i]);
			mpz_clear(ae_spdz_mac[i]);
		}

		mpz_clear(ppp);

		delete[] d;
		delete[] d_int;
		delete[] e;
		delete[] e_int;
		delete[] db;
		delete[] de;
		delete[] ae;
		delete[] db_spdz_mac;
		delete[] de_spdz_mac;
		delete[] ae_spdz_mac;
	}

	// derivative of x^2
	void x2_derivative(mpz_t *x, mpz_t *mac_x_spdz, int l1, int l3)
	{
		mpz_t ppp;
		mpz_t const1;
		mpz_init(ppp);
		mpz_init(const1);
		mpz_set_ui(ppp, 2305843009213693951UL);
		mpz_set_ui(const1, 2);

		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < l1 * l3; i++)
		{
			mpz_mul(x[i], x[i], const1);
			mpz_mul(mac_x_spdz[i], mac_x_spdz[i], const1);
			mpz_mod(x[i], x[i], ppp);
			mpz_mod(mac_x_spdz[i], mac_x_spdz[i], ppp);
		}

		mpz_clear(const1);
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

	void Gradient_Descent_DNN(mpz_t *x, mpz_t *w1, mpz_t *b1, mpz_t *w2, mpz_t *b2, mpz_t *y, mpz_t *mac_x, mpz_t *mac_w1, mpz_t *mac_b1, mpz_t *mac_w2, mpz_t *mac_b2, mpz_t *mac_y, int num_samples, int key_length)
	{
		// init
		mpz_t ppp;
		mpz_init(ppp);
		mpz_set_ui(ppp, 2305843009213693951UL);

		get_key(key_length);

		mpz_t *batch_input = new mpz_t[batch_size * inputnodes];
		mpz_t *batch_input_t = new mpz_t[batch_size * inputnodes];
		mpz_t *batch_label = new mpz_t[batch_size * outputnodes];
		mpz_t *batch_z1 = new mpz_t[batch_size * hiddennodes];
		mpz_t *batch_z1_t = new mpz_t[batch_size * hiddennodes];
		mpz_t *batch_a1 = new mpz_t[batch_size * hiddennodes];
		mpz_t *batch_a1_t = new mpz_t[batch_size * hiddennodes];
		mpz_t *batch_z2 = new mpz_t[batch_size * outputnodes];
		mpz_t *batch_a2 = new mpz_t[batch_size * outputnodes];
		mpz_t *dz1 = new mpz_t[batch_size * hiddennodes];
		mpz_t *dz1_t = new mpz_t[batch_size * hiddennodes];
		mpz_t *dz2 = new mpz_t[batch_size * outputnodes];
		mpz_t *dz2_t = new mpz_t[batch_size * outputnodes];
		mpz_t *middle1 = new mpz_t[batch_size * hiddennodes];
		mpz_t *middle2 = new mpz_t[hiddennodes * outputnodes];
		mpz_t *middle3 = new mpz_t[outputnodes];
		mpz_t *middle4 = new mpz_t[inputnodes * hiddennodes];
		mpz_t *middle5 = new mpz_t[hiddennodes];

		mpz_t *batch_input_mac = new mpz_t[batch_size];
		mpz_t *batch_input_t_mac = new mpz_t[inputnodes];
		mpz_t *batch_label_mac = new mpz_t[batch_size];
		mpz_t *batch_z1_mac = new mpz_t[batch_size];
		mpz_t *batch_z1_t_mac = new mpz_t[hiddennodes];
		mpz_t *batch_a1_mac = new mpz_t[batch_size];
		mpz_t *batch_a1_t_mac = new mpz_t[hiddennodes];
		mpz_t *batch_z2_mac = new mpz_t[batch_size];
		mpz_t *batch_a2_mac = new mpz_t[batch_size];
		mpz_t *dz1_mac = new mpz_t[hiddennodes];
		mpz_t *dz1_t_mac = new mpz_t[batch_size];
		mpz_t *dz2_mac = new mpz_t[batch_size];
		mpz_t *dz2_t_mac = new mpz_t[outputnodes];
		mpz_t *middle1_mac = new mpz_t[hiddennodes];
		mpz_t *middle2_mac = new mpz_t[hiddennodes];
		mpz_t middle3_mac;
		mpz_t *middle4_mac = new mpz_t[inputnodes];
		mpz_t middle5_mac;
		mpz_t *mac_b1_spdz = new mpz_t[hiddennodes];
		mpz_t *mac_b2_spdz = new mpz_t[outputnodes];
		mpz_t *batch_z1_mac_spdz = new mpz_t[batch_size * hiddennodes];
		mpz_t *batch_z1_t_mac_spdz = new mpz_t[batch_size * hiddennodes];
		mpz_t *batch_a1_mac_spdz = new mpz_t[batch_size * hiddennodes];
		mpz_t *batch_z2_mac_spdz = new mpz_t[batch_size * outputnodes];
		mpz_t *batch_a2_mac_spdz = new mpz_t[batch_size * outputnodes];
		mpz_t *batch_label_mac_spdz = new mpz_t[batch_size * outputnodes];
		mpz_t *dz1_mac_spdz = new mpz_t[batch_size * hiddennodes];
		mpz_t *dz1_t_mac_spdz = new mpz_t[batch_size * hiddennodes];
		mpz_t *dz2_mac_spdz = new mpz_t[batch_size * outputnodes];
		mpz_t *middle1_mac_spdz = new mpz_t[batch_size * hiddennodes];
		mpz_t *middle2_mac_spdz = new mpz_t[hiddennodes * outputnodes];
		mpz_t *middle3_mac_spdz = new mpz_t[outputnodes];
		mpz_t *middle4_mac_spdz = new mpz_t[inputnodes * hiddennodes];
		mpz_t *middle5_mac_spdz = new mpz_t[hiddennodes];

		mpz_init(middle3_mac);
		mpz_set_ui(middle3_mac, 0);
		mpz_init(middle5_mac);
		mpz_set_ui(middle5_mac, 0);
		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < batch_size * inputnodes; i++)
		{
			mpz_init(batch_input[i]);
			mpz_init(batch_input_t[i]);
		}
		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < batch_size * hiddennodes; i++)
		{
			mpz_init(batch_z1[i]);
			mpz_init(batch_z1_t[i]);
			mpz_init(batch_a1[i]);
			mpz_init(batch_a1_t[i]);
			mpz_init(dz1[i]);
			mpz_init(dz1_t[i]);
			mpz_init(middle1[i]);
			mpz_init(batch_z1_mac_spdz[i]);
			mpz_init(batch_z1_t_mac_spdz[i]);
			mpz_init(batch_a1_mac_spdz[i]);
			mpz_init(dz1_mac_spdz[i]);
			mpz_init(dz1_t_mac_spdz[i]);
			mpz_init(middle1_mac_spdz[i]);
		}
		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < batch_size * outputnodes; i++)
		{
			mpz_init(batch_label[i]);
			mpz_init(batch_z2[i]);
			mpz_init(batch_a2[i]);
			mpz_init(dz2[i]);
			mpz_init(dz2_t[i]);
			mpz_init(batch_z2_mac_spdz[i]);
			mpz_init(batch_a2_mac_spdz[i]);
			mpz_init(batch_label_mac_spdz[i]);
			mpz_init(dz2_mac_spdz[i]);
		}
		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < hiddennodes * outputnodes; i++)
		{
			mpz_init(middle2[i]);
			mpz_init(middle2_mac_spdz[i]);
		}
		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < hiddennodes * inputnodes; i++)
		{
			mpz_init(middle4[i]);
			mpz_init(middle4_mac_spdz[i]);
		}
		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < batch_size; i++)
		{
			mpz_init(batch_input_mac[i]);
			mpz_init(batch_label_mac[i]);
			mpz_init(batch_z1_mac[i]);
			mpz_init(batch_a1_mac[i]);
			mpz_init(batch_z2_mac[i]);
			mpz_init(batch_a2_mac[i]);
			mpz_init(dz1_t_mac[i]);
			mpz_init(dz2_mac[i]);
		}
		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < inputnodes; i++)
		{
			mpz_init(batch_input_t_mac[i]);
			mpz_init(middle4_mac[i]);
		}
		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < hiddennodes; i++)
		{
			mpz_init_set_ui(middle5[i], 0);
			mpz_init_set_ui(mac_b1_spdz[i], 0); // because the shares of b1 are initialized to all zeros
			mpz_init(dz1_mac[i]);
			mpz_init(middle1_mac[i]);
			mpz_init(middle2_mac[i]);
			mpz_init(batch_z1_t_mac[i]);
			mpz_init(batch_a1_t_mac[i]);
		}
		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < outputnodes; i++)
		{
			mpz_init_set_ui(mac_b2_spdz[i], 0); // the same as mac_b2_spdz
			mpz_init(dz2_t_mac[i]);
			mpz_init_set_ui(middle3[i], 0);
			mpz_init_set_ui(middle3_mac_spdz[i], 0);
		}

		// train epoch times
		for (int e = 0; e < epoch; e++)
		{
			cout << "第" << (e + 1) << "轮次开始执行" << endl;

			// set batch_input & batch_label
			for (int i = 0; i < num_samples; i += batch_size)
			{
				omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
				for (int j1 = 0; j1 < batch_size * inputnodes; j1++)
				{
					mpz_set(batch_input[j1], x[i * inputnodes + j1]);
				}
				omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
				for (int j2 = 0; j2 < batch_size * outputnodes; j2++)
				{
					mpz_set(batch_label[j2], x[i * outputnodes + j2]);
				}
				omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
				for (int j3 = 0; j3 < batch_size; j3++)
				{
					mpz_set(batch_input_mac[j3], mac_x[i + j3]);
					mpz_set(batch_label_mac[j3], mac_y[i + j3]);
				}

				// generate batch_label_mac_spdz
				transform_vectormac_to_spdzmac(batch_label, batch_label_mac, batch_label_mac_spdz, batch_size, outputnodes);

				// forward propagation
				forward_propagation(batch_input, w1, b1, batch_z1, batch_input_mac, mac_w1, mac_b1_spdz, batch_z1_mac, batch_z1_mac_spdz, batch_size, inputnodes, hiddennodes, key_length);
				Online_bitwise_mul(batch_z1, batch_z1, batch_z1_mac_spdz, batch_z1_mac_spdz, batch_a1, batch_a1_mac_spdz, batch_size, hiddennodes, key_length, 1);
				truncate(batch_a1, batch_a1_mac_spdz, accuracy, batch_size, hiddennodes, key_length);
				transform_spdzmac_to_vectormac(batch_a1_mac, batch_a1_mac_spdz, batch_size, hiddennodes);

				forward_propagation(batch_a1, w2, b2, batch_z2, batch_a1_mac, mac_w2, mac_b2_spdz, batch_z2_mac, batch_z2_mac_spdz, batch_size, hiddennodes, outputnodes, key_length);
				Online_bitwise_mul(batch_z2, batch_z2, batch_z2_mac_spdz, batch_z2_mac_spdz, batch_a2, batch_a2_mac_spdz, batch_size, outputnodes, key_length, 1);
				truncate(batch_a2, batch_a2_mac_spdz, accuracy, batch_size, outputnodes, key_length);
				// transform_spdzmac_to_vectormac(batch_a2_mac, batch_a2_mac_spdz, batch_size, outputnodes);

				// back propagation
				// generate dz2 & dz2_mac_spdz
				x2_derivative(batch_z2, batch_z2_mac_spdz, batch_size, outputnodes);
				omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
				for (int j4 = 0; j4 < batch_size * outputnodes; j4++)
				{
					mpz_sub(batch_a2[j4], batch_a2[j4], batch_label[j4]);
					mpz_mod(batch_a2[j4], batch_a2[j4], ppp);
					mpz_sub(batch_a2_mac_spdz[j4], batch_a2_mac_spdz[j4], batch_label_mac_spdz[j4]);
					mpz_mod(batch_a2_mac_spdz[j4], dz2_mac_spdz[j4], ppp);
				}
				Online_bitwise_mul(batch_a2, batch_z2, batch_a2_mac_spdz, batch_z2_mac_spdz, dz2, dz2_mac_spdz, batch_size, outputnodes, key_length, 0);
				truncate(dz2, dz2_mac_spdz, accuracy, batch_size, outputnodes, key_length);

				// generate dz1 & dz1_mac_spdz
				// compute z1_t & 2 * z1_t
				transform_spdzmac_to_vectormac(batch_z1_mac, batch_z1_mac_spdz, batch_size, hiddennodes);
				transpose(batch_z1, batch_z1_mac, batch_z1_t, batch_z1_t_mac, batch_size, hiddennodes);
				transform_vectormac_to_spdzmac(batch_z1_t, batch_z1_t_mac, batch_z1_t_mac_spdz, hiddennodes, batch_size);
				x2_derivative(batch_z1_t, batch_z1_t_mac_spdz, hiddennodes, batch_size);
				// compute dz2_t
				transform_spdzmac_to_vectormac(dz2_mac, dz2_mac_spdz, batch_size, outputnodes);
				transpose(dz2, dz2_mac, dz2_t, dz2_t_mac, batch_size, outputnodes);
				Online_nonsqure_mul(w2, dz2_t, mac_w2, dz2_t_mac, middle1, middle1_mac, hiddennodes, outputnodes, batch_size, key_length);
				transform_vectormac_to_spdzmac(middle1, middle1_mac, middle1_mac_spdz, hiddennodes, batch_size);
				// compute dz1
				Online_bitwise_mul(middle1, batch_z1_t, middle1_mac_spdz, batch_z1_t_mac_spdz, dz1, dz1_mac_spdz, hiddennodes, batch_size, key_length, 0);

				// update w2
				// compute a1_t
				transpose(batch_a1, batch_a1_mac, batch_a1_t, batch_a1_t_mac, batch_size, hiddennodes);
				// compute a1_t * dz2
				Online_nonsqure_mul(batch_a1_t, dz2, batch_a1_t_mac, dz2_mac, middle2, middle2_mac, hiddennodes, batch_size, outputnodes, key_length);
				// update
				transform_vectormac_to_spdzmac(middle2, middle2_mac, middle2_mac_spdz, hiddennodes, outputnodes);
				truncate(middle2, middle2_mac_spdz, accuracy + batch_size_bit + learningrate_DNN_bit, batch_size, outputnodes, key_length);
				transform_spdzmac_to_vectormac(middle2_mac, middle2_mac_spdz, hiddennodes, outputnodes);
				omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
				for (int j5 = 0; j5 < hiddennodes * outputnodes; j5++)
				{
					mpz_sub(w2[j5], w2[j5], middle2[j5]);
					mpz_mod(middle2[j5], middle2[j5], ppp);
				}
				omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
				for (int j6 = 0; j6 < hiddennodes; j6++)
				{
					mpz_sub(mac_w2[j6], mac_w2[j6], middle2_mac[j6]);
					mpz_mod(mac_w2[j6], mac_w2[j6], ppp);
				}

				// update b2
				// compute sum
				omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
				for (int j7 = 0; j7 < outputnodes; j7++)
				{
					for (int k = 0; k < batch_size; k++)
					{
						mpz_add(middle3[j7], middle3[j7], dz2_t[j7 * batch_size + k]);
					}
					mpz_mod(middle3[j7], middle3[j7], ppp);
				}
				omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
				for (int j8 = 0; j8 < outputnodes; j8++)
				{
					for (int k = 0; k < batch_size; b++)
					{
						mpz_add(middle3_mac_spdz[j8], middle3_mac_spdz[j8], dz2_mac_spdz[k * outputnodes + j8]);
					}
					mpz_mod(middle3_mac_spdz[j8], middle3_mac_spdz[j8], ppp);
				}
				truncate(middle3, middle3_mac_spdz, batch_size_bit + learningrate_DNN_bit, 1, outputnodes, key_length);
				// update
				transform_spdzmac_to_vectormac(middle3_mac, middle3_mac_spdz, 1, outputnodes);
				omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
				for (int j9 = 0; j9 < outputnodes; j9++)
				{
					mpz_sub(b2[j9], b2[j9], middle3[j9]);
					mpz_mod(b2[j9], b2[j9], ppp);
				}
				mpz_sub(mac_b1[0], mac_b1[0], middle3_mac);
				mpz_mod(mac_b2[0], mac_b2[0], ppp);

				// update w1
				// compute batch_input_t
				transpose(batch_input, batch_input_mac, batch_input_t, batch_input_t_mac, batch_size, inputnodes);
				// compute dz1_t
				transpose(dz1, dz1_mac, dz1_t, dz1_t_mac, hiddennodes, batch_size);
				// compute batch_input_t * dz1_t
				Online_nonsqure_mul(batch_input_t, dz1_t, batch_input_t_mac, dz1_t_mac, middle4, middle4_mac, inputnodes, batch_size, hiddennodes, key_length);
				transform_vectormac_to_spdzmac(middle4, middle4_mac, middle4_mac_spdz, inputnodes, hiddennodes);
				truncate(middle4, middle4_mac_spdz, accuracy + batch_size_bit + learningrate_DNN_bit, inputnodes, hiddennodes, key_length);
				transform_spdzmac_to_vectormac(middle4_mac, middle4_mac_spdz, inputnodes, hiddennodes);
				// update
				omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
				for (int j10 = 0; j10 < inputnodes * hiddennodes; j10++)
				{
					mpz_sub(w1[j10], w1[j10], middle4[j10]);
					mpz_mod(w1[j10], w1[j10], ppp);
				}
				omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
				for (int j10 = 0; j10 < inputnodes; j10++)
				{
					mpz_sub(mac_w1[j10], mac_w1[j10], middle4_mac[j10]);
					mpz_mod(mac_w1[j10], mac_w1[j10], ppp);
				}

				// update b1
				omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
				for (int j11 = 0; j11 < hiddennodes; j11++)
				{
					for (int k = 0; k < batch_size; k++)
					{
						mpz_add(middle5[j11], middle5[j11], dz1[j11 * batch_size + k]);
					}
					mpz_mod(middle5[j11], middle5[j11], ppp);
				}
				transform_vectormac_to_spdzmac(dz1_t, dz1_t_mac, dz1_t_mac_spdz, batch_size, hiddennodes);
				omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
				for (int j12 = 0; j12 < hiddennodes; j12++)
				{
					for (int k = 0; k < batch_size; k++)
					{
						mpz_add(middle5_mac_spdz[j12], middle5_mac_spdz[j12], dz1_t_mac_spdz[j12 + k * hiddennodes]);
					}
					mpz_mod(middle5_mac_spdz[j12], middle5_mac_spdz[j12], ppp);
				}
				truncate(middle5, middle5_mac_spdz, batch_size_bit + learningrate_DNN_bit, 1, hiddennodes, key_length);
				transform_spdzmac_to_vectormac(middle5_mac, middle5_mac_spdz, 1, hiddennodes);
				omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
				for (int j13 = 0; j13 < hiddennodes; j13++)
				{
					mpz_sub(b1[j13], b1[j13], middle5[j13]);
					mpz_mod(b1[j13], b1[j13], ppp);
				}
				mpz_sub(mac_b1[0], mac_b1[0], middle5_mac);
				mpz_mod(mac_b1[0], mac_b1[0], ppp);
			}
		}

		// check
		mac_check(w1, mac_w1, inputnodes, hiddennodes);
		mac_check(w2, mac_w2, hiddennodes, outputnodes);

		// print parameter
		// print w1
		int w1_x = inputnodes;
		int w1_y = hiddennodes;
		uint64_t *w1_int = new uint64_t[inputnodes * hiddennodes];
		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < inputnodes * hiddennodes; i++)
		{
			w1_int[i] = mpz_get_ui(w1[i]);
		}

		if (party != 1)
		{
			io->send_data(1, w1_int, inputnodes * hiddennodes * sizeof(uint64_t));
			io->flush(1);
			io->recv_data(1, w1_int, inputnodes * hiddennodes * sizeof(uint64_t));
			omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
			for (int i = 0; i < inputnodes * hiddennodes; i++)
			{
				mpz_set_ui(w1[i], w1_int[i]);
			}
		}
		else
		{
			uint64_t *tmp[nP + 1];
			for (int i = 1; i <= nP; ++i)
			{
				tmp[i] = new uint64_t[inputnodes * hiddennodes];
			}
			vector<future<void>> res;
			int party2 = 2;
			res.push_back(pool->enqueue([this, tmp, party2, w1_x, w1_y]()
										{ io->recv_data(party2, tmp[party2], inputnodes * hiddennodes * sizeof(uint64_t)); }));
			joinNclean(res);
			mpz_t *tmp_i = new mpz_t[inputnodes * hiddennodes];
			omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
			for (int i = 0; i < inputnodes * hiddennodes; ++i)
			{
				mpz_init(tmp_i[i]);
				mpz_set_ui(tmp_i[i], tmp[2][i]);
				mpz_add(w1[i], w1[i], tmp_i[i]);
				mpz_mod(w1[i], w1[i], ppp);
				w1_int[i] = mpz_get_ui(w1[i]);
			}
			res.push_back(pool->enqueue([this, w1_int, party2, w1_x, w1_y]()
										{
					io->send_data(party2, w1_int, inputnodes * hiddennodes * sizeof(uint64_t));
					io->flush(party2); }));
			joinNclean(res);
			for (int i = 1; i <= nP; ++i)
			{
				delete[] tmp[i];
			}
		}
		for (int i = 0; i < 5; i++)
		{
			cout << "Value of w1[";
			cout << i;
			gmp_printf("]: %Zd\n", w1[i]);
		}

		// delete
		mpz_clear(middle3_mac);
		mpz_clear(middle5_mac);
		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < batch_size * inputnodes; i++)
		{
			mpz_clear(batch_input[i]);
			mpz_clear(batch_input_t[i]);
		}
		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < batch_size * hiddennodes; i++)
		{
			mpz_clear(batch_z1[i]);
			mpz_clear(batch_z1_t[i]);
			mpz_clear(batch_a1[i]);
			mpz_clear(batch_a1_t[i]);
			mpz_clear(dz1[i]);
			mpz_clear(dz1_t[i]);
			mpz_clear(middle1[i]);
			mpz_clear(batch_z1_mac_spdz[i]);
			mpz_clear(batch_z1_t_mac_spdz[i]);
			mpz_clear(batch_a1_mac_spdz[i]);
			mpz_clear(dz1_mac_spdz[i]);
			mpz_clear(dz1_t_mac_spdz[i]);
			mpz_clear(middle1_mac_spdz[i]);
		}
		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < batch_size * outputnodes; i++)
		{
			mpz_clear(batch_label[i]);
			mpz_clear(batch_z2[i]);
			mpz_clear(batch_a2[i]);
			mpz_clear(dz2[i]);
			mpz_clear(dz2_t[i]);
			mpz_clear(batch_z2_mac_spdz[i]);
			mpz_clear(batch_a2_mac_spdz[i]);
			mpz_clear(batch_label_mac_spdz[i]);
			mpz_clear(dz2_mac_spdz[i]);
		}
		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < hiddennodes * outputnodes; i++)
		{
			mpz_clear(middle2[i]);
			mpz_clear(middle2_mac_spdz[i]);
		}
		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < hiddennodes * inputnodes; i++)
		{
			mpz_clear(middle4[i]);
			mpz_clear(middle4_mac_spdz[i]);
		}
		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < batch_size; i++)
		{
			mpz_clear(batch_input_mac[i]);
			mpz_clear(batch_label_mac[i]);
			mpz_clear(batch_z1_mac[i]);
			mpz_clear(batch_a1_mac[i]);
			mpz_clear(batch_z2_mac[i]);
			mpz_clear(batch_a2_mac[i]);
			mpz_clear(dz1_t_mac[i]);
			mpz_clear(dz2_mac[i]);
		}
		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < inputnodes; i++)
		{
			mpz_clear(batch_input_t_mac[i]);
			mpz_clear(middle4_mac[i]);
		}
		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < hiddennodes; i++)
		{
			mpz_clear(middle5[i]);
			mpz_clear(mac_b1_spdz[i]); // because the shares of b1 are initialized to all zeros
			mpz_clear(dz1_mac[i]);
			mpz_clear(middle1_mac[i]);
			mpz_clear(middle2_mac[i]);
			mpz_clear(batch_z1_t_mac[i]);
			mpz_clear(batch_a1_t_mac[i]);
		}
		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < outputnodes; i++)
		{
			mpz_clear(mac_b2_spdz[i]); // the same as mac_b2_spdz
			mpz_clear(dz2_t_mac[i]);
			mpz_clear(middle3[i]);
			mpz_clear(middle3_mac_spdz[i]);
		}

		delete [] batch_input;
		delete [] batch_input_t;
		delete [] batch_z1;
		delete [] batch_z1_t;
		delete [] batch_a1;
		delete [] batch_a1_t;
		delete [] batch_z2;
		delete [] batch_a2;
		delete [] dz1;
		delete [] dz1_t;
		delete [] dz2;
		delete [] dz2_t;
		delete [] middle1;
		delete [] middle2;
		delete [] middle3;
		delete [] middle4;
		delete [] middle5;

		delete [] batch_input_mac;
		delete [] batch_input_t_mac;
		delete [] batch_label_mac;
		delete [] batch_z1_mac;
		delete [] batch_z1_t_mac;
		delete [] batch_a1_mac;
		delete [] batch_a1_t_mac;
		delete [] batch_z2_mac;
		delete [] batch_a2_mac;
		delete [] dz1_mac;
		delete [] dz1_t_mac;
		delete [] dz2_mac;
		delete [] dz2_t_mac;
		delete [] middle1_mac;
		delete [] middle2_mac;
		delete [] middle4_mac;
		delete [] mac_b1_spdz;
		delete [] mac_b2_spdz;
		delete [] batch_z1_mac_spdz;
		delete [] batch_z1_t_mac_spdz;
		delete [] batch_a1_mac_spdz;
		delete [] batch_z2_mac_spdz;
		delete [] batch_a2_mac_spdz;
		delete [] batch_label_mac_spdz;
		delete [] dz1_mac_spdz;
		delete [] dz1_t_mac_spdz;
		delete [] dz2_mac_spdz;
		delete [] middle1_mac_spdz;
		delete [] middle2_mac_spdz;
		delete [] middle3_mac_spdz;
		delete [] middle4_mac_spdz;
		delete [] middle5_mac_spdz;
	}

	void transform_vectormac_to_spdzmac(mpz_t *x, mpz_t *x_mac, mpz_t *x_mac_spdz, int l1, int l3)
	{
		mpz_t ppp;
		mpz_init(ppp);
		mpz_set_ui(ppp, 2305843009213693951UL);
		get_transform_tuple(l1, l3);

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

	void transform_spdzmac_to_vectormac(mpz_t *x_mac, mpz_t *x_mac_spdz, int l1, int l2)
	{
		mpz_t ppp;
		mpz_init(ppp);
		mpz_set_ui(ppp, 2305843009213693951UL);

		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < l1; i++)
		{
			mpz_set_ui(x_mac[i], 0);
			for (int j = 0; j < l2; j++)
			{
				mpz_add(x_mac[i], x_mac[i], x_mac_spdz[i * l2 + j]);
			}
			mpz_mod(x_mac[i], x_mac[i], ppp);
		}
	}

	void transpose(mpz_t *x, mpz_t *x_mac, mpz_t *x_t, mpz_t *x_t_mac, int l1, int l2)
	{
		mpz_t ppp;
		mpz_t *d = new mpz_t[l1 * l2];
		mpz_t *derivative_mac = new mpz_t[l2];
		uint64_t *d_int = new uint64_t[l1 * l2];
		mpz_init(ppp);
		mpz_set_ui(ppp, 2305843009213693951UL);

		get_transpose(l1, l2);

		// open
		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < l1 * l2; i++)
		{
			mpz_init(d[i]);
			mpz_sub(d[i], x[i], a[i]);
			mpz_mod(d[i], d[i], ppp);
			d_int[i] = mpz_get_ui(d[i]);
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
			}
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

		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < l2; i++)
		{
			mpz_init(derivative_mac[i]);
		}

		// compute x_t and mac_x_t
		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < l2; i++)
		{
			mpz_set(x_t_mac[i], mac_a_t[i]);
			for (int j = 0; j < l1; j++)
			{
				mpz_set(x_t[i * l1 + j], x[j * l2 + i]);
				mpz_mul(derivative_mac[i], d[j * l2 + i], key[j]);
				mpz_add(x_t_mac[i], x_t_mac[i], derivative_mac[i]);
			}
			mpz_mod(x_t_mac[i], x_t_mac[i], ppp);
		}

		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < l1 * l2; i++)
		{
			mpz_clear(d[i]);
		}
		omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
		for (int i = 0; i < l2; i++)
		{
			mpz_clear(derivative_mac[i]);
		}

		delete[] d;
		delete[] d_int;
		delete[] derivative_mac;
	}

	void truncate(mpz_t *x, mpz_t *x_mac_spdz, int accuracy, int l1, int l3, int key_length)
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

		get_truncate_tuple(l1, l3);

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

	void forward_propagation(mpz_t *x, mpz_t *w, mpz_t *bias, mpz_t *output, mpz_t *mac_x, mpz_t *mac_w, mpz_t *mac_bias_spdz, mpz_t *mac_output, mpz_t *mac_output_spdz, int l1, int l2, int l3, int key_length)
	{
		// the mac_output here computed is not consist with mac_output_spdz(the spdz version is correct)
		mpz_t ppp;
		mpz_init(ppp);
		mpz_set_ui(ppp, 2305843009213693951UL);

		// compute xw
		Online_nonsqure_mul(x, w, mac_x, mac_w, output, mac_output, l1, l2, l3, key_length);

		// truncate
		transform_vectormac_to_spdzmac(output, mac_output, mac_output_spdz, l1, l3);
		truncate(output, mac_output_spdz, accuracy, l1, l3, key_length);
		// transform_spdzmac_to_vectormac(mac_output, mac_output_spdz, l1, l3);

		// compute xw+bias and mac
		for (int i = 0; i < l1; i++)
		{
			for (int j = 0; j < l3; j++)
			{
				mpz_add(output[i * l3 + j], output[i * l3 + j], bias[j]);
				mpz_mod(output[i * l3 + j], output[i * l3 + j], ppp);
				mpz_add(mac_output_spdz[i * l3 + j], mac_output_spdz[i * l3 + j], mac_bias_spdz[j]);
				mpz_mod(mac_output_spdz[i * l3 + j], mac_output_spdz[i * l3 + j], ppp);
			}
		}
	}
};
#endif // MSPDZ_H__