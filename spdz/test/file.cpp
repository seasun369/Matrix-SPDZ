//just test file stream
#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include <sstream>
struct triple{
    std::string party;
    std::vector<uint64_t> value;
};
int main( int argc, char** argv )
{
std::ifstream fin;
fin.open("/home/seasun/spdz/pre_data/triple_a.txt",std::ios::in);
if(!fin.is_open())
{
    std::cerr<<"cannot open the file";
}
char line[1024]={0};
std::vector<triple> Triple;
//从文件中提取“行”
while(fin.getline(line,sizeof(line)))
{
    //定义局部变量
    triple t;
    //从“行”中提取“单词”
    std::stringstream word(line);
    word>>t.party;
    uint64_t num;
    while(word>>num)
        t.value.push_back(num);
    Triple.push_back(t);
}
 std::cout<<stoi(Triple[0].party)<<"'s value and mac are:"<< Triple[0].value[0] <<" "<< Triple[0].value[1] << std::endl;  
 std::cout<<Triple[1].party<<"'s value and mac are:"<< Triple[1].value[0] <<" "<< Triple[0].value[1] << std::endl;
 std::cout<<Triple[2].party<<"'s value and mac are:"<< Triple[2].value[0] <<" "<< Triple[0].value[1] << std::endl;  
}