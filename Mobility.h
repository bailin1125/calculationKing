#pragma once
#ifndef MOBILITY
#define MOBILITY
#include "calculation_tools.h"
#include <string>
#include <string.h>
#include <istream>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <sstream>
#include <iomanip>
#include <vector>
using namespace std;
const int direction = 0;
const string result_name = "monility_res_energy";
//定义一些其他的常用常数

const double e_elec = 1.602177e-19;
const double kulun = e_elec;
const double plunk_const = 6.6206e-34;
//const double plunk_const = 4.1356e-15;//电子福特的单位
const double pi = 3.1415926;
const double kb = 1.38065e-23;//单位是JK
const double T_temp = 303.15;//25摄氏度
const double plunk_ba = plunk_const / 2 / pi;
const double m0 = 9.10956e-31;





//int run_num_mobi = 0;
//mobility的整天调用
void mobility_client(string config_name, string path);
//获得迁移率的相关结果
void get_mobility_data(string path);


//在制定路径下准备好迁移率的文件
bool prepare_mobility(string file_name, string path, string config);

//在当前路径下准备好单个的输入,relax
void prepare_mobility_singel_relax(double factor,cell cell_a,string name);

//在当前路径下准备好单个的输入,scf
void prepare_mobility_singel_scf(double factor,string name,cell& cell_a,string path);

//准备好连续的应变系数
vector<double> get_factor(double start, double increase,int num);

//string的conventer函数
string conventer(string input);

//获得迁移率其他参数
double mobility_const(void);

//获得迁移率
double get_mobility(const double C, const double m, const double E);
#endif MOBILITY