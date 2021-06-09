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
//����һЩ�����ĳ��ó���

const double e_elec = 1.602177e-19;
const double kulun = e_elec;
const double plunk_const = 6.6206e-34;
//const double plunk_const = 4.1356e-15;//���Ӹ��صĵ�λ
const double pi = 3.1415926;
const double kb = 1.38065e-23;//��λ��JK
const double T_temp = 303.15;//25���϶�
const double plunk_ba = plunk_const / 2 / pi;
const double m0 = 9.10956e-31;





//int run_num_mobi = 0;
//mobility���������
void mobility_client(string config_name, string path);
//���Ǩ���ʵ���ؽ��
void get_mobility_data(string path);


//���ƶ�·����׼����Ǩ���ʵ��ļ�
bool prepare_mobility(string file_name, string path, string config);

//�ڵ�ǰ·����׼���õ���������,relax
void prepare_mobility_singel_relax(double factor,cell cell_a,string name);

//�ڵ�ǰ·����׼���õ���������,scf
void prepare_mobility_singel_scf(double factor,string name,cell& cell_a,string path);

//׼����������Ӧ��ϵ��
vector<double> get_factor(double start, double increase,int num);

//string��conventer����
string conventer(string input);

//���Ǩ������������
double mobility_const(void);

//���Ǩ����
double get_mobility(const double C, const double m, const double E);
#endif MOBILITY