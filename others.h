#pragma once
#ifndef OTHERS
#define OTHERS
#include <fstream>
#include <iostream>
#include <math.h>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <vector>


/**
   其他辅助函数，临时性函数相关
*/
using namespace std;



double get_length(double* a);


//逆矩阵相关
double get_dot_product(double*a, double*b, int n);//求解点乘
double* get_cross_product(double*a, double*b, int n);//求解叉乘
double det(double **a, int n);
void inv(int n, double **a, double **b);
void trans(int ax, int ay, double **a, double **b);
void cal_rp(double **d, double x, double y, double z, double *rp);
double dot_product(double ***a, double ***b, int nx, int ny, int nz);
double dot_product(double ***a_r, double ***a_i, double ***b_r, double ***b_i, int nx, int ny, int nz);


//增加的补充，或者测试用函数
void plus_stress_mask(string path);//增加stress_mask
void plus_etot_error(string path, string etot_file, string report_file = "REPORT_scf_test");//增加error
void plus_init_magnetic(string path,int* type_list,int num,int value=1);//增加初始磁矩

bool check_relax(string path);//检查是不是relax好了
bool check_now_exist_file(string name);//检查有没有这个文件
vector<string> get_file_name(string file_name);//获得name列表
void plusRelaxQdiv(string path, string file_name);//将relax的qdiv结果作为初始的scf的magnetic

//vector持久化
//template <typename T>
void vector2file(vector<vector<double>>&data, ofstream &fout);
//按照空格分隔字符串
vector<string> cut_string(string  input);
#endif // !OTHERS
