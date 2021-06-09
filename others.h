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
   ����������������ʱ�Ժ������
*/
using namespace std;



double get_length(double* a);


//��������
double get_dot_product(double*a, double*b, int n);//�����
double* get_cross_product(double*a, double*b, int n);//�����
double det(double **a, int n);
void inv(int n, double **a, double **b);
void trans(int ax, int ay, double **a, double **b);
void cal_rp(double **d, double x, double y, double z, double *rp);
double dot_product(double ***a, double ***b, int nx, int ny, int nz);
double dot_product(double ***a_r, double ***a_i, double ***b_r, double ***b_i, int nx, int ny, int nz);


//���ӵĲ��䣬���߲����ú���
void plus_stress_mask(string path);//����stress_mask
void plus_etot_error(string path, string etot_file, string report_file = "REPORT_scf_test");//����error
void plus_init_magnetic(string path,int* type_list,int num,int value=1);//���ӳ�ʼ�ž�

bool check_relax(string path);//����ǲ���relax����
bool check_now_exist_file(string name);//�����û������ļ�
vector<string> get_file_name(string file_name);//���name�б�
void plusRelaxQdiv(string path, string file_name);//��relax��qdiv�����Ϊ��ʼ��scf��magnetic

//vector�־û�
//template <typename T>
void vector2file(vector<vector<double>>&data, ofstream &fout);
//���տո�ָ��ַ���
vector<string> cut_string(string  input);
#endif // !OTHERS
