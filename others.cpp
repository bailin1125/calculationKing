#include "others.h"
#include "config.h"
#include <cstdlib>
#include <math.h>
#include <unistd.h>
#include<iomanip>
#pragma warning(disable : 4996)

using namespace std;











//�����ĸ����õĺ���
double det(double **a, int n)
{
	int i, j, k;
	int ii, jj;
	int a_i, a_j;
	a_i = -1;
	a_j = -1;
	double m;
	m = 0;
	double **b;
	int f;
	b = (double **)malloc((n - 1) * sizeof(double*));
	for (i = 0; i < n - 1; i++)
	{
		b[i] = (double *)malloc((n - 1) * sizeof(double));
	}
	if (n == 1)
	{
		m = a[0][0];
	}
	else
	{
		j = 0;
		for (i = 0; i < n; i++)
		{
			f = 1;
			for (k = 0; k < i + j; k++)
			{
				f = f * (-1);
			}
			for (ii = 0; ii < n; ii++)
			{
				for (jj = 0; jj < n; jj++)
				{
					if (ii < i)
					{
						a_i = ii;
					}
					if (ii > i)
					{
						a_i = ii - 1;
					}
					if (jj < j)
					{
						a_j = jj;
					}
					if (jj > j)
					{
						a_j = jj - 1;
					}
					if (ii != i && jj != j)
					{
						b[a_i][a_j] = a[ii][jj];
					}
				}
			}
			m = m + a[i][j] * det(b, n - 1)*f;
		}
	}
	for (i = 0; i < n - 1; i++)
	{
		free(b[i]);
	}
	free(b);
	return m;
}
void inv(int n, double **a, double **b)
{
	int i, j, k;
	int ii, jj;
	int a_i, a_j;
	a_i = -1;
	a_j = -1;
	double det_a;
	double **c;
	int f;

	c = (double **)malloc((n - 1) * sizeof(double*));
	for (i = 0; i < n - 1; i++)
	{
		c[i] = (double *)malloc((n - 1) * sizeof(double));
	}

	det_a = det(a, n);
	if (fabs(det_a) < 0.000000001)
	{
		for (i = 0; i < n; i++)
		{
			for (j = 0; j < n; j++)
			{
				b[i][j] = 0;
			}
		}
	}
	else
	{
		for (i = 0; i < n; i++)
		{
			for (j = 0; j < n; j++)
			{
				f = 1;
				for (k = 0; k < i + j; k++)
				{
					f = f * (-1);
				}
				for (ii = 0; ii < n; ii++)
				{
					for (jj = 0; jj < n; jj++)
					{
						if (ii < i)
						{
							a_i = ii;
						}
						if (ii > i)
						{
							a_i = ii - 1;
						}
						if (jj < j)
						{
							a_j = jj;
						}
						if (jj > j)
						{
							a_j = jj - 1;
						}
						if (ii != i && jj != j)
						{
							c[a_i][a_j] = a[ii][jj];
						}
					}
				}
				b[j][i] = det(c, n - 1) / det_a * f;
			}
		}
	}
	for (i = 0; i < n - 1; i++)
	{
		free(c[i]);
	}
	free(c);
	return;
}

void trans(int ax, int ay, double **a, double **b)
{
	int i, j;
	for (i = 0; i < ay; i++)
	{
		for (j = 0; j < ax; j++)
		{
			b[j][i] = a[i][j];
		}
	}
	return;
}

void cal_rp(double **d, double x, double y, double z, double *rp)
{
	rp[0] = x * d[0][0] + y * d[1][0] + z * d[2][0];
	rp[1] = x * d[0][1] + y * d[1][1] + z * d[2][1];
	rp[2] = x * d[0][2] + y * d[1][2] + z * d[2][2];
	return;
}

double dot_product(double ***a, double ***b, int nx, int ny, int nz)
{
	int i, j, k;
	double sum = 0;
	for (i = 0; i < nx; i++)
	{
		for (j = 0; j < ny; j++)
		{
			for (k = 0; k < nz; k++)
			{
				sum = sum + a[i][j][k] * b[i][j][k];
			}
		}
	}
	return sum;
}
double dot_product(double ***a_r, double ***a_i, double ***b_r, double ***b_i, int nx, int ny, int nz)
{
	int i, j, k;
	double sum = 0;
	for (i = 0; i < nx; i++)
	{
		for (j = 0; j < ny; j++)
		{
			for (k = 0; k < nz; k++)
			{
				sum = sum + a_r[i][j][k] * b_r[i][j][k] + a_i[i][j][k] * b_i[i][j][k];
			}
		}
	}
	return sum;
}

//�������ֵ
double get_dot_product(double*a, double*b, int n) {
	double sum = 0;
	for (int i = 0; i < n; i++) {
		sum += a[i] * b[i];
	}
	return sum;
}

//�������ֵ
double* get_cross_product(double*a, double*b, int n) {
	double*res = new double[n];
	if (n != 3) {
		cout << "unsupported n value!" << endl;
		cin.get();
		return NULL;
	}
	res[0] = a[1] * b[2] - a[2] * b[1];
	res[1] = a[2] * b[0] - a[0] * b[2];
	res[2] = a[0] * b[1] - a[1] * b[0];
	return res;
}
void plus_etot_error(string path, string etot_file, string report_file) {
	if (IF_ETOT_PLUS_EROR == true) {
		chdir(path.c_str());
		FILE* in;
		string com = "grep 'ERROR' " + etot_file + "|wc -l";
		in = popen(com.c_str(), "r");
		char num[20];
		fgets(num, sizeof(num) - 1, in);
		pclose(in);
		if (num[0] > '0')
		{
			//cout << path << "has have stress!" << endl;
			//cin.get();
			return;
		}
		//���濪ʼ��ʽд��
		char errror_save[150][4];
		double error[4] = { 0 };
		ifstream fin;
		fin.open(path + "/" + report_file, ios::in);
		if (!fin.is_open()) {
			cout << " i can not find the file:" << path + report_file << endl;
			cin.get();
		}
		char buff[150];
		//ע�����������⣬�����Ǹ���������Ҫ�����ҽ����ж�
		while (fin.good() && fin.peek() != EOF) {
			fin.getline(buff, 149);
			if (strstr(buff, "WG_ERROR") != NULL)
				break;
			memset(buff, 0, 150);
		}

		int index = 0;
		for (int j = 0; j < 150; j++) {
			if (buff[j] >= '0' && buff[j] <= '9')
			{
				index = j;
				break;
			}
		}
		error[0] = atof(buff + index);
		error[0] *= 1e-4;

		//�����ǽ�string��ȡ��
		for (int i = 1; i < 4; i++) {
			fin.getline(errror_save[i], 150);
			//cout << errror_save[i] << endl;
			//cin.get();
			int index = 0;
			for (int j = 0; j < 150; j++) {
				if (errror_save[i][j] >= '0' && errror_save[i][j] <= '9')
				{
					index = j;
					break;
				}
			}
			error[i] = atof(errror_save[i] + index);
			error[i] *= 1e-4;
			//cout << error[i];
		}
		fin.close();
		//Ȼ����etot���濪ʼд�µĶ���
		fstream fplus;
		fplus.open(path + "/" + etot_file, ios::app);
		if (!fplus.is_open()) {
			cout << " i can not find the file:" << path + "/" + etot_file << endl;
			cin.get();
		}
		fplus << "    WG_ERROR  =   " << std::scientific << error[0] << endl;
		fplus << "    E_ERROR   =   " << std::scientific << error[1] << endl;
		fplus << "    RHO_ERROR =   " << std::scientific << error[2] << endl;
		fplus << "    RHO_RELATIVE_ERROR      =  " << std::scientific << error[3] << endl;
		fplus.close();

	}
	return;
}

double get_length(double* a) {
	return pow(pow(a[0], 2) + pow(a[1], 2) + pow(a[2], 2), 0.5);
}

void plus_stress_mask(string path)
{
	if (IF_STRESS == true)
	{
		chdir(path.c_str());
		FILE* in;
		in = popen("grep 'STRESS' atom.config|wc -l", "r");
		char num[20];
		fgets(num, sizeof(num) - 1, in);
		pclose(in);
		if (num[0] > '0')
		{
			//cout << path << "has have stress!" << endl;
			//cin.get();
			return;
		}

		fstream fin;
		fin.open(path, ios::app);
		if (!fin.is_open())
		{
			cout << "i cann not find the file" << path << endl;
			cin.get();
		}
		if (CALCULATION_MODE == 1)
		{
			fin << "STRESS_MASK" << endl;
			fin << "  1 0 0" << endl;
			fin << "  0 0 0" << endl;
			fin << "  0 0 0" << endl;
			fin.close();
		}
		else if (CALCULATION_MODE == 2)
		{
			fin << "STRESS_MASK" << endl;
			fin << "  1 0 0" << endl;
			fin << "  0 1 0" << endl;
			fin << "  0 0 0" << endl;
			fin.close();
		}
		else
		{
			cout << "warging ! stress mask do not work!" << endl;

		}
		return;
	}
}

bool check_relax(string path)
{
	FILE* fi;
	string com = "cat RELAXSTEPS|wc -l";
	char num[10];
	chdir(path.c_str());
	fi = popen(com.c_str(), "r");
	fgets(num, 9, fi);
	int line = atoi(num) + 1;
	pclose(fi);
	return line < RE_STEP;
}

bool check_now_exist_file(string name)
{
	//��鵱ǰĿ¼���Ƿ����ĳ���ļ�
	string com;
	com = "ls " + name;
	FILE* fi;
	fi = popen(com.c_str(), "r");
	char buff[20];
	fgets(buff, 19, fi);
	pclose(fi);
	//cout << buff << endl;
	//cin.get();
	if (buff[0] == name[0])
	{
		return true;
	}
	else return false;
}

vector<string> get_file_name(string file_name)
{
	ifstream fin;
	fin.open(file_name, ios::in);
	vector<string> name;
	if (!fin.is_open())
	{
		cout << "i can not find the file:" << file_name << endl;
		cin.get();
	}
	char buff[100];
	while (fin.good() && fin.peek() != EOF)
	{
		fin.getline(buff, 100);
		name.push_back(string(buff));
		memset(buff, 0, 100);
	}
	fin.close();
	//system(("rm " + file_name).c_str());
	return name;

}

//template <typename T>
void vector2file(vector<vector<double>>&data, ofstream &fout) {
	//ʹ�÷��ͽ���������ó���
	if (data.size() == 0)
		return;
	fout.precision(9);
	for (int i = 0; i < data.size(); i++) {
		for (int j = 0; j < data[0].size(); j++) {
			 fout<< data[i][j] << "\t";
		}
		fout << endl;
	}
	return;
}
void plus_init_magnetic(string path,int*type_list,int type_num,int value) {
	//���ӳ�ʼ�ž�,��Ҫ�Ĳ����Ǽ���Ķ�Ӧ·��
	if (IF_PULS_MARNETIC == true) {
		chdir(path.c_str());
		FILE* in;
		in = popen("grep 'MAGNETIC' atom.config|wc -l", "r");
		char num[20];
		fgets(num, sizeof(num) - 1, in);
		pclose(in);
		if (num[0] > '0')
		{
			//cout << path << "has have stress!" << endl;
			//cin.get();
			return;
		}
		fstream fin;
		fin.open("atom.config", ios::app);
		if (!fin.is_open())
		{
			cout << "i cann not find the file" << path << endl;
			cin.get();
		}
		//���濪ʼ��ʽд��žأ�
		for (int i = 0; i < type_num; i++) {
			fin << type_list[i] << "  " << value << endl;
		}
		fin.close();
	}
	return;
}

vector<string> cut_string(string word) {
	//���ڴ�ŷָ����ַ��� 
	vector<string> res;
	string result;
	//���ַ�������input�� 
	stringstream input(word);
	//���������result�У�������res�� 
	while (input >> result)
		res.push_back(result);
	return res;

}

void plusRelaxQdiv(string path, string file_name) {
	//��relax��qdiv�����Ϊ��ʼ��scf��magnetic
	if (IF_ACCERATION_SCF) {
		chdir(path.c_str());
		system("awk  '{print $1\"\t\"$2}' OUT.QDIV_relax >tmp");
		system("echo \"MAGNETIC\" >>atom.config");
		system("sed -i \"1d\" tmp");
		system("cat tmp >> atom.config");
		system("rm tmp");
		cout << "read relax step Magnetic and output to scf finished!" << endl;
		return;
	}
	
}