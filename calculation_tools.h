#pragma once

#ifndef CALCULATION_TOO
#define CALCULATION_TOO
#include <string>
#include <string.h>
#include <istream>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <sstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <map>
#include "stack"
#include<numeric>
#include "config.h"

const int Meatal_xuhao[]= { 3, 4, 11, 12, 13, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 37, 38, 39, 40, 41, 42, 43, 44, 45,46, 47, 48, 49, 50,  51,55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83,84, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111 };
const int Agroup_matal[]= {3,4,11,12,19,20,37,38,55,56,87,88,13,31,32,49,50,51,81,82,83,84};
const int main_group_element[16] = { 1, 5,6,7,8,9,14,15,16,17,33,34,35,52,53,84 };//at为金属元素
const int rare_gas[6] = { 2,10,18,36,54,86 };
const int metal_num = 89;
const int main_groupnum = 16;
const int rare_gasnum = 6;
static char atom_name[120][3] = { " ","H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg" };

using namespace std;


void singel_test();//单元测试的code
int  get_information(string config_name, string data_path);//统计相关信息




//有关read_report的部分
struct E_E {
	int spin;
	int kpt;
	int band;
	double energy;
};


/**
	有关读取REPORT，获得带隙有效质量等信息的基类
**/
class report {
public:	
	int  read_report(char *name_report,int spin);
	int n_x, n_y, n_z;
	double orig_fermi = 0;//根据能带获得的费米能级
	int mp_x, mp_y, mp_z;
	struct E_E *eig;
	double ***eig_kns;
	int num_band;
	int num_kpt;
	double **kp_p;//这个是倒空间的，写在report里面的k点坐标
	vector<vector<double>>base_kpp;//这个是分数坐标的k点坐标
	int spin;
	int num_electron;
	string soucrce;//数据的来源
	//在基础上进行扩充
	string report_name;
	double diff = 0;
	double energy = 0;
	


	int cbm_band = 0;
	int vbm_band = 0;
	//补充有关vbm和cbm的获取
	double* vbm_kppt;
	double*cbm_kppt;
	double vbm_energy;
	double cbm_energy;
	int vbm_kpt;
	int cbm_kpt;
	//获得倒空间的精确度
	static double get_accuracy(string file_name = "REPORT");
	//将vasp的in。kpt转化为pwmat的in。kpt
	static void vaspKpt2pwmat(string file_name, string out_path,int cut_num=19);
	//上面的批量化
	static void allvaspKpt2pwmat(string filename, string inpath, string outpath);
	//获得有效质量
	void get_valid_m(int mode);//输入是计算的模式
	double get_gap(string path,vector<string>& step);//输入是有结构的路径
	double get_kpoint_bandgap(int kpt);//输出指定k点对应的带隙大小

	//判断这个结构的指定的方向下的带隙情况
	int check_direction_gap(vector<string>flag,string file_name,string cell_name,string file_name2="OUT.KPT");
	static void check_direction_gap_multi(vector<string>flag, string name_file, string path, string kpt_name = "OUT.KPT");
	//输出report两种spin的gap情况
	vector<double>get_half_situation(string report_name);
	double get_singel_number(string flag);
	~report() {
		free(eig);
		for (int i = 0; i < num_kpt; i++) {
			for (int j = 0; j < num_band; j++) {
				free(eig_kns[i][j]);
			}
			free(eig_kns[i]);
			free(kp_p[i]);
		}
		free(eig_kns);
		free(kp_p);
	}
private:	
	vector<vector<double>>kp_fractional;
	void report_k2_fractional(string file_name);//读取晶胞结构文件，将report的真实坐标还原为分数坐标
	vector<double> get_diff(string file_name);//获得diff的相关信息，就三个数值
	vector<int> base_kpt_flag;//标记这个k点的所属方向，0是原点，1是x，2是y，3是z，4是其他
	int judge_base_kpt_singel(vector<double>&input);//输入k点三个坐标，输出这个k点所属方向是什么
	vector<double>get_allk_energy(int spin, int band);//获取指定spin和band的所有的k点的列表
	//有关能带提取相关方法	
	vector<vector<double>> get_band_value_onkpt(vector<int>&klist,string file_name);//获取指定k点列表对应的所有的本征能量值
	vector<int>get_all_kptlist(string flag);//获取指定标记下对应k点的序号
	bool if_band_cross(vector<double>&eigen);//判断指定的一系列的本征能量值是否穿过了能带
	bool read_base_kpt(string file_name);//读取指定路径下面的分数坐标的k点
	void inrich_source_information(string config_name);//补充数据来源的信息
};
int Read_N123(char *name, const char *tag, int *n1, int *n2, int *n3);
int Read_double(char *name, const char *tag, double *out);
int Read_int(char *name, const char *tag, int *out);
int Read_kpt_num(char *name, const char *tag, int *out);
int Read_eig(char *name, struct E_E *eig, int num_band, int num_kpt, int spin, double ***eig_kns);
int Read_kpp(char *name, double **kp_p, int num_kpt);


/*从report继承而来，用来*/
class Bandstru :public report {	
public:
	//static bool Fermi_in_orig;//标记费米能级的获取，true为自己计算获得
	static int generate_all_bandstu(string name, string path);
	int k_num = 0;
	int band_num = 0;
	double fermi = 0;
	int orig_start_band = 0;//选定的画band的基准band数
	vector<double>k_step;
	vector<vector<double>> band_eig;
	/*构造函数，读取band文件*/
	Bandstru(string file_name,int flag=1);
	//核心方法，输出可以直接画图的band信息
	int drawBandstru(int start, int expand);
	//补充画dos的结果
	int drawDos(string path,string out_name);

	static void output_band(string out_name,vector<Bandstru*> bandlist);//输出band文件
	int real_start = 0;//设定画的band的边界值
	int real_end = 0;
private:
	void init_start_band_num();//设定画的band的基准数值
	
};


//有关元素信息的储存类
class element
{
public:
	int atomic_num; //atomic number
	double vdw_radius_min;
	double vdw_radius_max;
	//van der waals radius,min and max
	double cov_radius; //covalence radius
	int num_metal_radius;
	double *metal_radius;

	double posi_ionic_ridus;
	double nega_ionic_ridus;
};

//关于atom.config部分
class atom {
private:
	//char atom_name[120][3] = { " ","H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg" };
public:
	int num_atom;
	double **d;
	double **rd;
	double **p;
	double **rp;
	int *atom_type;
	void read_atom(char *name);
	void print_name(int i, FILE *in);
	int type_num;
	int type_name[30];
	int flag_read;
	void freeatom();
};


//有关结构及相关输出的工具类
class cell
{
	static string cif_path;
	const static double SuperO_rule;
public:
	bool Valid = true;
	cell(const char *jiegou_name, int flag=0);
    ~cell();
	int num = 0;
	double **letice;
	double **p;
	double ***p_real;
	double ***real_position;
	double **ridus;
	int *type;
	int type_num;
	int *type_save;
	vector<vector<int>> type_count;
	
	//检查输出文件中部分的结果
	double check_diff(string path, string file_name = "REPORT");
	static double check_energy(string path, string file_name = "REPORT");
	bool check_contain_transition_metal();
	int check_main_transition_situation();//确定结构含有元素的情况，1是只有主族，2是只过渡金属，3是二者都有，4是其他（错误
	bool check_contain_water();
	bool check_contain_supoerO();//判断结构是否含有超氧结构

	//用来变形结构或者创造新的结构
	
	void supper_cell(int x, int y, int z,string path);
	void supper_cell_md(string path);//读取路径下面的结构，为md进行扩胞
	void plus_vacum(vector<double>& vacum, string path,string file_name);
	void output_config(string path,string file_name="atom.config");
	void trans_form(string path,string file_name,vector<vector<int>>&input);

	
	//获得结构中某些性质
	string  get_formula_name();
	void get_center_enviro();
	void get_center_valence(string file_name = "OUT.QDIV_scf");
	void more_k_toget_band(string path);
	vector<double> get_onemore_point(vector<double> a, vector<double>b, double step);
	void output(string file_name="atom.config");
	double get_positive_distance();
	string get_cif_spacegroup(string name);//获得体材料的空间群
	static string get_cif_name(string name);//输入低维结构的名字，输出体材料的名字


	//其他辅助性的私有变量
	int center_num = 0;
	int envoroment_num = 0;
	int* center = NULL;
	int* envoroment = NULL;
	double *center_valence = NULL;
	int *my_classify;
	
	
private:
	int *positive;	
	map<int, int>type_and_num;
	
	
};

/*
	从cell继承而来主要是用来获得脚本中相关能量之类的
*/
class Cal_cell :public cell {
	//每个是2D，对应绑定有bulk的相关信息
public:
	bool Valid = true;
	string name = "";
	double total_energy=0;
	//这个是普通的per/atom的结合能
	double formation_energy_peratom=0;
	//正常的per square的结合能，2D是除以面积，对于0D或者1D是使用等效面积法
	double formation_energy_square = 0;

	double suqre = 0;//优化后的2D面积，在其他维情况时，对应的是等效面积
	double org_square = 0;//原始bulk对应的2D面积
	double ridus = 0;//在0D和1D情况下的等效半径
	double h_length = 0;//1D情况下那个h

	Cal_cell(string name, char*jiegou, int flag) :cell(jiegou, flag), name(name) {

	};
	void get_org_cell(string path);
	void get_fromation_energy_peratom();
	void get_formation_energy_persquare();
	Cal_cell* bulk;
	cell * cell_org;
	vector<double>  get_all_magnetic();//获得该结构下所有自旋的情况可能
	int static get_all_magnetic_for_all(string name, string path);//批量化的产生所有电子的所有情况


	//输出所有感兴趣的信息
	void output_inforamtion(ofstream& fout);//打印相关信息

};


//从cell继承的专门用于扩展模板化结构的
class Pattern_cell :public cell{	
public:
	string pattern_name = "";
	Pattern_cell(string name,char* jiegou,int flag):cell(jiegou,flag),pattern_name(name){
		//初始化列表显式初始化子类		
	};	
	bool get_pattern(string& name, vector<int>element,cell &cell_pattern,double factor=1,string flag="1T");//根据输入的结构进行结构化产生新结构

};
void generate_pattern_configs(string filename,string path,string flag,element* e);



//从cell继承过来用来进行waniner计算的类
class Wannier_cell :public cell {	
public:
	static string pbe_path;
	string source_name = "";
	string metal_name = "";
	string main_group_name = "";
	string real_name = "";
	map<string, double> lambda;
	vector<double> valid_cb;//输出有效质量vb，格式，格点x距离差，三个能量；格点y距离差，三个能量
	vector<double>valid_vb;//上同
	vector<double> ordered_lambda_cb;
	vector<double> ordered_lambda_vb;
	bool lmbda_valid = true;
	Wannier_cell(const char*jiegou, int flag) :cell(jiegou, flag) {
		//获得金属和非金属的名字
		for (int i = 0; i < this->num; i++) {
			if (this->my_classify[i] == 1) {
				metal_name = atom_name[this->type[i]];
			}
			else if (this->my_classify[i] == 2) {
				main_group_name = atom_name[this->type[i]];
			}
		}
	};
	static void get_source_dependent();
	static int find_close_kppt(double* kpoint, report& rp);
	static int find_close_band(double energy, int kpoint,report&rp);
	void get_source_name();
	void get_wain90(string path,string source_path);//创造wan90文件
	static int get_start_line(string flag = "begin unit_cell_cart");
	static int get_start_all_line();
	void get_wannier_param();//获得in.wannier.param文件
	void run_wannier(string now_path, string source);
	double get_singel_ilmbda(string path);//产生单个的ilmbda
	void get_ilmbda(string path);//产生所有的ilmbda

	void output_ilmbda(string file_name);
	
	vector<int> get_left_and_right(int start, int end, int num, int band, string flag);

	//添加运行wkm的函数
	void run_wkm(string now_path, string source);
	void prepare_singel_file(string dir_path,string base_path,string source_path,string flag,string x);//为一个文件夹准备好所有的输入文件

	void run_last_scf(string now_path, string source); //进行最后的计算的函数
	void get_valid_m(string now_path,string name);//老方法获取wkm的有效质量
	void get_valid_m_new(string now_path, string name);//新的方法获取有效质量
private:
	double dict2tofind_double(map<string, double> mat, string key);//从字典中获得数值
	

};
void generate_winf90(string source_path, string dirpath, string name);



//基于Lobister的相关计算类
class Lobister {
public:
	bool sucessful = true;//标志这个计算成功的标志
	Lobister(string name,string path) {
		this->name = name;
		this->path = path;
	}
	string name = "";
	string path = "";//每个结构计算的路径
	string comm;//关于每个cohp的相关说明
	int time = 0;//ridus的个数
	void getCohp();//获得每个结构的cohp的数值
	void runCohp();//执行每个结构下面的cohp的计算
	void getRidus();//获得每个结构的键长信息统计
	void getElectronNum();//获取结构的外层电子数目的统计
	void getEferAndStep();//获取费米能级以及步长信息
	void getStartRange();//获取初始情况下的能量下限数值
	//对每个计算执行完整计算流程
	bool runProcess() {
		bool res = false;
		getElementList();
		getStartRange();
		//runCohp();
		getCohp();		
		getRidus();
		res = true;
		return res;
	}
	string static sample_path;//样本输入文件的路径
	static string COHP_GET;//获得cohp的方法，是简单读取还是自己分能级累加
	static void runAllProcesses(string file_name, string path, string out_name);//执行全部的计算流程并且输出相关文件
	vector<pair<string, int>> element_info;
	static string getLineString(int line, string file_name);//获取指定文件的指定第几行
	
private:
	vector<double> cohp;//对于每个结构的,此时每个结构均可能有不同的数值了
	vector<double> ridus;//每个结构的键长信息的统计
	vector<double> min_cohp;
	vector<double> max_cohp;
	string getBaseFunctions(string name,bool org_flag);//提供不同的手段，可以从当前文件夹读取或者之类的
	void  getElementList();//获得元素列表的信息
	vector<vector<string>> ridusElment;//进行运算的相应的元素对
	int getNoForElementPair(vector<string>&ele_info);//对于指定的元素连接，得到是NO.第几个的对应数据

	//0605讨论之后确定的方案
	double ferimi = 0;//费米能级
	double step = 0;//步长
	double electron_num = 0;//价电子数，这个需要读取元素信息然后累加
	vector<int> interactions;//获取读取cohp的N
	double startRange = 0;//进行分能级计算之后的能量下限

	vector<double> up_differ ;//自旋向上的cohp的插值
	vector<double> down_differ;//自旋向下的cohp的插值

};







int classify_metal_maingroup(int &atomic_number);//0 rare gas,1 metal ,2 main group
bool map_com(pair<int, int>a, pair<int, int>b);

const double ridus_plus_factor = 1.2;
const double val_radius_factor = 1.2;
double dis(double *p1, double *p2);
void read_element(element *e, string &file_element_r,string & ioinci);

void detect_spin_eualone(string & path, string & file_name, string & out_name, int **flag,double rule= SPIN_check_rule);
int detect_free(const string &choose);//查看有多少空闲节点
int generate_config(string name="BiI3");
int if_finish(string path,int flag=0);
int generate_jos(string name, string data_path,string dir_path, int num, string flag,int* yanshi_flag,int **task_flag,ofstream & fout);
int generate_jos_temp(string name, string data_path, string dir_path, int num, string flag, int* yanshi_flag, int **task_flag, ofstream & fout);
void pbs_got(string name, string flag, string dir_path,string lines);
int if_no_yanshi(string path, int *yanshi_flag);

int  sole_generate_jobs(string name, string dir_path,  string flag,report* rp=NULL);//用来产生单个的脚本文件
int make_efm(report& rep,string path, int KPT_X, int KPT_Y, int KPT_Z);//输入是迭代变量和路径，输出是新的IN.KPT


void give_mp123(string name,int *mp);//输入是结构名字，输出是填充mp
int check_now_erfen_pk(string path,int last_flag,int kpt_num=1);
//整理相关结果的函数，重载版本，前者是用于生成相关文件，后者用于快速打印结果
void get_last_result(string path,double& gap,double &ey,double& ex,double& vy,double& vz,int kpt_num,double** inter,vector<string>&step,int flag=1);
void get_last_result(string path, vector<string>&step, vector<vector<double>>&table,int flag);

void generate_last_file(string outname, int** flag, int total, string name);
int  backup_status(string config_name,string &name1, string& name2, int**flag,int num,int** mp,int linenumber);
void full_fill_yanshi(string yanshi_name, string data_path);
void get_all_output(string &out_name, string& chazhi_name, string & calculation, string *config, int**flag, int num);





vector<int> calculate_n123(double a, double b, double c, int ECUT2);

double get_vecror_mo(double*a, int n);

void fullfill_inter_real_location(double**inter,int kp,cell &cell_a,string cell_path,int fullfille_idex,
	int band,int spin,string flag);

void puls_magnetic(string path,string flag="F");//添加初始磁矩
void bothTestAfsituation(string path,string name);//测试反铁磁的情况,输入是计算路径以及该结构的名字

void test(element *e);//其他的相关测试

#endif CALCULATION_TOO