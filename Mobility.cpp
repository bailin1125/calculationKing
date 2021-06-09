#include "Mobility.h"
#include "calculation_tools.h"
#include "others.h"
#include <cstdlib>
#include <math.h>
#include <unistd.h>
#include<iomanip>
#include <map>
#pragma warning(disable : 4996)

using namespace std;

//mobility的整天调用
void mobility_client(string config_name,string calculation_path) {
	int index = 0;
	string config_path = "/share2/wangz/0721_nanotube/config/";
	string input_file = "mobil_input";
	vector<string> name = get_file_name(config_name);
	//依次为C，两个有效质量，两个E，以及两个迁移率
	vector<vector<double>> total_result(name.size());
	for (int i = 0; i < name.size(); i++)
	{
		total_result[i].resize(7);
	}

	//for (auto x : name) {
	//	cout << "start to " + string(x) << endl;
	//	chdir(calculation_path.c_str());
	//	//prepare_mobility(x, calculation_path, config_path + x);
	//	get_mobility_data(calculation_path + x);
	//	//cin.get();
	//	index++;
	//}
	//cout << "mobility module finished!" << endl;

	//读取有效质量和其他三个数值的文件，输出相关的结果
	ifstream fin;
	fin.open(input_file,ios::in);
	int i = 0;
	string temp;
	char buff[200];
	while (fin.good() && fin.peek() != EOF) {
		fin >> temp;
		cout << temp << endl;
		for (int j = 0; j < 5; j++) {
			fin >> total_result[i][j];
		}
		fin.getline(buff,199);
		i++;
	}
	fin.close();
	for (int j = 0; j < i; j++) {
		total_result[j][5] = get_mobility(total_result[j][0], total_result[j][1], total_result[j][3]);
		total_result[j][6] = get_mobility(total_result[j][0], total_result[j][2], total_result[j][4]);
	}
	//最后输出最后两个结果
	ofstream fout;
	fout.open("mobil_total_res", ios::out);
	for (int i = 0; i < name.size(); i++) {
		fout << name[i] << "  " << total_result[i][5] << "  " << total_result[i][6] << endl;
	}
	fout.close();
	cout << "total work done!" << endl;
	cin.get();
	return;
}


string conventer(string input) {
	int index = input.find_first_of("g")+2;
	return input.substr(index);
}


//在制定路径下准备好迁移率的文件
bool prepare_mobility(string file_name, string path,string config){
	//1 确定拉伸的系数
	vector<double> factor = get_factor(-0.01, 0.0025, 9);
	cell cell_org(config.c_str());
	
	//2 创造一系列的输入文件	
	string temp_path;	
	string com = "mkdir " + file_name;	
	system(com.c_str());	
	for (double fac : factor) {
		cout << "start to " +file_name+ to_string(fac)<< endl;
		chdir((path + file_name).c_str());
		temp_path= path + file_name + "/"+to_string(fac);
		com = "mkdir " + temp_path;
		system(com.c_str());
		chdir(temp_path.c_str());
		//或者开始计算relax,后来我们发现根本不用重新relax
		//prepare_mobility_singel_relax(fac, cell_org, file_name);

		//或者最后计算scf
		com = "cp " + config + " ./atom.config";
		system(com.c_str());
		prepare_mobility_singel_scf(fac, file_name, cell_org,temp_path);
	}
	return true;

}

//在当前路径下准备好单个的输入,准备好realx的
void prepare_mobility_singel_relax(double factor,cell cell_a,string name) {

	//1 进行晶格常数拉伸
	for (int i = 0; i < 3; i++) {
		cell_a.letice[direction][i] *= (1 + factor);
	}
	// 2 准备好其他文件
	if (sole_generate_jobs("atom.config", "./", "relax") != 0) {
		cout << "relax encounter error!" << endl;
		cin.get();
	}
	
	// 3 提交上文件
	pbs_got(name, "relax", "./", "gpu2");
	system("qsub aaa.pbs");
	return;
}

//在当前路径下准备好单个的输入,scf
void prepare_mobility_singel_scf(double factor,string name,cell & cell_a,string path) {
	cout << "prepare " + name + " scf calculation!" << endl;
	//1 进行晶格常数拉伸
	for (int i = 0; i < 3; i++) {
		cell_a.letice[direction][i] *= (1 + factor);
	}
	cell_a.output();
	//再返回去
	for (int i = 0; i < 3; i++) {
		cell_a.letice[direction][i] /= (1 + factor);
	}
	// 2 准备好其他文件
	if (sole_generate_jobs("atom.config", "./", "scf") != 0) {
		cout << "scf encounter error!" << endl;
		cin.get();
	}

	//准备好延时文件
	//复制堰势文件到指定目录下
	for (int i = 0; i < cell_a.type_num; i++)
	{
		string com;
		com = "cp /share/home/wangz/NCPP-SG15-PBE/" + string(atom_name[cell_a.type_save[i]]) + ".SG15.PBE.UPF" + " ./";
		//cout << com << endl;
		system(com.c_str());

	}

	// 3 提交上文件
	/*if (run_num_mobi< 12) {
		pbs_got(name, "scf", path, "gpu3");
		run_num_mobi++;
	}
	else
	{*/
		pbs_got(name, "scf", path, "gpu2");
	/*}*/
	
	system("qsub aaa.pbs");
	return;
}


//准备好连续的应变系数
vector<double> get_factor(double start, double increase,int num) {
	vector<double> res;
	for (int i = 0; i < num; i++) {
		res.push_back(start + i * increase);
	}
	return res;
}


//将相关结果打印到当前路径下
void get_mobility_data(string path) {
	chdir(path.c_str());
	string com;
	com = "ls >name";
	system(com.c_str());
	vector<string> name = get_file_name("name");
	
	//建立结果的map,这个是总能量的
	map<string, double> result;
	string convert;
	for (string x : name) {
		if (x == "name")
			continue;
		convert = conventer(x);
		double temp_energy=cell::check_energy(path +"/"+ x);
		result.insert(pair<string,double>(convert,temp_energy));
	}

	//最后输出到相关文件里面去
	chdir(path.c_str());
	ofstream fout;
	fout.open(result_name, ios::app);
	if (!fout.is_open()) {
		cout << "i can not find the name:" << result_name << endl;
		cin.get();
	}
	map<string, double> ::iterator it = result.begin();
	fout.setf(ios::fixed);
	for (it = result.begin(); it != result.end(); it++) {
		fout <<it->first << "    " <<setprecision(9)<< it->second << endl;
	}
	fout.close();


	//然后输出vbm和cbm的能量，这个是针对vbm和cbm的	
	map<string, double> vbm;
	map<string, double> cbm;
	for (string x : name) {
		if (x == "name"||strstr(x.c_str(), "monil")!=NULL)
			continue;
		report rp;
		rp.read_report(const_cast<char*>((path + "/" + x+"/REPORT").c_str()),1);
		convert = conventer(x);
		vbm.insert(pair<string, double>(convert, rp.vbm_energy));
		cbm.insert(pair<string, double>(convert, rp.cbm_energy));
	}
	//分别进行输出
	fout.open(result_name + "_vb", ios::out);
	it = vbm.begin();
	while (it != vbm.end()) {
		fout << it->first << "    " << setprecision(9) << it->second << endl;
		it++;
	}
	fout.close();

	fout.open(result_name + "_cb", ios::out);
	it = cbm.begin();
	while (it != cbm.end()) {
		fout << it->first << "    " << setprecision(9) << it->second << endl;
		it++;
	}
	fout.close();

	return;

}

double mobility_const(void) {	
	double fenzi = e_elec*kulun* plunk_ba*plunk_ba;
	double fenmu = pow(2 * pi*kb*T_temp, 0.5);
	return fenzi / fenmu;
}

//获得迁移率,注意这里返回的结果是*e3 cm2/V S
double get_mobility(const double C, const double m, const double E) {
	double const_number = mobility_const();
	double fenmu = pow(m*m0, 1.5)*pow(E*e_elec, 2);
	return const_number * C*1e13 / fenmu*10;
}