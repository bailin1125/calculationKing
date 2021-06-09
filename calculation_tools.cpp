
#include "calculation_tools.h"
#include "others.h"
#include "config.h"
#include <cstdlib>
#include <math.h>
#include <unistd.h>
#include<iomanip>
#include <algorithm>
#pragma warning(disable : 4996)

using namespace std;
char a[112][3] = { " ","H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lw", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg" };


void singel_test() {
	//单元测试的code
	report rep;
	string test_path = "/share2/wangz/0527_got1T/test/calculation/GaAs.config/REPORT";
	double near_location[3];
	int cbm_kpt_xuhao[] = { 0,1,5,36,180,6,30 };
	int vbm_kpt_xuhao[] = { 0,1,5,36,180,6,30 };
	//打印一下需要的东西
	rep.read_report(const_cast<char*>(test_path.c_str()), 1);
	//CBM附近的
	//能量
	cout << "start cbm!" << endl;
	cout << rep.cbm_kpt << endl;
	cout << rep.eig_kns[rep.cbm_kpt][rep.cbm_band - 1][0] / HARTREE << endl;
	for (int x : cbm_kpt_xuhao) {
		cout << "kpt: " << x << endl;
		cout << "energy:" << rep.eig_kns[x][rep.cbm_band - 1][0] / HARTREE << endl;
		cout << "location: " << rep.kp_p[x][0] << "," << rep.kp_p[x][1] << "," << rep.kp_p[x][2] << endl;
	}


	//坐标
	cout << rep.cbm_kppt[0] << "," << rep.cbm_kppt[1] << "," << rep.cbm_kppt[2] << endl;
	//cout << dis(near_location, rep.cbm_kppt);


	//VBM附近的
	//能量
	cout << "start vbm!" << endl;
	cout << rep.vbm_kpt << endl;
	cout << rep.eig_kns[rep.vbm_kpt][rep.vbm_band - 1][0] / HARTREE << endl;
	for (int x : vbm_kpt_xuhao) {
		cout << "kpt: " << x << endl;
		cout << "energy:" << rep.eig_kns[x][rep.vbm_band - 1][0] / HARTREE << endl;
		cout << "location: " << rep.kp_p[x][0] << "," << rep.kp_p[x][1] << "," << rep.kp_p[x][2] << endl;
	}
	//坐标
	cout << rep.vbm_kppt[0] << "," << rep.vbm_kppt[1] << "," << rep.vbm_kppt[2] << endl;
	//cout << dis(near_location, rep.vbm_kppt);
	cin.get();

}

int  get_information(string config_name, string data_path) {
	string unsual = "check_time";
	vector<string> unsual_name = get_file_name(config_name);
	/*vector<vector<double>> data_save(unsual_name.size());
	for (int i = 0; i < data_save.size(); i++) {
		data_save[i].resize(4);
	}*/

	ofstream tt;
	tt.open("formula_stat", ios::out);
	tt << "name" << " " << "formula" << " " << "center_element" << " " << "enviro_element" << " " << "average_center_val";
	tt << " " << "space_group";
	tt << endl;

	/*tt << "name" << " " << "formula" << " " << "if_contain_transmetal" << " " << "if_contain_water" << " "
		<<"energy" << " " << "diff" << " " << endl;*/
	ostringstream room;
	for (string x : unsual_name)
	{
		cout << "check " << x << endl;
		cell* cell_a = new cell(const_cast<char*>((data_path + x).c_str()), 0);
		//cell_a->get_center_enviro();
		//cell_a->get_center_valence(calculation_path+x+"/OUT.QDIV_scf");
		tt << x << " ";
		tt << cell_a->get_formula_name() << " ";
		tt << cell_a->num << " ";
		room << cell_a->check_contain_transition_metal();
		tt << string(room.str()) << " ";
		room.str("");
		room << cell_a->check_contain_water();
		tt << string(room.str()) << " ";
		room.str("");
		cout << x << " " << cell_a->get_cif_spacegroup(x) << endl;
		tt << cell_a->get_cif_spacegroup(x);
		//下面是输出能量和diff数值
		//tt << setprecision(12) <<cell_a->check_energy(calculation_path + x, "REPORT_scf_test")<<" ";
		//cin.get();
		//tt << cell_a->check_diff(calculation_path + x, "REPORT_scf_test") << " ";
		//tt << setprecision(12) << cell_a->check_energy(calculation_path + x, "REPORT_spin2_test") << " ";
		//tt << setprecision(12) << cell_a->check_diff(calculation_path +x, "REPORT_scf_test")<<" ";
		//tt << setprecision(12) << cell_a->check_energy(calculation_path + x, "REPORT") - cell_a->check_energy(calculation_path + x, "REPORT_spin2_test");
		//输出centr和周围元素
		/*for (int i = 0; i < cell_a->center_num; i++) {
			tt << aa[cell_a->center[i]];
			tt << ",";
		}
		tt << " ";
		for (int i = 0; i < cell_a->envoroment_num; i++) {
			tt << aa[cell_a->envoroment[i]];
			tt << ",";
		}
		tt << " ";
		for (int i = 0; i < cell_a->center_num; i++) {
			tt << setprecision(6)<<cell_a->center_valence[i];
			tt << ",";
		}*/
		tt << endl;


	}
	tt.close();
	cout << "output analyse complete!" << endl;
	cin.get();


	//有关测试分析数据的变量
	string need_to_check_path = "/share/home/wangz/2d_search/ridus_cut/result/result_ionicandval/0103/three_data/";
	string check_name = "check_water";
	int check_num = 9626;
	string check_out_name = "contain_water";
	ifstream fcheck;
	ofstream fcheck_out;
	fcheck.open(check_name, ios::in);
	fcheck_out.open(check_out_name, ios::out);
	string temp_check;
	for (int i = 0; i < check_num&& fcheck.good(); i++)
	{
		fcheck >> temp_check;
		if (i % 100 == 0)
		{
			cout << "has gone " << i << endl;
		}
		cell cell_a(const_cast<char*>((need_to_check_path + temp_check).c_str()));
		if (cell_a.check_contain_water() == true)
		{
			fcheck_out << temp_check << endl;
		}
	}
	fcheck.close();
	fcheck_out.close();
	cout << "检查含有水的2D结构的整理！" << endl;
	return 0;
	//先补充堰势文件
	string yanshi_name;
	string new_yanshi_path;
	full_fill_yanshi(yanshi_name, new_yanshi_path);
	cout << "fullfille yanshi compleste!" << endl;
	cin.get();
}

//这个是日常杂事需要处理的
void test(element*e ) {
	
	//绘制DOS结果
	Bandstru* bandstru = new Bandstru("/share2/wangz/0329_superO/0419_testAF/calculation/23163.cif.config_nou/pbe_dos", 2);
	bandstru->drawDos("/share2/wangz/0329_superO/0419_testAF/calculation/23163.cif.config_nou/pbe_dos", "DOS_save");
	Bandstru* bandstru1 = new Bandstru("/share2/wangz/0329_superO/0419_testAF/calculation/23163.cif.config_nou/hse_dos", 2);
	bandstru1->drawDos("/share2/wangz/0329_superO/0419_testAF/calculation/23163.cif.config_nou/hse_dos", "DOS_save");
	

	//批量进行输出能带结构
	Bandstru::generate_all_bandstu("/share2/wangz/0329_superO/0410_dirac/calculation/soc/name", 
		"/share2/wangz/0329_superO/0410_dirac/calculation/soc/");

	//批量产生IN,KPT
	/*report::allvaspKpt2pwmat("/share2/wangz/0329_superO/0404_plusconfig/kpath/narrow_path/name","/share2/wangz/0329_superO/0404_plusconfig/kpath/narrow_path/","/share2/wangz/0329_superO/0404_plusconfig/kpath/narrow_path/inkpt/");
	cin.get();*/

	//批量产生不同自旋情况的电子
	//Cal_cell::get_all_magnetic_for_all("/share2/wangz/1216_half_metal/0102/name", "/share2/wangz/1216_half_metal/0102/config/");

	//批量产生所有结构的半金属信息
	//string cal_path = "/share2/wangz/1216_half_metal/0114/strain_control/calculation_more/";
	//vector<string>half_name = get_file_name("/share2/wangz/1216_half_metal/0114/strain_control/calculation_more/name");
	//for (auto name : half_name) {
	//	report rt;
	//	cout << name << " ";
	//	rt.read_report(const_cast<char*>((cal_path+name+"/REPORT").c_str()),2);
	//	rt.get_half_situation(const_cast<char*>((cal_path + name+"/REPORT").c_str()));
	//}
	//cout << "finish!" << endl;
	//cin.get();
	////单个测试某个结构是不是半金属
	/*report rt;
	rt.read_report("/share2/wangz/0329_superO/0331_conductor/calculation/hse_dos/23163.config/REPORT_scf",2);
	cout << rt.orig_fermi << endl;
	vector<double>res = rt.get_half_situation("/share2/wangz/1216_half_metal/0114/hse_band_v2/calculation/432515.config/REPORT_scf");
	cout << res[0] << "," << res[1] << "," << res[2];
	cin.get();*/



	////针对某个结构进行变形
	//cell cell_a("/share2/wangz/1216_half_metal/0114/all_spin_situation/config/432515.config", 0);
	//cell_a.supper_cell(2, 2, 1,"/share2/wangz/1216_half_metal/0114/all_spin_situation/config/");
	//cin.get();

	//ofstream fout;
	//
	//vector<double> vacum = { -1,12,12};//重新定义真空层
	//vector<string>file_name = get_file_name("/share/home/wangz/high/1118_1d");
	//string config_path = "/share2/wangz/1007_formation_energy/singel_tes/1118_vacum_test/config/base/1d_after_relax/";
	//string out_path = "/share2/wangz/1007_formation_energy/singel_tes/1118_vacum_test/config/1d_change/";
	//for (auto name : file_name) {
	//	cell cell_a(const_cast<char*>((config_path + name ).c_str()), 0);
	//	cell_a.plus_vacum(vacum, out_path,name+"_12");
	//	cout << "transformed to :" << name << endl;
	//}
	//cout << "transform completed!" << endl;
	//cin.get();
	//vector<string> step(4, "");
	//vector<string> name = get_file_name(name_path+pattern_name);
	////每个结构对应一个二维数组，存储一定K点坐标及相关能量
	////用来储存相关的k点和对应的能量
	//vector<vector<double>> table(Electro_CAL_NUM, vector<double>(4));
	//for (string single_name : name) {
	//	cout << single_name << endl;
	//	//fout.open(name_path + "0908_"+single_name, ios::out);
	//	//cell cell_a((base_path  + single_name).c_str());
	//	//cell_a.trans_form(base_path + "config" + "/",single_name+"_trans");
	//	//fout.precision(9);
	//	//fout << single_name << "  " << cell_a.get_formula_name() << " " << cell_a.get_positive_distance() << "   "<< cell_a.check_energy(base_path +  single_name, "REPORT_relax")<<endl;
	//	//cell_a.plus_vacum(vacum, out_path, single_name);
	//	//get_last_result(test_path + single_name, step, table,1);
	//	//vector2file(table, fout);
	//	//cout << "trans form completed!" << endl;
	//	//fout.close();
	//}
	//




	//get_information("/share2/wangz/1007_formation_energy/space_group/lowd_name", "/share2/wangz/1007_formation_energy/space_group/config/");
	////这次的需求是，打印出diff，结合能，以及结构元素的判断情况
	//希望把这部分抽象出来，不要重复的代码
	


	string twod_name = "/share/home/wangz/high/1122_1d_formation_square";
	string bulk_cal_path = "/share2/wangz/1007_formation_energy/1d_new_bulk/calculation/";//存放体结构的计算路径
	string twod_cal_path_a = "/share2/wangz/1007_formation_energy/1d_new/calculation/";//存放2d部分计算路径
	string twod_cal_path_b = "/share2/wangz/1007_formation_energy/2d/calculation/";//存放2d部分计算路径
	string twod_config_path = "/share/home/wangz/2d_search/ridus_cut/result/result_ionicandval/0103/my_1d/config/";//原始的切下来的结构，用来计算一些相关数据


	vector<string> name = get_file_name(twod_name);
	ofstream fout;
	fout.open("formtion_1d", ios::out);
	double diff = 0;//diff数值
	int element_judge = 0;//元素判断情况
	for (auto singel : name) {
		cout << "start the :" << singel << endl;
		try {
			Cal_cell cal(singel, const_cast<char*>((twod_cal_path_a + singel+"-0_1d.config/atom.config").c_str()), 0);
			if (cal.Valid == false)
				continue;
			//首先获得元素判断情况
			element_judge = cal.check_main_transition_situation();
			/*bool super_flag = cal.check_contain_supoerO();
			if (super_flag == false)
				continue;*/
			fout << singel;
			cal.total_energy = cell::check_energy(twod_cal_path_a + singel + "-0_1d.config", "REPORT");
			cout << cal.total_energy << endl;
			diff = cal.check_diff(twod_cal_path_a + singel + "-0_1d.config", "REPORT");
			if (cal.total_energy == -100 || diff == -100) {
				cout << "???" << endl;
				cal.total_energy = cell::check_energy(twod_cal_path_b + singel + "-0_1d.config", "REPORT");
				diff = cal.check_diff(twod_cal_path_b + singel + "-0_1d.config", "REPORT");
			}
			if (cal.total_energy == -100) {
				cout << singel << "what the hell!怎么都没找到！！！" << endl;
				cout << "可能是没有yanshi或者其他问题，先不管了" << endl;
				cal.formation_energy_peratom = -100;

			}

			//然后开始绑定体结构
			Cal_cell* bulk = NULL;
			bulk = new Cal_cell(singel + "bulk", const_cast<char*>((bulk_cal_path + singel + ".cif.config" + "/atom.config").c_str()), 0);
			bulk->total_energy = cell::check_energy(bulk_cal_path + singel + ".cif.config", "REPORT");
			cal.bulk = bulk;
			cout << cal.bulk->total_energy << endl;
			//绑定初始结构
			cal.get_org_cell(twod_config_path + singel + "-0_1d.config");
			double res = 0;
			if (cal.formation_energy_peratom != -100) {
				cout << "计算结合能~" << endl;
				cal.get_formation_energy_persquare();
			}
			fout.precision(12);
			cal.output_inforamtion(fout);
			fout << endl;
			cout << "finished:" << singel << endl;
		}
		catch (...) {
			cout << "may encounter some wrong!" << endl;
		}



	}
	fout.close();
	//有关形成能的相关结果
	cout << "all total work done!" << endl;
	cin.get();

}














double report::get_accuracy(string file_name) {
	//获取该路径下的精确度
	//先只写1D情况的

	if (CALCULATION_MODE == 1) {
		double a, b;
		ifstream fin;
		fin.open(file_name, ios::in);
		if (!fin.is_open()) {
			cout << "file :" << file_name << " :doensnt exist!" << endl;
			cin.get();
			return -1;
		}
		char buffer[300];
		while (fin.good() && fin.peek() != EOF) {
			fin.getline(buffer, 300);
			if (strstr(buffer, "total number of K-point") != NULL) {
				//说明开始读到了
				fin >> a;
				fin.getline(buffer, 300);
				fin >> b;
				//cout << "res is:" << abs(a - b) << endl;
				//cin.get();
				fin.close();
				return abs(abs(a) - abs(b));

			}
		}
		fin.close();

	}
	else
		return -1;
}














void puls_magnetic(string path,string flag)
{
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
		cell cell_a(const_cast<char*>((path + "/atom.config").c_str()), 0);

		fstream fin;
		fin.open(path + "/atom.config", ios::app);
		if (!fin.is_open())
		{
			cout << "i cann not find the file" << path << endl;
			cin.get();
		}
		//然后开始书写初始磁矩
		fin << "MAGNETIC" << endl;
		if (flag == "F") {
			for (int i = 0; i < cell_a.num; i++) {
				//重新改变增加磁矩的策略
				//金属的全部默认为3			
				if (classify_metal_maingroup(cell_a.type[i]) == 1) {
					fin << "	" << cell_a.type[i] << " " << "2.0" << endl;
				}
				else
				{
					fin << "	" << cell_a.type[i] << " " << "0.0" << endl;
				}


			}
		}
		else if (flag == "AF") {
			bool positive = true;
			for (int i = 0; i < cell_a.num; i++) {
				//重新改变增加磁矩的策略
				//金属的全部默认为3			
				if (classify_metal_maingroup(cell_a.type[i]) == 1) {
					if(positive)
						fin << "	" << cell_a.type[i] << " " << "2.0" << endl;
					else
						fin << "	" << cell_a.type[i] << " " << "-2.0" << endl;
					positive = !positive;
				}
				else
				{
					fin << "	" << cell_a.type[i] << " " << "0.0" << endl;
				}


			}
		}
		
		fin.close();
		return;
	}
}

//给结构创造真空层
void cell::plus_vacum(vector<double>& vacum, string path,string file_name) {
	//输入情况每个真空层，每个数值小于1e-2视为不改变
	vector<double> now_vacum(3, 0);
	//第一个是最小值，第二个是最大值,记录的是原子的序号
	vector<double>a_save = { 100,-100 };
	vector<double>b_save = { 100,-100 };
	vector<double>c_save = { 100,-100 };
	vector<int> abc_xuhao_save(6, 0);
	vector<double>length_axle(3, -100);
	for (int i = 0; i < 3; i++) {
		length_axle[i] = pow(letice[i][0] * letice[i][0] + letice[i][1] * letice[i][1] + letice[i][2] * letice[i][2], 0.5);
	}
	for (int i = 0; i < this->num; i++) {
		if (this->p[i][0] < a_save[0]) {
			a_save[0] = this->p[i][0];
			abc_xuhao_save[0] = i;
		}
		if (this->p[i][0] > a_save[1]) {
			a_save[1] = this->p[i][0];
			abc_xuhao_save[1] = i;
		}

		if (this->p[i][1] < b_save[0]) {
			b_save[0] = this->p[i][1];
			abc_xuhao_save[2] = i;
		}
		if (this->p[i][1] > b_save[1]) {
			b_save[1] = this->p[i][1];
			abc_xuhao_save[3] = i;
		}

		if (this->p[i][2] < c_save[0]) {
			c_save[0] = this->p[i][2];
			abc_xuhao_save[4] = i;
		}
		if (this->p[i][2] > c_save[1]) {
			c_save[1] = this->p[i][2];
			abc_xuhao_save[5] = i;
		}


	}
	//首先获得当前的真空层是多少
	//需要修改真空层的定义，这里是投影到那个轴的距离
	for (int i = 0; i < 3; i++) {
		if (i == 0)
			now_vacum[i] = (1 + a_save[0] - a_save[1])*length_axle[i];
		else if (i == 1)
			now_vacum[i] = (1 + b_save[0] - b_save[1])*length_axle[i];
		else if (i == 2)
			now_vacum[i] = (1 + c_save[0] - c_save[1])*length_axle[i];
	}
	//打印当前的真空层大小
	
	for (double x : now_vacum)
		cout <<"now_cacum:"<< x << ",";
	
	//首先看看现在的letice是什么
	cout << "origin: letice:" << endl;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			cout << this->letice[i][j] << " ";
		}
		cout << endl;
	}
	vector<double> factor(3, 0);
	for (int i = 0; i < 3; i++) {
		if (vacum[i] < 1e-2)
			continue;
		factor[i] = (get_length(this->letice[i]) - now_vacum[i] + vacum[i]) / get_length(this->letice[i]);
		//cout << factor[i] << endl;
		for (int j = 0; j < 3; j++) {
			this->letice[i][j] = factor[i] * letice[i][j];
		}
	}
	//然后对于分数反变换回去
	double**letice_r = new double*[3];
	for (int i = 0; i < 3; i++) {
		letice_r[i] = new double[3];
	}

	inv(3, this->letice, letice_r);

	for (int i = 0; i < this->num; i++)
	{
		for (int m = 0; m < 3; m++)
		{
			this->p[i][m] = this->real_position[13][i][0] * letice_r[0][m] + this->real_position[13][i][1] * letice_r[1][m] + this->real_position[13][i][2] * letice_r[2][m];
			//cout << fenshu_plus[i][m] << ",";
		}
	}
	//看看现在的letice是什么
	/*cout << "origin: letice:" << endl;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			cout << this->letice[i][j] << " ";
		}
		cout << endl;
	}*/


	//最后输出结果文件
	output_config(path,file_name);
	//cin.get();
	return;

}
double cell::get_positive_distance() {
	//输出结构中的相邻晶胞的阳离子之间的距离
	double res = 1e5;
	int center = (yanshen - 1) / 2;
	double temp_dis = 0;
	//首先检测有没有阳离子
	bool posi_flag = false;
	int index = -1;//标记要记录的阳离子元素是什么
	for (int i = 0; i < num; i++) {
		if (my_classify[i] == 1) {
			posi_flag = true;
			index = type[i];
			true;
		}			
	}
	if(!posi_flag)
	{
		//如果没有阳离子，那么指定一个元素来检查距离
		index = type_count[0][1];

	}

	for (int i = 0; i < num; i++) {
		if (my_classify[i] == 1||type[i]==index) {
			for (int j = 0; j < yanshen*num; j++) {
				if (type[j%num] == type[i]) {
					temp_dis = dis(this->real_position[center][i], this->real_position[j / num][j%num]);
					if (abs(temp_dis) < 1e-2)
						continue;
					res = min(res, temp_dis);
				}

			}
		}
	}
	return res;
}
double cell::check_energy(string path, string file_name) {
	//从report文件中，读取能量信息
	string comm = "grep 'E_tot(eV)' " + file_name + " |tail -1";

	string check_comm = "ls |grep '" + file_name + "' |wc -l";

	char buf[200];
	char result[1024];
	FILE* fp = NULL;
	int i = 0;

	chdir(path.c_str());
	//system("pwd");
	//首先需要确认有这个文件
	fp = popen(check_comm.c_str(), "r");
	fgets(buf, sizeof(buf) - 1, fp);
	if (buf[0] == '0')
	{
		pclose(fp);
		return -100;
	}
	pclose(fp);


	fp = popen(comm.c_str(), "r");
	fgets(result, sizeof(result) - 1, fp);
	/*while (fgets(buf, sizeof(buf) - 1, fp) != NULL)
	{
		strcat(result, buf);
	}*/
	pclose(fp);
	fp = NULL;
	//cout << result << endl;
	//然后是从字符串中读出能量
	int start = 0;
	int end = 0;
	for (int i = 0; i < 200; i++) {
		if (result[i] == '-') {
			start = i;
			break;
		}
	}
	for (int i = start; i < 200; i++) {
		if (result[i] == ' ' || result[i] == '	')
		{
			end = i;
			break;
		}
	}
	if (start == 0 || end == 0) {
		cout << "identify the energy failed!" << endl;
		cout << "caiclulation may encounter error!" << endl;
		//cin.get();
		return -100;
	}
	//cout << start << "," << end << endl;
	string res = string(result);
	//cout << result << endl;

	string sub = res.substr(start, end);
	//cout << sub << endl;
	double cha = atof(sub.c_str());
	//cout << cha << endl;
	//cin.get();	
	return cha;
}

double cell::check_diff(string path, string file_name)
{
	char buf[100];
	char result[4096];
	FILE* fp = NULL;
	int i = 0;

	chdir(path.c_str());
	//首先需要确认有这个文件
	if (check_now_exist_file(file_name) == false) {
		file_name = file_name + "_scf";
	}
	string comm = "grep 'diff' " + file_name + " |tail -1";

	string check_comm = "ls |grep '" + file_name + "' |wc -l";
	fp = popen(check_comm.c_str(), "r");
	fgets(buf, sizeof(buf) - 1, fp);
	if (buf[0] == '0')
	{
		pclose(fp);
		return -100;
	}
	pclose(fp);


	fp = popen(comm.c_str(), "r");
	while (fgets(buf, sizeof(buf) - 1, fp) != NULL)
	{
		strcat(result, buf);
	}
	pclose(fp);
	fp = NULL;
	//然后是从字符串中读出最后一个数字

	string res = string(result);
	//cout << result << endl;
	int loc = res.find_last_of(" ");
	string sub = res.substr(loc, res.size() - 1);
	//cout << sub << endl;
	double cha = atof(sub.c_str());
	//cout << cha << endl;
	//cin.get();	
	return cha;
}
int cell::check_main_transition_situation() {
	bool contain_trans = false;
	bool contain_main = false;
	contain_trans = this->check_contain_transition_metal();
	//所谓的主族元素就是A的金属+之前认为的负价元素
	vector<int>main_element(Agroup_matal, Agroup_matal + 22);
	for (int i = 0; i < main_groupnum; i++) {
		main_element.push_back(main_group_element[i]);
	}
	sort(main_element.begin(), main_element.end());
	//下面判断是不是含有主族元素
	for (int i = 0; i < this->num; i++) {
		for (int j = 0; j < main_element.size(); j++) {
			if (this->type[i] == main_element[j])
			{
				contain_main = true;
				break;
			}
			if (main_element[j] > this->type[i])
				break;
		}
		if (contain_main == true)
			break;
	}
	//然后是判断结果
	if (contain_main && !contain_trans) {
		return 1;
	}
	else if (!contain_main&& contain_trans) {
		return 2;
	}
	else if (contain_main&& contain_trans) {
		return 3;
	}
	else {
		cout << "结构元素判断有误，请检查！" << endl;
		cin.get();
		return 4;
	}
}

bool cell::check_contain_transition_metal()
{
	//判断该结构是不是含有过渡金属
	vector<int>meatal_xuhao(Meatal_xuhao, Meatal_xuhao + 89);
	vector<int>A_matal(Agroup_matal, Agroup_matal + 22);
	bool trans_flag = true;
	bool meatal_flag = false;
	for (int i = 0; i < num; i++)
	{
		for (int j = 0; j < meatal_xuhao.size(); j++)
		{
			if (this->type[i] == meatal_xuhao[j])
			{
				meatal_flag = true;
				//cout << a[this->type[i]] << " :metal " << endl;
				for (int k = 0; k < A_matal.size(); k++)
				{
					if (this->type[i] == A_matal[k])
					{
						trans_flag = false;
						//cout << a[this->type[i]] << " :pure metal " << endl;
					}
				}
				if (meatal_flag && trans_flag)
				{
					return true;
				}
				break;
			}
		}
		meatal_flag = false;
		trans_flag = true;
	}
	return false;

}












void generate_winf90(string source_path, string dirpath, string name) {
	//输出是原路径，用来拷贝相关的文件，目的路径，以及输入的结构文件名字
	vector<string> config_name = get_file_name(name);
	Wannier_cell::get_source_dependent();//建立路径依赖
	string super_cell_path = "/share2/wangz/0527_got1T/test/super_cell/";
	string output_path = "/share/home/wangz/high/";
	string out_name = "wkm_valid_mass";
	
	for (string x : config_name) {
		const char *jiegou = (dirpath + x + "/atom.config").c_str();
		//cout << (dirpath + x + "/atom.config" )<< endl;
		cout << "start to " + x << endl;

		Wannier_cell wan((dirpath + x + "/atom.config").c_str(), 0);
		cout << wan.metal_name << "," << wan.main_group_name << endl;
		/*wan.supper_cell(3, 3, 3, super_cell_path);
		system(("mv " + super_cell_path + "atom_super.config " + super_cell_path + x).c_str());
		cout << "got super cell conplete!" << endl;*/
		wan.real_name = x;
		wan.get_source_name();
		cout << wan.source_name << endl;
		wan.get_wain90(dirpath + x, source_path);
		wan.run_wannier(dirpath + x, source_path);
		cout << "has generated " << x << endl;
		/*wan.run_wkm(dirpath + x, source_path);
		cout << "wkm calculation has procedd! " << x << endl;*/
		/*wan.get_ilmbda(dirpath + x);
		wan.output_ilmbda(output_path + wan.real_name + "_ilmbda");*/
		//cout << "ilmbda has got! " << x << endl;
		/*if (wan.lmbda_valid == true)
		{
			wan.run_last_scf(dirpath + x, source_path);
			cout << "last wkm scf has procedded!" << endl;
		}
		else {
			cout << "config:" << wan.real_name << " :scf encounter error!please check!" << endl;
		}*/
		/*wan.get_valid_m_new(dirpath + x,out_name);*/
		//wan.more_k_toget_band(dirpath + x);

	}

	return;
}

void Wannier_cell::run_wannier(string now_path, string source) {
	//进入这个路径进行相关脚本运行
	chdir(now_path.c_str());
	string com = "source /share/home/wangz/pwmat.intel.cuda.env";
	execv(com.c_str(), NULL);
	com = "source /share/home/wangz/PWMAT/20200902/pwmat.env";
	execv(com.c_str(), NULL);
	get_source_dependent();
	cout << "get environment dependent succeed!" << endl;

	com = "cp " + source + "wannier90.chk  ./";
	system(com.c_str());
	com = "cp " + source + "wannier90.pw2wan  ./";
	system(com.c_str());
	//产生cb的
	system("mkdir CB");
	system("mv wannier90.win_cb wannier90.win");
	cout << "run wannier90 command!" << endl;
	system("wannier90.x -pp wannier90");
	cout << "run pw2 command " << endl;
	//system("pw2wannier90.x < wannier90.pw2wan >check_cb");
	com = "cp " + source + "pwrun_cb ./";
	system(com.c_str());
	system("./pwrun_cb");
	cout << "run wannier90 last!" << endl;
	system("wannier90.x wannier90");
	cout << "generate the xsf file completed!" << endl;
	cout << "try to generate the IN.WANNIER  files!" << endl;
	system("ls |grep 'xsf' >need");
	int start = 1;
	vector<string> xsf_file = get_file_name("need");
	for (string xsf_name : xsf_file) {
		com = "convert_xsf_inwannier.x " + xsf_name + " CB/IN.WANNIER_000" + to_string(start++);
		system(com.c_str());
		cout << "got IN.WANNIER" + to_string(start - 1) << endl;
	}
	xsf_file.clear();
	//然后把所有的都移动到CB文件下
	system("mv UNK* *.xsf need wannier90.win wannier90.eig wannier90.amn wannier90.mmn wannier90.nnkp prefix.wf* CB/");
	cout << "CB got wannier function successfuly!" << endl;

	system("mkdir VB");
	system("mv wannier90.win_vb wannier90.win");
	cout << "run wannier90 command" << endl;
	system("wannier90.x -pp wannier90");
	cout << "run pw2 command!" << endl;
	com = "cp " + source + "pwrun_vb ./";
	system(com.c_str());
	system("./pwrun_vb");
	//system("pw2wannier90.x < wannier90.pw2wan >check_vb");
	cout << "run wannier90 last!" << endl;
	system("wannier90.x wannier90");
	cout << "generate the xsf file completed!" << endl;
	cout << "try to generate IN.WANNIER files!" << endl;
	system("ls |grep 'xsf' >need");
	start = 1;
	xsf_file = get_file_name("need");
	for (string xsf_name : xsf_file) {
		com = "convert_xsf_inwannier.x " + xsf_name + " VB/IN.WANNIER_000" + to_string(start++);
		system(com.c_str());
		cout << "got IN.WANNIER" + to_string(start - 1) << endl;
	}
	//然后把所有的都移动到CB文件下
	system("mv UNK* *.xsf need wannier90.win wannier90.eig wannier90.amn wannier90.mmn wannier90.nnkp prefix.wf* VB/");
	system("rm check_cb check_vb");
	cout << "VB got wannier function successfuly!" << endl;
}

vector<int> Wannier_cell::get_left_and_right(int start, int end, int num, int band, string flag) {
	vector<int> res(4, 0);
	res[0] = start;
	res[3] = end;
	if (flag == "CBM") {
		//如果是CBM的话需要从cbm开始往上
		res[1] = band - 1;
		res[2] = band + num;
		if (res[2] >= end) {
			res[2] = -100;
		}
	}
	else if (flag == "VBM") {
		//VBM的话需要往下数
		res[1] = band - num;
		if (res[1] <= start) {
			res[1] = -100;
		}
		res[2] = band + 1;
	}
	else
	{
		cout << "wrong flag!please check！" << endl;
	}
	return res;
}
void Wannier_cell::get_wain90(string path, string source_path) {
	//输入参数一个是结构所在路径，一个是制作好的文件路径
	chdir(path.c_str());

	//system("sleep 2");
	//先替换超胞参数
	string com;
	ofstream fout;
	//将堰势拷贝进去
	com = "cp *.UPF prefix.save/";
	system(com.c_str());


	//首先构造cb的,需要写进阳离子的位置
	com = "cp " + source_path + this->source_name + "_cb" + " wannier90.win_cb";
	system(com.c_str());
	//对于cb需要考虑report读取获得的相关信息
	fout.open("wannier90.win_cb", ios::app);
	report rp;
	rp.read_report("REPORT", 1);
	double efer = rp.get_singel_number("E_Fermi(eV)=");
	double cbm = rp.get_singel_number("CBM");

	//首先打印excludes_band
	vector<int>excludes = get_left_and_right(1, 80, 30, rp.cbm_band, "CBM");
	if (excludes[2] == -100) {
		fout << " exclude_bands = " << excludes[0] << "-" << excludes[1] << endl;
	}
	else
		fout << " exclude_bands = " << excludes[0] << "-" << excludes[1] << "," << excludes[2] << "-" << excludes[3] << endl;
	fout << " wannier_plot = T" << endl;
	fout << " wannier_plot_supercell = 4 1 1" << endl;


	fout << endl;
	fout << " dis_win_max = 100" << endl;
	fout << setprecision(3) << " dis_win_min = " << efer - 0.3 << endl;
	fout << setprecision(3) << " dis_froz_max = " << cbm + 0.05 << endl;
	fout << setprecision(3) << " dis_froz_min = " << efer - 0.3 << endl;
	fout << " dis_mix_ratio = 1" << endl;
	fout << endl;
	fout << " !spin=up" << endl;

	fout << "search_shells=24" << endl;
	fout << endl;
	fout << " begin projections" << endl;

	//区分三种情况
	if (this->source_name.find("SnBr2") != string::npos) {
		//单独给SnBr2开一条专线
		fout << " " << this->metal_name << ":l=1" << endl;
		fout << " " << this->main_group_name << ":l=1" << endl;
	}
	else if (this->source_name == "Si")
	{
		fout << " " << this->main_group_name << ":l=1" << endl;
	}
	else if (this->source_name == "GaAs") {
		fout << " " << this->metal_name << ":l=0;l=1" << endl;
		fout << " " << this->main_group_name << ":l=1" << endl;
	}
	else if (this->source_name.find("Ti") != string::npos)
	{
		fout << " " << this->metal_name << ":l=0;l=2" << endl;
	}
	else if (this->source_name.find("Mg") != string::npos)
	{
		fout << " " << this->metal_name << ":l=0" << endl;
		fout << " " << this->main_group_name << ":l=0;l=1" << endl;
	}
	else
	{
		fout << " " << this->metal_name << ":l=0" << endl;
		fout << " " << this->main_group_name << ":l=0;l=1" << endl;

	}
	fout << "  end projections" << endl;
	fout << endl;
	//书写晶胞信息
	fout << "  begin unit_cell_cart" << endl;
	fout << "  Ang" << endl;
	fout.setf(ios::fixed);
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++) {
			fout << setprecision(6) << this->letice[i][j] << "   ";
		}
		fout << endl;
	}
	fout << "  end unit_cell_cart" << endl;
	fout << endl;
	//书写晶格分量和格点设置
	fout << "  begin atoms_frac" << endl;
	for (int i = 0; i < this->num; i++) {
		fout << a[this->type[i]] << " ";
		for (int j = 0; j < 3; j++) {
			fout << setprecision(6) << this->p[i][j] << "   ";
		}
		fout << endl;
	}
	fout << "  end atoms_frac" << endl;
	fout << endl;
	fout << "  mp_grid= 50 1 1" << endl;
	fout << endl;
	fout << "  begin kpoints" << endl;
	fout.close();



	//接着是产生vb的
	com = "cp " + source_path + this->source_name + "_vb" + " wannier90.win_vb";
	system(com.c_str());
	fout.open("wannier90.win_vb", ios::app);
	excludes = get_left_and_right(1, 80, 14, rp.vbm_band, "VBM");
	if (excludes[1] == -100) {
		fout << " exclude_bands = " << excludes[2] << "-" << excludes[3] << endl;
	}
	else
	{
		fout << " exclude_bands = " << excludes[0] << "-" << excludes[1] << "," << excludes[2] << "-" << excludes[3] << endl;
	}

	fout << " wannier_plot = T" << endl;
	fout << " wannier_plot_supercell = 4 1 1" << endl;
	//打印上不是很需要的
	fout << endl;
	fout << " !dis_win_max =" << endl;
	fout << " !dis_win_min =" << endl;
	fout << " !dis_froz_max =" << endl;
	fout << " !dis_froz_min =" << endl;
	fout << " !dis_mix_ratio =" << endl;
	fout << endl;
	fout << "  !spin=up" << endl;
	fout << "search_shells=24" << endl;
	fout << endl;
	fout << " begin projections" << endl;

	if (this->source_name.find("SnBr2") != string::npos) {
		//单独给SnBr2开一条专线
		fout << " " << this->metal_name << ":l=0" << endl;
		fout << " " << this->main_group_name << ":l=1" << endl;
	}
	else if (this->source_name == "Si")
	{
		fout << " " << this->main_group_name << ":l=1" << endl;
	}
	else if (this->source_name == "GaAs") {
		fout << " " << this->metal_name << ":l=1" << endl;
		fout << " " << this->main_group_name << ":l=1" << endl;
	}
	else if (this->source_name.find("Zn") != string::npos)
	{
		for (int i = 0; i < this->num; i++) {
			if (this->my_classify[i] == 2) {
				fout << " " << a[this->type[i]] << ":l=1" << endl;
				break;
			}
		}
	}
	else
	{
		for (int i = 0; i < this->num; i++) {
			if (this->my_classify[i] == 2) {
				fout << " " << a[this->type[i]] << ":l=0;l=1" << endl;
				break;
			}
		}
	}

	fout << "  end projections" << endl;
	fout << endl;
	//书写晶胞信息
	fout << "  begin unit_cell_cart" << endl;
	fout << "  Ang" << endl;
	fout.setf(ios::fixed);
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++) {
			fout << setprecision(6) << this->letice[i][j] << "   ";
		}
		fout << endl;
	}
	fout << "  end unit_cell_cart" << endl;
	fout << endl;
	//书写晶格分量和格点设置
	fout << "  begin atoms_frac" << endl;
	for (int i = 0; i < this->num; i++) {
		fout << a[this->type[i]] << " ";
		for (int j = 0; j < 3; j++) {
			fout << setprecision(6) << this->p[i][j] << "   ";
		}
		fout << endl;
	}
	fout << "  end atoms_frac" << endl;
	fout << endl;
	fout << "  mp_grid= 50 1 1" << endl;
	fout << endl;
	fout << "  begin kpoints" << endl;
	fout.close();


	//添加上最后的k点坐标
	//int line_num= Wannier_cell::get_start_line();//返回的是带着标记的那一行行号
	int all_line = Wannier_cell::get_start_all_line();
	com = "cat OUT.KPT|tail -" + to_string(all_line - 2) + " >temp";
	system(com.c_str());
	com = "cat temp |cut -d ' ' -f 1-20 >>wannier90.win_vb";
	system(com.c_str());
	com = "cat temp |cut -d ' ' -f 1-20 >>wannier90.win_cb";
	system(com.c_str());
	system("echo '  End kpoints' >>wannier90.win_vb");
	system("echo '  End kpoints' >>wannier90.win_cb");
	system("rm temp");
	return;
}
int Wannier_cell::get_start_all_line() {
	//获得全部的行数
	string com = "cat OUT.KPT |wc -l >temp";
	system(com.c_str());
	ifstream fin;
	fin.open("temp", ios::in);
	int res = 0;
	fin >> res;
	fin.close();
	system("rm temp");
	//cout << res << endl;
	return res;
}

void Wannier_cell::output_ilmbda(string file_name) {
	//给定参数是文件名，输出是将获得的ilmbda输出到文件里面

	ofstream fout;
	fout.open(file_name, ios::out);
	if (fout.is_open()) {
		//将map持久化到文件
		map<string, double>::iterator it = this->lambda.begin();
		while (it != this->lambda.end()) {
			fout << it->first << " " << it->second << endl;
			it++;
		}
	}
	else
	{
		cout << "generate outname:" << file_name << " wrong!" << endl;
		cin.get();
	}
	fout.close();
	return;
}
double  Wannier_cell::get_singel_ilmbda(string path) {
	//读取一定目录下面的三个report，返回lambda
	double e[3] = { 0 };
	vector<string> name = { "0.0","0.5","1.0" };
	for (int i = 0; i < 3; i++) {
		if (if_finish(path + "/" + name[i]) != 1) {
			cout << "path " << name[i] << " unfinashed!" << endl;
			return -1000;
			cin.get();
		}
		else
		{
			e[i] = cell::check_energy(path + "/" + name[i]);
		}
	}
	return e[2] - e[0] - 4 * (e[1] - e[0]) + e[2] - e[0];

}
void Wannier_cell::get_ilmbda(string path) {
	//通过读取相关的文件,输入参数是读取路径
	chdir(path.c_str());

	//首先是获取cb下面的
	//需要注意将获得的ilmbda填充到数组中
	chdir((path + "/CB_wkm").c_str());
	string com;
	com = "ls >name";
	system(com.c_str());
	vector<string> name = get_file_name("name");
	for (string x : name) {
		if (x == "name")
			continue;
		double ilmbda = get_singel_ilmbda(path + "/CB_wkm" + "/" + x);
		if (abs(ilmbda + 1000) < 1e-2)
			this->lmbda_valid = false;
		this->lambda["CB " + x] = ilmbda;
		cout << "cb " << x << ": " << ilmbda << endl;
	}
	chdir((path + "/CB_wkm").c_str());
	com = "rm name";
	system(com.c_str());

	//接着读取vb下面的
	chdir((path + "/VB_wkm").c_str());
	com = "ls >name";
	system(com.c_str());
	name.clear();
	name = get_file_name("name");
	for (string x : name) {
		if (x == "name")
			continue;
		double ilmbda = get_singel_ilmbda(path + "/VB_wkm" + "/" + x);
		if (abs(ilmbda + 1000) < 1e-2)
			this->lmbda_valid = false;
		this->lambda["VB " + x] = ilmbda;
		cout << "vb " << x << ": " << ilmbda << endl;
	}
	chdir((path + "/VB_wkm").c_str());
	com = "rm name";
	system(com.c_str());

	//将字典变成vector储存
	if (this->source_name.find("Ti") != string::npos) {
		//首先是cb的
		ordered_lambda_cb.push_back(dict2tofind_double(lambda, "CB " + metal_name + "_s"));
		ordered_lambda_cb.push_back(dict2tofind_double(lambda, "CB " + metal_name + "_dz2"));
		ordered_lambda_cb.push_back(dict2tofind_double(lambda, "CB " + metal_name + "_dxz"));
		ordered_lambda_cb.push_back(dict2tofind_double(lambda, "CB " + metal_name + "_dxz"));
		ordered_lambda_cb.push_back(dict2tofind_double(lambda, "CB " + metal_name + "_dx2-y2"));
		ordered_lambda_cb.push_back(dict2tofind_double(lambda, "CB " + metal_name + "_dx2-y2"));

		//然后是VB的
		vector<double >temp_vb;
		temp_vb.push_back(dict2tofind_double(lambda, "VB " + main_group_name + "_s"));
		temp_vb.push_back(dict2tofind_double(lambda, "VB " + main_group_name + "_p"));
		temp_vb.push_back(dict2tofind_double(lambda, "VB " + main_group_name + "_p"));
		temp_vb.push_back(dict2tofind_double(lambda, "VB " + main_group_name + "_p"));
		ordered_lambda_vb = temp_vb;
		ordered_lambda_vb.insert(ordered_lambda_vb.end(), temp_vb.begin(), temp_vb.end());

	}
	else if (this->source_name.find("Mg") != string::npos) {
		//首先是cb的
		vector<double >temp_cb;
		temp_cb.push_back(dict2tofind_double(lambda, "CB " + main_group_name + "_s"));
		temp_cb.push_back(dict2tofind_double(lambda, "CB " + main_group_name + "_p"));
		temp_cb.push_back(dict2tofind_double(lambda, "CB " + main_group_name + "_p"));
		temp_cb.push_back(dict2tofind_double(lambda, "CB " + main_group_name + "_p"));
		ordered_lambda_cb.push_back(dict2tofind_double(lambda, "CB " + metal_name + "_s"));
		ordered_lambda_cb.insert(ordered_lambda_cb.end(), temp_cb.begin(), temp_cb.end());
		ordered_lambda_cb.insert(ordered_lambda_cb.end(), temp_cb.begin(), temp_cb.end());

		//然后是vb的
		vector<double >temp_vb;
		temp_vb.push_back(dict2tofind_double(lambda, "VB " + main_group_name + "_s"));
		temp_vb.push_back(dict2tofind_double(lambda, "VB " + main_group_name + "_p"));
		temp_vb.push_back(dict2tofind_double(lambda, "VB " + main_group_name + "_p"));
		temp_vb.push_back(dict2tofind_double(lambda, "VB " + main_group_name + "_p"));
		ordered_lambda_vb = temp_vb;
		ordered_lambda_vb.insert(ordered_lambda_vb.end(), temp_vb.begin(), temp_vb.end());
	}
	else
	{
		//最后是含Zn的这一类
		//首先是cb的
		vector<double >temp_cb;
		temp_cb.push_back(dict2tofind_double(lambda, "CB " + main_group_name + "_s"));
		temp_cb.push_back(dict2tofind_double(lambda, "CB " + main_group_name + "_p"));
		temp_cb.push_back(dict2tofind_double(lambda, "CB " + main_group_name + "_p"));
		temp_cb.push_back(dict2tofind_double(lambda, "CB " + main_group_name + "_p"));
		ordered_lambda_cb.push_back(dict2tofind_double(lambda, "CB " + metal_name + "_s"));
		ordered_lambda_cb.insert(ordered_lambda_cb.end(), temp_cb.begin(), temp_cb.end());
		ordered_lambda_cb.insert(ordered_lambda_cb.end(), temp_cb.begin(), temp_cb.end());

		//然后是vb的
		double p = dict2tofind_double(lambda, "VB " + main_group_name + "_p");
		for (int i = 0; i < 6; i++) {
			ordered_lambda_vb.push_back(p);
		}

	}
	return;
}
double Wannier_cell::dict2tofind_double(map<string, double> mat, string key)
{
	//获得相关的double数值，从map中
	bool falg = false;
	map<string, double>::iterator it = mat.begin();
	for (it = mat.begin(); it != mat.end(); it++) {
		if (it->first.find(key) != string::npos) {
			falg = true;
			return it->second;
		}
	}
	if (falg == false) {
		return -1000;
	}
}
void Wannier_cell::get_wannier_param() {
	ofstream fout;
	fout.setf(ios::fixed);

	fout.open("IN.WANNIER_PARAM", ios::out);
	//注意先打印vb，再打印cb
	report* rp = new report();
	rp->read_report("REPORT", 1);//这里存疑注意一下
	fout << " " << 4 * rp->n_x << " " << 4 * rp->n_y << " " << rp->n_z << endl;
	fout << " " << this->ordered_lambda_vb.size() + this->ordered_lambda_cb.size() << " " << "1" << endl;
	fout << " " << "------------------" << endl;

	for (int i = 0; i < 2; i++) {
		//注意要写两次一模一样的
		//首先产生vb的
		for (int i = 0; i < this->ordered_lambda_vb.size(); i++) {
			fout << " " << i + 1 << " " << setprecision(6) << ordered_lambda_vb[i] << endl;
		}
		//然后是cb的
		for (int j = this->ordered_lambda_vb.size() + 1; j < (this->ordered_lambda_vb.size() + this->ordered_lambda_cb.size() + 1); j++)
		{
			fout << " " << j << " " << setprecision(6) << ordered_lambda_cb[j - this->ordered_lambda_vb.size() - 1] << endl;
		}
		fout << " " << "------------------" << endl;
	}
	fout << " " << "1  0  0  0" << endl;
	fout << " " << "------------------" << endl;
	fout.close();
	return;
}

void Wannier_cell::run_last_scf(string now_path, string source) {
	//进行最后的计算的函数
	if (this->lmbda_valid == true) {
		chdir(now_path.c_str());
		string com;
		get_wannier_param();
		cout << "got wannier param completed!" << endl;
		//需要准备in.wannier的文件，注意格式
		//还是先把vb的拿过来
		for (int i = 0; i < ordered_lambda_vb.size(); i++) {
			if ((i + 1) >= 10)
				com = "cp VB/IN.WANNIER_000" + to_string(i + 1) + " IN.WANNIER_000" + to_string(i + 1) + ".u";
			else
				com = "cp VB/IN.WANNIER_000" + to_string(i + 1) + " IN.WANNIER_0000" + to_string(i + 1) + ".u";
			system(com.c_str());
			if ((i + 1) >= 10)
				com = "cp VB/IN.WANNIER_000" + to_string(i + 1) + " IN.WANNIER_000" + to_string(i + 1) + ".d";
			else
				com = "cp VB/IN.WANNIER_000" + to_string(i + 1) + " IN.WANNIER_0000" + to_string(i + 1) + ".d";
			system(com.c_str());
		}
		//然后是CB的拿过来
		for (int i = ordered_lambda_vb.size(); i < ordered_lambda_vb.size() + ordered_lambda_cb.size(); i++) {

			if ((i + 1) >= 10)
				com = "cp CB/IN.WANNIER_000" + to_string(i - ordered_lambda_vb.size() + 1) + " IN.WANNIER_000" + to_string(i + 1) + ".u";
			else
				com = "cp CB/IN.WANNIER_000" + to_string(i - ordered_lambda_vb.size() + 1) + " IN.WANNIER_0000" + to_string(i + 1) + ".u";
			system(com.c_str());
			if ((i + 1) >= 10)
				com = "cp CB/IN.WANNIER_000" + to_string(i - ordered_lambda_vb.size() + 1) + " IN.WANNIER_000" + to_string(i + 1) + ".d";
			else
				com = "cp CB/IN.WANNIER_000" + to_string(i - ordered_lambda_vb.size() + 1) + " IN.WANNIER_0000" + to_string(i + 1) + ".d";
			system(com.c_str());
		}
		cout << this->real_name << " :prepare inwannier completed!" << endl;
		//下面开始是准备计算文件以及脚本

		com = "mv aaa.pbs aaa_scf.pbs";
		system(com.c_str());
		com = "mv etot.input etot_scf.input";
		system(com.c_str());
		com = "mv REPORT REPORT_scf";
		system(com.c_str());
		com = "mv output output_scf";
		system(com.c_str());


		report* rp = new report();
		rp->read_report("REPORT_scf", 1);
		sole_generate_jobs("atom.config", "./", "WKM_LAST", rp);
		cout << "got etot for last scf completed！" << endl;
		string default_name = "gpu3";
		//cout << real_name << endl;
		com = real_name;
		int detect_flag = detect_free("gpu2");
		if (detect_flag == 1)
			pbs_got(com, "scf", now_path, "gpu2");//这里存疑，需要根据
		else if (detect_flag == 2)
			pbs_got(com, "scf", now_path, "gpu3");//这里存疑，需要根据
		else
			pbs_got(com, "scf", now_path, default_name);
		chdir(now_path.c_str());
		system("qsub aaa.pbs");
		return;
	}


}


void Wannier_cell::get_valid_m_new(string now_path, string name) {
	//新的方法获取有效质量
	//主要是进行对于etot.input进行修改和重跑来获得的
	//需要准备IN.KPT等文件
	string com;

	chdir(now_path.c_str());
	//首先准备IN.KPT
	int all_line = Wannier_cell::get_start_all_line();
	system("cp OUT.KPT BACKUP_KPT");
	system("mv OUT.KPT IN.KPT");
	ifstream fin;
	fin.open("IN.KPT", ios::in);
	int kpt_num = 0;
	fin >> kpt_num;
	cout << kpt_num << endl;
	cin.get();
	com = "cat IN.KPT|tail -" + to_string(all_line - 2) + " >temp";
	system(com.c_str());
	int inter_num = 1;//在每一个方向上插入多少个点
	kpt_num += 8 * inter_num * 2;
	ofstream fout;
	fout.open("IN.KPT", ios::out);
	fout << "            " << kpt_num << endl;
	fout << "            2           0" << endl;
	fout.close();
	system("cat temp>>IN.KPT");
	report rp;
	rp.read_report(const_cast<char*>((now_path + "/REPORT").c_str()), 1);
	//注意这里的数值
	double inter_x_length = 0, inter_y_length = 0;
	if (this->real_name == "BaI2.config" || this->real_name == "CdI2.config" || this->real_name == "HgI2.config" || this->real_name == "SrI2.config") {
		inter_x_length = (1.0 / 12.0) / (inter_num + 1);
		inter_y_length = (1.0 / 12.0) / (inter_num + 1);
	}
	else
	{
		inter_x_length = (1.0 / 20.0) / (inter_num + 1);
		inter_y_length = (1.0 / 20.0) / (inter_num + 1);
	}
	int new_num = 2 * inter_num + 1;
	fout.open("IN.KPT", ios::app);
	fout.setf(ios::fixed);
	for (int i = 0; i < new_num; i++) {
		for (int j = 0; j < new_num; j++) {
			fout << setprecision(8) << "     " << rp.cbm_kppt[0] + i * inter_x_length << "     " << rp.cbm_kppt[1] + j * inter_y_length << "      " << "0.00000000" << "      " << "0.00000000" << endl;
		}
	}
	//同样产生针对VB的	
	for (int i = 0; i < new_num; i++) {
		for (int j = 0; j < new_num; j++) {
			fout << setprecision(8) << "     " << rp.vbm_kppt[0] + i * inter_x_length << "     " << rp.vbm_kppt[1] + j * inter_y_length << "      " << "0.00000000" << "      " << "0.00000000" << endl;
		}
	}
	fout.close();
	system("rm temp");
	fin.close();

	//然后是针对etot.input
	system("cp etot.input  etot_last");
	system("sed -i '7d' etot.input");
	system("echo '    IN.KPT=T' >>etot.input");

	//最后是将report保存一下
	system("cp REPORT REPORT_last");
	system("qsub aaa.pbs");
	return;


}
void Wannier_cell::get_valid_m(string now_path, string out_name) {
	//获得wkm下的有效质量
	//具体做法是比较pbe下的找到的cbm和vbm的位置，相对应的在wkm的report中寻找
	//需要注意几个k点生成是12 12 1 的
	//添加输出功能，注意输出的格式，需要输出相关的用来插值的数值
	int colom = 20;
	if (this->real_name.find("BaI2") != string::npos || this->real_name.find("CdI2") != string::npos || real_name.find("HgI2") != string::npos || real_name.find("SrI2") != string::npos) {
		colom = 12;
	}
	chdir(now_path.c_str());
	report rp;
	rp.read_report(const_cast<char*>((Wannier_cell::pbe_path + real_name + "/REPORT").c_str()), 1);;
	//相当于获取了位置
	report rp_now;
	rp_now.read_report("REPORT", 2);
	//需要获得第几个k点和第几条band和它是最接近的
	//先是vb的
	int kk = find_close_kppt(rp.vbm_kppt, rp_now);
	int band = find_close_band(rp.vbm_energy, kk, rp_now);
	int save_band = band;
	string base_path = "/share/home/wangz/high/";
	ofstream fout, fout2;
	fout.open(base_path + out_name, ios::app);
	fout2.open(base_path + out_name + "_data", ios::app);
	fout << real_name << " " << endl;
	//fout << band << " ";
	//先计算x方向的有效质量
	double a1 = 0, a2 = 0, a3 = 0, b = 0;
	double a4 = 0, a5 = 0, a6 = 0, a7 = 0;//分别对应的左上右上左下右下的能量
	a2 = rp_now.eig_kns[kk][band][rp_now.spin - 1] / HARTREE;

	if (kk % colom == 0)
	{
		a1 = a3 = rp_now.eig_kns[kk + 1][band][rp_now.spin - 1] / HARTREE;
		b = dis(rp_now.kp_p[kk], rp_now.kp_p[kk + 1]);
		if (kk / colom == 0) {
			a4 = a5 = a6 = a7 = rp_now.eig_kns[kk + 1 + colom][band][rp_now.spin - 1] / HARTREE;
		}
		else if (kk / colom == colom - 1) {
			a4 = a5 = a6 = a7 = rp_now.eig_kns[kk + 1 - colom][band][rp_now.spin - 1] / HARTREE;
		}
		else
		{
			a4 = a5 = rp_now.eig_kns[kk + 1 - colom][band][rp_now.spin - 1] / HARTREE;
			a6 = a7 = rp_now.eig_kns[kk + 1 + colom][band][rp_now.spin - 1] / HARTREE;
		}
	}
	else if (kk % colom == colom - 1)
	{
		a1 = a3 = rp_now.eig_kns[kk - 1][band][rp_now.spin - 1] / HARTREE;
		b = dis(rp_now.kp_p[kk], rp_now.kp_p[kk - 1]);
		if (kk / colom == 0) {
			a4 = a5 = a6 = a7 = rp_now.eig_kns[kk - 1 + colom][band][rp_now.spin - 1] / HARTREE;
		}
		else if (kk / colom == colom - 1) {
			a4 = a5 = a6 = a7 = rp_now.eig_kns[kk - 1 - colom][band][rp_now.spin - 1] / HARTREE;
		}
		else
		{
			a4 = a5 = rp_now.eig_kns[kk - 1 - colom][band][rp_now.spin - 1] / HARTREE;
			a6 = a7 = rp_now.eig_kns[kk - 1 + colom][band][rp_now.spin - 1] / HARTREE;
		}
	}
	else
	{
		a1 = rp_now.eig_kns[kk - 1][band][rp_now.spin - 1] / HARTREE;
		a3 = rp_now.eig_kns[kk + 1][band][rp_now.spin - 1] / HARTREE;
		b = dis(rp_now.kp_p[kk], rp_now.kp_p[kk + 1]);

		if (kk / colom == 0) {
			a4 = a6 = rp_now.eig_kns[kk - 1 + colom][band][rp_now.spin - 1] / HARTREE;
			a5 = a7 = rp_now.eig_kns[kk + 1 + colom][band][rp_now.spin - 1] / HARTREE;
		}
		else if (kk / colom == colom - 1) {
			a4 = a6 = rp_now.eig_kns[kk - 1 - colom][band][rp_now.spin - 1] / HARTREE;
			a5 = a7 = rp_now.eig_kns[kk + 1 - colom][band][rp_now.spin - 1] / HARTREE;
		}
		else
		{
			a4 = rp_now.eig_kns[kk - 1 - colom][band][rp_now.spin - 1] / HARTREE;
			a5 = rp_now.eig_kns[kk + 1 - colom][band][rp_now.spin - 1] / HARTREE;
			a6 = rp_now.eig_kns[kk - 1 + colom][band][rp_now.spin - 1] / HARTREE;
			a7 = rp_now.eig_kns[kk + 1 + colom][band][rp_now.spin - 1] / HARTREE;
		}
	}

	cout << real_name << " :x_vacum " << 1 / ((a1 + a3 - 2 * a2) / b / b) << endl;
	fout << 1 / ((a1 + a3 - 2 * a2) / b / b) << " ";
	valid_vb.push_back(b);
	valid_vb.push_back(a1);
	valid_vb.push_back(a2);
	valid_vb.push_back(a3);
	/*if (!(kk % 20 == 0 || kk % 20 == 19)) {

	}
	else
	{

		cout << "kx point wrong number!" << endl;
	}*/

	a2 = rp_now.eig_kns[kk][band][rp_now.spin - 1] / HARTREE;
	if (kk / colom == 0)
	{
		a3 = a1 = rp_now.eig_kns[kk + colom][band][rp_now.spin - 1] / HARTREE;
		b = dis(rp_now.kp_p[kk], rp_now.kp_p[kk + colom]);
	}
	else if (kk / colom == colom - 1)
	{
		a3 = a1 = rp_now.eig_kns[kk - colom][band][rp_now.spin - 1] / HARTREE;
		b = dis(rp_now.kp_p[kk], rp_now.kp_p[kk - colom]);
	}
	else
	{
		a1 = rp_now.eig_kns[kk - colom][band][rp_now.spin - 1] / HARTREE;
		a3 = rp_now.eig_kns[kk + colom][band][rp_now.spin - 1] / HARTREE;
		b = dis(rp_now.kp_p[kk], rp_now.kp_p[kk + colom]);
	}


	cout << real_name << " :y_vacum " << 1 / ((a1 + a3 - 2 * a2) / b / b) << endl;
	fout << 1 / ((a1 + a3 - 2 * a2) / b / b) << " ";
	valid_vb.push_back(b);
	valid_vb.push_back(a1);
	valid_vb.push_back(a2);
	valid_vb.push_back(a3);
	//别忘了添加上其他四个的数值




	valid_vb.push_back(a4);
	valid_vb.push_back(a5);
	valid_vb.push_back(a6);
	valid_vb.push_back(a7);
	//然后是CB的
	kk = find_close_kppt(rp.cbm_kppt, rp_now);
	band = find_close_band(rp.cbm_energy, kk, rp_now);
	if (band - save_band > 1)
	{
		cout << "unsula band !" << endl;
		cin.get();
	}
	if (band == save_band)
		band += 1;
	//fout << band << " ";
	//先计算x方向的有效质量
	if (kk % colom == 0)
	{
		a1 = a3 = rp_now.eig_kns[kk + 1][band][rp_now.spin - 1] / HARTREE;
		b = dis(rp_now.kp_p[kk], rp_now.kp_p[kk + 1]);
		if (kk / colom == 0) {
			a4 = a5 = a6 = a7 = rp_now.eig_kns[kk + 1 + colom][band][rp_now.spin - 1] / HARTREE;
		}
		else if (kk / colom == colom - 1) {
			a4 = a5 = a6 = a7 = rp_now.eig_kns[kk + 1 - colom][band][rp_now.spin - 1] / HARTREE;
		}
		else
		{
			a4 = a5 = rp_now.eig_kns[kk + 1 - colom][band][rp_now.spin - 1] / HARTREE;
			a6 = a7 = rp_now.eig_kns[kk + 1 + colom][band][rp_now.spin - 1] / HARTREE;
		}
	}
	else if (kk % colom == colom - 1)
	{
		a1 = a3 = rp_now.eig_kns[kk - 1][band][rp_now.spin - 1] / HARTREE;
		b = dis(rp_now.kp_p[kk], rp_now.kp_p[kk - 1]);
		if (kk / colom == 0) {
			a4 = a5 = a6 = a7 = rp_now.eig_kns[kk - 1 + colom][band][rp_now.spin - 1] / HARTREE;
		}
		else if (kk / colom == colom - 1) {
			a4 = a5 = a6 = a7 = rp_now.eig_kns[kk - 1 - colom][band][rp_now.spin - 1] / HARTREE;
		}
		else
		{
			a4 = a5 = rp_now.eig_kns[kk - 1 - colom][band][rp_now.spin - 1] / HARTREE;
			a6 = a7 = rp_now.eig_kns[kk - 1 + colom][band][rp_now.spin - 1] / HARTREE;
		}
	}
	else
	{
		a1 = rp_now.eig_kns[kk - 1][band][rp_now.spin - 1] / HARTREE;
		a3 = rp_now.eig_kns[kk + 1][band][rp_now.spin - 1] / HARTREE;
		b = dis(rp_now.kp_p[kk], rp_now.kp_p[kk + 1]);

		if (kk / colom == 0) {
			a4 = a6 = rp_now.eig_kns[kk - 1 + colom][band][rp_now.spin - 1] / HARTREE;
			a5 = a7 = rp_now.eig_kns[kk + 1 + colom][band][rp_now.spin - 1] / HARTREE;
		}
		else if (kk / colom == colom - 1) {
			a4 = a6 = rp_now.eig_kns[kk - 1 - colom][band][rp_now.spin - 1] / HARTREE;
			a5 = a7 = rp_now.eig_kns[kk + 1 - colom][band][rp_now.spin - 1] / HARTREE;
		}
		else
		{
			a4 = rp_now.eig_kns[kk - 1 - colom][band][rp_now.spin - 1] / HARTREE;
			a5 = rp_now.eig_kns[kk + 1 - colom][band][rp_now.spin - 1] / HARTREE;
			a6 = rp_now.eig_kns[kk - 1 + colom][band][rp_now.spin - 1] / HARTREE;
			a7 = rp_now.eig_kns[kk + 1 + colom][band][rp_now.spin - 1] / HARTREE;
		}
	}
	cout << real_name << " :x_elec " << 1 / ((a1 + a3 - 2 * a2) / b / b) << endl;
	fout << 1 / ((a1 + a3 - 2 * a2) / b / b) << " ";
	valid_cb.push_back(b);
	valid_cb.push_back(a1);
	valid_cb.push_back(a2);
	valid_cb.push_back(a3);
	/*if (!(kk % 20 == 0 || kk % 20 == 19)) {

	}
	else
	{
		cout << "kx point wrong number!" << endl;
	}*/
	a2 = rp_now.eig_kns[kk][band][rp_now.spin - 1] / HARTREE;
	if (kk / colom == 0)
	{
		a3 = a1 = rp_now.eig_kns[kk + colom][band][rp_now.spin - 1] / HARTREE;
		b = dis(rp_now.kp_p[kk], rp_now.kp_p[kk + colom]);
	}
	else if (kk / colom == colom - 1)
	{
		a3 = a1 = rp_now.eig_kns[kk - colom][band][rp_now.spin - 1] / HARTREE;
		b = dis(rp_now.kp_p[kk], rp_now.kp_p[kk - colom]);
	}
	else
	{
		a1 = rp_now.eig_kns[kk - colom][band][rp_now.spin - 1] / HARTREE;
		a3 = rp_now.eig_kns[kk + colom][band][rp_now.spin - 1] / HARTREE;
		b = dis(rp_now.kp_p[kk], rp_now.kp_p[kk + colom]);
	}
	cout << real_name << " :y_elec " << 1 / ((a1 + a3 - 2 * a2) / b / b) << endl;
	fout << 1 / ((a1 + a3 - 2 * a2) / b / b) << endl;
	valid_cb.push_back(b);
	valid_cb.push_back(a1);
	valid_cb.push_back(a2);
	valid_cb.push_back(a3);

	valid_cb.push_back(a4);
	valid_cb.push_back(a5);
	valid_cb.push_back(a6);
	valid_cb.push_back(a7);
	/*if (!(kk / 20 == 0 || kk / 20 == 19))
	{
	}
	else
	{
		cout << "ky point wrong number!" << endl;
	}*/
	for (double x : valid_vb)
		fout2 << x << " ";
	fout2 << endl;
	for (double x : valid_cb)
		fout2 << x << " ";
	fout2 << endl;
	fout.close();
	fout2.close();
	return;
}



void Wannier_cell::prepare_singel_file(string dir_path, string base_path, string source_path, string flag, string x) {
	//为1个文件夹准备所有的输入文件
	string default_name = "gpu3";
	chdir(dir_path.c_str());
	system(("cp " + base_path + "/atom_super.config " + " atom.config").c_str());//获得atom.config
	system(("cp " + base_path + "/*.UPF " + "./").c_str());//获得堰势文件
	system(("cp " + base_path + "/REPORT " + "./REPORT_stepa").c_str());//获得堰势文件
	report* rp = new report();
	rp->read_report("REPORT_stepa", 1);
	if (flag == "wkm_vb")
		sole_generate_jobs("atom.config", "./", "wkm_vb", rp);//获得etot.input
	else if (flag == "wkm_cb")
		sole_generate_jobs("atom.config", "./", "wkm_cb", rp);//获得etot.input
	int detect_flag = detect_free("gpu2");
	if (detect_flag == 1)
		pbs_got(real_name + x, "scf", dir_path, "gpu2");//这里存疑，需要根据
	else if (detect_flag == 2)
		pbs_got(real_name + x, "scf", dir_path, "gpu3");//这里存疑，需要根据
	else
		pbs_got(real_name + x, "scf", dir_path, default_name);
	//获得aaa.pbs
	chdir(dir_path.c_str());
	system(("cp " + source_path + "IN.S_WKM_" + x + " " + "IN.S_WKM").c_str());//获得in。swkm
	//最重要要获得IN,WANNIER文件
	return;
}
void Wannier_cell::run_wkm(string nowpath, string source_path) {
	//直接运行
	string com;
	vector<string> test = { "0.0","0.5","1.0" };
	chdir(nowpath.c_str());
	this->supper_cell(3, 3, 3, nowpath + "/");//产生超胞文件,产生的文件名是atom_super.config
	system("mkdir CB_wkm VB_wkm");
	string temp_despath;

	//先处理VB的东西
	if (this->source_name.find("Si") != string::npos) {
		//添加上相关的有关测试的
		chdir(nowpath.c_str());
		com = "mkdir -p VB_wkm/" + main_group_name + "_1s";
		system(com.c_str());
		temp_despath = nowpath + "/VB_wkm/" + main_group_name + "_1s";
		for (string x : test) {
			chdir((temp_despath).c_str());
			system(("mkdir " + x).c_str());
			prepare_singel_file(temp_despath + "/" + x, nowpath, source_path, "wkm_vb", x);
			system(("cp " + nowpath + "/VB/IN.WANNIER_1 " + temp_despath + "/" + x + "/IN.WANNIER_00001.u").c_str());
			system(("qsub " + temp_despath + "/" + x + "/aaa.pbs").c_str());
			cout << "VB:" << this->real_name << ":" << x << " :1_s start!" << endl;
		}

		chdir(nowpath.c_str());
		com = "mkdir -p VB_wkm/" + main_group_name + "1_p";
		system(com.c_str());
		temp_despath = nowpath + "/VB_wkm/" + main_group_name + "1_p";
		for (string x : test) {
			chdir((temp_despath).c_str());
			system(("mkdir " + x).c_str());
			prepare_singel_file(temp_despath + "/" + x, nowpath, source_path, "wkm_vb", x);
			system(("cp " + nowpath + "/VB/IN.WANNIER_2 " + temp_despath + "/" + x + "/IN.WANNIER_00001.u").c_str());
			system(("qsub " + temp_despath + "/" + x + "/aaa.pbs").c_str());
			cout << "VB:" << this->real_name << ":" << x << " :1_p start!" << endl;
		}

		chdir(nowpath.c_str());
		com = "mkdir -p VB_wkm/" + main_group_name + "2_s";
		system(com.c_str());
		temp_despath = nowpath + "/VB_wkm/" + main_group_name + "2_s";
		for (string x : test) {
			chdir((temp_despath).c_str());
			system(("mkdir " + x).c_str());
			prepare_singel_file(temp_despath + "/" + x, nowpath, source_path, "wkm_vb", x);
			system(("cp " + nowpath + "/VB/IN.WANNIER_5 " + temp_despath + "/" + x + "/IN.WANNIER_00001.u").c_str());
			system(("qsub " + temp_despath + "/" + x + "/aaa.pbs").c_str());
			cout << "VB:" << this->real_name << ":" << x << " :2_s start!" << endl;
		}

		chdir(nowpath.c_str());
		com = "mkdir -p VB_wkm/" + main_group_name + "2_p";
		system(com.c_str());
		temp_despath = nowpath + "/VB_wkm/" + main_group_name + "2_p";
		for (string x : test) {
			chdir((temp_despath).c_str());
			system(("mkdir " + x).c_str());
			prepare_singel_file(temp_despath + "/" + x, nowpath, source_path, "wkm_vb", x);
			system(("cp " + nowpath + "/VB/IN.WANNIER_6 " + temp_despath + "/" + x + "/IN.WANNIER_00001.u").c_str());
			system(("qsub " + temp_despath + "/" + x + "/aaa.pbs").c_str());
			cout << "VB:" << this->real_name << ":" << x << " :2_p start!" << endl;
		}

	}
	else if (this->source_name.find("GaAs") != string::npos) {
		//添加上相关的有关测试的
		chdir(nowpath.c_str());
		com = "mkdir -p VB_wkm/" + main_group_name + "_p";
		system(com.c_str());
		temp_despath = nowpath + "/VB_wkm/" + main_group_name + "_p";
		for (string x : test) {
			chdir((temp_despath).c_str());
			system(("mkdir " + x).c_str());
			prepare_singel_file(temp_despath + "/" + x, nowpath, source_path, "wkm_vb", x);
			system(("cp " + nowpath + "/VB/IN.WANNIER_13 " + temp_despath + "/" + x + "/IN.WANNIER_00001.u").c_str());
			system(("qsub " + temp_despath + "/" + x + "/aaa.pbs").c_str());
			cout << "VB:" << this->real_name << ":" << x << " :_p start!" << endl;
		}

		chdir(nowpath.c_str());
		com = "mkdir -p VB_wkm/" + metal_name + "_p";
		system(com.c_str());
		temp_despath = nowpath + "/VB_wkm/" + metal_name + "_p";
		for (string x : test) {
			chdir((temp_despath).c_str());
			system(("mkdir " + x).c_str());
			prepare_singel_file(temp_despath + "/" + x, nowpath, source_path, "wkm_vb", x);
			system(("cp " + nowpath + "/VB/IN.WANNIER_1 " + temp_despath + "/" + x + "/IN.WANNIER_00001.u").c_str());
			system(("qsub " + temp_despath + "/" + x + "/aaa.pbs").c_str());
			cout << "VB:" << this->real_name << ":" << x << " :_p start!" << endl;
		}

	}
	else if (this->source_name.find("Zn") != string::npos) {
		com = "mkdir -p VB_wkm/" + main_group_name + "_p";
		system(com.c_str());
		temp_despath = nowpath + "/VB_wkm/" + main_group_name + "_p";
		for (string x : test) {
			chdir((temp_despath).c_str());
			system(("mkdir " + x).c_str());
			prepare_singel_file(temp_despath + "/" + x, nowpath, source_path, "wkm_vb", x);
			system(("cp " + nowpath + "/VB/IN.WANNIER_0001 " + temp_despath + "/" + x + "/IN.WANNIER_00001.u").c_str());
			system(("qsub " + temp_despath + "/" + x + "/aaa.pbs").c_str());
			cout << "VB:" << this->real_name << ":" << x << " :p start!" << endl;
		}
	}
	else
	{
		com = "mkdir -p VB_wkm/" + main_group_name + "_p " + "VB_wkm/" + main_group_name + "_s";
		system(com.c_str());
		temp_despath = nowpath + "/VB_wkm/" + main_group_name + "_s";
		for (string x : test) {
			chdir((temp_despath).c_str());
			system(("mkdir " + x).c_str());
			prepare_singel_file(temp_despath + "/" + x, nowpath, source_path, "wkm_vb", x);
			system(("cp " + nowpath + "/VB/IN.WANNIER_0001 " + temp_despath + "/" + x + "/IN.WANNIER_00001.u").c_str());
			system(("qsub " + temp_despath + "/" + x + "/aaa.pbs").c_str());
			cout << "VB:" << this->real_name << ":" << x << " :s start!" << endl;
		}
		temp_despath = nowpath + "/VB_wkm/" + main_group_name + "_p";
		for (string x : test) {
			chdir((temp_despath).c_str());
			system(("mkdir " + x).c_str());
			prepare_singel_file(temp_despath + "/" + x, nowpath, source_path, "wkm_vb", x);
			system(("cp " + nowpath + "/VB/IN.WANNIER_0002 " + temp_despath + "/" + x + "/IN.WANNIER_00001.u").c_str());
			system(("qsub " + temp_despath + "/" + x + "/aaa.pbs").c_str());
			cout << "VB:" << this->real_name << ":" << x << " :p start!" << endl;
		}
	}








	//然后处理CB的东西	
	chdir(nowpath.c_str());
	if (source_name.find("Si") != string::npos) {
		chdir(nowpath.c_str());
		com = "mkdir -p CB_wkm/" + main_group_name + "_1s";
		system(com.c_str());
		temp_despath = nowpath + "/CB_wkm/" + main_group_name + "_1s";

		for (string x : test) {
			chdir((temp_despath).c_str());
			system(("mkdir " + x).c_str());
			prepare_singel_file(temp_despath + "/" + x, nowpath, source_path, "wkm_cb", x);
			system(("cp " + nowpath + "/CB/IN.WANNIER_1 " + temp_despath + "/" + x + "/IN.WANNIER_00001.u").c_str());
			system(("qsub " + temp_despath + "/" + x + "/aaa.pbs").c_str());
			cout << "CB:" << this->real_name << ":" << x << " :1p start!" << endl;
		}

		chdir(nowpath.c_str());
		com = "mkdir -p CB_wkm/" + main_group_name + "_2s";
		system(com.c_str());
		temp_despath = nowpath + "/CB_wkm/" + main_group_name + "_2s";

		for (string x : test) {
			chdir((temp_despath).c_str());
			system(("mkdir " + x).c_str());
			prepare_singel_file(temp_despath + "/" + x, nowpath, source_path, "wkm_cb", x);
			system(("cp " + nowpath + "/CB/IN.WANNIER_5 " + temp_despath + "/" + x + "/IN.WANNIER_00001.u").c_str());
			system(("qsub " + temp_despath + "/" + x + "/aaa.pbs").c_str());
			cout << "CB:" << this->real_name << ":" << x << " :2p start!" << endl;
		}

		chdir(nowpath.c_str());
		com = "mkdir -p CB_wkm/" + main_group_name + "_1p";
		system(com.c_str());
		temp_despath = nowpath + "/CB_wkm/" + main_group_name + "_1p";

		for (string x : test) {
			chdir((temp_despath).c_str());
			system(("mkdir " + x).c_str());
			prepare_singel_file(temp_despath + "/" + x, nowpath, source_path, "wkm_cb", x);
			system(("cp " + nowpath + "/CB/IN.WANNIER_2 " + temp_despath + "/" + x + "/IN.WANNIER_00001.u").c_str());
			system(("qsub " + temp_despath + "/" + x + "/aaa.pbs").c_str());
			cout << "CB:" << this->real_name << ":" << x << " :1p start!" << endl;
		}

		chdir(nowpath.c_str());
		com = "mkdir -p CB_wkm/" + main_group_name + "_2p";
		system(com.c_str());
		temp_despath = nowpath + "/CB_wkm/" + main_group_name + "_2p";

		for (string x : test) {
			chdir((temp_despath).c_str());
			system(("mkdir " + x).c_str());
			prepare_singel_file(temp_despath + "/" + x, nowpath, source_path, "wkm_cb", x);
			system(("cp " + nowpath + "/CB/IN.WANNIER_6 " + temp_despath + "/" + x + "/IN.WANNIER_00001.u").c_str());
			system(("qsub " + temp_despath + "/" + x + "/aaa.pbs").c_str());
			cout << "CB:" << this->real_name << ":" << x << " :2p start!" << endl;
		}
	}
	else if (source_name.find("GaAs") != string::npos) {
		chdir(nowpath.c_str());
		com = "mkdir -p CB_wkm/" + main_group_name + "_p";
		system(com.c_str());
		temp_despath = nowpath + "/CB_wkm/" + main_group_name + "_p";
		for (string x : test) {
			chdir((temp_despath).c_str());
			system(("mkdir " + x).c_str());
			prepare_singel_file(temp_despath + "/" + x, nowpath, source_path, "wkm_cb", x);
			system(("cp " + nowpath + "/CB/IN.WANNIER_17 " + temp_despath + "/" + x + "/IN.WANNIER_00001.u").c_str());
			system(("qsub " + temp_despath + "/" + x + "/aaa.pbs").c_str());
			cout << "CB:" << this->real_name << ":" << x << " :p start!" << endl;
		}


		chdir(nowpath.c_str());
		com = "mkdir -p CB_wkm/" + metal_name + "_s ";
		system(com.c_str());
		temp_despath = nowpath + "/CB_wkm/" + metal_name + "_s";

		for (string x : test) {
			chdir((temp_despath).c_str());
			system(("mkdir " + x).c_str());
			prepare_singel_file(temp_despath + "/" + x, nowpath, source_path, "wkm_cb", x);
			system(("cp " + nowpath + "/CB/IN.WANNIER_1 " + temp_despath + "/" + x + "/IN.WANNIER_00001.u").c_str());
			system(("qsub " + temp_despath + "/" + x + "/aaa.pbs").c_str());
			cout << "CB:" << this->real_name << ":" << x << " :p start!" << endl;
		}

		chdir(nowpath.c_str());
		com = "mkdir -p CB_wkm/" + metal_name + "_p ";
		system(com.c_str());
		temp_despath = nowpath + "/CB_wkm/" + metal_name + "_p";

		for (string x : test) {
			chdir((temp_despath).c_str());
			system(("mkdir " + x).c_str());
			prepare_singel_file(temp_despath + "/" + x, nowpath, source_path, "wkm_cb", x);
			system(("cp " + nowpath + "/CB/IN.WANNIER_2 " + temp_despath + "/" + x + "/IN.WANNIER_00001.u").c_str());
			system(("qsub " + temp_despath + "/" + x + "/aaa.pbs").c_str());
			cout << "CB:" << this->real_name << ":" << x << " :p start!" << endl;
		}
	}





	//首先是公共的东西,在s轨道相同
	/*chdir(nowpath.c_str());
	com = "mkdir -p CB_wkm/" + metal_name + "_s ";
	system(com.c_str());
	temp_despath = nowpath + "/CB_wkm/"  + metal_name+"_s";

	for (string x : test) {
		chdir((temp_despath).c_str());
		system(("mkdir " + x).c_str());
		prepare_singel_file(temp_despath + "/" + x, nowpath, source_path, "wkm_cb", x);
		system(("cp " + nowpath + "/CB/IN.WANNIER_0001 " + temp_despath + "/" + x + "/IN.WANNIER_00001.u").c_str());
		system(("qsub " + temp_despath + "/" + x + "/aaa.pbs").c_str());
		cout << "CB:" << this->real_name << ":" << x << " :s start!" << endl;
	}




	chdir(nowpath.c_str());
	if (source_name.find("Ti") == string::npos) {
		chdir(nowpath.c_str());
		com = "mkdir -p CB_wkm/" + main_group_name + "_s ";
		system(com.c_str());
		temp_despath = nowpath + "/CB_wkm/" + main_group_name + "_s";

		for (string x : test) {
			chdir((temp_despath).c_str());
			system(("mkdir " + x).c_str());
			prepare_singel_file(temp_despath + "/" + x, nowpath, source_path, "wkm_cb", x);
			system(("cp " + nowpath + "/CB/IN.WANNIER_0002 " + temp_despath + "/" + x + "/IN.WANNIER_00001.u").c_str());
			system(("qsub " + temp_despath + "/" + x + "/aaa.pbs").c_str());
			cout << "CB:" << this->real_name << ":" << x << " :s start!" << endl;
		}
		chdir(nowpath.c_str());
		com = "mkdir -p CB_wkm/" + main_group_name + "_p ";
		system(com.c_str());
		temp_despath = nowpath + "/CB_wkm/" + main_group_name + "_p";

		for (string x : test) {
			chdir((temp_despath).c_str());
			system(("mkdir " + x).c_str());
			prepare_singel_file(temp_despath + "/" + x, nowpath, source_path, "wkm_cb", x);
			system(("cp " + nowpath + "/CB/IN.WANNIER_0003 " + temp_despath + "/" + x + "/IN.WANNIER_00001.u").c_str());
			system(("qsub " + temp_despath + "/" + x + "/aaa.pbs").c_str());
			cout << "CB:" << this->real_name << ":" << x << " :p start!" << endl;
		}

	}
	else
	{
		com = "mkdir -p CB_wkm/"  + metal_name + "_dz2 "+"CB_wkm/"+metal_name+"_dxz " + "CB_wkm/" + metal_name + "_dx2-y2";
		system(com.c_str());

		temp_despath = nowpath + "/CB_wkm/" + metal_name + "_dz2";

		for (string x : test) {
			chdir(temp_despath.c_str());
			system(("mkdir " + x).c_str());
			prepare_singel_file(temp_despath + "/" + x, nowpath, source_path, "wkm_cb", x);
			system(("cp " + nowpath + "/CB/IN.WANNIER_0002 " + temp_despath + "/" + x + "/IN.WANNIER_00001.u").c_str());
			system(("qsub " + temp_despath + "/" + x + "/aaa.pbs").c_str());
			cout << "CB:" << this->real_name << ":" << x << " :d start!" << endl;

		}
		temp_despath = nowpath + "/CB_wkm/" + metal_name + "_dxz";
		for (string x : test) {
			chdir(temp_despath.c_str());
			system(("mkdir " + x).c_str());
			prepare_singel_file(temp_despath + "/" + x, nowpath, source_path, "wkm_cb", x);
			system(("cp " + nowpath + "/CB/IN.WANNIER_0003 " + temp_despath + "/" +x + "/IN.WANNIER_00001.u").c_str());
			system(("qsub " + temp_despath + "/" + x + "/aaa.pbs").c_str());
			cout << "CB:" << this->real_name << ":" << x << " :d start!" << endl;

		}

		temp_despath = nowpath + "/CB_wkm/" + metal_name + "_dx2-y2";
		for (string x : test) {
			chdir(temp_despath.c_str());
			system(("mkdir " + x).c_str());
			prepare_singel_file(temp_despath + "/" + x, nowpath, source_path, "wkm_cb", x);
			system(("cp " + nowpath + "/CB/IN.WANNIER_0005 " + temp_despath + "/" + x + "/IN.WANNIER_00001.u").c_str());
			system(("qsub " + temp_despath + "/" + x + "/aaa.pbs").c_str());
			cout << "CB:" << this->real_name << ":" << x << " :d start!" << endl;

		}

	}*/
	return;


}
int Wannier_cell::find_close_kppt(double* kpoint, report& rp)
{
	//根据输入的k点坐标判断最接近的k点的序号
	double cha = 1000;
	cout << "search " << kpoint[0] << " " << kpoint[1] << " " << kpoint[2] << endl;
	int k = 1000;
	for (int i = 0; i < rp.num_kpt; i++) {
		//cout << rp.kp_p[i][0] << " " << rp.kp_p[i][1] << " " << rp.kp_p[i][2] << endl;
		double temp = abs(rp.kp_p[i][0] - kpoint[0]) + abs(rp.kp_p[i][1] - kpoint[1]) + abs(rp.kp_p[i][2] - kpoint[2]);
		if (temp < cha) {
			cha = temp;
			k = i;
		}
	}
	cout << "most close k  " << rp.kp_p[k][0] << " " << rp.kp_p[k][1] << " " << rp.kp_p[k][2] << endl;
	return k;

}
int Wannier_cell::find_close_band(double energy, int kppint, report&rp)
{
	//根据输入的能量，在指定k点上寻找对应的band数目
	cout << "search energy:" << energy << endl;
	double cha = 1000;
	double temp = 0;
	int band = 1000;
	double org = 0;
	for (int i = 0; i < rp.num_band; i++) {

		temp = abs(rp.eig_kns[kppint][i][rp.spin - 1] - energy);
		//cout << rp.eig_kns[kppint][i][rp.spin-1] << endl;
		if (temp < cha) {
			cha = temp;
			band = i;
			org = rp.eig_kns[kppint][i][rp.spin - 1];
		}
	}
	cout << "org " << org << endl;
	return band;
}
void Wannier_cell::get_source_dependent() {
	//获得依赖
	string com = "source /share/home/wangz/pwmat.intel.cuda.env";
	system(com.c_str());
	com = "source /share/home/wangz/PWMAT/20200414/pwmat.env";
	system(com.c_str());
	com = "source ~/.bashrc";
	system(com.c_str());
	return;
}
int Wannier_cell::get_start_line(string flag) {
	//读取一个文件
	ifstream fin;
	int start = 0;
	char buffer[200];
	fin.open("OUT.KPT", ios::in);
	if (!fin.is_open()) {
		cout << "i can not find the file :wan90!" << endl;
		cin.get();
	}
	while (fin.good() && fin.peek() != EOF) {
		fin.getline(buffer, 200);
		start++;
		if (strstr(buffer, flag.c_str()) != NULL) {
			break;
		}
	}
	fin.close();
	return start;
}
void Wannier_cell::get_source_name() {
	//读取元素列表，确定当前的来源是什么

	this->source_name = "SnBr2";
	return;

	bool falg_s = false, flag_f = false, flag_mg = false;
	bool flag_si = false, flag_as = false;
	bool flag_snbr = false;
	for (int i = 0; i < this->type_num; i++) {

		//补充上测试的两个结构
		if (this->type[i] == 50 || this->type[i] == 32)
			flag_snbr = true;
		else if (this->type[i] == 14)
			flag_si = true;
		else if (this->type[i] == 31)
			flag_as = true;
		else if (this->type_save[i] == 16 || this->type_save[i] == 34 || this->type_save[i] == 52)
			falg_s = true;
		else if (this->type_save[i] == 9 || this->type_save[i] == 17 || this->type_save[i] == 35 || this->type_save[i] == 53)
			flag_f = true;
		else if (this->type_save[i] == 12 || this->type_save[i] == 20 || this->type_save[i] == 38 || this->type_save[i] == 56)
			flag_mg = true;


	}
	//然后开始断言
	if (flag_snbr == true)
	{
		this->source_name = "SnBr2";
		return;
	}
	else if (flag_si == true) {
		this->source_name = "Si";
		return;
	}
	else if (flag_as == true) {
		this->source_name = "GaAs";
		return;
	}
	else if (falg_s == true) {
		this->source_name = "TiS2";
		return;
	}
	else if (flag_f == true && flag_mg == true) {
		this->source_name = "MgF";
		return;
	}
	else
	{
		this->source_name = "ZnF";
		return;
	}
}





void detect_spin_eualone(string & path, string & file_name, string & out_name, int** flag, double rule)
{
	string comm = "grep 'diff' REPORT_test |tail -1";
	string check_comm = "ls |grep 'REPORT_test' |wc -l";
	char temp[100];
	char buf[100];
	char result[4096];
	FILE* fp = NULL;
	ifstream fin;
	ofstream fout;
	fout.open(out_name, ios::out);
	fin.open(file_name, ios::in);
	if (!fin.is_open())
	{
		cout << "cann not find the file " << file_name << endl;
		cin.get();
	}
	int i = 0;
	while (fin.good() && fin.peek() != EOF)
	{
		fin.getline(temp, 100);
		//cout << temp << endl;
		if (temp == "" || temp[0]<'0' || temp[0]>'9')
			break;
		if (flag[i++][0] != 99)
			continue;

		chdir((path + string(temp)).c_str());
		//首先需要确认有这个文件
		fp = popen(check_comm.c_str(), "r");
		fgets(buf, sizeof(buf) - 1, fp);
		if (buf[0] == '0')
		{
			continue;
			pclose(fp);
		}
		pclose(fp);


		fp = popen(comm.c_str(), "r");
		while (fgets(buf, sizeof(buf) - 1, fp) != NULL)
		{
			strcat(result, buf);
		}
		pclose(fp);
		fp = NULL;
		//然后是从字符串中读出最后一个数字

		string res = string(result);
		//cout << result << endl;
		int loc = res.find_last_of(" ");
		string sub = res.substr(loc, res.size() - 1);
		//cout << sub << endl;
		double cha = atof(sub.c_str());
		//cout << cha << endl;
		//cin.get();

		if (cha < rule)
		{
			fout << temp << endl;
			//cout << temp <<endl;
			flag[i - 1][0] = 0;
		}
		else
		{
			flag[i - 1][0] = 87;
		}
		memset(temp, 0, sizeof(temp));
		memset(result, 0, sizeof(result));
		memset(buf, 0, sizeof(buf));

	}
	fin.close();
	fout.close();
	return;
}








bool cell::check_contain_water()
{
	//输入是一个结构，输出是判断该结构有无存在水分子
	double rule = 1.426;
	double distance = 0;
	int bond_num = 0;
	for (int i = 0; i < num; i++)
	{
		if (this->type[i] == 8)
		{
			bond_num = 0;
			for (int j = 0; j < yanshen*num; j++)
			{
				if (type[j%num] == 1)
				{
					distance = dis(real_position[13][i], real_position[j / num][j%num]);
					if (distance<1.426 && distance>1e-5)
					{
						bond_num++;
					}
				}
			}
			if (bond_num == 2)
			{
				return true;
			}
		}
	}
	return false;
}

void full_fill_yanshi(string yanshi_name, string data_path)
{
	int i = 0, j = 0;
	//读取相应的堰势文件,缺少的从path中寻找补充进来
	ifstream fin;
	fin.open(yanshi_name, ios::in);
	if (!fin.is_open())
	{
		cout << "i can not find the file " << yanshi_name << endl;
		cin.get();
	}
	string temp;
	int yanshi[120] = { 0 };
	i = 0;
	while (fin.good() && fin.peek() != EOF)
	{
		fin >> temp;
		for (i = 0; i < 120; i++)
		{
			if (strcmp(a[i], temp.c_str()) == 0)
			{
				yanshi[i] = 1;
				cout << temp << "  " << a[i] << endl;
				//cin.get();
				break;
			}

		}
		if (i == 120)
		{
			cout << "ukonwn element!" << temp << endl;
			cin.get();
		}
	}
	fin.close();
	string com;
	for (i = 1; i < 112; i++)
	{
		if (yanshi[i] == 0 && a[i][0] != ' ')
		{
			//说明没有，需要从path中复制
			com = "cp " + data_path + a[i] + "-sp.PD04.PBE.UPF " + "/share/home/wangz/NCPP-SG15-PBE/" + a[i] + ".SG15.PBE.UPF";
			//cout << com << endl;
			system(com.c_str());
			com.clear();

			com = "cp " + data_path + a[i] + "-d.PD04.PBE.UPF " + "/share/home/wangz/NCPP-SG15-PBE/" + a[i] + ".SG15.PBE.UPF";
			//cout << com << endl;
			system(com.c_str());
			com.clear();
		}
	}
	return;
}


int  backup_status(string config_name, string &name1, string& name2, int**flag, int TASK_NUM, int** mp, int line_number)
{
	//建立一个重载版本
	int i = 0, j = 0;
	string com;
	ifstream fin;
	//第一步读取名字
	vector<string> all_name = get_file_name(config_name);


	//第二步读取status文件
	fin.open(name1, ios::in);
	if (!fin.is_open())
	{
		cout << "i can not find the backup file!" << endl;
		cin.get();
		return 1;
	}
	j = 0;
	char temp[100];
	for (i = 0; i < TASK_NUM && fin.peek() != EOF; i++)
	{
		fin.get(temp, 89, '\t');
		fin >> flag[i][0] >> flag[i][1] >> flag[i][2];
		for (j = 0; j < 3; j++)
		{
			fin >> mp[i][j];
		}
		fin.getline(temp, 20);

	}
	fin.close();
	/*cout << temp << " " << flag[i-1][0] << flag[i-1][1] << flag[i-1][2] << mp[i-1][0] << mp[i-1][1] << mp[i-1][2] << endl;
	cin.get();*/
	//最后，根据log_all进行补充
	if (line_number == -1)
		return 0;
	fin.open(name2, ios::in);
	if (!fin.is_open())
	{
		cout << "i can not find the backup file!" << endl;
		cin.get();
		return 1;

	}
	char line[100];
	char temp_name[30];

	//需要注意我们只需要读取最后几行的信息就行
	for (i = 0; i < line_number - 1; i++)
	{
		fin.getline(line, 99);
	}

	int push = 0;
	int index = 0;
	while (fin.good() && fin.peek() != EOF)
	{
		fin.getline(line, 90);
		cout << line << endl;
		cin.get();
		for (j = 0; j < strlen(line); j++)
		{
			if (line[j] == ':' || line[j] == ' ')
			{
				temp_name[push] = '\0';
				break;
			}
			temp_name[push++] = line[j];
		}
		push = 0;
		//然后根据不同的情况给出flag值
		for (i = 0; i < TASK_NUM; i++)
		{
			if (all_name[i] == temp_name)
			{
				break;
			}
		}
		/*cout << temp_name << endl;
		cout << all_name[i] << endl;
		cin.get();*/
		if (i == TASK_NUM)
		{
			cout << "can not find the name:" << temp_name << endl;
			cin.get();
		}
		else
		{
			//这里注意，71和61的二义性，需要进行检测
			if (strstr(line, "relax finished!") != NULL)
			{
				flag[i][0] = 2;
				//cout << line << endl;
				//cin.get();
			}
			else if (strstr(line, "nonscf do not have band gap") != NULL)
			{
				flag[i][0] = 20;
				//cout << line << endl;
				//cout << flag[i][0] << endl;
				//cin.get();
			}
			else if (strstr(line, "scf do not have band gap") != NULL)
			{
				flag[i][0] = 47;
				//cout << line << endl;
				//cout << flag[i][0] << endl;
				//cin.get();
			}
			else if (strstr(line, "nonscf 1st finished") != NULL)
			{
				flag[i][0] = 51;
				/*cout << line << endl;
				cout << flag[i][0] << endl;
				cin.get();*/
			}

			else if (strstr(line, "nonscf 1st error") != NULL)
			{
				flag[i][0] = 57;
				/*cout << line << endl;
				cout << flag[i][0] << endl;
				cin.get();*/
			}
			else if (strstr(line, "nonscf wait!") != NULL)
			{
				flag[i][0] = 61;
				flag[i][2]++;
				mp[i][0] = 2 * mp[i][0];
				mp[i][1] = 2 * mp[i][1];
				mp[i][2] = 1;
				//cout << line[strlen(line) - 1] << endl;
				//cin.get();

			}
			else if (strstr(line, "1stnonscf start!") != NULL)
			{
				flag[i][0] = 5;
				//cout << line[strlen(line) - 1] << endl;
				//cin.get();

			}

			else if (strstr(line, "nonscf  finished") != NULL)
			{
				flag[i][0] = 71;
				/*int max = (mp[i][0]> mp[i][1] ? mp[i][0] : mp[i][1]);
				if (max >= RULE_NUM)
				{
					flag[i][0] = 99;
				}*/
				/*cout << line << endl;
				cout << flag[i][0] << endl;
				cin.get();*/
			}
			else if (strstr(line, "nonscf start!") != NULL)
			{
				flag[i][0] = 71;
				flag[i][2]++;
				//每次读取这个mp需要更新
				mp[i][0] = 2 * mp[i][0];
				mp[i][1] = 2 * mp[i][1];
				mp[i][2] = 1;

			}

			else if (strstr(line, "two_cut method") != NULL)
			{
				flag[i][0] = 30;
				/*cout << line << endl;
				cout << flag[i][0] << endl;
				cin.get();*/
			}

			else if (strstr(line, "nonscf  error") != NULL)
			{
				flag[i][0] = 67;
				/*cout << line << endl;
				cout << flag[i][0] << endl;
				cin.get();*/
			}
			else if (strstr(line, "scf start!") != NULL)
			{
				flag[i][0] = 3;
			}
			else if (strstr(line, "scf finished!") != NULL)
			{
				flag[i][0] = 4;
			}
			else if (strstr(line, "scf error") != NULL)
			{
				flag[i][0] = 37;
			}
			else if (strstr(line, "succeed") != NULL)
			{
				flag[i][0] = 99;
			}
			else
			{
				cout << "unknown the flag :" << line << endl;
				cin.get();
			}
		}
		/*cout << line << " " << flag[i][0] << endl;*/
	}




	return 0;
}





int classify_metal_maingroup(int &atomic_number)//0 rare gas,1 metal ,2 main group
{

	int j = 0, k = 0, m = 0;
	while (j < main_groupnum || k < metal_num || m < rare_gasnum)
	{
		if (rare_gas[m] == atomic_number)
		{
			return 0;
		}
		if (Meatal_xuhao[k] == atomic_number)
		{
			return 1;
		}
		if (main_group_element[j] == atomic_number)
		{
			return 2;
		}


		if (Meatal_xuhao[k] < atomic_number)
		{
			k++;
			continue;
		}
		if (rare_gas[m] < atomic_number)
		{
			m++;
			continue;
		}
		if (main_group_element[j] < atomic_number)
		{
			j++;
			continue;
		}
	}
	cout << "what's wrong!" << atomic_number << endl;
	cin.get();
	return -1;

}

string  cell::get_formula_name()
{
	int i = 0;
	int flag[120] = { 0 };
	for (i = 0; i < num; i++)
	{
		flag[type[i]]++;
	}
	string name;
	for (i = 0; i < 120; i++)
	{
		if (flag[i] != 0)
		{
			name += a[i];
			name += to_string(flag[i]);
		}

	}
	//cout << name << endl;
	//cin.get();
	return name;
}






//void cell::more_k_toget_band(string path)
//



void cell::more_k_toget_band(string path) {
	//产生更多的高对称点来画出良好的能带图
	chdir(path.c_str());
	//首先得到cbm和vbm的那条线的5个点k坐标
	//因为偏差的都是-45度，所以相差6就完事啦
	//K	0.5052386279	0.8750989735	0.0000000000
	//K'	-0.5052386279	-0.8750989735	-0.0000000000
	//	L	0.7578579418	0.4375494867	0.1306455197
	//	L'	-0.7578579418	-0.4375494867	-0.1306455197
	//	M	0.7578579418	0.4375494867	0.0000000000
	vector<vector<double>> vbm, cbm;
	vector<int> vbm_k, cbm_k;
	report rp;
	rp.read_report("REPORT", 1);
	//填充相关的序号，先是vbm,然后才是vbm
	for (int i = 0; i < 25; i++) {
		if (abs(rp.vbm_kpt - i) % 6 == 0) {
			//cout << i << endl;
			vbm_k.push_back(i);
		}
	}
	for (int i = 25; i < 50; i++) {
		if (abs(rp.cbm_kpt - i) % 6 == 0) {
			//cout << i << endl;
			cbm_k.push_back(i);
		}
	}

	ifstream fin;
	char buff[200];
	fin.open("IN.KPT", ios::in);
	fin.getline(buff, 199);
	fin.getline(buff, 199);
	int xuhao = 0;
	vector<vector<double>>all_k;
	while (fin.good() && fin.peek() != EOF) {
		vector<double> temp(3);
		fin >> temp[0];
		fin >> temp[1];
		fin >> temp[2];
		fin.getline(buff, 199);
		all_k.push_back(temp);
	}
	fin.close();
	//然后就是根据序号进行扩展了
	for (int xuhao : vbm_k)
		vbm.push_back(all_k[xuhao]);
	for (int xuhao : cbm_k)
		cbm.push_back(all_k[xuhao]);



	vector<vector<double>> vbm_plus, cbm_plus;
	double rule = 1.1;//标准是1.0244这里长了一点
	double every_step = 0.1;
	int cycle_num = 1;
	double now_rule = pow(vbm[0][0] * vbm[0][0] + vbm[0][1] * vbm[0][1], 0.5);
	while (now_rule < rule) {
		vector<double> plus = get_onemore_point(vbm[0], vbm[1], every_step*cycle_num);
		vbm_plus.push_back(plus);
		cycle_num++;
		now_rule = pow(plus[0] * plus[0] + plus[1] * plus[1], 0.5);
	}
	cycle_num = 1;
	now_rule = pow(cbm[0][0] * cbm[0][0] + cbm[0][1] * cbm[0][1], 0.5);
	while (now_rule < rule) {
		vector<double> plus = get_onemore_point(cbm[0], cbm[1], every_step*cycle_num);
		cbm_plus.push_back(plus);
		cycle_num++;
		now_rule = pow(plus[0] * plus[0] + plus[1] * plus[1], 0.5);
	}

	//然后开始构建IN.KPT文件
	system("mv IN.KPT IN.KPT_old");
	system("mv REPORT REPORT_last_non");
	ofstream fout;
	fout.open("IN.KPT", ios::out);
	int all_num = vbm.size() + cbm.size() + vbm_plus.size() + cbm_plus.size();
	fout << "     " << all_num << endl;
	fout << "    2    0" << endl;
	fout.setf(ios::fixed);
	for (vector<double> temp : vbm) {
		fout << setprecision(6) << "    " << temp[0] << "  " << temp[1] << "  " << temp[2] << "  " << "1.000000" << endl;
	}
	for (vector<double> temp : cbm) {
		fout << setprecision(6) << "    " << temp[0] << "  " << temp[1] << "  " << temp[2] << "  " << "1.000000" << endl;
	}
	for (vector<double> temp : vbm_plus) {
		fout << setprecision(6) << "    " << temp[0] << "  " << temp[1] << "  " << temp[2] << "  " << "0.000000" << endl;
	}
	for (vector<double> temp : cbm_plus) {
		fout << setprecision(6) << "    " << temp[0] << "  " << temp[1] << "  " << temp[2] << "  " << "0.000000" << endl;
	}

	fout.close();
	system("qsub aaa.pbs");

	return;
}

vector<double>  cell::get_onemore_point(vector<double> a, vector<double>b, double step) {
	//根据输入的步长和两个点，产生一个新的共线的点
	double x = abs(a[0] - b[0]);
	double y = abs(a[1] - b[1]);
	vector<double> res(3);
	//我们要向大的增
	if (a[0] > b[0]) {
		res[0] = a[0] + step * x;
		res[1] = a[1] + step * y;
		res[2] = 0.0;
	}
	else {
		res[0] = b[0] + step * x;
		res[1] = b[1] + step * y;
		res[2] = 0.0;
	}
	return res;

}


void cell::output(string file_name) {
	//简单的输出文件
	ofstream fout;
	fout.open(file_name, ios::out);
	fout << "  " << this->num << endl;
	fout << "  " << "Lattice vector" << endl;
	fout.setf(ios::fixed);
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			fout << setprecision(9) << this->letice[i][j] << "    ";
		}
		fout << endl;
	}
	fout << "  " << "Position" << endl;
	for (int i = 0; i < this->num; i++) {
		fout << "  " << this->type[i] << "    " << setprecision(12) << this->p[i][0] << "    " << this->p[i][1] << "    " << this->p[i][2] << "    " << "1    1    1" << endl;
	}

	fout.close();
	return;
}



void cell::output_config(string path, string file_name) {
	//用来输出文件的
	ofstream fout;
	fout.open(path + file_name, ios::out);
	fout << "    " << this->num << endl;
	fout << " Lattice vector\n";

	int temp = 0;
	fout.setf(ios::fixed);
	//下面开始输出
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			fout << "    " << this->letice[i][j];
		}
		fout << endl;
	}
	fout << " Position\n";
	for (int i = 0; i < this->num; i++) {
		fout << "    " << this->type[i] << "    " << this->p[i][0] << "     " << this->p[i][1] << "    " << this->p[i][2] << "    " << "1    1    1" << endl;
	}
	fout.close();
	return;
}



void cell::supper_cell(int x, int y, int z, string path) {
	//对晶胞进行扩胞，产生结果的config文件
	ofstream fout;
	fout.open(path + "atom_super.config", ios::out);
	fout << "    " << this->num*x*y*z << endl;
	fout << " LATTICE\n";

	int temp = 0;
	fout.setf(ios::fixed);
	for (int i = 0; i < 3; i++)
	{
		if (i == 0)
			temp = x;
		else if (i == 1)
			temp = y;
		else
			temp = z;
		for (int j = 0; j < 3; j++)
		{
			fout << "   " << setprecision(9) << letice[i][j] * temp << "    ";
		}
		fout << endl;
	}
	fout << " POSITION\n";



	for (int i = 0; i < num; i++)
	{
		for (int ii = 0; ii < x; ii++)
		{
			for (int jj = 0; jj < y; jj++)
			{
				for (int kk = 0; kk < z; kk++)
				{
					fout << "    " << type[i] << "    " << p[i][0] / (double)x + ii * 1.0 / (double)x << "    " << p[i][1] / (double)y + jj * 1.0 / (double)y << "    " << p[i][2] / (double)z + kk * 1.0 / (double)z << "    1  1  1" << endl;
				}
			}
		}
	}
	fout.close();
	return;
}

void cell::supper_cell_md(string path) {
	//读取路径下面的结构，为md进行扩胞
	if (CALCULATION_MODE == 1) {
		double realLength = pow(pow(this->letice[0][0], 2) + pow(this->letice[0][1], 2) + pow(this->letice[0][2], 2), 0.5);
		int times = ceil(MD_LENGTH / realLength);
		//在之前修改源文件
		chdir(path.c_str());
		system("mv atom.config atom_org.config");
		supper_cell(times,1,1,path);
		//修改后改回atom.config
		system("mv atom_super.config atom.config");
		return;
	}
	else {
		cout << "unsupported calculation mode!" << endl;
		cin.get();
		return;
	}
}
void pbs_got(string name, string flag, string dir_path, string lines)
{
	//将对应的aaa.pbs复制到对应路径
	//根据输入的falg和队列名字确定aaa.pbs
	//脚本同目录应事先准备relax.pbs,scf.pbs,nonscf.pbs
	string comm;
	chdir("/share/home/wangz/high");

	if (flag == "relax")
	{
		comm = "cp relax.pbs " + string("aaa.pbs");
		system(comm.c_str());
		comm = "sed -i s/CaI2/" + string("relax_") + name + "/g aaa.pbs";
		system(comm.c_str());

	}
	else if (flag == "scf")
	{
		if (IF_WKM == false) {
			if (IF_MD == true) {
				comm = "cp md.pbs " + string("aaa.pbs");
			}
			else if (IF_HSE == true)
			{
				comm = "cp scf_hse.pbs " + string("aaa.pbs");
			}
			else
				comm = "cp scf.pbs " + string("aaa.pbs");
		}
		else
			comm = "cp scf_wkm.pbs " + string("aaa.pbs");

		system(comm.c_str());
		if(IF_MD==false)
			comm = "sed -i s/CaI2/" + string("scf_") + name + "/g aaa.pbs";
		else
			comm = "sed -i s/CaI2/" + string("md_") + name + "/g aaa.pbs";		
		system(comm.c_str());
	}
	else if (flag == "nonscf")
	{
		if (IF_HSE == true)
		{
			comm = "cp nonscf.pbs " + string("aaa.pbs");
		}
		else
			comm = "cp nonscf.pbs " + string("aaa.pbs");
		system(comm.c_str());
		comm = "sed -i s/CaI2/" + string("nonscf_") + name + "/g aaa.pbs";
		system(comm.c_str());
	}
	else if (flag == "dos")
	{
		if (IF_HSE == true)
		{
			comm = "cp dos.pbs " + string("aaa.pbs");
		}
		else
			comm = "cp dos.pbs " + string("aaa.pbs");
		system(comm.c_str());
		comm = "sed -i s/CaI2/" + string("dos_") + name + "/g aaa.pbs";
		system(comm.c_str());
	}
	else
	{
		cout << "unkonwn flags!please check!" << endl;
		cin.get();
		return;
	}
	comm = "sed -i s/gpu2/" + lines + "/g aaa.pbs";
	system(comm.c_str());
	comm = "mv aaa.pbs " + dir_path + "/" + string("aaa.pbs");
	system(comm.c_str());
	chdir(dir_path.c_str());
	return;
}

void give_mp123(string temp_name, int *mp)
{
	char str2[500];
	FILE *atom;
	double d[3][3];
	double dd[3];
	int *type;
	int num;
	int i, j;
	int t[20];
	int kp[3];
	atom = fopen((temp_name + "/atom.config").c_str(), "r");
	if (atom == NULL)
	{
		printf("ERROR opening :%s\n", (temp_name + "/atom.config").c_str());
		cin.get();
		return;
	}
	fscanf(atom, "%d", &num);
	fgets(str2, 500, atom);
	fgets(str2, 200, atom);
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			fscanf(atom, "%lf", &d[i][j]);
		}
		fgets(str2, 200, atom);
	}
	//fgets(str2, 200, atom);
	fgets(str2, 200, atom);
	type = (int *)malloc(num * sizeof(int));
	for (i = 0; i < num; i++)
	{
		fscanf(atom, "%d", &type[i]);
		fgets(str2, 200, atom);
	}
	for (i = 0; i < 20; i++)
	{
		t[i] = 0;
	}

	//然后产生mp123
	for (i = 0; i < 3; i++)
	{
		dd[i] = 0;
		for (j = 0; j < 3; j++)
		{
			dd[i] = d[i][j] * d[i][j] + dd[i];
		}
		dd[i] = sqrt(dd[i]);
	}
	if (CALCULATION_MODE == 0)
	{
		kp[0] = kp[1] = kp[2] = 1;
	}
	else if (CALCULATION_MODE == 2)
	{
		kp[0] = pow(double(K_rule) / num * dd[0] * dd[1], 1.0 / 2.0) / dd[0] + 0.5;
		kp[1] = pow(double(K_rule) / num * dd[0] * dd[1], 1.0 / 2.0) / dd[1] + 0.5;
		kp[2] = 1;
	}
	else if (CALCULATION_MODE == 1)
	{
		kp[0] = pow(double(K_rule) / num * dd[0] * dd[1], 1.0 / 2.0) / dd[0] + 0.5;
		kp[1] = kp[2] = 1;
	}
	else if (CALCULATION_MODE == 3) {
		kp[0] = pow(double(K_rule) / num * dd[0] * dd[1], 1.0 / 2.0) / dd[0] + 0.5;
		kp[1] = pow(double(K_rule) / num * dd[0] * dd[1], 1.0 / 2.0) / dd[1] + 0.5;
		kp[2] = pow(1000.0 / num * dd[0] * dd[1] * dd[2], 1.0 / 3.0) / dd[2] + 0.5;
	}
	else {
		cout << "unknwon calculation mode number!" << endl;
		cin.get();
	}
	for (i = 0; i < 3; i++)
	{
		if (kp[i] == 0)
		{
			kp[i] = 1;
		}
	}
	for (i = 0; i < 3; i++)
	{
		mp[i] = kp[i];
	}

	fclose(atom);
	return;
}
int generate_jos(string name, string data_path, string dir_path, int NUM, string flag, int*yanshi_flag, int** task_flag, ofstream & fout)
{
	//用来得到计算的脚本
	//输入是结构名字，以及脚本目录下面的对应结构文件
	//输出是在指定路径下面，产生对应文件夹，并且根据flag产生对应的脚本文件
	//目前针对的是2d结构
	//循环读入名字，进行一些列操作

	char str[200];
	char str1[200];
	char str2[500];
	FILE *in;
	FILE *atom;
	double d[3][3];
	double dd[3];
	int *type;
	int num;
	int i, j;
	int t[20];
	int n_type;
	int flag_type;
	FILE *etot;
	int kp[3];
	int count = 0;
	ifstream fin;
	string temp_name;
	fin.open(name, ios::in);
	if (!fin.is_open())
	{
		cout << "i can not find the file!" << name << endl;
		cin.get();
		return 0;
	}
	string com;
	while (fin.good() && fin.peek() != EOF && (count++ < NUM))
	{
		//建立文件夹以及atom.config
		fin >> temp_name;
		cout << temp_name << endl;
		com = "mkdir " + dir_path + temp_name;
		system(com.c_str());
		com.clear();
		com = "cp " + data_path + temp_name + " " + dir_path + temp_name + "/atom.config";
		system(com.c_str());
		//开始读取结构文件，产生脚本
		//atom = fopen((temp_name+".config").c_str(), "r");
		atom = fopen((data_path + temp_name).c_str(), "r");
		if (atom == NULL)
		{
			printf("ERROR opening %s\n", str1);
			return 1;
		}
		fscanf(atom, "%d", &num);
		fgets(str2, 500, atom);
		fgets(str2, 200, atom);
		for (i = 0; i < 3; i++)
		{
			for (j = 0; j < 3; j++)
			{
				fscanf(atom, "%lf", &d[i][j]);
			}
			fgets(str2, 200, atom);
		}
		//fgets(str2, 200, atom);
		fgets(str2, 200, atom);
		type = (int *)malloc(num * sizeof(int));
		for (i = 0; i < num; i++)
		{
			fscanf(atom, "%d", &type[i]);
			fgets(str2, 200, atom);
		}
		for (i = 0; i < 20; i++)
		{
			t[i] = 0;
		}
		n_type = 0;
		for (i = 0; i < num; i++)
		{
			flag_type = 0;
			for (j = 0; j < n_type; j++)
			{
				if (type[i] == t[j])
				{
					flag_type = 1;
				}
			}
			if (flag_type == 0)
			{
				t[n_type] = type[i];
				n_type++;
			}
		}
		//printf("n_type=%d\n", n_type);


		//复制堰势文件到指定目录下
		if (Soc_on == false) {
			for (i = 0; i < n_type; i++)
			{
				if (yanshi_flag[t[i]] != 1)
				{
					fout << temp_name << ": lack yanshi file!" << endl;
					cout << temp_name << ": lack yanshi file!" << endl;
					//cin.get();
					task_flag[count - 1][0] = 7;
					//cin.get();
					break;
				}
				com = "cp /share/home/wangz/NCPP-SG15-PBE/" + string(a[t[i]]) + ".SG15.PBE.UPF" + " " + dir_path + temp_name;
				//cout << com << endl;
				system(com.c_str());

			}
		}
		else {
			for (i = 0; i < n_type; i++)
			{
				if (yanshi_flag[t[i]] != 1)
				{
					fout << temp_name << ": lack yanshi file!" << endl;
					cout << temp_name << ": lack yanshi file!" << endl;
					//cin.get();
					task_flag[count - 1][0] = 7;
					//cin.get();
					break;
				}
				com = "cp /share/home/wangz/SG15.SOC/NCPP-SG15-PBE-SOC/SG15.PBE.SOC/" + string(a[t[i]]) + ".SG15.PBE.SOC.UPF" + " " + dir_path + temp_name;
				//cout << com << endl;
				system(com.c_str());

			}
		}
		//cin.get();
		//开始写etot.input

		fclose(atom);

		//下面根据不同的需要填写不同的内容
		if (flag == "relax")
		{
			fout << "    MP_N123  =  " << kp[0] << " " << kp[1] << " " << kp[2] << " " << "0 0 0" << endl;
			fout << "    XCFUNCTIONAL = PBE\n";
			fout << "    JOB = RELAX\n";
			fout << "    RELAX_DETAIL = 1 100 0.02 1 0.1\n";
			fout << "    ECUT = 50\n";
			fout << "    ECUT2 = 200\n";
			fout << "    SPIN = 2\n";
			fout << "    VDW= DFT-D2\n";


		}
		else if (flag == "scf")
		{
			fout << "    MP_N123  =  " << kp[0] << " " << kp[1] << " " << kp[2] << " " << "0 0 0" << endl;
			fout << "    XCFUNCTIONAL = PBE\n";
			fout << "    JOB = SCF\n";
			//fout << "    RELAX_DETAIL = 1 100 0.02 1 0.1\n";
			fout << "    ECUT = 50\n";
			fout << "    ECUT2 = 200\n";
			fout << "    SPIN = 2\n";

		}
		else if (flag == "nonscf")
		{
			//fout << "    MP_N123  =  " << kp[0] << " " << kp[1] << " " << kp[2] << " " << "0 0 0" << endl;
			fout << "    XCFUNCTIONAL = PBE\n";
			fout << "    JOB = NONSCF\n";
			//fout << "    RELAX_DETAIL = 1 100 0.02 1 0.1\n";
			fout << "    ECUT = 50\n";
			fout << "    ECUT2 = 200\n";
			fout << "    SPIN = 2\n";

			fout << "    CHARGE_DECOMP = T\n";
			fout << "    IN.VR=T\n";
			fout << "    IN.KPT=T\n";
		}
		else if (flag == "generate")
		{
			cout << "just generate the file!!" << endl;

		}
		else
		{
			cout << "unknown clacualtion name!please chekc!" << endl;
			cin.get();
			return 0;
		}
		fout.close();

	}
	fin.close();
	return 0;

}


int if_no_yanshi(string path, int *yanshi_flag)
{
	FILE* atom;
	char str[200];
	char str1[200];
	char str2[500];
	double d[3][3];
	double dd[3];
	int *type;
	int num;
	int i, j;
	int t[20];
	int n_type;
	int flag_type;
	int kp[3];
	int flag = 0;
	atom = fopen((path + "/atom.config").c_str(), "r");
	if (atom == NULL)
	{
		printf("ERROR opening %s\n", (path + "/atom.config").c_str());
		cin.get();
		return 1;
	}
	fscanf(atom, "%d", &num);
	fgets(str2, 500, atom);
	fgets(str2, 200, atom);
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			fscanf(atom, "%lf", &d[i][j]);
		}
		fgets(str2, 200, atom);
	}
	//fgets(str2, 200, atom);
	fgets(str2, 200, atom);
	type = (int *)malloc(num * sizeof(int));
	for (i = 0; i < num; i++)
	{
		fscanf(atom, "%d", &type[i]);
		fgets(str2, 200, atom);
	}
	for (i = 0; i < 20; i++)
	{
		t[i] = 0;
	}
	n_type = 0;
	for (i = 0; i < num; i++)
	{
		flag_type = 0;
		for (j = 0; j < n_type; j++)
		{
			if (type[i] == t[j])
			{
				flag_type = 1;
			}
		}
		if (flag_type == 0)
		{
			t[n_type] = type[i];
			n_type++;
		}
	}
	//printf("n_type=%d\n", n_type);


	//复制堰势文件到指定目录下
	for (i = 0; i < n_type; i++)
	{
		if (yanshi_flag[t[i]] != 1)
		{

			flag = 1;
			break;
		}

	}
	fclose(atom);
	return flag;

}




int generate_jos_temp(string name, string data_path, string dir_path, int NUM, string flag, int*yanshi_flag, int** task_flag, ofstream & fout)
{
	//用来得到计算的脚本
	//输入是结构名字，以及脚本目录下面的对应结构文件
	//输出是在指定路径下面，产生对应文件夹，并且根据flag产生对应的脚本文件
	//目前针对的是2d结构
	//循环读入名字，进行一些列操作

	char str[200];
	char str1[200];
	char str2[500];
	FILE *in;
	FILE *atom;
	double d[3][3];
	double dd[3];
	int *type;
	int num;
	int i, j;
	int t[20];
	int n_type;
	int flag_type;
	FILE *etot;
	int kp[3];
	int count = 0;
	ifstream fin;
	string temp_name;
	fin.open(name, ios::in);
	if (!fin.is_open())
	{
		cout << "i can not find the file!" << name << endl;
		cin.get();
		return 0;
	}

	int cycle_num = 0;
	string com;
	while (fin.good() && fin.peek() != EOF && (count++ < NUM))
	{
		//建立文件夹以及atom.config
		fin >> temp_name;
		cout << temp_name << endl;
		/*com = "mkdir " + dir_path + temp_name;
		system(com.c_str());
		com.clear();*/
		com = "cp " + data_path + temp_name + " " + dir_path + temp_name + "/atom.config";
		system(com.c_str());
		//开始读取结构文件，产生脚本
		//atom = fopen((temp_name+".config").c_str(), "r");
		atom = fopen((data_path + temp_name).c_str(), "r");
		if (atom == NULL)
		{
			printf("ERROR opening %s\n", str1);
			return 1;
		}
		fscanf(atom, "%d", &num);
		fgets(str2, 500, atom);
		fgets(str2, 200, atom);
		for (i = 0; i < 3; i++)
		{
			for (j = 0; j < 3; j++)
			{
				fscanf(atom, "%lf", &d[i][j]);
			}
			fgets(str2, 200, atom);
		}
		//fgets(str2, 200, atom);
		fgets(str2, 200, atom);
		type = (int *)malloc(num * sizeof(int));
		for (i = 0; i < num; i++)
		{
			fscanf(atom, "%d", &type[i]);
			fgets(str2, 200, atom);
		}
		for (i = 0; i < 20; i++)
		{
			t[i] = 0;
		}
		n_type = 0;
		for (i = 0; i < num; i++)
		{
			flag_type = 0;
			for (j = 0; j < n_type; j++)
			{
				if (type[i] == t[j])
				{
					flag_type = 1;
				}
			}
			if (flag_type == 0)
			{
				t[n_type] = type[i];
				n_type++;
			}
		}
		//printf("n_type=%d\n", n_type);


		//复制堰势文件到指定目录下
		for (i = 0; i < n_type; i++)
		{
			if (yanshi_flag[t[i]] != 1)
			{
				fout << temp_name << ": lack yanshi file!" << endl;
				cout << temp_name << ": lack yanshi file!" << endl;
				//cin.get();
				task_flag[cycle_num][0] = 7;
				break;
			}
			com = "cp /share/home/wangz/NCPP-SG15-PBE/" + string(a[t[i]]) + ".SG15.PBE.UPF" + " " + dir_path + temp_name;
			//cout << com << endl;
			system(com.c_str());

		}
		//cin.get();
		//开始写etot.input
		string etot_file = dir_path + temp_name + "/etot.input";
		ofstream fout;
		fout.open(etot_file, ios::out);
		//fout<<"    4       1\n";
		if (num > 40)
		{
			fout << "    4       1\n";
		}
		if (num <= 40 && num > 20)
		{
			fout << "    2       2\n";
		}
		if (num <= 20)
		{
			fout << "    1       4\n";
		}
		fout << "    IN.ATOM  =  atom.config\n";
		for (i = 0; i < n_type; i++)
		{
			fout << "    IN.PSP" << i + 1 << "  =  " + string(a[t[i]]) + ".SG15.PBE.UPF\n";
		}
		/*for (i=0;i<n_type;i++)
		{
			if (t[i]==28)
			{
				fprintf(etot, "    LDAU_PSP%d = 2 6.0\n", i+1);
			}
			if (t[i]==23)
			{
				fprintf(etot, "    LDAU_PSP%d = 2 3.1\n", i+1);
			}
			if (t[i]==24)
			{
				fprintf(etot, "    LDAU_PSP%d = 2 3.5\n", i+1);
			}
			if (t[i]==25)
			{
				fprintf(etot, "    LDAU_PSP%d = 2 3.9\n", i+1);
			}
			if (t[i]==26)
			{
				fprintf(etot, "    LDAU_PSP%d = 2 4.0\n", i+1);
			}
			if (t[i]==27)
			{
				fprintf(etot, "    LDAU_PSP%d = 2 3.4\n", i+1);
			}
			if (t[i]==29)
			{
				fprintf(etot, "    LDAU_PSP%d = 2 4.0\n", i+1);
			}
			if (t[i]==42)
			{
				fprintf(etot, "    LDAU_PSP%d = 2 3.5\n", i+1);
			}
		}*/

		//然后产生mp123
		for (i = 0; i < 3; i++)
		{
			dd[i] = 0;
			for (j = 0; j < 3; j++)
			{
				dd[i] = d[i][j] * d[i][j] + dd[i];
			}
			dd[i] = sqrt(dd[i]);
		}

		kp[0] = pow(double(K_rule) / num * dd[0] * dd[1], 1.0 / 2.0) / dd[0] + 0.5;
		kp[1] = pow(double(K_rule) / num * dd[0] * dd[1], 1.0 / 2.0) / dd[1] + 0.5;
		//        kp[2]=pow(1000.0/num*dd[0]*dd[1]*dd[2],1.0/3.0)/dd[2]+0.5;
		for (i = 0; i < 3; i++)
		{
			if (kp[i] == 0)
			{
				kp[i] = 1;
			}
		}
		//kp[0]=1;
		//kp[1]=1;
		kp[2] = 1;
		fclose(atom);

		//下面根据不同的需要填写不同的内容
		if (flag == "relax")
		{
			fout << "    MP_N123  =  " << kp[0] << " " << kp[1] << " " << kp[2] << " " << "0 0 0" << endl;
			fout << "    XCFUNCTIONAL = PBE\n";
			fout << "    JOB = RELAX\n";
			fout << "    RELAX_DETAIL = 1 100 0.02 1 0.1\n";
			fout << "    ECUT = 50\n";
			fout << "    ECUT2 = 200\n";
			fout << "    SPIN = 2\n";
			fout << "    VDW= DFT-D2\n";


		}
		else if (flag == "scf")
		{
			fout << "    MP_N123  =  " << kp[0] << " " << kp[1] << " " << kp[2] << " " << "0 0 0" << endl;
			fout << "    XCFUNCTIONAL = PBE\n";
			fout << "    JOB = SCF\n";
			//fout << "    RELAX_DETAIL = 1 100 0.02 1 0.1\n";
			fout << "    ECUT = 50\n";
			fout << "    ECUT2 = 200\n";
			fout << "    SPIN = 2\n";

		}
		else if (flag == "nonscf")
		{
			//fout << "    MP_N123  =  " << kp[0] << " " << kp[1] << " " << kp[2] << " " << "0 0 0" << endl;
			fout << "    XCFUNCTIONAL = PBE\n";
			fout << "    JOB = NONSCF\n";
			//fout << "    RELAX_DETAIL = 1 100 0.02 1 0.1\n";
			fout << "    ECUT = 50\n";
			fout << "    ECUT2 = 200\n";
			fout << "    SPIN = 2\n";

			fout << "    CHARGE_DECOMP = T\n";
			fout << "    IN.VR=T\n";
			fout << "    IN.KPT=T\n";
		}
		else if (flag == "generate")
		{
			cout << "just generate the file!!" << endl;

		}
		else
		{
			cout << "unknown clacualtion name!please chekc!" << endl;
			cin.get();
			return 0;
		}
		fout.close();
		cycle_num++;
	}
	fin.close();
	return 0;

}


int  sole_generate_jobs(string name, string dir_path, string flag, report* rp)
{
	//根据一个文件单个产生脚本，便于管理
	//产生etot.input
	char str[200];
	char str1[200];
	char str2[500];
	FILE *in;
	FILE *atom;
	double d[3][3];
	double dd[3];
	int *type;
	int num;
	int i, j;
	int t[20];
	int n_type;
	int flag_type;
	FILE *etot;
	int kp[3];


	//开始读取结构文件，产生脚本,此时结构文件就在目录里面
	atom = fopen(name.c_str(), "r");
	if (atom == NULL)
	{
		printf("ERROR opening %s\n", str1);
		return 1;
	}
	fscanf(atom, "%d", &num);
	fgets(str2, 500, atom);
	fgets(str2, 200, atom);
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			fscanf(atom, "%lf", &d[i][j]);
		}
		fgets(str2, 200, atom);
	}
	//fgets(str2, 200, atom);
	fgets(str2, 200, atom);
	type = (int *)malloc(num * sizeof(int));
	for (i = 0; i < num; i++)
	{
		fscanf(atom, "%d", &type[i]);
		fgets(str2, 200, atom);
	}
	for (i = 0; i < 20; i++)
	{
		t[i] = 0;
	}
	n_type = 0;
	for (i = 0; i < num; i++)
	{
		flag_type = 0;
		for (j = 0; j < n_type; j++)
		{
			if (type[i] == t[j])
			{
				flag_type = 1;
			}
		}
		if (flag_type == 0)
		{
			t[n_type] = type[i];
			n_type++;
		}
	}
	//printf("n_type=%d\n", n_type);


	//开始写etot.input
	string etot_file;
	if (flag.find("wkm") != string::npos) {
		etot_file = "etot.input";
	}
	else
	{
		etot_file = dir_path + "/etot.input";
	}

	ofstream fout;
	fout.open(etot_file, ios::out);
	if (num > 40 || IF_HSE)
	{
		fout << "    4       1\n";
	}
	else if (num <= 40 && num > 20)
	{
		fout << "    2       2\n";
	}
	else
	{
		fout << "    1       4\n";
	}
	fout << "    IN.ATOM  =  atom.config\n";
	if (Soc_on == false) {
		for (i = 0; i < n_type; i++)
		{
			fout << "    IN.PSP" << i + 1 << "  =  " + string(a[t[i]]) + ".SG15.PBE.UPF\n";
		}
	}
	else {
		for (i = 0; i < n_type; i++)
		{
			fout << "    IN.PSP" << i + 1 << "  =  " + string(a[t[i]]) + ".SG15.PBE.SOC.UPF\n";
		}
	}
	if (IF_PLUS_U == true)
	{
		for (i = 0; i < n_type; i++)
		{

			if (t[i] == 28)
			{
				fout << "    LDAU_PSP" << i + 1 << "= 2 6.0" << endl;
			}
			if (t[i] == 23)
			{
				fout << "    LDAU_PSP" << i + 1 << "= 2 3.1" << endl;
			}
			if (t[i] == 24)
			{
				fout << "    LDAU_PSP" << i + 1 << "= 2 3.5" << endl;
			}
			if (t[i] == 25)
			{
				fout << "    LDAU_PSP" << i + 1 << "= 2 3.9" << endl;
			}
			if (t[i] == 26)
			{
				fout << "    LDAU_PSP" << i + 1 << "=  2 4.0" << endl;
			}
			if (t[i] == 27)
			{
				fout << "    LDAU_PSP" << i + 1 << "=  2 3.4" << endl;
			}
			if (t[i] == 29)
			{
				fout << "    LDAU_PSP" << i + 1 << "= 2 4.0" << endl;
			}
			if (t[i] == 42)
			{
				fout << "    LDAU_PSP" << i + 1 << "= 2 3.5" << endl;
			}
			if (t[i] == 22)
			{
				fout << "    LDAU_PSP" << i + 1 << "= 2 4.2" << endl;
			}
			if (t[i] == 40)
			{
				fout << "    LDAU_PSP" << i + 1 << "= 2 4.2" << endl;
			}
			/*if (t[i] >= 23 && t[i] <= 29)
			{
				cout << "plus u config" << name << endl;
				cin.get();
			}*/

		}
	}


	//然后产生mp123
	for (i = 0; i < 3; i++)
	{
		dd[i] = 0;
		for (j = 0; j < 3; j++)
		{
			dd[i] = d[i][j] * d[i][j] + dd[i];
		}
		dd[i] = sqrt(dd[i]);
	}
	double kk = 1000;
	if(flag!="relax")
		kk = K_rule;
	if (IF_K_GENERATE_DIMENTIONAL == true)
		kk = K_rule_for_towD;
	//如果是relax的话，不受任何约束，就是直接1000
	if (CALCULATION_MODE == 0)
	{
		kp[0] = kp[1] = kp[2] = 1;
	}
	else if (CALCULATION_MODE == 2)
	{		
		kp[0] = pow(double(kk) / num * dd[0] * dd[1], 1.0 / 2.0) / dd[0] + 0.5;
		kp[1] = pow(double(kk) / num * dd[0] * dd[1], 1.0 / 2.0) / dd[1] + 0.5;				
		if (IF_ITERATION == false&&flag=="scf") {
			kp[0] = K_point_time * kp[0];
			kp[1] = K_point_time * kp[1];
		}
		kp[2] = 1;
	}
	else if (CALCULATION_MODE == 1)
	{
		kp[0] = pow(double(kk) / num * dd[0] * dd[1], 1.0 / 2.0) / dd[0] + 0.5;
		if (IF_ITERATION == false && flag == "scf") {
			kp[0] = K_point_time * kp[0];
		}
		kp[1] = kp[2] = 1;
	}
	else if (CALCULATION_MODE == 3) {
		kp[0] = pow(double(kk) / num * dd[0] * dd[1], 1.0 / 2.0) / dd[0] + 0.5;
		kp[1] = pow(double(kk) / num * dd[0] * dd[1], 1.0 / 2.0) / dd[1] + 0.5;
		kp[2] = pow(double(kk) / num * dd[0] * dd[1] * dd[2], 1.0 / 3.0) / dd[2] + 0.5;		
	}
	else {
		cout << "unknwon calculation mode number!" << endl;
		cin.get();
	}
	for (i = 0; i < 3; i++)
	{
		if (kp[i] == 0)
		{
			kp[i] = 1;
		}
	}

	fclose(atom);

	//下面根据不同的需要填写不同的内容
	if (flag == "relax")
	{
		
		if (IF_SURE_KPOINT == true && CALCULATION_MODE == 2) {
			if(num<10)
				fout << "    MP_N123  =  " << 8 << " " << 8 << " " << 1 << " " << "0 0 0 2" << endl;
			else if(num>50)
				fout << "    MP_N123  =  " << 2 << " " << 2 << " " << 1 << " " << "0 0 0 2" << endl;
			else
				fout << "    MP_N123  =  " << 6 << " " << 6 << " " << 1 << " " << "0 0 0 2" << endl;
		}
		else if (IF_SURE_KPOINT == true && CALCULATION_MODE == 3) {
			if (num < 10)
				fout << "    MP_N123  =  " << 6 << " " << 6 << " " << 6 << " " << "0 0 0" << endl;
			else if (num > 50)
				fout << "    MP_N123  =  " << 1 << " " << 1 << " " << 1 << " " << "0 0 0" << endl;
			else
				fout << "    MP_N123  =  " << 2<< " " << 2<< " " << 2<< " " << "0 0 0" << endl;
		}
		else			
			fout << "    MP_N123  =  " << kp[0] << " " << kp[1] << " " << kp[2] << " " << "0 0 0" << endl;

		fout << "    XCFUNCTIONAL = PBE\n";
		fout << "    JOB = RELAX\n";
		//fout << "    RELAX_DETAIL = 1 100 0.02 1 0.1\n";
		//如果是0D的话就不再考虑什么晶格的收敛了
		if(CALCULATION_MODE==0)
			fout << "    RELAX_DETAIL = 1 500 0.01\n";
		else if (IF_STRESS==true)
			fout << "    RELAX_DETAIL = 1 500 0.01 1 0.1\n";
		else
			fout << "    RELAX_DETAIL = 1 500 0.01\n";
		/*fout << "    ECUT = 50\n";
		fout << "    ECUT2 = 200\n";*/
		fout << "    ECUT = 70\n";
		fout << "    ECUT2 = 280\n";
		if (spin_check == 1)
		{
			fout << "    SPIN = 1\n";
		}
		else {
			if (Soc_on == true) {
				fout << "    SPIN = 22\n";
			}
			else
			{
				fout << "    SPIN = 2\n";
			}
			
		}
			
		fout << "    CHARGE_DECOMP = T\n";
		if(IF_VDW_CHECK==true)
			fout << "    VDW= DFT-D2\n";


	}
	else if (flag == "WKM_LAST") {


		fout << "    XCFUNCTIONAL = PBEWKM\n";
		fout << "    NUM_BAND = 60\n";
		fout << "    MP_N123  =  " << 20 << " " << 20 << " " << 1 << " " << "0 0 0 2" << endl;
		fout << "    ECUT = 50\n";
		fout << "    ECUT2 = 100\n";
		fout << "    JOB = SCF\n";
		fout << "    SPIN = 2\n";
		//最后需要添加上原始的n123
		fout << "    N123 = " << rp->n_x << " " << rp->n_y << " " << rp->n_z << endl;
	}
	else if (flag == "scf")
	{
		if (Soc_on == true) {
			//计算soc自旋轨道耦合的情况，把k点确定的情况纳入进来
			if (IF_SURE_KPOINT == true && CALCULATION_MODE == 2) {
				if (num < 10)
					fout << "    MP_N123  =  " << 8 << " " << 8 << " " << 1 << " " << "0 0 0 2" << endl;
				else if (num > 50)
					fout << "    MP_N123  =  " << 2 << " " << 2 << " " << 1 << " " << "0 0 0 2" << endl;
				else
					fout << "    MP_N123  =  " << 6 << " " << 6 << " " << 1 << " " << "0 0 0 2" << endl;
			}
			else if (IF_SURE_KPOINT == true && CALCULATION_MODE == 3) {
				if (num < 10)
					fout << "    MP_N123  =  " << 8 << " " << 8 << " " << 8 << " " << "0 0 0" << endl;
				else if (num > 50)
					fout << "    MP_N123  =  " << 2 << " " << 2 << " " << 2 << " " << "0 0 0" << endl;
				else
					fout << "    MP_N123  =  " << 4 << " " << 4 << " " << 4 << " " << "0 0 0" << endl;
			}
			else
				fout << "    MP_N123  =  " << kp[0] << " " << kp[1] << " " << kp[2] << " " << "0 0 0" << endl;


			if (IF_HSE == true)
			{
				fout << "    XCFUNCTIONAL = HSE\n";
			}
			else
			{
				fout << "    XCFUNCTIONAL = PBE\n";
			}
			fout << "    CHARGE_DECOMP = T\n";
			fout << "    ECUT = 60\n";
			fout << "    ECUT2 = 240\n";			
			fout << "    SPIN = 22\n";
		}
		else if (IF_WKM == true) {
			if (IF_SURE_KPOINT == true)
				fout << "    MP_N123  =  " << default_kpoint[0] << " " << default_kpoint[1] << " " << default_kpoint[2] << " " << "0 0 0 2" << endl;
			else
				fout << "    MP_N123  =  " << kp[0] << " " << kp[1] << " " << kp[2] << " " << "0 0 0 2" << endl;
			//fout << "    MP_N123  =  " <<20 << " " << 20 << " " << 1 << " " << "0 0 0 2" << endl;
			fout << "    XCFUNCTIONAL = PBE\n";
			if (IF_HAND_N123 == true) {
				vector<int> n123 = calculate_n123(dd[0], dd[1], dd[2], 200);
				fout << "    N123 = " << n123[0] << " " << n123[1] << " " << n123[2] << endl;
			}
			fout << "    NUM_BAND = 80\n";
			fout << "    PWSCF_OUTPUT = T\n";
			fout << "    ECUT = 50\n";
			fout << "    ECUT2 = 200\n";
			fout << "    SPIN = 1\n";//如果是跑wkm的话，一定只能设1
		}
		else if (IF_MD == true) {
			//有关进行md的计算
			fout << "    MP_N123  =  " << 2 << " " << 1 << " " << 1 << " " << "0 0 0 2" << endl;
			fout << "    MD_DETAIL = 6 3000 2 300 300" << endl;
			fout << "    IN.ATOM= atom.config" << endl;
			fout << "    ECUT = 70\n";
			fout << "    ECUT2 = 280\n";
			fout << "    CHARGE_DECOMP = T\n";
			fout << "    XCFUNCTIONAL = PBE\n";
			fout << "    IN.MDOPT = T\n";
		}
		else {
			//正常的情况
			if (IF_SURE_KPOINT == true && CALCULATION_MODE == 2) {
				if (num < 10)
					fout << "    MP_N123  =  " << 8 << " " << 8 << " " << 1 << " " << "0 0 0 2" << endl;
				else if (num > 50)
					fout << "    MP_N123  =  " << 2 << " " << 2 << " " << 1 << " " << "0 0 0 2" << endl;
				else
					fout << "    MP_N123  =  " << 6 << " " << 6 << " " << 1 << " " << "0 0 0 2" << endl;
			}
			else if (IF_SURE_KPOINT == true && CALCULATION_MODE == 3) {
				if (num < 10)
					fout << "    MP_N123  =  " << 8 << " " << 8 << " " << 8 << " " << "0 0 0" << endl;
				else if (num > 50)
					fout << "    MP_N123  =  " << 2 << " " << 2 << " " << 2 << " " << "0 0 0" << endl;
				else
					fout << "    MP_N123  =  " << 4 << " " << 4 << " " << 4 << " " << "0 0 0" << endl;
			}
			else
				if(IF_HSE==true)
					fout << "    MP_N123  =  " <<ceil(HSE_minaval* kp[0]) << " " << ceil(HSE_minaval*kp[1]) << " " << ceil(HSE_minaval*kp[2]) << " " << "0 0 0" << endl;
				else
					fout << "    MP_N123  =  " << kp[0] << " " << kp[1] << " " << kp[2] << " " << "0 0 0" << endl;


			if (IF_HSE == true)
			{
				fout << "    XCFUNCTIONAL = HSE\n";
			}
			else
			{
				fout << "    XCFUNCTIONAL = PBE\n";
			}
			fout << "    CHARGE_DECOMP = T\n";
			if (IF_DOS == false) {
				fout << "    ECUT = 70\n";
				fout << "    ECUT2 = 280\n";
			}
			else
			{
				fout << "    ECUT = 60\n";
				fout << "    ECUT2 = 240\n";
			}
			if (spin_check == 1)
			{
				fout << "    SPIN = 1\n";
			}
			else
				fout << "    SPIN = 2\n";
		}




		//下面通用的设置
		if(IF_MD==true)
			fout << "    JOB = MD\n";
		else
			fout << "    JOB = SCF\n";
		if(IF_VDW_CHECK)
			fout << "    VDW= DFT-D2\n";

		
	}
	else if (flag == "nonscf")
	{
		if (IF_HSE == true)
		{
			fout << "    XCFUNCTIONAL = HSE\n";
		}
		else
		{
			fout << "    XCFUNCTIONAL = PBE\n";
		}
		fout << "    JOB = NONSCF\n";
		fout << "    ECUT = 60\n";		
		fout << "    ECUT2 = 240\n";
		
		if (IF_HAND_N123 == true && IF_WKM == true) {
			vector<int> n123 = calculate_n123(dd[0], dd[1], dd[2], 100);
			fout << "    N123 = " << n123[0] << " " << n123[1] << " " << n123[2] << endl;
		}
		//fout << "    NUM_BAND = 60\n";
		if (spin_check == 1)
		{
			fout << "    SPIN = 1\n";
		}
		else {
			if(Soc_on==true)
				fout << "    SPIN = 22\n";
			else
				fout << "    SPIN = 2\n";
		}	

		fout << "    CHARGE_DECOMP = T\n";
		fout << "    IN.VR=T\n";
		fout << "    IN.KPT=T\n";
	}
	else if (flag.find("wkm") != string::npos) {
		//补充上wkm的计算过程
		fout << "    JOB=WKM" << endl;
		fout << "    CHARGE_DECOMP = T\n";
		fout << "    SPIN = 2\n";
		fout << "    ECUT = 50\n";
		fout << "    ECUT2 = 100\n";
		//输出三个精度
		fout << "    E_error = 1E-5\n";
		fout << "    Rho_error = 1E-6\n";
		fout << "    WG_error = 0\n";

		fout << "    N123 = " << 3 * rp->n_x << " " << 3 * rp->n_y << " " << 3 * rp->n_z << endl;
		fout << "    MP_N123 = 2 2 2 0 0 0 2" << endl;
		fout << "    XCFUNCTIONAL = PBE\n";
		if (flag == "wkm_cb") {
			fout << "    NUM_ELECTRON =" << 27 * rp->num_electron << endl;
		}
		else if (flag == "wkm_vb") {
			fout << "    NUM_ELECTRON =" << 27 * rp->num_electron - 1 << endl;
		}
		else
		{
			cout << "unkonwn wkm calculation_flag!" << endl;
			cin.get();
		}



	}
	else if (flag == "dos") {
		//进行dos计算
		if(IF_HSE==true)
			fout << "    MP_N123  =  " << ceil(HSE_minaval* kp[0]) << " " << ceil(HSE_minaval*kp[1]) << " " << ceil(HSE_minaval*kp[2]) << " " << "0 0 0" << endl;
		else
			fout << "    MP_N123  =  " << kp[0] << " " << kp[1] << " " << kp[2] << " " << "0 0 0" << endl;
		if(IF_HSE==false)
			fout << "    XCFUNCTIONAL = PBE\n";
		else
			fout << "    XCFUNCTIONAL = HSE\n";
		fout << "    JOB = DOS\n";
		fout << "    ECUT = 60\n";
		fout << "    ECUT2 = 240\n";
		if (spin_check == 1)
		{
			fout << "    SPIN = 1\n";
		}
		else
			fout << "    SPIN = 2\n";
		fout << "    IN.WG = T\n";

	}
	else
	{
		cout << "unknown clacualtion name!please chekc!" << endl;
		cin.get();
		return 0;
	}
	fout.close();

	return 0;

}



int generate_config(string name)
{
	//根据输入的结构产生对应的结构文件
	//根据关键词产生对应批量结构
	if (name == "BiI3")
	{
		vector<int> ele1 = { 9,14,15,16,17,35,53 };
		vector<int>ele2 = { 13,21,24,26,27,31,39,44,45,49,71,77,79,81 };

	}
	else
	{
		cout << "unknown flag name!" << endl;
		cin.get();
		return 0;
	}

}

int if_finish(string path, int flag)
{
	//输入是一个路径
	//输出是计算是否结束，0表示没有report没结束，1表示正常结束 ，2 表示非正常结束
	//首先判断有没有report文件
	FILE* in;
	int line = 0;
	string comm = "ls " + path + "|grep \"^REPORT$\" |wc -l ";
	in = popen(comm.c_str(), "r");
	int num = 0;
	fscanf(in, "%d", &num);
	pclose(in);
	if (num == 0)
	{
		//cout << "unifinshed!" << endl;
		return 0;
	}
	else
	{
		//如果发现有REPORT文件，还要注意给时间让他复制完
		//cout << "found REPORT,but sleep 3s wait to copy finished!" << endl;
		//system("sleep 2");
		comm = "cat " + path + "/REPORT | grep \"total computation time\" | wc -l";
		in = popen(comm.c_str(), "r");
		int res = 0;
		fscanf(in, "%d", &res);
		pclose(in);

		if (res == 0)
		{
			cout << "task :" << path << " :encounter errors!please check!" << endl;
			//cin.get();
			return 2;
		}

		//这里还要监测对于relax的情况
		if (flag == 1)
		{
			ifstream fin;
			fin.open(path + "/RELAXSTEPS", ios::in);
			if (!fin.is_open())
			{
				cout << "can not find the file :RELAXSTEPS！" << endl;
				cin.get();
			}

			char temp[200];
			while (fin.good() && fin.peek() != EOF)
			{
				fin.getline(temp, 170);
				line++;
			}
			fin.close();

			if (line < 100)
			{
				//cout << "right report!" << endl;
				//cin.get();
				return 1;
			}
			else
			{
				//cout << "task :" << path << " :encounter errors!please check!" << endl;
				//cin.get();
				return 2;
			}

		}
		else
		{
			if (res == 1)
				return 1;
		}



	}

}


void read_element(element *e, string &file_element_r, string& ionic_ridus)
{
	int i, j, k;
	char atom_name[120][3] = { "D", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg" };

	//先初始化信息
	for (i = 0; i < 120; i++)
	{
		//strcpy(e[i].name, atom_name[i]);
		e[i].atomic_num = i;
		//e[i].num_common_val = 0;
		e[i].num_metal_radius = 0;
		//e[i].num_unusual_val = 0;
	}
	i = 0, j = 0, k = 0;
	int line_in;
	ifstream fin;
	int atom_num_in;
	double vdw_r_min_in, vdw_r_max_in;
	double cov_r;
	int num_metal_r;
	double *temp;
	char str[5];
	char str_temp[500];
	int num_common_in;
	int num_unsual_in;
	//读取共价键范德华键和金属键

	fin.open(file_element_r.c_str(), ios::in);
	if (!fin.is_open())
	{
		cout << "i can not find the file !" << file_element_r << endl;
		cin.get();
	}
	while (fin.peek() != EOF && fin.good())
	{
		fin >> atom_num_in;
		fin >> vdw_r_min_in;
		fin >> vdw_r_max_in;
		fin >> cov_r;
		fin >> num_metal_r;

		e[atom_num_in].metal_radius = new double[num_metal_r];
		double *temp = new double[num_metal_r];
		for (i = 0; i < num_metal_r; i++)
		{
			fin >> temp[i];
		}
		e[atom_num_in].atomic_num = atom_num_in;
		e[atom_num_in].vdw_radius_max = vdw_r_max_in;
		e[atom_num_in].vdw_radius_min = vdw_r_min_in;
		e[atom_num_in].cov_radius = cov_r;
		e[atom_num_in].num_metal_radius = num_metal_r;
		for (i = 0; i < num_metal_r; i++)
		{
			e[atom_num_in].metal_radius[i] = temp[i];
		}
		free(temp);
	}
	fin.close();
	//然后读取离子键的最大值和最小值
	fin.open(ionic_ridus, ios::in);
	if (!fin.is_open()) {
		cout << "i can not find the file:" << ionic_ridus << endl;
		cin.get();
		return;
	}
	i = 0;
	while (fin.good() && fin.peek() != EOF) {
		fin >> e[i].posi_ionic_ridus;
		fin >> e[i].nega_ionic_ridus;
		i++;
	}
	fin.close();
	return;
}


double dis(double *p1, double *p2)
{
	return sqrt(pow(p1[0] - p2[0], 2) + pow(p1[1] - p2[1], 2) + pow(p1[2] - p2[2], 2));
}

int detect_free(const string& choose)
{
	//检测有多少空余节点，返回空余节点的个数
	//输入是监测哪个队列
	//返回0表示没有空闲，1表示gpu2空闲，2表示gpu3空闲
	vector<int> nouse;
	vector<int> wrong = { 32 };
	if (IF_USE_ALL_QUEUE == false)
	{
		if (choose == "gpu2")
		{
			for (int i = 0; i < 25; i++)
			{
				nouse.push_back(i);
			}
		}
		else if (choose == "gpu3")
		{
			if (IF_QUEUE_CHECK == false) {
				for (int i = 31; i < 61; i++)
				{
					nouse.push_back(i);
				}
			}
			else {
				for (int i = 39; i < 61; i++)
				{
					nouse.push_back(i);
				}
			}
			
		}
		else
		{
			cout << "wrong line name!" << endl;
			cin.get();
			return 0;
		}
	}
	else
	{
		for (int i = 0; i < 25; i++)
		{
			nouse.push_back(i);
		}
		if (IF_QUEUE_CHECK == false) {
			for (int i = 31; i < 61; i++)
			{
				nouse.push_back(i);
			}
		}
		else {
			for (int i = 39; i < 61; i++)
			{
				nouse.push_back(i);
			}
		}
	}
	//删掉监控的错误节点
	vector<int>::iterator it;
	for (auto i : wrong)
	{
		for (it = nouse.begin(); it != nouse.end(); it++)
		{
			if (*it == i)
				nouse.erase(it);
		}
	}
	/*for (auto i : nouse)
	{
		cout << i << endl;
	}
	cin.get();*/
	int i = 0, num = 0, j = 0;
	char line[200];
	string str = "pbsnodes |grep \"state = \" > free_nums";
	system(str.c_str());
	ifstream fin;
	fin.open("free_nums", ios::in);
	if (!fin.is_open())
	{
		cout << "generate the free_nums falied!" << endl;
		cin.get();
		return 0;
	}
	while (fin.good() && fin.peek() != EOF)
	{
		fin.getline(line, 100);
		i++;//此时检测的是第i-1个节点
		if (strstr(line, "free") != NULL)
		{
			for (j = 0; j < nouse.size(); j++)
			{
				if (nouse[j] == (i - 1))
				{
					cout << i-1 << endl;
					num++;
					//补充上这一个快速版本
					fin.close();
					system("rm free_nums");
					//cout << i - 1 << "free,gpu2" << endl;
					//cin.get();
					if (i - 1 < 25)
					{
						//cout << i - 1 << "free,gpu2" << endl;
						//cin.get();
						return 1;
					}

					else if ((i - 1) >= 31)
					{
						//cout << i - 1 << "free,gpu3" << endl;
						//cin.get();
						return 2;
					}
					else
						return 0;


					break;

				}
			}
		}
	}
	fin.close();
	system("rm free_nums");
	return (num || IGNORE_WAIT_LINE);
}


cell::cell(const char *name, int flag)
{
	int i, j, k;
	//cout << "expand the :" << cengshu << "layer" << endl;
	char temp[300];
	double x_pian = 0.0;
	double y_pian = 0.0;
	double z_pian = 0.0;
	//strcpy(wenjian, "atom1.config");
	FILE *in;
	in = fopen(name, "rt");
	//system("pause");
	if (in == NULL)
	{
		printf("error of rading atom.config!\n");
		printf("the filename is :%s\n", name);
		this->Valid = false;
		//cin.get();
		return;
	}
	this->Valid = true;
	fscanf(in, "%d", &num);
	if (flag == 1)
		positive = new int[num];
	type = (int *)malloc(num * sizeof(int));
	letice = (double **)malloc(3 * sizeof(double *));
	for (i = 0; i < 3; i++)
	{
		letice[i] = (double *)malloc(3 * sizeof(double));
	}
	p = (double **)malloc(num * sizeof(double *));
	for (i = 0; i < num; i++)
	{
		p[i] = (double *)malloc(3 * sizeof(double));
	}
	real_position = (double ***)malloc(yanshen * sizeof(double **));
	for (i = 0; i < yanshen; i++)
	{
		real_position[i] = (double **)malloc(num * sizeof(double *));
		for (k = 0; k < num; k++)
			real_position[i][k] = (double *)malloc(3 * sizeof(double));
	}

	p_real = new double **[yanshen];
	for (i = 0; i < yanshen; i++)
	{
		p_real[i] = new double *[num];
		for (k = 0; k < num; k++)
		{
			p_real[i][k] = new double[3];
		}
	}
	if (flag == 1)
	{
		ridus = new double*[num];
		for (i = 0; i < num; i++)
		{
			ridus[i] = new double[2];
		}
	}

	while (fgets(temp, 300, in) != NULL)
	{
		if (strstr(temp, "VECTOR") != NULL || strstr(temp, "vector") != NULL || strstr(temp, "LATTICE") != NULL)
			break;
	}
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			fscanf(in, "%lf", &letice[i][j]);
		}
		fgets(temp, 300, in);
	}

	//fgets(temp, 300, in);
	//cout << temp << endl;
	fgets(temp, 300, in);
	//cout << temp << endl;
	char line[10];
	for (i = 0; i < num; i++)
	{

		fscanf(in, "%d", &type[i]);
		//cout << type[i] << endl;
		fscanf(in, "%lf", &p[i][0]);
		fscanf(in, "%lf", &p[i][1]);
		fscanf(in, "%lf", &p[i][2]);
		if (flag == 1)
		{
			fscanf(in, "%s", line);
			if (line[0] == 'p')
			{
				positive[i] = 1;
				fscanf(in, "%lf", &ridus[i][0]);
			}
			else if (line[0] == 's')//读到s的时候说明是单质
			{
				fscanf(in, "%s", line);
				fscanf(in, "%lf", &ridus[i][0]);
			}
			else//如果是negative的话，就读两个ridus值
			{
				positive[i] = -1;
				fscanf(in, "%lf", &ridus[i][0]);//第一个是负价态对应的半径
				fscanf(in, "%lf", &ridus[i][1]);//对应正价态的半径
			}
		}
		fgets(temp, 300, in);
	}
	//int x_xishu = 0;
	//int y_xishu = 0;
	//int z_zishu = 0;

	//增加记录type的记录
	int temp_save[120] = { 0 };
	type_num = 0;
	for (i = 0; i < num; i++)
	{
		temp_save[type[i]]++;
	}
	for (i = 0; i < 120; i++)
	{
		if (temp_save[i] != 0)
		{
			type_num++;
			type_count.push_back(vector<int>{temp_save[i],i});
		}			
	}
	std::sort(type_count.begin(), type_count.end());
	type_save = new int[type_num];
	
	j = 0;
	for (i = 0; i < 120; i++)
	{
		if (temp_save[i] != 0)
			type_save[j++] = i;
	}

	for (i = 0; i < yanshen; i++)
	{

		for (j = 0; j < num; j++)
		{
			//x_xishu = i/3;
			//y_xishu = i / 3;
			//z_zishu = (i % 9) / 3;

			real_position[i][j][0] = letice[0][0] * p[j][0] + letice[1][0] * p[j][1] + letice[2][0] * p[j][2] + (i % cengshu - ((cengshu - 1) / 2)) * letice[0][0] + (i % (cengshu * cengshu) / cengshu - ((cengshu - 1) / 2)) * letice[1][0] - (i / (cengshu * cengshu) - ((cengshu - 1) / 2)) * letice[2][0];
			real_position[i][j][1] = letice[0][1] * p[j][0] + letice[1][1] * p[j][1] + letice[2][1] * p[j][2] + (i % cengshu - ((cengshu - 1) / 2)) * letice[0][1] + (i % (cengshu * cengshu) / cengshu - ((cengshu - 1) / 2)) * letice[1][1] - (i / (cengshu * cengshu) - ((cengshu - 1) / 2)) * letice[2][1];
			real_position[i][j][2] = letice[0][2] * p[j][0] + letice[1][2] * p[j][1] + letice[2][2] * p[j][2] + (i % cengshu - ((cengshu - 1) / 2)) * letice[0][2] + (i % (cengshu * cengshu) / cengshu - ((cengshu - 1) / 2)) * letice[1][2] - (i / (cengshu * cengshu) - ((cengshu - 1) / 2)) * letice[2][2];
		}
	}
	for (i = 0; i < yanshen; i++)
	{
		for (j = 0; j < num; j++)
		{
			p_real[i][j][0] = (i % cengshu - ((cengshu - 1) / 2)) + p[j][0];
			p_real[i][j][1] = (i % (cengshu * cengshu) / cengshu - ((cengshu - 1) / 2)) + p[j][1];
			p_real[i][j][2] = -(i / (cengshu * cengshu) - ((cengshu - 1) / 2)) + p[j][2];
		}
	}
	my_classify = new int[num];
	//开始按照之前的结果进行标定
	for (i = 0; i < num; i++)
	{
		my_classify[i] = classify_metal_maingroup(type[i]);

	}
	fclose(in);
}
cell::~cell() {
	int i = 0, j = 0, k = 0;
	//free(type);	
	if (this->Valid == true) {
		for (i = 0; i < 3; i++)
		{
			free(letice[i]);
		}
		free(letice);
		for (i = 0; i < num; i++)
		{
			free(p[i]);
		}
		free(p);

		for (i = 0; i < yanshen; i++)
		{
			for (k = 0; k < num; k++)
				free(real_position[i][k]);
			free(real_position[i]);
		}
		free(real_position);


		for (i = 0; i < yanshen; i++)
		{
			for (k = 0; k < num; k++)
			{
				delete[]p_real[i][k];
			}
			delete[]p_real[i];
		}
		delete[]p_real;
	}
	
}
void cell::trans_form(string path, string file_name,vector<vector<int>>&input) {
	//需要自己手写变形的策略
	/*for (int i = 0; i < num; i++) {
		if (my_classify[i] == 1) {
			this->type[i] = 29;
		}
		else if (my_classify[i] == 2) {
			this->type[i] = 35;
		}
	}*/

	//读取元素替换策略，读取两个vector，分别指明顺序的第i个原子应该为什么元素
	for (auto data : input) {
		this->type[data[0]] = data[1];
	}
	output_config(path, file_name);
	return;
}
void cell::get_center_enviro() {
	for (int i = 0; i < num; i++) {
		type_and_num[type[i]]++;
	}
	vector<pair<int, int>> vec;
	for (map<int, int>::iterator it = type_and_num.begin(); it != type_and_num.end(); it++) {
		vec.push_back(pair<int, int>(it->first, it->second));
	}
	//按照原子数从小到大排序
	sort(vec.begin(), vec.end(), map_com);
	//cout << vec[0].first << "," << vec[0].second;
	//cout << vec[1].first << "," << vec[1].second;
	//用来获得元素的中间和环境元素
	center = new int();
	envoroment = new int();
	if (this->type_num == 2) {

		center[0] = vec[0].first;
		envoroment[0] = vec[1].first;
		center_num = envoroment_num = 1;
	}
	else if (this->type_num == 1)
	{
		center[0] = envoroment[0] = type[0];
		center_num = envoroment_num = 1;
	}
	else
	{
		//需要针对是否是金属元素进行标定
		center = new int[type_num];
		envoroment = new int[type_num];
		for (int i = 0; i < type_num; i++) {
			if (classify_metal_maingroup(type_save[i]) == 1)
			{
				center[center_num++] = type_save[i];
			}
			else
			{
				envoroment[envoroment_num++] = type_save[i];
			}
		}

	}
	//检查一下当前的获取元素
	/*cout << "center ";
	for (int i = 0; i < center_num; i++) {
		cout << a[center[i]] << ",";
	}
	cout << "enviroment ";
	for (int i = 0; i < envoroment_num; i++) {
		cout << a[envoroment[i]] << ",";
	}

	cout << "please check!" << endl;
	cin.get();*/
	return;
}
void cell::get_center_valence(string file_name) {
	//输入参数是那个指定的包括路径的文件名称
	ifstream fin;
	fin.open(file_name, ios::in);
	if (!fin.is_open()) {
		cout << "i can not find the file:" << file_name << endl;
		cin.get();
	}
	vector<vector<double>> center_val(center_num);
	for (int i = 0; i < center_num; i++) {
		center_val[i].resize(1);
	}
	int  which_center = 0;
	bool flag = false;
	char buffer[200];
	int temp_type;
	double temp_val;
	fin.getline(buffer, 200);
	while (fin.good() && fin.peek() != EOF) {
		fin >> temp_type;
		fin >> temp_val;
		//cout << temp_type << "," << temp_val << endl;
		for (int i = 0; i < center_num; i++) {
			if (temp_type == center[i]) {
				which_center = i;
				flag = true;
				break;
			}
		}
		if (flag)
		{
			center_val[which_center].push_back(temp_val);
			//cout << "got:" << temp_val << endl;
		}
		fin.getline(buffer, 200);
		flag = false;
	}
	fin.close();
	//下面正式开始计算价态的平均值

	center_valence = new double(center_num);
	for (int i = 0; i < center_num; i++) {
		double temp_sum = 0;
		for (double temp : center_val[i]) {
			temp_sum += temp;
		}
		center_valence[i] = temp_sum / (center_val[i].size() - 1);
		//cout << "center:" << a[center[i]] << ":" << center_valence[i] << endl;
	}
	//cin.get();
	return;
}
string cell::get_cif_name(string name) {
	//输入一个低维的名字，输出是体材料的名字，只输出数字
	size_t start = name.find_first_of('-');
	if (start == string::npos) {
		cout << "wrong low dimensioanl name!" << name << endl;
		cin.get();
		return "";
	}
	string bulk = name.substr(0, start - 1);
	return bulk;
}
string cell::get_cif_spacegroup(string name) {	
		//根据文件名获得体材料的空间群信息
		string space_group = "";
		string bulk = cell::get_cif_name(name);
		ifstream fin;
		fin.open(cif_path + bulk+".cif", ios::in);
		if (!fin.is_open()) {
			cout << "i can not find the file:" << cell::cif_path + bulk << endl;
			cin.get();
			return "";
		}
		char buffer[200];
		while (fin.good()&&fin.peek()!=EOF)
		{
			fin.getline(buffer, 199);
			if (strstr(buffer, "symmetry_space_group_name_H-M") != NULL)
				break;
		}
		fin.close();
		string test = string(buffer);
		int start = test.find_first_of('\'');
		if (start == string::npos) {
			cout << "no such spacegroup data!" << endl;
			return "";
		}
		space_group = test.substr(start+1);
		space_group[space_group.size() - 1] = '\0';
		return space_group;

}
void atom::read_atom(char *name)
{
	FILE *in;
	int i, j;
	char str[500];
	double det;
	int flag;
	in = fopen(name, "r");
	if (in == NULL)
	{
		flag_read = 0;
		printf("Do not have file %s\n", name);
		return;
	}
	flag_read = 1;
	fscanf(in, "%d", &num_atom);
	d = (double **)malloc(3 * sizeof(double *));
	rd = (double **)malloc(3 * sizeof(double *));
	for (i = 0; i < 3; i++)
	{
		d[i] = (double *)malloc(3 * sizeof(double));
		rd[i] = (double *)malloc(3 * sizeof(double));
	}
	atom_type = (int *)malloc(num_atom * sizeof(int));
	p = (double **)malloc(num_atom * sizeof(double *));
	for (i = 0; i < num_atom; i++)
	{
		p[i] = (double *)malloc(3 * sizeof(double));
	}
	rp = (double **)malloc(num_atom * sizeof(double *));
	for (i = 0; i < num_atom; i++)
	{
		rp[i] = (double *)malloc(3 * sizeof(double));
	}

	fgets(str, 500, in);
	fgets(str, 500, in);
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			fscanf(in, "%lf", &d[i][j]);
		}
		fgets(str, 500, in);
	}
	fgets(str, 500, in);
	//fgets(str, 500, in);
	for (i = 0; i < num_atom; i++)
	{
		fscanf(in, "%d", &atom_type[i]);
		//printf("%d\n", atom_type[i]);
		for (j = 0; j < 3; j++)
		{
			fscanf(in, "%lf", &p[i][j]);
		}
		fgets(str, 500, in);
	}
	type_num = 0;

	for (i = 0; i < num_atom; i++)
	{
		flag = 0;
		for (j = 0; j < type_num; j++)
		{
			if (atom_type[i] == type_name[j])
			{
				flag = 1;
				break;
			}
		}
		if (flag == 0)
		{
			type_name[type_num] = atom_type[i];
			type_num++;
		}
	}

	for (i = 0; i < num_atom; i++)
	{
		for (j = 0; j < 3; j++)
		{
			rp[i][j] = p[i][0] * d[0][j] + p[i][1] * d[1][j] + p[i][2] * d[2][j];
		}
	}



	det = d[0][0] * d[1][1] * d[2][2] + d[0][1] * d[1][2] * d[2][0] + d[0][2] * d[1][0] * d[2][1] - d[0][0] * d[1][2] * d[2][1] - d[0][1] * d[1][0] * d[2][2] - d[0][2] * d[1][1] * d[2][0];
	rd[0][0] = (d[1][1] * d[2][2] - d[1][2] * d[2][1]) / det * 2 * PI*BOHR;
	rd[0][1] = (d[1][2] * d[2][0] - d[1][0] * d[2][2]) / det * 2 * PI*BOHR;
	rd[0][2] = (d[1][0] * d[2][1] - d[2][0] * d[1][1]) / det * 2 * PI*BOHR;

	rd[1][0] = (d[2][1] * d[0][2] - d[2][2] * d[0][1]) / det * 2 * PI*BOHR;
	rd[1][1] = (d[2][2] * d[0][0] - d[2][0] * d[0][2]) / det * 2 * PI*BOHR;
	rd[1][2] = (d[2][0] * d[0][1] - d[2][1] * d[0][0]) / det * 2 * PI*BOHR;

	rd[2][0] = (d[0][1] * d[1][2] - d[1][1] * d[0][2]) / det * 2 * PI*BOHR;
	rd[2][1] = (d[0][2] * d[1][0] - d[1][2] * d[0][0]) / det * 2 * PI*BOHR;
	rd[2][2] = (d[0][0] * d[1][1] - d[1][0] * d[0][1]) / det * 2 * PI*BOHR;

	/*rd[0][0] = (d[1][0] * d[2][0] - d[1][2] * d[2][1]) / det * 2 * PI*BOHR;
	rd[0][1] = (d[1][2] * d[2][0] - d[1][0] * d[2][2]) / det * 2 * PI*BOHR;
	rd[0][2] = (d[1][0] * d[2][1] - d[2][0] * d[1][1]) / det * 2 * PI*BOHR;

	rd[1][0] = (d[2][1] * d[0][2] - d[2][2] * d[0][1]) / det * 2 * PI*BOHR;
	rd[1][1] = (d[2][2] * d[0][0] - d[2][0] * d[0][2]) / det * 2 * PI*BOHR;
	rd[1][2] = (d[2][0] * d[0][1] - d[2][1] * d[0][0]) / det * 2 * PI*BOHR;

	rd[2][0] = (d[0][1] * d[1][2] - d[1][1] * d[0][2]) / det * 2 * PI*BOHR;
	rd[2][1] = (d[0][2] * d[1][0] - d[1][2] * d[0][0]) / det * 2 * PI*BOHR;
	rd[2][2] = (d[0][0] * d[1][1] - d[1][0] * d[0][1]) / det * 2 * PI*BOHR;
	*/fclose(in);
	return;
}

void atom::print_name(int i, FILE *in)
{
	fprintf(in, "%s", atom_name[i]);
	return;
}

void atom::freeatom() {
	int i;
	if (flag_read == 1)
	{
		for (i = 0; i < 3; i++)
		{
			free(d[i]);
			free(rd[i]);
		}
		free(d);
		free(rd);
		free(atom_type);
		for (i = 0; i < num_atom; i++)
		{
			free(p[i]);
		}
		free(p);
		for (i = 0; i < num_atom; i++)
		{
			free(rp[i]);
		}
		free(rp);
	}
}


int Read_eig(char *name, struct E_E *eig, int num_band, int num_kpt, int spin, double ***eig_kns)
{
	FILE *in;
	char str[300];
	int flag = 1;
	int i, j, k;
	in = fopen(name, "r");
	if (in == NULL)
	{
		printf("No inputfile: %s\n", name);
		return 1;
	}
	for (i = 0; i < spin; i++)
	{
		for (j = 0; j < num_kpt; j++)
		{
			flag = 1;
			while (fgets(str, 300, in) != NULL)
			{
				if (strstr(str, "eigen energies, in eV") != NULL)
				{
					for (k = 0; k < num_band; k++)
					{
						fscanf(in, "%lf", &eig[k*num_kpt*spin + i * num_kpt + j].energy);
						//printf("%lf\n", eig[k*num_kpt*spin + i * num_kpt + j].energy);
						eig_kns[j][k][i] = eig[k*num_kpt*spin + i * num_kpt + j].energy;
						eig[k*num_kpt*spin + i * num_kpt + j].kpt = j;
						eig[k*num_kpt*spin + i * num_kpt + j].spin = i;
						eig[k*num_kpt*spin + i * num_kpt + j].band = k;
					}
					flag = 0;
					break;
				}
			}
		}
	}
	fclose(in);
	if (flag == 0)
	{
		return 0;
	}
	else
	{
		printf("Not enough %s in input file %s\n", "eigen energies", name);
		return 1;
	}
}



int Read_kpp(char *name, double **kp_p, int num_kpt)
{
	FILE *in;
	in = fopen(name, "r");
	double temp;
	int i;
	char str[200];

	while (fgets(str, 200, in) != NULL)
	{
		if (strstr(str, "total number of K-point:") != NULL)
		{
			/* for (i=0;i<strlen(str);i++)
			 {
				 if (str[i]==':')
				 {
					 break;
				 }
			 }
			 fseek(in, 2 * (-strlen(str) + i + 1), SEEK_CUR);*/
			break;
		}
	}
	//fgets(str, 200, in);
	for (i = 0; i < num_kpt; i++)
	{
		fscanf(in, "%lf", &kp_p[i][0]);
		fscanf(in, "%lf", &kp_p[i][1]);
		fscanf(in, "%lf", &kp_p[i][2]);
		fscanf(in, "%lf", &temp);
	}

	fclose(in);

	return 1;
}

int Comp(const void *p1, const void *p2)
{
	if (((struct E_E *)p1)->energy > ((struct E_E *)p2)->energy)
	{
		return 1;
	}
	if (((struct E_E *)p1)->energy < ((struct E_E *)p2)->energy)
	{
		return -1;
	}
	return 0;
}
bool map_com(pair<int, int>a, pair<int, int>b) {
	return a.second < b.second;
}
int  report::read_report(char *name_report, int spin)
{
	this->report_name = string(name_report);
	int i, j;
	if (Read_N123(name_report, "N123", &n_x, &n_y, &n_z) == 1)
	{
		//cout << name_report << " get N123 encounter wrong!" << endl;
		//cin.get();
		return 1;
	}

	if (Read_N123(name_report, "MP_N123", &mp_x, &mp_y, &mp_z) == 1)
	{
		//return;
	}
	if (Read_int(name_report, "NUM_BAND", &num_band) == 1)
	{
		cout << name_report << " get NUM BAND encounter wrong!" << endl;
		cin.get();
		return 1;
	}
	if (Read_int(name_report, "NUM_ELECTRON", &num_electron) == 1)
	{
		//cout << name_report << " get NUM_ELECTRION encounter wrong!" << endl;
		//cin.get();
		return 1;
	}
	// if (Read_int(name_report, "SPIN", &spin)==1)
	// {
	 //    return;
	// }
	this->spin = spin;
	if (Soc_on == true)
		this->spin = 1;
	//	printf("spin :%d\n",spin);
	if (Read_kpt_num(name_report, "total number of K-point", &num_kpt) == 1)
	{
		//cout << name_report << " get K-POINT encounter wrong!" << endl;
		//cin.get();
		return 1;
	}
	eig_kns = (double ***)malloc(num_kpt * sizeof(double **));
	for (i = 0; i < num_kpt; i++)
	{
		eig_kns[i] = (double **)malloc(num_band * sizeof(double *));
	}
	for (i = 0; i < num_kpt; i++)
	{
		for (j = 0; j < num_band; j++)
		{
			eig_kns[i][j] = (double *)malloc(this->spin * sizeof(double));
		}
	}

	eig = (struct E_E *)malloc(this->spin*num_kpt*num_band * sizeof(struct E_E));
	if (Read_eig(name_report, eig, num_band, num_kpt, this->spin, eig_kns) == 1)
	{
		return 1;
	}
	kp_p = (double **)malloc(num_kpt * sizeof(double *));
	for (i = 0; i < num_kpt; i++)
	{
		kp_p[i] = (double *)malloc(3 * sizeof(double));
	}
	Read_kpp(name_report, kp_p, num_kpt);
	qsort((void *)&eig[0], (size_t)num_band*num_kpt*this->spin, sizeof(struct E_E), Comp);
	/*for (i = 0; i < num_band; i++)
	{
		for (j = 0; j < spin*num_kpt; j++)
		{
			printf("%lf,", eig[i*spin*num_kpt + j].energy);
		}
		printf("%\n");
	}*/
	//修改费米能级的数值
	//this->orig_fermi = (eig[num_kpt * num_electron * this->spin / 2 - 1].energy + eig[num_kpt * num_electron * this->spin / 2].energy) / 2;//获取根据二者得到的费米能级
	this->orig_fermi = eig[num_kpt * num_electron * this->spin / 2 - 1].energy;
	cbm_band = eig[num_kpt * num_electron * this->spin / 2].band + 1;
	vbm_band = eig[num_kpt * num_electron * this->spin / 2 - 1].band + 1;

	vbm_kppt = kp_p[eig[num_kpt * num_electron * this->spin / 2 - 1].kpt];
	cbm_kppt = kp_p[eig[num_kpt * num_electron * this->spin / 2].kpt];
	vbm_energy = eig[num_kpt * num_electron * this->spin / 2 - 1].energy;
	cbm_energy = eig[num_kpt * num_electron *this->spin / 2].energy;

	cbm_kpt = eig[num_kpt*num_electron*this->spin / 2].kpt;//应该是12
	vbm_kpt = eig[num_kpt*num_electron*this->spin / 2 - 1].kpt;//应该是37
	return 0;
}

int Read_double(char *name, const char *tag, double *out)
{
	FILE *in;
	char str[300];
	int i;
	int flag = 1;
	int j = 0;
	in = fopen(name, "r");
	if (in == NULL)
	{
		printf("No inputfile: %s\n", name);
		return 1;
	}
	while (fgets(str, 300, in) != NULL)
	{
		if (strstr(str, tag) != NULL)
		{
			flag = 0;
			for (i = 0; i < strlen(str); i++)
			{
				if (str[i] == '=')
				{
					break;
				}
			}
			char aa[100];
			//	printf("%s\n", tag);
			for (j = i + 1; j < strlen(str) - 1; j++)
			{
				aa[j - i - 1] = str[j];
			}
			aa[strlen(str) - 1] = '\0';
			*out = atof(aa);
			/*fseek(in, 2 * (-strlen(str) + i + 1), SEEK_CUR);
			fscanf(in, "%lf", out);*/
			break;
		}
	}
	fclose(in);
	if (flag == 0)
	{
		return 0;
	}
	else
	{
		printf("No %s in input file %s\n", tag, name);
		return 1;
	}
}

int Read_int(char *name, const char *tag, int *out)
{
	FILE *in;
	char str[300];
	int i;
	int j;
	int flag = 1;
	in = fopen(name, "r");
	if (in == NULL)
	{
		printf("No inputfile: %s\n", name);
		return 1;
	}
	while (fgets(str, 300, in) != NULL)
	{
		if (strstr(str, tag) != NULL)
		{
			flag = 0;
			for (i = 0; i < strlen(str); i++)
			{
				if (str[i] == '=')
				{
					break;
				}
			}
			char aa[100];
			//	printf("%s\n", tag);
			int mm = 0;
			for (j = i + 1; j < strlen(str) - 1; j++)
			{
				if (str[j] == ' ' || str[j] == '\t')
					continue;
				aa[mm++] = str[j];
			}
			aa[mm] = '\0';
			*out = atoi(aa);
			/*fseek(in, 2 * (int(strlen(str)) - i -1), SEEK_CUR);
			printf("%d\n", 2 * (-int(strlen(str)) + i + 1));
			fscanf(in, "%s\n", aa);
			fscanf(in, "%d", out);*/
			//	printf("%d\n", *out);
			break;
		}
	}
	fclose(in);
	if (flag == 0)
	{
		return 0;
	}
	else
	{
		printf("No %s in input file %s\n", tag, name);
		return 1;
	}
}

int Read_N123(char *name, const char *tag, int *n1, int *n2, int *n3)
{
	FILE *in;
	char str[300];
	int i;
	int flag = 1;
	in = fopen(name, "r");
	if (in == NULL)
	{
		printf("No inputfile: %s\n", name);
		return 1;
	}
	while (fgets(str, 300, in) != NULL)
	{
		if (strstr(str, tag) != NULL)
		{
			flag = 0;
			for (i = 0; i < strlen(str); i++)
			{
				if (str[i] == '=')
				{
					break;
				}
			}
			//printf("%d\n", strlen(str));			
			fseek(in, (-int(strlen(str)) + i + 1), SEEK_CUR);
			//fseek(in, -7, SEEK_CUR);
			//char strr[200];
			//fscanf(in,"%s\n",strr);
			fscanf(in, "%d", n1);
			fscanf(in, "%d", n2);
			fscanf(in, "%d", n3);
			//	printf("%s\n", tag);
			//	printf("%d,%d,%d\n", *n1,*n2,*n3);
			break;
		}
	}
	fclose(in);
	if (flag == 0)
	{
		return 0;
	}
	else
	{
		//printf("No %s in input file %s\n", tag, name);
		return 1;
	}
}

int Read_kpt_num(char *name, const char *tag, int *out)
{
	FILE *in;
	char str[300];
	int i;
	int flag = 1;
	int j = 0;
	in = fopen(name, "r");
	if (in == NULL)
	{
		printf("No inputfile: %s\n", name);
		return 1;
	}
	while (fgets(str, 300, in) != NULL)
	{
		if (strstr(str, tag) != NULL)
		{
			flag = 0;
			for (i = 0; i < strlen(str); i++)
			{
				if (str[i] == ':')
				{
					break;
				}
			}
			char aa[100];
			int mm = 0;
			//	printf("%s\n", tag);
			for (j = i + 1; j < strlen(str) - 1; j++)
			{
				if (str[j] == ' ' || str[j] == '\t')
					continue;
				aa[mm++] = str[j];
			}
			aa[mm] = '\0';
			*out = atoi(aa);
			//	printf("%d\n", *out);
			break;
		}
	}
	fclose(in);
	if (flag == 0)
	{
		return 0;
	}
	else
	{
		printf("No %s in input file %s\n", tag, name);
		return 1;
	}
}

double report::get_gap(string path,vector<string>&step)
{
	//report rep;
	//rep.read_report(const_cast<char *>((path+"/REPORT").c_str()));
	atom a;
	a.read_atom(const_cast<char *>((path + "/atom.config").c_str()));
	double **d_temp;
	int i, j;
	int flag = 0;
	double *max, *min;
	double f;
	double kp_temp[3];
	int kp_int[3];
	FILE *out;
	double gap = 0;
	max = (double *)malloc(num_band * spin * sizeof(double));
	min = (double *)malloc(num_band * spin * sizeof(double));
	d_temp = (double **)malloc(3 * sizeof(double *));
	for (i = 0; i < 3; i++)
	{
		d_temp[i] = (double *)malloc(3 * sizeof(double));
	}
	for (i = 0; i < num_band * spin; i++)
	{
		max[i] = -9999;
		min[i] = 9999;
	}

	f = (eig[num_kpt * num_electron * spin / 2 - 1].energy + eig[num_kpt * num_electron * spin / 2].energy) / 2;
	this->orig_fermi = f;
	//printf("f=%lf\n", f);
	for (i = 0; i < num_band * num_kpt * spin; i++)
	{
		if (max[eig[i].band + num_band * eig[i].spin] < eig[i].energy)
		{
			max[eig[i].band + num_band * eig[i].spin] = eig[i].energy;
		}
		if (min[eig[i].band + num_band * eig[i].spin] > eig[i].energy)
		{
			min[eig[i].band + num_band * eig[i].spin] = eig[i].energy;
		}
	}
	flag = 0;
	for (i = 0; i < num_band * spin; i++)
	{
		//printf ("%d\t%lf\t%lf\n",i, max[i]-f, min[i]-f);
		if ((max[i] - f) * (min[i] - f) < 0)
		{

			if (flag == 0)
			{
				printf("There is no band gap:\n");
			}
			flag = 1;
			//printf("Spin %d Band %d\n", i / rep.num_band, i % rep.num_band);
		}
	}
	if (flag == 1)
	{
		printf("Go cross fermi level\n");
		//printf("Notes: All the number counts from 0\n");
		free(max);
		free(min);
		return -1;
	}

	if (flag == 0)
	{
		//printf("1/52 is paramter.\n");
		//cout << KPT_X << KPT_Y << KPT_Z << endl;

		gap = eig[num_kpt * num_electron * spin / 2].energy - eig[num_kpt * num_electron * spin / 2 - 1].energy;
		printf("Band Gap= %lf\n", gap);
		printf("near 13 %d\tKPT %d\tSPIN %d\n", eig[num_kpt * num_electron * spin / 2].band, eig[num_kpt * num_electron * spin / 2].kpt + 1, eig[num_kpt * num_electron * spin / 2].spin);
		printf("near 38 %d\tKPT %d\tSPIN %d\n", eig[num_kpt * num_electron * spin / 2 - 1].band, eig[num_kpt * num_electron * spin / 2 - 1].kpt + 1, eig[num_kpt * num_electron * spin / 2 - 1].spin);
		//注意要记录相应的步长了
		step[0] = to_string(kp_p[eig[num_kpt*num_electron*spin / 2].kpt][0]) + ","+to_string(kp_p[eig[num_kpt*num_electron*spin / 2].kpt][1]) + ","+to_string(kp_p[eig[num_kpt*num_electron*spin / 2].kpt][2]);
		step[1] = to_string(kp_p[eig[num_kpt*num_electron*spin / 2-1].kpt][0]) + "," + to_string(kp_p[eig[num_kpt*num_electron*spin / 2-1].kpt][1]) + "," + to_string(kp_p[eig[num_kpt*num_electron*spin / 2-1].kpt][2]);

		printf("CBM step original is:%lf,%lf,%lf\n", kp_p[eig[num_kpt*num_electron*spin / 2].kpt][0], kp_p[eig[num_kpt*num_electron*spin / 2].kpt][1], kp_p[eig[num_kpt*num_electron*spin / 2].kpt][2]);
		printf("VBM step original is:%lf,%lf,%lf\n", kp_p[eig[num_kpt*num_electron*spin / 2 - 1].kpt][0], kp_p[eig[num_kpt*num_electron*spin / 2 - 1].kpt][1], kp_p[eig[num_kpt*num_electron*spin / 2 - 1].kpt][2]);
	
	
	}
	return gap;

}
double report::get_singel_number(string flag) {
	//输入是一个单一的
	if (flag == "CBM") {
		return this->eig[num_kpt * num_electron * spin / 2].energy;
	}
	else if (flag == "VBM") {
		return this->eig[num_kpt * num_electron * spin / 2 - 1].energy;
	}
	ifstream fin;
	fin.open(this->report_name, ios::in);
	if (!fin.is_open()) {
		cout << "i can not find the file" << report_name << endl;
		cin.get();
		return 0;
	}
	char buffer[150];
	while (fin.good() && fin.peek() != EOF) {
		fin.getline(buffer, 149);
		if (strstr(buffer, flag.c_str()) != NULL)
			break;
	}
	double res = 0;
	string a = string(buffer);
	int start = a.find_first_of(flag[flag.size() - 1]);
	res = atof(a.substr(start + 1, a.size() - 1).c_str());
	cout << "got " << flag << ":" << res << endl;
	//cin.get();
	fin.close();
	return res;
}

int make_efm(report& rep, string path, int KPT_X, int KPT_Y, int KPT_Z)
{
	//report rep;
	//rep.read_report(const_cast<char *>((path+"/REPORT").c_str()));
	atom a;
	a.read_atom(const_cast<char *>((path + "/atom.config").c_str()));
	double **d_temp;
	int i, j;
	int flag = 0;
	double *max, *min;
	double f;
	double kp_temp[3];
	int kp_int[3];
	FILE *out;

	max = (double *)malloc(rep.num_band * rep.spin * sizeof(double));
	min = (double *)malloc(rep.num_band * rep.spin * sizeof(double));
	d_temp = (double **)malloc(3 * sizeof(double *));
	for (i = 0; i < 3; i++)
	{
		d_temp[i] = (double *)malloc(3 * sizeof(double));
	}
	for (i = 0; i < rep.num_band * rep.spin; i++)
	{
		max[i] = -9999;
		min[i] = 9999;
	}

	f = (rep.eig[rep.num_kpt * rep.num_electron * rep.spin / 2 - 1].energy + rep.eig[rep.num_kpt * rep.num_electron * rep.spin / 2].energy) / 2;
	//printf("f=%lf\n", f);
	for (i = 0; i < rep.num_band * rep.num_kpt * rep.spin; i++)
	{
		if (max[rep.eig[i].band + rep.num_band * rep.eig[i].spin] < rep.eig[i].energy)
		{
			max[rep.eig[i].band + rep.num_band * rep.eig[i].spin] = rep.eig[i].energy;
		}
		if (min[rep.eig[i].band + rep.num_band * rep.eig[i].spin] > rep.eig[i].energy)
		{
			min[rep.eig[i].band + rep.num_band * rep.eig[i].spin] = rep.eig[i].energy;
		}
	}
	flag = 0;
	for (i = 0; i < rep.num_band * rep.spin; i++)
	{
		//printf ("%d\t%lf\t%lf\n",i, max[i]-f, min[i]-f);
		if ((max[i] - f) * (min[i] - f) < 0)
		{

			if (flag == 0)
			{
				printf("There is no band gap:\n");
			}
			flag = 1;
			//printf("Spin %d Band %d\n", i / rep.num_band, i % rep.num_band);
		}
	}
	if (flag == 1)
	{
		printf("Go cross fermi level\n");
		//printf("Notes: All the number counts from 0\n");
		free(max);
		free(min);
		return 1;
	}

	if (flag == 0)
	{
		//printf("1/52 is paramter.\n");
		//cout << KPT_X << KPT_Y << KPT_Z << endl;
		printf("Band Gap= %lf\n", rep.eig[rep.num_kpt * rep.num_electron * rep.spin / 2].energy - rep.eig[rep.num_kpt * rep.num_electron * rep.spin / 2 - 1].energy);
		printf("near 13 %d\tKPT %d\tSPIN %d\n", rep.eig[rep.num_kpt * rep.num_electron * rep.spin / 2].band, rep.eig[rep.num_kpt * rep.num_electron * rep.spin / 2].kpt + 1, rep.eig[rep.num_kpt * rep.num_electron * rep.spin / 2].spin);
		printf("near 38 %d\tKPT %d\tSPIN %d\n", rep.eig[rep.num_kpt * rep.num_electron * rep.spin / 2 - 1].band, rep.eig[rep.num_kpt * rep.num_electron * rep.spin / 2 - 1].kpt + 1, rep.eig[rep.num_kpt * rep.num_electron * rep.spin / 2 - 1].spin);
		printf("CBM step original is:%lf,%lf,%lf\n", rep.kp_p[rep.eig[rep.num_kpt*rep.num_electron*rep.spin / 2].kpt][0], rep.kp_p[rep.eig[rep.num_kpt*rep.num_electron*rep.spin / 2].kpt][1], rep.kp_p[rep.eig[rep.num_kpt*rep.num_electron*rep.spin / 2].kpt][2]);
		printf("VBM step original is:%lf,%lf,%lf\n", rep.kp_p[rep.eig[rep.num_kpt*rep.num_electron*rep.spin / 2 - 1].kpt][0], rep.kp_p[rep.eig[rep.num_kpt*rep.num_electron*rep.spin / 2 - 1].kpt][1], rep.kp_p[rep.eig[rep.num_kpt*rep.num_electron*rep.spin / 2 - 1].kpt][2]);

		//printf("Notes: All the number counts from 1\n");

		inv(3, a.rd, d_temp);
		/*for (i = 0; i < 3; i++)
		{
			for (j = 0; j < 3; j++)
			{
				printf("%lf\t", d_temp[i][j]);
			}
			printf("\n");
		}*/

		out = fopen((path + "/IN.KPT").c_str(), "wb");
		if (CALCULATION_MODE == 2)
		{
			fprintf(out, "    50\n");
		}
		else if (CALCULATION_MODE == 1)
		{
			fprintf(out, "    10\n");
		}
		fprintf(out, "    2        0\n");
		kp_temp[0] = rep.kp_p[rep.eig[rep.num_kpt * rep.num_electron * rep.spin / 2].kpt][0] * d_temp[0][0] + rep.kp_p[rep.eig[rep.num_kpt * rep.num_electron * rep.spin / 2].kpt][1] * d_temp[1][0] + rep.kp_p[rep.eig[rep.num_kpt * rep.num_electron * rep.spin / 2].kpt][2] * d_temp[2][0];
		kp_temp[1] = rep.kp_p[rep.eig[rep.num_kpt * rep.num_electron * rep.spin / 2].kpt][0] * d_temp[0][1] + rep.kp_p[rep.eig[rep.num_kpt * rep.num_electron * rep.spin / 2].kpt][1] * d_temp[1][1] + rep.kp_p[rep.eig[rep.num_kpt * rep.num_electron * rep.spin / 2].kpt][2] * d_temp[2][1];
		kp_temp[2] = rep.kp_p[rep.eig[rep.num_kpt * rep.num_electron * rep.spin / 2].kpt][0] * d_temp[0][2] + rep.kp_p[rep.eig[rep.num_kpt * rep.num_electron * rep.spin / 2].kpt][1] * d_temp[1][2] + rep.kp_p[rep.eig[rep.num_kpt * rep.num_electron * rep.spin / 2].kpt][2] * d_temp[2][2];
		//区分正负
		for (i = 0; i < 3; i++)
		{
			if (i == 0)
			{
				if (kp_temp[i] < 0.0)
				{
					kp_int[i] = int(kp_temp[i] * KPT_X - 0.5);
				}
				else
				{
					kp_int[i] = int(kp_temp[i] * KPT_X + 0.5);
				}
				continue;
			}
			if (i == 1)
			{
				if (kp_temp[i] < 0.0)
				{
					kp_int[i] = int(kp_temp[i] * KPT_Y - 0.5);
				}
				else
				{
					kp_int[i] = int(kp_temp[i] * KPT_Y + 0.5);
				}
				continue;
			}
			if (i == 2)
			{
				if (kp_temp[i] < 0.0)
				{
					kp_int[i] = int(kp_temp[i] * KPT_Z - 0.5);
				}
				else
				{
					kp_int[i] = int(kp_temp[i] * KPT_Z + 0.5);
				}
				continue;
			}
		}
		//	printf("%lf,%lf,%lf\n",kp_temp[0],kp_temp[1],kp_temp[2]);
		//printf("%d,%d,%d\n", kp_int[0], kp_int[1], kp_int[2]);
		//printf("%lf,%lf,%lf\n", (double)kp_int[0] / KPT_X, (double)kp_int[1] / KPT_Y, (double)kp_int[2] / KPT_Z);
		//找到了1/52单位的最低值，然后产生1/104的新区间
		if (CALCULATION_MODE == 2)
		{
			for (i = -2; i < 3; i++)
			{
				for (j = -2; j < 3; j++)
				{
					fprintf(out, "%lf\t%lf\t%lf\t1.00000\n", ((double)kp_int[0] + i) / KPT_X, ((double)kp_int[1] + j) / KPT_Y, ((double)kp_int[2]) / KPT_Z);
				}
			}
		}
		else if (CALCULATION_MODE == 1)
		{
			for (i = -2; i < 3; i++)
			{
				fprintf(out, "%lf\t%lf\t%lf\t1.00000\n", ((double)kp_int[0] + i) / KPT_X, ((double)kp_int[1]) / KPT_Y, ((double)kp_int[2]) / KPT_Z);
			}
		}



		//然后是cbm
		kp_temp[0] = rep.kp_p[rep.eig[rep.num_kpt * rep.num_electron * rep.spin / 2 - 1].kpt][0] * d_temp[0][0] + rep.kp_p[rep.eig[rep.num_kpt * rep.num_electron * rep.spin / 2 - 1].kpt][1] * d_temp[1][0] + rep.kp_p[rep.eig[rep.num_kpt * rep.num_electron * rep.spin / 2 - 1].kpt][2] * d_temp[2][0];
		kp_temp[1] = rep.kp_p[rep.eig[rep.num_kpt * rep.num_electron * rep.spin / 2 - 1].kpt][0] * d_temp[0][1] + rep.kp_p[rep.eig[rep.num_kpt * rep.num_electron * rep.spin / 2 - 1].kpt][1] * d_temp[1][1] + rep.kp_p[rep.eig[rep.num_kpt * rep.num_electron * rep.spin / 2 - 1].kpt][2] * d_temp[2][1];
		kp_temp[2] = rep.kp_p[rep.eig[rep.num_kpt * rep.num_electron * rep.spin / 2 - 1].kpt][0] * d_temp[0][2] + rep.kp_p[rep.eig[rep.num_kpt * rep.num_electron * rep.spin / 2 - 1].kpt][1] * d_temp[1][2] + rep.kp_p[rep.eig[rep.num_kpt * rep.num_electron * rep.spin / 2 - 1].kpt][2] * d_temp[2][2];

		//区分正负
		for (i = 0; i < 3; i++)
		{
			if (i == 0)
			{
				if (kp_temp[i] < 0.0)
				{
					kp_int[i] = int(kp_temp[i] * KPT_X - 0.5);
				}
				else
				{
					kp_int[i] = int(kp_temp[i] * KPT_X + 0.5);
				}
				continue;
			}
			if (i == 1)
			{
				if (kp_temp[i] < 0.0)
				{
					kp_int[i] = int(kp_temp[i] * KPT_Y - 0.5);
				}
				else
				{
					kp_int[i] = int(kp_temp[i] * KPT_Y + 0.5);
				}
				continue;
			}
			if (i == 2)
			{
				if (kp_temp[i] < 0.0)
				{
					kp_int[i] = int(kp_temp[i] * KPT_Z - 0.5);
				}
				else
				{
					kp_int[i] = int(kp_temp[i] * KPT_Z + 0.5);
				}
				continue;
			}
		}
		if (CALCULATION_MODE == 2)
		{
			for (i = -2; i < 3; i++)
			{
				for (j = -2; j < 3; j++)
				{
					//printf("%lf\t%lf\t%lf\n", (double)kp_int[0]+0.5*i, (double)kp_int[1]+0.5*j, (double)kp_int[2]);
					fprintf(out, "%lf\t%lf\t%lf\t1.00000\n", ((double)kp_int[0] + i) / KPT_X, ((double)kp_int[1] + j) / KPT_Y, ((double)kp_int[2]) / KPT_Z);
				}
			}
		}
		else if (CALCULATION_MODE == 1)
		{
			for (i = -2; i < 3; i++)
			{
				//printf("%lf\t%lf\t%lf\n", (double)kp_int[0]+0.5*i, (double)kp_int[1]+0.5*j, (double)kp_int[2]);
				fprintf(out, "%lf\t%lf\t%lf\t1.00000\n", ((double)kp_int[0] + i) / KPT_X, ((double)kp_int[1]) / KPT_Y, ((double)kp_int[2]) / KPT_Z);

			}
		}


		fclose(out);

		/*
		for (i=0;i<rep.num_kpt;i++)
		{
			kp_temp[0]=rep.kp_p[i][0]*d_temp[0][0]+rep.kp_p[i][1]*d_temp[1][0]+rep.kp_p[i][2]*d_temp[2][0];
			kp_temp[1]=rep.kp_p[i][0]*d_temp[0][1]+rep.kp_p[i][1]*d_temp[1][1]+rep.kp_p[i][2]*d_temp[2][1];
			kp_temp[2]=rep.kp_p[i][0]*d_temp[0][2]+rep.kp_p[i][1]*d_temp[1][2]+rep.kp_p[i][2]*d_temp[2][2];
			printf("%d\t%d\t%d\n", int(kp_temp[0]*26+0.01), int(kp_temp[1]*26+0.01), int(kp_temp[2]*1+0.01));
		}
		*/
	}

	return 0;
}









int check_now_erfen_pk(string path, int last_flag, int kpt_flag)
{
	//作用是检查当前的vbm和cbm是否在中间的位置上，是返回1，不是返回0
	//需要区分产生k点交叠的情况
	//设置两个情况，一个是之前的读取，一个是最后的读取
	int i, j, k;
	double a1, a2, a3;
	double b;
	report rep;
	atom a;
	rep.read_report(const_cast<char*>((path + "/REPORT").c_str()), 1);

	a.read_atom(const_cast<char*>((path + "/atom.config").c_str()));
	if (CALCULATION_MODE == 2)
	{
		if (rep.eig[rep.num_kpt*rep.num_electron*rep.spin / 2].kpt == 12 && rep.eig[rep.num_kpt*rep.num_electron*rep.spin / 2 - 1].kpt == 37)
			return 1;
	}
	else if (CALCULATION_MODE == 1)
	{
		if (rep.eig[rep.num_kpt*rep.num_electron*rep.spin / 2].kpt == 2 && rep.eig[rep.num_kpt*rep.num_electron*rep.spin / 2 - 1].kpt == 7)
			return 1;
	}
	else
	{
		cout << "calculaiton do not work on check now erfen!" << endl;
		cin.get();
	}
	//注意这里我们人工已经调整了精度，存在序号不对的情况，分情况去除
	int vbm = 5, cbm = 5;
	vbm = rep.eig[rep.num_kpt*rep.num_electron*rep.spin / 2].kpt;//应该是12
	cbm = rep.eig[rep.num_kpt*rep.num_electron*rep.spin / 2 - 1].kpt;//应该是37
	//在这里我们需要核实是不是二分法正确的
	ifstream fin;
	if (last_flag == 1)
	{
		fin.open((path + "/IN.KPT").c_str(), ios::in);
		if (!fin.is_open())
		{
			fin.open((path + "/IN.KPT_" + to_string(kpt_flag)).c_str(), ios::in);
		}
	}
	else if (last_flag == 0)
	{
		fin.open((path + "/IN.KPT_" + to_string(kpt_flag)).c_str(), ios::in);
	}
	if (!fin.is_open())
	{
		cout << "no report!" << path << endl;
		cin.get();
	}
	//在这里我们尝试使用map类键值对
	map<int, string> sure;
	char buffer[100];
	for (int i = 1; i < 52 && fin.peek() != EOF; i++)
	{
		fin.getline(buffer, 99);
		//cout << buffer << endl;
		sure[i] = string(buffer);
	}
	fin.close();
	vector<int>  vbm_xuhao;
	vector<int> cbm_xuhao;
	if (CALCULATION_MODE == 2)
	{
		vbm_xuhao = { 7,8,9,12,13,14,17,18,19 };
		cbm_xuhao = { 32,33,34,37,38,39,42,43,44 };
	}
	else if (CALCULATION_MODE == 1)
	{
		vbm_xuhao = { 2,3,4 };
		cbm_xuhao = { 7,8,9 };
	}
	//cout << vbm << ":" << cbm << endl;
	bool vbm_ok = false, cbm_ok = false;
	for (auto x : vbm_xuhao)
	{
		if (sure[x + 2] == sure[vbm + 3])
		{
			vbm_ok = true;
			//cout << x + 2 << endl;
			//cout << "vbm_ok:"<<sure[x + 2] << endl;
		}
	}

	for (auto x : cbm_xuhao)
	{
		if (sure[x + 2] == sure[cbm + 3])
		{
			cbm_ok = true;
			//cout << x + 2 << endl;
			//cout << "cbm_ok:" << sure[x + 2] << endl;
		}
	}
	if (vbm_ok&& cbm_ok)
	{
		return 1;
	}
	else
	{
		//printf("Not right VBM or CBM\t %d %d\n", rep.eig[rep.num_kpt*rep.num_electron*rep.spin / 2].kpt, rep.eig[rep.num_kpt*rep.num_electron*rep.spin / 2 - 1].kpt);
		return 0;
	}




}
double report::get_kpoint_bandgap(int kpt) {
	//输出指定k点对应的带隙大小
	double energy_1 = this->eig_kns[kpt][this->num_electron/2][0];
	double energy_2 = this->eig_kns[kpt][this->num_electron/2-1][0];
	return energy_1 - energy_2;
}
void report::get_valid_m(int calculation_mode) {
	//输出有效质量,注意区分1D，2D还有3D的情况
	int band_vbm = eig[num_kpt*num_electron*spin / 2].band;
	int band_cbm = eig[num_kpt*num_electron*spin / 2 - 1].band;
	int spin_vbm = eig[num_kpt*num_electron*spin / 2].spin;
	int spin_cbm = eig[num_kpt*num_electron*spin / 2 - 1].spin;
	int vbm = 5, cbm = 5;
	vbm = eig[num_kpt*num_electron*spin / 2].kpt;//应该是12
	cbm = eig[num_kpt*num_electron*spin / 2 - 1].kpt;//应该是37
	int top = 0, botom = 0, left = 0, right = 0, left_top = 0, right_top = 0, left_end = 0, right_end = 0;
	cout << "vbm: " << vbm << "," << "cbm: " << cbm << endl;
	double ex, ey, vx, vy;
	int**inter = new int*[6];
	for (int i = 0; i < 6; i++) {
		inter[i] = new int[9];
	}
	double a1 = 0, a2 = 0, a3 = 0, b = 0;
	double averge_x = 0, averge_y = 0;


	if (calculation_mode == 2)
	{
		//下面要区分这两种情况

		top = vbm - 5;
		botom = vbm + 5;
		left = vbm - 1;
		right = vbm + 1;
		left_top = top - 1;
		right_top = top + 1;
		left_end = botom - 1;
		right_end = botom + 1;
		int vbm_location[9] = { left_end,botom,right_end,left,vbm,right,left_top,top,right_top };

		//填充进每一个inter的相关，便于后面输出
		for (int i = 0; i < 9; i++) {
			inter[0][i] = kp_p[vbm_location[i]][0];
			inter[1][i] = kp_p[vbm_location[i]][1];
		}

		//下面是快速获得x和y方向的
		inter[2][0] = eig_kns[left_end][band_vbm][spin_vbm] / HARTREE;
		inter[2][1] = eig_kns[botom][band_vbm][spin_vbm] / HARTREE;
		inter[2][2] = eig_kns[right_end][band_vbm][spin_vbm] / HARTREE;

		inter[2][3] = eig_kns[left][band_vbm][spin_vbm] / HARTREE;
		inter[2][4] = eig_kns[vbm][band_vbm][spin_vbm] / HARTREE;
		inter[2][5] = eig_kns[right][band_vbm][spin_vbm] / HARTREE;

		inter[2][6] = eig_kns[left_top][band_vbm][spin_vbm] / HARTREE;
		inter[2][7] = eig_kns[top][band_vbm][spin_vbm] / HARTREE;
		inter[2][8] = eig_kns[right_top][band_vbm][spin_vbm] / HARTREE;

		//注意这里我们根据新的修改后的方法，利用真实坐标距离求解计算
		/*a1 0.2 0.2 0->实际坐标
			a2 0.2 0.3 0
			a3 0.2 0.1 0

			a1 0.2 0.1 0
			a2 0.2 0.2 0
			a3 0.2 0.3 0
			b为转化为实际坐标的距离*/
			//同时输出变换之后的坐标位置

			//然后是x的距离，y的距离
		averge_x = (dis(kp_p[vbm], kp_p[left]) + dis(kp_p[vbm], kp_p[right])) / 2;
		averge_y = (dis(kp_p[vbm], kp_p[top]) + dis(kp_p[vbm], kp_p[botom])) / 2;
		//注意价带中是空穴，导带中是电子
		//首先是电子的有效质量，分为两个方向，先y后x
		//这里认为是同一个band

		b = averge_y;
		a1 = eig_kns[top][band_vbm][spin_vbm] / HARTREE;
		a2 = eig_kns[vbm][band_vbm][spin_vbm] / HARTREE;
		a3 = eig_kns[botom][band_vbm][spin_vbm] / HARTREE;
		//printf("%lf\t", 1 / ((a1 + a3 - 2 * a2) / b / b));


		ey = 1 / ((a1 + a3 - 2 * a2) / b / b);
		cout << "ey: " << ey << endl;
		b = averge_x;
		a1 = eig_kns[left][band_vbm][spin_vbm] / HARTREE;
		a2 = eig_kns[vbm][band_vbm][spin_vbm] / HARTREE;
		a3 = eig_kns[right][band_vbm][spin_vbm] / HARTREE;
		//printf("%lf\t", 1 / ((a1 + a3 - 2 * a2) / b / b));
		ex = 1 / ((a1 + a3 - 2 * a2) / b / b);
		cout << "ex: " << ex << endl;

		//然后同理是针对cbm的
		top = cbm - 5;
		botom = cbm + 5;
		left = cbm - 1;
		right = cbm + 1;
		left_top = top - 1;
		right_top = top + 1;
		left_end = botom - 1;
		right_end = botom + 1;

		int cbm_location[9] = { left_end,botom,right_end,left,cbm,right,left_top,top,right_top };
		for (int i = 0; i < 9; i++) {
			inter[3][i] = kp_p[cbm_location[i]][0];
			inter[4][i] = kp_p[cbm_location[i]][1];
		}


		//然后填充inter		
		inter[5][0] = eig_kns[left_end][band_cbm][spin_cbm] / HARTREE;
		inter[5][1] = eig_kns[botom][band_cbm][spin_cbm] / HARTREE;
		inter[5][2] = eig_kns[right_end][band_cbm][spin_cbm] / HARTREE;

		inter[5][3] = eig_kns[left][band_cbm][spin_cbm] / HARTREE;
		inter[5][4] = eig_kns[cbm][band_cbm][spin_cbm] / HARTREE;
		inter[5][5] = eig_kns[right][band_cbm][spin_cbm] / HARTREE;

		inter[5][6] = eig_kns[left_top][band_cbm][spin_cbm] / HARTREE;
		inter[5][7] = eig_kns[top][band_cbm][spin_cbm] / HARTREE;
		inter[5][8] = eig_kns[right_top][band_cbm][spin_cbm] / HARTREE;


		//然后是x的距离，y的距离
		//注意这里的x和y的距离注意退化到实空间去做



		//投影后的x和y

		averge_x = (dis(kp_p[cbm], kp_p[left]) + dis(kp_p[cbm], kp_p[right])) / 2;
		averge_y = (dis(kp_p[cbm], kp_p[top]) + dis(kp_p[cbm], kp_p[botom])) / 2;
		b = averge_y;
		a1 = eig_kns[top][band_cbm][spin_cbm] / HARTREE;
		a2 = eig_kns[cbm][band_cbm][spin_cbm] / HARTREE;
		a3 = eig_kns[botom][band_cbm][spin_cbm] / HARTREE;
		//printf("%lf\t", 1 / ((a1 + a3 - 2 * a2) / b / b));
		vy = 1 / ((a1 + a3 - 2 * a2) / b / b);
		cout << "vy: " << vy << endl;
		b = averge_x;
		a1 = eig_kns[left][band_cbm][spin_cbm] / HARTREE;
		a2 = eig_kns[cbm][band_cbm][spin_cbm] / HARTREE;
		a3 = eig_kns[right][band_cbm][spin_cbm] / HARTREE;
		//printf("%lf\n", 1 / ((a1 + a3 - 2 * a2) / b / b));
		vx = 1 / ((a1 + a3 - 2 * a2) / b / b);
		cout << "vx" << vx << endl;
		//cin.get();

		return;
	}
	else if (calculation_mode == 1)
	{
		//注意是1d的情况时候没有插值，就是两个质量

		left = vbm - 1;
		right = vbm + 1;
		a2 = eig_kns[vbm][band_vbm][spin_vbm] / HARTREE;


		if (left == -1)
		{
			left = num_kpt - 1;
			a3 = eig_kns[right][band_vbm][spin_vbm] / HARTREE;
			a1 = a3;
		}

		else if (right == num_kpt)
		{
			right = 0;
			a1 = eig_kns[left][band_vbm][spin_vbm] / HARTREE;
			a3 = a1;
		}
		else
		{
			a1 = eig_kns[left][band_vbm][spin_vbm] / HARTREE;
			a3 = eig_kns[right][band_vbm][spin_vbm] / HARTREE;
		}
		/*a1 = eig_kns[left][band_vbm][spin_vbm] / HARTREE;
		a3 = eig_kns[right][band_vbm][spin_vbm] / HARTREE;*/

		inter[0][3] = eig_kns[left][band_vbm][spin_vbm] / HARTREE;
		inter[0][4] = eig_kns[vbm][band_vbm][spin_vbm] / HARTREE;
		inter[0][5] = eig_kns[right][band_vbm][spin_vbm] / HARTREE;
		inter[0][9] = dis(kp_p[vbm], kp_p[right]);
		cout << kp_p[left][0] << "," << kp_p[left][1] << "," << kp_p[left][2] << endl;
		cout << kp_p[vbm][0] << "," << kp_p[vbm][1] << "," << kp_p[vbm][2] << endl;

		b = min(dis(kp_p[left], kp_p[vbm]), dis(kp_p[vbm], kp_p[right]));
		//printf("%lf\t", 1 / ((a1 + a3 - 2 * a2) / b / b));
		cout << a1 << "," << a2 << "," << a3 << "," << b << endl;
		ex = 1 / ((a1 + a3 - 2 * a2) / b / b);
		cout << "ex: " << ex << endl;
		left = cbm - 1;
		right = cbm + 1;
		a2 = eig_kns[cbm][band_cbm][spin_cbm] / HARTREE;
		if (left == -1)
			left = num_kpt - 1;
		if (right == num_kpt)
			right = 0;
		a1 = eig_kns[left][band_cbm][spin_cbm] / HARTREE;
		a3 = eig_kns[right][band_cbm][spin_cbm] / HARTREE;

		inter[1][3] = eig_kns[left][band_cbm][spin_cbm] / HARTREE;
		inter[1][4] = eig_kns[cbm][band_cbm][spin_cbm] / HARTREE;
		inter[1][5] = eig_kns[right][band_cbm][spin_cbm] / HARTREE;
		inter[1][9] = dis(kp_p[cbm], kp_p[right]);
		b = max(dis(kp_p[left], kp_p[cbm]), dis(kp_p[cbm], kp_p[right]));

		//printf("%lf\n", 1 / ((a1 + a3 - 2 * a2) / b / b));
		//cout << a1 << "," << a2 << "," << a3 << "," << b << endl;
		vx = 1 / ((a1 + a3 - 2 * a2) / b / b);
		cout << "vx: " << vx << endl;
		//cin.get();
		return;
	}
	else if (calculation_mode == 3) {
		//针对三维情况下的获取有效质量
		//首先是vbm的
		int row = mp_x;
		int z_top = vbm + 1;
		int z_botom = vbm - 1;
		int x_plus = vbm + row * row;
		int x_minux = vbm - row * row;
		int y_plus = vbm + row;
		int y_minus = vbm - row;
		double e_zbotom = 0, e_xminux = 0, e_yminus = 0;
		double b_x = 0, b_y = 0, b_z = 0;
		double e_vbm = eig_kns[vbm][band_vbm][spin_vbm] / HARTREE;
		double e_ztop = eig_kns[z_top][band_vbm][spin_vbm] / HARTREE;
		double e_xplus = eig_kns[x_plus][band_vbm][spin_vbm] / HARTREE;
		double e_yplus = eig_kns[y_plus][band_vbm][spin_vbm] / HARTREE;
		if (z_botom < 0)
		{
			e_zbotom = eig_kns[z_top][band_vbm][spin_vbm] / HARTREE;
			b_z = dis(kp_p[vbm], kp_p[z_top]);
		}
		else
		{
			e_zbotom = eig_kns[z_botom][band_vbm][spin_vbm] / HARTREE;
			b_z = dis(kp_p[vbm], kp_p[z_botom]);
		}

		if (x_minux < 0)
		{
			e_xminux = eig_kns[x_plus][band_vbm][spin_vbm] / HARTREE;
			b_x = dis(kp_p[vbm], kp_p[x_plus]);
		}
		else
		{
			e_xminux = eig_kns[x_minux][band_vbm][spin_vbm] / HARTREE;
			b_x = dis(kp_p[vbm], kp_p[x_minux]);
		}

		if (y_minus < 0) {
			e_yminus = eig_kns[y_plus][band_vbm][spin_vbm] / HARTREE;
			b_y = dis(kp_p[vbm], kp_p[y_plus]);
		}
		else
		{
			e_yminus = eig_kns[y_minus][band_vbm][spin_vbm] / HARTREE;
			b_y = dis(kp_p[vbm], kp_p[y_minus]);
		}

		double valid_x = 0, valid_y = 0, valid_z = 0;
		valid_x = 1 / ((e_xplus + e_xminux - 2 * e_vbm) / b_x / b_x);
		valid_y = 1 / ((e_yplus + e_yminus - 2 * e_vbm) / b_y / b_y);
		valid_z = 1 / ((e_ztop + e_zbotom - 2 * e_vbm) / b_z / b_z);
		cout << "vb:  x:" << valid_x << " y:" << valid_y << " z:" << valid_z << endl;

		//同理开始是cbm的
		z_top = cbm + 1;
		z_botom = cbm - 1;
		x_plus = cbm + row * row;
		x_minux = cbm - row * row;
		y_plus = cbm + row;
		y_minus = cbm - row;
		e_zbotom = 0, e_xminux = 0, e_yminus = 0;
		b_x = 0, b_y = 0, b_z = 0;
		double e_cbm = eig_kns[cbm][band_cbm][spin_cbm] / HARTREE;
		e_ztop = eig_kns[z_top][band_cbm][spin_cbm] / HARTREE;
		e_xplus = eig_kns[x_plus][band_cbm][spin_cbm] / HARTREE;
		e_yplus = eig_kns[y_plus][band_cbm][spin_cbm] / HARTREE;
		if (z_botom < 0)
		{
			e_zbotom = eig_kns[z_top][band_cbm][spin_cbm] / HARTREE;
			b_z = dis(kp_p[cbm], kp_p[z_top]);
		}
		else
		{
			e_zbotom = eig_kns[z_botom][band_cbm][spin_cbm] / HARTREE;
			b_z = dis(kp_p[cbm], kp_p[z_botom]);
		}

		if (x_minux < 0)
		{
			e_xminux = eig_kns[x_plus][band_cbm][spin_cbm] / HARTREE;
			b_x = dis(kp_p[cbm], kp_p[x_plus]);
		}
		else
		{
			e_xminux = eig_kns[x_minux][band_cbm][spin_cbm] / HARTREE;
			b_x = dis(kp_p[cbm], kp_p[x_minux]);
		}

		if (y_minus < 0) {
			e_yminus = eig_kns[y_plus][band_cbm][spin_cbm] / HARTREE;
			b_y = dis(kp_p[cbm], kp_p[y_plus]);
		}
		else
		{
			e_yminus = eig_kns[y_minus][band_cbm][spin_cbm] / HARTREE;
			b_y = dis(kp_p[cbm], kp_p[y_minus]);
		}

		valid_x = 0, valid_y = 0, valid_z = 0;
		valid_x = 1 / ((e_xplus + e_xminux - 2 * e_cbm) / b_x / b_x);
		valid_y = 1 / ((e_yplus + e_yminus - 2 * e_cbm) / b_y / b_y);
		valid_z = 1 / ((e_ztop + e_zbotom - 2 * e_cbm) / b_z / b_z);
		cout << "cb:  x:" << valid_x << " y:" << valid_y << " z:" << valid_z << endl;



	}

}


//重载的版本，输入路径，但是只打印出来相关结果,有重复代码了，需要注意了…………
void get_last_result(string path, vector<string>&step, vector<vector<double>>&table,int flag)
{
	//中间里面的是步长的记录

	int i, j, k;
	double a1, a2, a3;
	double b;
	report rep;
	cell cell_a(const_cast<char*>((path + "/atom.config").c_str()));
	double gap = 0;
	double ex = 0, vx = 0;
	if (flag == 1)
	{
		if (rep.read_report(const_cast<char*>((path + "/REPORT").c_str()), spin_check) == 1)
		{
			//说明没这个文件
			cout << "no report!" << endl;
			return;
		}

	}
	else if (flag == 2)
	{
		//读2的的话，视为30，只读取带隙
		if (rep.read_report(const_cast<char*>((path + "/REPORT").c_str()), spin_check) == 1)
		{
			//说明没这个文件
			cout << "no report!" << endl;
			cin.get();
			return;
		}
		double gap = rep.get_gap(path, step);
		cout << "gap: " << gap << endl;
		return;
	}
	else
	{
		cout << "unkonwn flag ,check!" << endl;
		cin.get();
	}
	gap = rep.get_gap(path, step);
	int vbm = 5, cbm = 5;
	vbm = rep.eig[rep.num_kpt*rep.num_electron*rep.spin / 2].kpt;//应该是12
	cbm = rep.eig[rep.num_kpt*rep.num_electron*rep.spin / 2 - 1].kpt;//应该是37
	//cout << vbm << "," << cbm << endl;
	//下面开始输出相关结果
	if (CALCULATION_MODE == 1) {
		int top = 0, botom = 0, left = 0, right = 0, left_top = 0, right_top = 0, left_end = 0, right_end = 0;
		
		//然后分别针对vbm和cbm进行输出
		int band_vbm = rep.eig[rep.num_kpt*rep.num_electron*rep.spin / 2].band;
		int band_cbm = rep.eig[rep.num_kpt*rep.num_electron*rep.spin / 2 - 1].band;
		int spin_vbm = rep.eig[rep.num_kpt*rep.num_electron*rep.spin / 2].spin;
		int spin_cbm = rep.eig[rep.num_kpt*rep.num_electron*rep.spin / 2 - 1].spin;
		//开始进行记录table
		for (int i = 0; i < Electro_CAL_NUM; i++) {
			int temp_index = 0;
			if (i < Electro_CAL_NUM / 2) {
				temp_index = vbm - Electro_CAL_NUM / 4 + i%(Electro_CAL_NUM / 2);
				if (temp_index < 0 || temp_index >= rep.num_kpt)
					temp_index = 2 * vbm - temp_index;
				for (int j = 0; j < 3; j++)
					table[i][j] = rep.kp_p[temp_index][j];
				table[i][3] = rep.eig_kns[temp_index][band_vbm][spin_vbm] / HARTREE;
			}
			else
			{
				temp_index = cbm - Electro_CAL_NUM / 4 + i % (Electro_CAL_NUM / 2);
				if (temp_index < 0 || temp_index >= rep.num_kpt)
					temp_index = 2 * cbm - temp_index;
				for (int j = 0; j < 3; j++)
					table[i][j] = rep.kp_p[temp_index][j];
				table[i][3] = rep.eig_kns[temp_index][band_cbm][spin_cbm] / HARTREE;
			}
			
		}




		left = vbm - 1;
		right = vbm + 1;
		a2 = rep.eig_kns[vbm][band_vbm][spin_vbm] / HARTREE;

		if (left == -1)
		{
			left = vbm + 1;
		}
		else if (right == rep.num_kpt)
		{
			right = vbm - 1;
		}
		a1 = rep.eig_kns[left][band_vbm][spin_vbm] / HARTREE;
		a3 = rep.eig_kns[right][band_vbm][spin_vbm] / HARTREE;

		b = min(dis(rep.kp_p[left], rep.kp_p[vbm]), dis(rep.kp_p[vbm], rep.kp_p[right]));
		cout << rep.kp_p[left][0] << "," << rep.kp_p[right][0] << endl;
		//printf("%lf\t", 1 / ((a1 + a3 - 2 * a2) / b / b));
		//cout << a1 << "," << a2 << "," << a3 << "," << b << endl;
		ex = 1 / ((a1 + a3 - 2 * a2) / b / b);

		left = cbm - 1;
		right = cbm + 1;
		a2 = rep.eig_kns[cbm][band_cbm][spin_cbm] / HARTREE;
		if (left == -1)
			left = cbm + 1;
		if (right == rep.num_kpt)
			right = cbm - 1;
		a1 = rep.eig_kns[left][band_cbm][spin_cbm] / HARTREE;
		a3 = rep.eig_kns[right][band_cbm][spin_cbm] / HARTREE;

		
		b = min(dis(rep.kp_p[left], rep.kp_p[cbm]), dis(rep.kp_p[cbm], rep.kp_p[right]));

		//printf("%lf\n", 1 / ((a1 + a3 - 2 * a2) / b / b));
		//cout << a1 << "," << a2 << "," << a3 << "," << b << endl;
		vx = 1 / ((a1 + a3 - 2 * a2) / b / b);
		//cin.get();
		
	}
	else
	{
		cout << "unsupported calculation mode num!" << endl;
		cin.get();
		return;
	}
	//最后打印出来
	cout << "gap:" << gap << ",ex:" << ex << ",vx:" << vx << endl;
	return;
}

void get_last_result(string path, double& gap, double &ey, double& ex, double& vy, double& vx, int kpt_flag, double** inter, vector<string>&step,int flag)
{
	int i, j, k;
	double a1, a2, a3;
	double b;
	report rep;
	cell cell_a(const_cast<char*>((path + "/atom.config").c_str()));
	if (flag == 1)
	{
		if (rep.read_report(const_cast<char*>((path + "/REPORT").c_str()), spin_check) == 1)
		{
			//说明没这个文件
			cout << "no report!" << endl;
			cin.get();
			gap = ex = ey = vy = vx = -100;
			return;
		}

	}
	else if (flag == 2)
	{
		//读2的的话，视为30，只读取带隙
		if (rep.read_report(const_cast<char*>((path + "/REPORT").c_str()), spin_check) == 1)
		{
			//说明没这个文件
			cout << "no report!" << endl;
			cin.get();
			gap = ex = ey = vy = vx = -100;
			return;
		}
		gap=rep.get_gap(path,step);
		return;
	}
	else
	{
		cout << "unkonwn flag ,check!" << endl;
		cin.get();
	}	
	gap = rep.get_gap(path,step);
	int vbm = 5, cbm = 5;
	vbm = rep.eig[rep.num_kpt*rep.num_electron*rep.spin / 2].kpt;//应该是12
	cbm = rep.eig[rep.num_kpt*rep.num_electron*rep.spin / 2 - 1].kpt;//应该是37
	cout << vbm << "," << cbm << endl;


	//在这里我们需要核实是不是二分法正确的
	//只有在非3D的情况下才会检视IN.KPT文件,并且迭代方法为true时候才检查
	bool vbm_ok = false, cbm_ok = false;
	if (CALCULATION_MODE != 3&&IF_ITERATION) {
		ifstream fin;
		fin.open((path + "/IN.KPT").c_str(), ios::in);
		if (!fin.is_open())
		{
			cout << "i can not find the file" << path + "/IN.KPT" << endl;
			fin.open((path + "/IN.KPT_" + to_string(kpt_flag)).c_str(), ios::in);
			if (!fin.is_open())
			{
				cout << "wrong!" << path << endl;
				cin.get();
			}

		}
		//在这里我们尝试使用map类键值对
		map<int, string> sure;
		char buffer[100];
		for (int i = 1; i < 52 && fin.peek() != EOF; i++)
		{
			fin.getline(buffer, 99);
			//cout << buffer << endl;
			sure[i] = string(buffer);
		}
		fin.close();
		vector<int>  vbm_xuhao;
		vector<int> cbm_xuhao;
		if (CALCULATION_MODE == 2)
		{
			vbm_xuhao = { 7,8,9,12,13,14,17,18,19 };
			cbm_xuhao = { 32,33,34,37,38,39,42,43,44 };
		}
		else if (CALCULATION_MODE == 1)
		{
			vbm_xuhao = { 2,3,4 };
			cbm_xuhao = { 7,8,9 };
		}
		//cout << vbm << ":" << cbm << endl;

		for (auto x : vbm_xuhao)
		{
			if (sure[x + 2] == sure[vbm + 3])
			{
				vbm_ok = true;
				vbm = x - 1;
				/*cout << x + 2 << endl;
				cout << "vbm_ok:"<<x-1 << endl;*/
			}
		}

		for (auto x : cbm_xuhao)
		{
			if (sure[x + 2] == sure[cbm + 3])
			{
				cbm_ok = true;
				cbm = x - 1;
				/*	cout << x + 2 << endl;
					cout << "cbm_ok:" << x-1<< endl;*/
			}
		}
	}

	//下面开始输出,输出的条件是，要么发现点合适，要么根本不用迭代
	if ((vbm_ok&& cbm_ok) || IF_ITERATION==false)
	{
		int top = 0, botom = 0, left = 0, right = 0, left_top = 0, right_top = 0, left_end = 0, right_end = 0;

		//然后分别针对vbm和cbm进行输出
		int band_vbm = rep.eig[rep.num_kpt*rep.num_electron*rep.spin / 2].band;
		int band_cbm = rep.eig[rep.num_kpt*rep.num_electron*rep.spin / 2 - 1].band;
		int spin_vbm = rep.eig[rep.num_kpt*rep.num_electron*rep.spin / 2].spin;
		int spin_cbm = rep.eig[rep.num_kpt*rep.num_electron*rep.spin / 2 - 1].spin;

		if (CALCULATION_MODE == 0) {
		
			//0d的话什么都不用干，gap上面已经读出来了

			return;
		}
		else if (CALCULATION_MODE == 2)
		{
			//下面要区分这两种情况
			double averge_x = 0, averge_y = 0;
			int hang = 0;
			int max_num = 0;
			if (cell_a.num < 10) {
				hang = 6;
				max_num = 36;
			}
			else if (cell_a.num >= 10) {
				hang = 2;
				max_num = 4;
			}
			top = vbm - hang;
			botom = vbm + hang;
			left = vbm - 1;
			right = vbm + 1;
			left_top = top - 1;
			right_top = top + 1;
			left_end = botom - 1;
			right_end = botom + 1;
			int vbm_location[9] = { left_end,botom,right_end,left,vbm,right,left_top,top,right_top };
			//对这些点进行映射，防止到达边界
			for (int i = 0; i < 9; i++) {
				if (vbm_location[i] < 0) {
					vbm_location[i] = 0;
				}
				else if (vbm_location[i] > (max_num - 1)) {
					vbm_location[i] = max_num - 1;
				}
			}
			//填充进每一个inter的相关，便于后面输出
			for (int i = 0; i < 9; i++) {
				inter[0][i] = rep.kp_p[vbm_location[i]][0];
				inter[1][i] = rep.kp_p[vbm_location[i]][1];
			}
			//下面是快速获得x和y方向的
			inter[2][0] = rep.eig_kns[vbm_location[0]][band_vbm][spin_vbm] / HARTREE;
			inter[2][1] = rep.eig_kns[vbm_location[1]][band_vbm][spin_vbm] / HARTREE;
			inter[2][2] = rep.eig_kns[vbm_location[2]][band_vbm][spin_vbm] / HARTREE;

			inter[2][3] = rep.eig_kns[vbm_location[3]][band_vbm][spin_vbm] / HARTREE;
			inter[2][4] = rep.eig_kns[vbm_location[4]][band_vbm][spin_vbm] / HARTREE;
			inter[2][5] = rep.eig_kns[vbm_location[5]][band_vbm][spin_vbm] / HARTREE;

			inter[2][6] = rep.eig_kns[vbm_location[6]][band_vbm][spin_vbm] / HARTREE;
			inter[2][7] = rep.eig_kns[vbm_location[7]][band_vbm][spin_vbm] / HARTREE;
			inter[2][8] = rep.eig_kns[vbm_location[8]][band_vbm][spin_vbm] / HARTREE;

				//然后是x的距离，y的距离
			averge_x = (dis(rep.kp_p[vbm], rep.kp_p[vbm_location[3]]) + dis(rep.kp_p[vbm], rep.kp_p[vbm_location[5]])) / 2;
			averge_y = (dis(rep.kp_p[vbm], rep.kp_p[vbm_location[7]]) + dis(rep.kp_p[vbm], rep.kp_p[vbm_location[1]])) / 2;
			//注意价带中是空穴，导带中是电子
			//首先是电子的有效质量，分为两个方向，先y后x
			//这里认为是同一个band

			b = averge_y;
			a1 = rep.eig_kns[vbm_location[7]][band_vbm][spin_vbm] / HARTREE;
			a2 = rep.eig_kns[vbm][band_vbm][spin_vbm] / HARTREE;
			a3 = rep.eig_kns[vbm_location[1]][band_vbm][spin_vbm] / HARTREE;
			//printf("%lf\t", 1 / ((a1 + a3 - 2 * a2) / b / b));


			ey = 1 / ((a1 + a3 - 2 * a2) / b / b);
			b = averge_x;
			a1 = rep.eig_kns[vbm_location[3]][band_vbm][spin_vbm] / HARTREE;
			a2 = rep.eig_kns[vbm][band_vbm][spin_vbm] / HARTREE;
			a3 = rep.eig_kns[vbm_location[5]][band_vbm][spin_vbm] / HARTREE;
			//printf("%lf\t", 1 / ((a1 + a3 - 2 * a2) / b / b));
			ex = 1 / ((a1 + a3 - 2 * a2) / b / b);


			//然后同理是针对cbm的
			top = cbm - hang;
			botom = cbm + hang;
			left = cbm - 1;
			right = cbm + 1;
			left_top = top - 1;
			right_top = top + 1;
			left_end = botom - 1;
			right_end = botom + 1;

			int cbm_location[9] = { left_end,botom,right_end,left,cbm,right,left_top,top,right_top };
			//对这些点进行映射，防止到达边界
			for (int i = 0; i < 9; i++) {
				if (cbm_location[i] < 0) {
					cbm_location[i] = 0;
				}
				else if (cbm_location[i] > (max_num - 1)) {
					cbm_location[i] = max_num - 1;
				}
			}
			for (int i = 0; i < 9; i++) {
				inter[3][i] = rep.kp_p[cbm_location[i]][0];
				inter[4][i] = rep.kp_p[cbm_location[i]][1];
			}
			
			//填充完整版的信息
			/*fullfill_inter_real_location(inter, left_end, cell_a, path, 0, band_cbm, spin_cbm, "cb");
			fullfill_inter_real_location(inter, botom, cell_a, path, 1, band_cbm, spin_cbm, "cb");
			fullfill_inter_real_location(inter, right_end, cell_a, path, 2, band_cbm, spin_cbm, "cb");
			fullfill_inter_real_location(inter, left, cell_a, path, 3, band_cbm, spin_cbm, "cb");
			fullfill_inter_real_location(inter, cbm, cell_a, path, 4, band_cbm, spin_cbm, "cb");
			fullfill_inter_real_location(inter, right, cell_a, path, 5, band_cbm, spin_cbm, "cb");
			fullfill_inter_real_location(inter, left_top, cell_a, path, 6, band_cbm, spin_cbm, "cb");
			fullfill_inter_real_location(inter, top, cell_a, path, 7, band_cbm, spin_cbm, "cb");
			fullfill_inter_real_location(inter, right_top, cell_a, path, 8, band_cbm, spin_cbm, "cb");*/
			//然后填充inter		
			inter[5][0] = rep.eig_kns[cbm_location[0]][band_cbm][spin_cbm] / HARTREE;
			inter[5][1] = rep.eig_kns[cbm_location[1]][band_cbm][spin_cbm] / HARTREE;
			inter[5][2] = rep.eig_kns[cbm_location[2]][band_cbm][spin_cbm] / HARTREE;

			inter[5][3] = rep.eig_kns[cbm_location[3]][band_cbm][spin_cbm] / HARTREE;
			inter[5][4] = rep.eig_kns[cbm_location[4]][band_cbm][spin_cbm] / HARTREE;
			inter[5][5] = rep.eig_kns[cbm_location[5]][band_cbm][spin_cbm] / HARTREE;

			inter[5][6] = rep.eig_kns[cbm_location[6]][band_cbm][spin_cbm] / HARTREE;
			inter[5][7] = rep.eig_kns[cbm_location[7]][band_cbm][spin_cbm] / HARTREE;
			inter[5][8] = rep.eig_kns[cbm_location[8]][band_cbm][spin_cbm] / HARTREE;





			//然后是x的距离，y的距离
			//注意这里的x和y的距离注意退化到实空间去做



			//投影后的x和y

			averge_x = (dis(rep.kp_p[cbm], rep.kp_p[cbm_location[3]]) + dis(rep.kp_p[cbm], rep.kp_p[cbm_location[5]])) / 2;
			averge_y = (dis(rep.kp_p[cbm], rep.kp_p[cbm_location[7]]) + dis(rep.kp_p[cbm], rep.kp_p[cbm_location[1]])) / 2;
			b = averge_y;
			a1 = rep.eig_kns[cbm_location[7]][band_cbm][spin_cbm] / HARTREE;
			a2 = rep.eig_kns[cbm][band_cbm][spin_cbm] / HARTREE;
			a3 = rep.eig_kns[cbm_location[1]][band_cbm][spin_cbm] / HARTREE;
			//printf("%lf\t", 1 / ((a1 + a3 - 2 * a2) / b / b));
			vy = 1 / ((a1 + a3 - 2 * a2) / b / b);
			b = averge_x;
			a1 = rep.eig_kns[cbm_location[3]][band_cbm][spin_cbm] / HARTREE;
			a2 = rep.eig_kns[cbm][band_cbm][spin_cbm] / HARTREE;
			a3 = rep.eig_kns[cbm_location[5]][band_cbm][spin_cbm] / HARTREE;
			//printf("%lf\n", 1 / ((a1 + a3 - 2 * a2) / b / b));
			vx = 1 / ((a1 + a3 - 2 * a2) / b / b);
			/*cout << gap << vx << endl;
			cin.get();*/

			return;
		}
		else if (CALCULATION_MODE == 1)
		{
			//注意是1d的情况时候没有插值，就是两个质量

			left = vbm - 1;
			right = vbm + 1;
			a2 = rep.eig_kns[vbm][band_vbm][spin_vbm] / HARTREE;

			if (left == -1)
			{
				left = vbm + 1;
			}
			else if (right == rep.num_kpt)
			{
				right = vbm - 1;
			}
			a1 = rep.eig_kns[left][band_vbm][spin_vbm] / HARTREE;
			a3 = rep.eig_kns[right][band_vbm][spin_vbm] / HARTREE;

			inter[0][3] = rep.eig_kns[left][band_vbm][spin_vbm] / HARTREE;
			inter[0][4] = rep.eig_kns[vbm][band_vbm][spin_vbm] / HARTREE;
			inter[0][5] = rep.eig_kns[right][band_vbm][spin_vbm] / HARTREE;
			inter[0][9] = dis(rep.kp_p[vbm], rep.kp_p[right]);

			b = min(dis(rep.kp_p[left], rep.kp_p[vbm]), dis(rep.kp_p[vbm], rep.kp_p[right]));
			//printf("%lf\t", 1 / ((a1 + a3 - 2 * a2) / b / b));
			//cout << a1 << "," << a2 << "," << a3 << "," << b << endl;
			ex = 1 / ((a1 + a3 - 2 * a2) / b / b);

			left = cbm - 1;
			right = cbm + 1;
			a2 = rep.eig_kns[cbm][band_cbm][spin_cbm] / HARTREE;
			if (left == -1)
				left = cbm + 1;
			if (right == rep.num_kpt)
				right = cbm - 1;
			a1 = rep.eig_kns[left][band_cbm][spin_cbm] / HARTREE;
			a3 = rep.eig_kns[right][band_cbm][spin_cbm] / HARTREE;

			inter[1][3] = rep.eig_kns[left][band_cbm][spin_cbm] / HARTREE;
			inter[1][4] = rep.eig_kns[cbm][band_cbm][spin_cbm] / HARTREE;
			inter[1][5] = rep.eig_kns[right][band_cbm][spin_cbm] / HARTREE;
			inter[1][9] = dis(rep.kp_p[cbm], rep.kp_p[right]);
			b = min(dis(rep.kp_p[left], rep.kp_p[cbm]), dis(rep.kp_p[cbm], rep.kp_p[right]));

			//printf("%lf\n", 1 / ((a1 + a3 - 2 * a2) / b / b));
			//cout << a1 << "," << a2 << "," << a3 << "," << b << endl;
			vx = 1 / ((a1 + a3 - 2 * a2) / b / b);
			//cin.get();
			return;
		}
		else if (CALCULATION_MODE == 3) {
			//针对三维情况下的获取有效质量
			//首先是vbm的
			return;
			int row = rep.mp_x;
			int z_top = vbm + 1;
			int z_botom = vbm - 1;
			int x_plus = vbm + row * row;
			int x_minux = vbm - row * row;
			int y_plus = vbm + row;
			int y_minus = vbm - row;
			double e_zbotom = 0, e_xminux = 0, e_yminus = 0;
			double b_x = 0, b_y = 0, b_z = 0;
			double e_vbm = rep.eig_kns[vbm][band_vbm][spin_vbm] / HARTREE;
			double e_ztop = rep.eig_kns[z_top][band_vbm][spin_vbm] / HARTREE;
			double e_xplus = rep.eig_kns[x_plus][band_vbm][spin_vbm] / HARTREE;
			double e_yplus = rep.eig_kns[y_plus][band_vbm][spin_vbm] / HARTREE;
			if (z_botom < 0)
			{
				e_zbotom = rep.eig_kns[z_top][band_vbm][spin_vbm] / HARTREE;
				b_z = dis(rep.kp_p[vbm], rep.kp_p[z_top]);
			}
			else
			{
				e_zbotom = rep.eig_kns[z_botom][band_vbm][spin_vbm] / HARTREE;
				b_z = dis(rep.kp_p[vbm], rep.kp_p[z_botom]);
			}

			if (x_minux < 0)
			{
				e_xminux = rep.eig_kns[x_plus][band_vbm][spin_vbm] / HARTREE;
				b_x = dis(rep.kp_p[vbm], rep.kp_p[x_plus]);
			}
			else
			{
				e_xminux = rep.eig_kns[x_minux][band_vbm][spin_vbm] / HARTREE;
				b_x = dis(rep.kp_p[vbm], rep.kp_p[x_minux]);
			}

			if (y_minus < 0) {
				e_yminus = rep.eig_kns[y_plus][band_vbm][spin_vbm] / HARTREE;
				b_y = dis(rep.kp_p[vbm], rep.kp_p[y_plus]);
			}
			else
			{
				e_yminus = rep.eig_kns[y_minus][band_vbm][spin_vbm] / HARTREE;
				b_y = dis(rep.kp_p[vbm], rep.kp_p[y_minus]);
			}

			double valid_x = 0, valid_y = 0, valid_z = 0;
			valid_x = 1 / ((e_xplus + e_xminux - 2 * e_vbm) / b_x / b_x);
			valid_y = 1 / ((e_yplus + e_yminus - 2 * e_vbm) / b_y / b_y);
			valid_z = 1 / ((e_ztop + e_zbotom - 2 * e_vbm) / b_z / b_z);
			cout << "vb:  x:" << valid_x << " y:" << valid_y << " z:" << valid_z << endl;

			//同理开始是cbm的
			z_top = cbm + 1;
			z_botom = cbm - 1;
			x_plus = cbm + row * row;
			x_minux = cbm - row * row;
			y_plus = cbm + row;
			y_minus = cbm - row;
			e_zbotom = 0, e_xminux = 0, e_yminus = 0;
			b_x = 0, b_y = 0, b_z = 0;
			double e_cbm = rep.eig_kns[cbm][band_cbm][spin_cbm] / HARTREE;
			e_ztop = rep.eig_kns[z_top][band_cbm][spin_cbm] / HARTREE;
			e_xplus = rep.eig_kns[x_plus][band_cbm][spin_cbm] / HARTREE;
			e_yplus = rep.eig_kns[y_plus][band_cbm][spin_cbm] / HARTREE;
			if (z_botom < 0)
			{
				e_zbotom = rep.eig_kns[z_top][band_cbm][spin_cbm] / HARTREE;
				b_z = dis(rep.kp_p[cbm], rep.kp_p[z_top]);
			}
			else
			{
				e_zbotom = rep.eig_kns[z_botom][band_cbm][spin_cbm] / HARTREE;
				b_z = dis(rep.kp_p[cbm], rep.kp_p[z_botom]);
			}

			if (x_minux < 0)
			{
				e_xminux = rep.eig_kns[x_plus][band_cbm][spin_cbm] / HARTREE;
				b_x = dis(rep.kp_p[cbm], rep.kp_p[x_plus]);
			}
			else
			{
				e_xminux = rep.eig_kns[x_minux][band_cbm][spin_cbm] / HARTREE;
				b_x = dis(rep.kp_p[cbm], rep.kp_p[x_minux]);
			}

			if (y_minus < 0) {
				e_yminus = rep.eig_kns[y_plus][band_cbm][spin_cbm] / HARTREE;
				b_y = dis(rep.kp_p[cbm], rep.kp_p[y_plus]);
			}
			else
			{
				e_yminus = rep.eig_kns[y_minus][band_cbm][spin_cbm] / HARTREE;
				b_y = dis(rep.kp_p[cbm], rep.kp_p[y_minus]);
			}

			valid_x = 0, valid_y = 0, valid_z = 0;
			valid_x = 1 / ((e_xplus + e_xminux - 2 * e_cbm) / b_x / b_x);
			valid_y = 1 / ((e_yplus + e_yminus - 2 * e_cbm) / b_y / b_y);
			valid_z = 1 / ((e_ztop + e_zbotom - 2 * e_cbm) / b_z / b_z);
			cout << "cb:  x:" << valid_x << " y:" << valid_y << " z:" << valid_z << endl;



		}

	}
	else
	{
		//printf("Not right VBM or CBM\t %d %d\n", rep.eig[rep.num_kpt*rep.num_electron*rep.spin / 2].kpt, rep.eig[rep.num_kpt*rep.num_electron*rep.spin / 2 - 1].kpt);
		gap = ex = ey = vy = vx = -100;
		cout << gap << vx << endl;
		cin.get();
		return;
	}

	//这里同样需要修改，因为考虑到临近数值的问题



}


void get_all_output(string &out_name, string& chazhi_name, string & calculation, string *config, int**flag, int num)
{
	//用来输出相关的结果，包括带隙有效质量等。

	int i = 0, j = 0;
	ofstream fout, fout2;
	double ***inter = new double**[num];
	double** resul = new double*[num];
	for (i = 0; i < num; i++)
	{
		resul[i] = new double[7];
		inter[i] = new double*[6];
		for (j = 0; j < 6; j++)
		{
			inter[i][j] = new double[9];
		}
	}
	//还有存储vbm和cbm的起点的内容
	vector<vector<string>> step(num, vector<string>(2, ""));
	fout2.open(chazhi_name + "_name", ios::out);
	fout.open(out_name, ios::out);
	if (!fout.is_open())
	{
		cout << out_name << "open error! " << endl;
		cin.get();
	}

	/*if(CALCULATION_MODE==1)
		fout << "org_name" << "\t" << "formula_name" << "\t" << "positive_dis" << "\t" << "comment"<<"\t"  <<"Gap" << "\t" << "Electro_x" << "\t"  << "Vacum_x" <<  "\t" << "diff" << "\t"<<"relax_energy"<<"\t"<<endl;
	else*/
	fout << "org_name" << "\t" << "formula_name" << "\t" << "positive_dis" << "\t" << "comment" << "\t" << "Gap" << "\t" << "Electro_y" << "\t" << "Electro_x" << "\t" << "Vacum_y" << "\t" << "Vacum_x" << "\t" << "Electry_got" << "\t" << "Vacum_got" << "\t" << "vbm" << "\t" << "cbm" << "\t" << "diff" << "\t" << "relax_energy" << "\t" << "gamma_gap" << "\t";
	if (IF_output_half == true) {
		fout << "up_gap" << "\t" << "down_gap" << "\t" << "half_flag" << "\t";
	}
	fout << endl;

	for (i = 0; i < num; i++)
	{		
		cell cell_a(const_cast<char*>((calculation + config[i] + "/atom.config").c_str()));
		cout << "has generate:" << config[i] << endl;
		fout << config[i] << "\t";//输出org_name
		fout << cell_a.get_formula_name() << "\t";//输出formula_name
		fout << cell_a.get_positive_distance() << "\t";//输出阳离子之间的距离
		//下面开始打印和commment，以及之后的数据

		
		if (flag[i][0] == 7) {
			fout << "No_yanshi_file!" << "\t";
		}
		else if (flag[i][0] == 17) {
			fout << "Relax_unfinished!" << "\t";			
		}
		else if (flag[i][0] == 37) {
			fout << "Scf_encounter_wrong!" << "\t";
			fout << cell_a.check_energy(calculation + config[i], "REPORT_relax") << "\t";
		}
		else if (flag[i][0] == 57 || flag[i][0] == 67) {
			fout << "Nonscf_encounter_wrong!" << "\t";				
			get_last_result(calculation + config[i], resul[i][0], resul[i][1], resul[i][2], resul[i][3], resul[i][4], flag[i][2], inter[i], step[i],2);
			fout << resul[i][0] << "\t";
			fout << cell_a.check_diff(calculation + config[i], "REPORT") << "\t";
			fout << cell_a.check_energy(calculation + config[i], "REPORT_relax") << "\t";
		}
		else if (flag[i][0] ==30) {
			fout << "Two_cut methond wrong!" << "\t";			
			get_last_result(calculation + config[i], resul[i][0], resul[i][1], resul[i][2], resul[i][3], resul[i][4], flag[i][2], inter[i], step[i],2);
			fout << resul[i][0] << "\t";
			fout << cell_a.check_diff(calculation + config[i], "REPORT") << "\t";
			fout << cell_a.check_energy(calculation + config[i], "REPORT_relax") << "\t";
		}
		else if (flag[i][0] == 20 || flag[i][0] == 47) {
			fout << "No band gap!" << "\t";			
			fout << 0 << "\t";
			fout << cell_a.check_diff(calculation + config[i], "REPORT") << endl;
			fout << cell_a.check_energy(calculation + config[i], "REPORT_relax") << "\t";
		}
		else if (flag[i][0] % 10 == 1 || flag[i][0] == 3) {
			fout << "Calculation_unfinished!wait" << "\t";
		}
		else if (flag[i][0] == 99) {
			fout << "succeed" << "\t";
			fout2 << config[i] << "\n";
			get_last_result(calculation + config[i], resul[i][0], resul[i][1], resul[i][2], resul[i][3], resul[i][4], flag[i][2], inter[i],step[i]);
			if (CALCULATION_MODE == 2)
			{
				resul[i][5] = 2 / (1 / (resul[i][1]) + 1 / (resul[i][2]));
				resul[i][6] = 2 / (1 / (resul[i][3]) + 1 / (resul[i][4]));
			}
			else if (CALCULATION_MODE == 1)
			{
				resul[i][5] = resul[i][6] = -100;
			}
			else if (CALCULATION_MODE == 0) {
				for (int j = 1; j < 7; j++) {
					resul[i][j] = -100;
				}
			}
			for (j = 0; j < 7; j++)
			{
				fout << resul[i][j] << "\t";
			}
			fout << step[i][0] << "\t";
			fout << step[i][1] << "\t";
			fout << cell_a.check_diff(calculation + config[i], "REPORT") << "\t";
			fout << cell_a.check_energy(calculation + config[i], "REPORT_relax") << "\t";
			report rt;
			rt.read_report(const_cast<char*>((calculation+config[i]+"/REPORT").c_str()),2);
			fout << rt.get_kpoint_bandgap(0) << "\t";
			
			//有关半金属的输出
			if (IF_output_half == true) {
				
				vector < double >res = rt.get_half_situation(calculation + config[i] + "/REPORT");
				fout << res[0] << "\t" << res[1] << "\t";
				if (res[0] * res[1] < 0 && res[2]>0) {
					fout << "1";
				}
				else if (res[0] * res[1] > 0 && res[2] > 0)
					fout << "0";
				else if (res[2] < 0)
					fout << "11";
				else {
					cout << "unkonwn the flag!please check!" << endl;
					cin.get();
				}
			}
			
		}
		else
		{
			cout << "unkonwn flag!please make sure!" << config[i] << endl;
			cin.get();
		}
		fout << endl;


	}
	fout.close();
	fout2.close();


	//正式输出得到的九个的数字的坐标
	//cout << "now generate the inter_place number:" << endl;
	//fout2.open(chazhi_name, ios::out);
	////设置定点6位有效数字
	//fout2.setf(ios::fixed);
	//fout2.precision(10);

	//for (i = 0; i < num; i++)
	//{
	//	if (flag[i][0] == 99)
	//	{
	//		cout << "generate the inter use number:" << config[i] << endl;
	//		//现在每个结构会有6行,每行是9个数字
	//		//首先是vbm的三行,分别是x的坐标，y的坐标和z的坐标
	//		//注意填充相关的数据
	//		for (int k = 0; k < 6; k++) {
	//			for (int j = 0; j < 9; j++) {
	//				fout2 << inter[i][k][j] << ",";
	//			}
	//			fout2 << endl;
	//		}
	//	}
	//}
	//fout2.close();
	/*for (i = 0; i <num; i++)
	{
		for (j = 0; j < 2; j++)
		{
			delete[]inter[i][j];
		}
		delete[]inter[i];
		delete[]resul[i];
	}
	delete[]inter;	*/
	delete[]resul;
	cout << "generate the result file completed!" << endl;
	return;
}


void generate_pattern_configs(string filename, string path, string flag, element* e) {
	//输入参数是待定的原始结构文件位置,指定路径和falg标识什么类型结构，对应不同的产生规则
	//输出是往指定路径输出相应结构
	char* target = const_cast<char*>(filename.c_str());
	vector<vector<int>> elements;
	Pattern_cell *pattern = new Pattern_cell(flag, target, 0);
	if (flag == "1T") {

		//构建元素列表		
		vector<int> center = { 22,40,50,72 };
		vector<int> env1 = { 16,34,52 };
		vector<int> center2 = { 12,20,38,56,30,48,80 };
		vector<int> env2 = { 9,17,35,53 };
		for (int i = 0; i < center.size(); i++) {
			for (int j = 0; j < env1.size(); j++) {
				elements.push_back(vector<int>{center[i], env1[j]});
			}
		}
		for (int i = 0; i < center2.size(); i++) {
			for (int j = 0; j < env2.size(); j++) {
				elements.push_back(vector<int>{center2[i], env2[j]});
			}
		}
		double now = e[16].nega_ionic_ridus + e[42].posi_ionic_ridus;
		//cout << now << endl;
		for (int i = 0; i < elements.size(); i++) {

			cell c_pattern(const_cast<char*>(filename.c_str()), 0);
			string dir_name = path + a[elements[i][0]] + a[elements[i][1]] + "2.config";
			double factor = 1;//用于坐标放缩的系数，按照阴离子半径进行放缩

			factor = (e[elements[i][1]].nega_ionic_ridus + e[elements[i][0]].posi_ionic_ridus) / now;
			if (factor < 0.8)
				factor = 1;
			pattern->get_pattern(dir_name, elements[i], c_pattern, factor);

		}


		cout << "generate done!please check!" << endl;
		cin.get();
		return;
	}
	else if (flag == "type6") {
		//这个有两类，分类来做,注意是三元的
		vector<int> Lusu = { 9,17,35,53 };
		vector<int> plus_three = { 33,51,83 };
		vector<int> S_type = { 16,34,52 };

		for (int i = 0; i < Lusu.size(); i++) {
			for (int j = 0; j < Lusu.size(); j++) {
				elements.push_back(vector<int>{50, Lusu[i], Lusu[j]});
			}
		}
		for (int i = 0; i < plus_three.size(); i++) {
			for (int j = 0; j < S_type.size(); j++) {
				for (int k = 0; k < Lusu.size(); k++) {
					elements.push_back(vector<int>{plus_three[i], S_type[j], Lusu[k]});
				}
			}
		}
		for (int i = 0; i < elements.size(); i++) {

			cell c_pattern(const_cast<char*>(filename.c_str()), 0);
			string dir_name = path + a[elements[i][0]] + a[elements[i][1]] + a[elements[i][2]] + "_type6.config";
			double factor = 1;//用于坐标放缩的系数，按照阴离子半径进行放缩
			//factor = (e[elements[i][1]].nega_ionic_ridus + e[elements[i][0]].posi_ionic_ridus) / now;
			if (factor < 0.8)
				factor = 1;
			pattern->get_pattern(dir_name, elements[i], c_pattern, factor, "type6");

		}


	}
	else if (flag == "type10") {
		//同样也是三元的，分成两类
		vector<int> S_type = { 16,34,52 };
		vector<int> plus_three_one = { 31,49,81 };
		vector<int> plus_three_two = { 33,51,83 };
		vector<int >plus_four = { 22,40,72 };
		vector<int> plus_two = { 22,50 };

		for (int i = 0; i < S_type.size(); i++) {
			for (int j = 0; j < plus_three_one.size(); j++) {
				for (int k = 0; k < plus_three_two.size(); k++) {
					elements.push_back(vector<int>{S_type[i], plus_three_one[j], plus_three_two[k]});
				}
			}
		}
		//然后是第二类
		for (int i = 0; i < S_type.size(); i++) {
			for (int j = 0; j < plus_four.size(); j++) {
				for (int k = 0; k < plus_two.size(); k++) {
					elements.push_back(vector<int>{S_type[i], plus_four[j], plus_two[k]});
				}
			}
		}

		for (int i = 0; i < elements.size(); i++) {

			cell c_pattern(const_cast<char*>(filename.c_str()), 0);
			string dir_name = path + a[elements[i][0]] + a[elements[i][1]] + a[elements[i][2]] + "_type10.config";
			double factor = 1;//用于坐标放缩的系数，按照阴离子半径进行放缩
			//factor = (e[elements[i][1]].nega_ionic_ridus + e[elements[i][0]].posi_ionic_ridus) / now;
			if (factor < 0.8)
				factor = 1;
			pattern->get_pattern(dir_name, elements[i], c_pattern, factor, "type10");

		}

	}
	else if (flag == "SnBr2_three") {
		//对应的是三元的情况
		double now = 0;
		vector<int> Lusu = { 9,17,35,53 };
		vector<int> center1 ={ 33,51,83,23,41,73 };
		vector<int> center2 = {16,34,52};
		for (int i = 0; i < center1.size(); i++) {
			for (int j = 0; j < Lusu.size(); j++) {
				for (int k = 0; k < center2.size(); k++) {
					elements.push_back(vector<int>{center1[i], Lusu[j],center2[k]});
				}
				
			}
		}
		for (int i = 0; i < elements.size(); i++) {
			string dir_name = path + flag + "_"+a[elements[i][0]] + a[elements[i][1]] + a[elements[i][2]] + ".config";
			double factor = 1;//用于坐标放缩的系数，按照阴离子半径进行放缩

			//factor = (e[elements[i][1]].nega_ionic_ridus + e[elements[i][0]].posi_ionic_ridus) / now;
			if (factor < 0.8)
				factor = 1;
			//第一个元素是中心元素，第二个是环境元素
			pattern->get_pattern(dir_name, elements[i], *pattern, factor);

		}

	}
	else if (flag == "1_4_A"||flag=="1_4_B"||flag=="2_1"||flag=="3_1"||flag=="SnBr2") {
		//对应结果是1:4的情况
		double now = 0;
		vector<int> Lusu = { 9,17,35,53 };
		vector<int> center;
		if (flag == "1_4_A" || flag == "1_4_B") {
			center = { 32,50,82,22,40,72 };
			now = e[53].nega_ionic_ridus + e[22].posi_ionic_ridus;
		}
		else if (flag == "2_1") {
			center = {21,39,71,31,49,81};
			now = e[9].nega_ionic_ridus + e[22].posi_ionic_ridus;
		}
		else if (flag == "3_1") {
			center = {33,51,83,23,41,73};
			now = e[9].nega_ionic_ridus + e[23].posi_ionic_ridus;
		}
		else if (flag == "SnBr2") {
			center = { 32,50,82,22,40,72 };
			now = e[35].nega_ionic_ridus + e[50].posi_ionic_ridus;
		}
		
		for (int i = 0; i < center.size(); i++) {
			for (int j = 0; j < Lusu.size(); j++) {
				elements.push_back(vector<int>{center[i], Lusu[j]});
			}
		}
		
		//cout << now << endl;
		for (int i = 0; i < elements.size(); i++) {			
			string dir_name = path + flag+a[elements[i][0]] + a[elements[i][1]] + ".config";
			double factor = 1;//用于坐标放缩的系数，按照阴离子半径进行放缩

			factor = (e[elements[i][1]].nega_ionic_ridus + e[elements[i][0]].posi_ionic_ridus) / now;
			if (factor < 0.8)
				factor = 1;
			//第一个元素是中心元素，第二个是环境元素
			pattern->get_pattern(dir_name, elements[i], *pattern, factor);

		}

	}
	else if (flag == "2_1") {
		//如果是第二类的
		vector<int> Lusu = { 9,17,35,53 };
		vector<int> center = { 32,50,82,22,40,72 };
		for (int i = 0; i < center.size(); i++) {
			for (int j = 0; j < Lusu.size(); j++) {
				elements.push_back(vector<int>{center[i], Lusu[j]});
			}
		}
		double now = e[53].nega_ionic_ridus + e[22].posi_ionic_ridus;
		//cout << now << endl;
		for (int i = 0; i < elements.size(); i++) {
			string dir_name = path + flag + a[elements[i][0]] + a[elements[i][1]] + ".config";
			double factor = 1;//用于坐标放缩的系数，按照阴离子半径进行放缩

			factor = (e[elements[i][1]].nega_ionic_ridus + e[elements[i][0]].posi_ionic_ridus) / now;
			if (factor < 0.8)
				factor = 1;
			//第一个元素是中心元素，第二个是环境元素
			pattern->get_pattern(dir_name, elements[i], *pattern, factor);

		}
	}
	else
		return;

}

bool Pattern_cell::get_pattern(string& name, vector<int>element, cell &cell_pattern, double factor, string flag) {
	//根据输入的结构进行结构化产生新结构
	ofstream fout;
	fout.open(name, ios::out);
	if (!fout.is_open()) {
		cout << "can not generate the file!" << endl;
		return false;
	}
	//开始正式写入结构文件
	//设置显示格式
	fout.setf(ios::fixed);
	fout << cell_pattern.num << endl;
	fout << " LATTICE" << endl;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			fout << "  " << setprecision(9) << cell_pattern.letice[i][j] * factor << "\t";
		}
		fout << endl;
	}
	fout << " POSITION" << endl;
	//从这里进行修改，直接替换元素的相关信息
	if (1) {
		//cout << "aa" << endl;
		for (int i = 0; i < cell_pattern.num; i++) {
			if (i < 4) {
				fout << "  " << element[0] << "    ";
			}
			else if (i < 8) {
				fout << "  " << element[1] << "    ";
			}
			else
				fout << "  " << element[2] << "    ";
			for (int j = 0; j < 3; j++) {
				fout << setprecision(9) << cell_pattern.p[i][j] << "    ";
			}
			fout << "1     1     1" << endl;
		}
	}
	fout.close();
	return true;

	for (int i = 0; i<cell_pattern.num; i++) {
		if (cell_pattern.my_classify[i] == 1) {
			fout << "  " << element[0] << "    ";
		}
		else
		{
			fout << "  " << element[1] << "    ";
		}			
		for (int j = 0; j < 3; j++) {
			fout << setprecision(9) << cell_pattern.p[i][j] << "    ";
		}
		fout << "1     1     1" << endl;
	}
	




	if (flag == "1T") {
		//这里我们规定，第一行是正的，第二行是负的
		for (int i = 0; i < 3; i++) {
			if (i == 0) {
				fout << "  " << element[0] << "    ";
				for (int j = 0; j < 3; j++) {
					fout << setprecision(9) << cell_pattern.p[0][j] << "    ";
				}
				fout << "1     1     1" << endl;
			}
			else
			{
				fout << "  " << element[1] << "    ";
				for (int j = 0; j < 3; j++) {
					fout << setprecision(9) << cell_pattern.p[i][j] << "    ";
				}
				fout << "1     1     1" << endl;
			}
		}
	}
	else if (flag == "type6") {
		for (int i = 0; i < this->num; i++) {
			int index = i / 2;
			switch (index) {
			case 0:
				fout << "  " << element[0] << "    ";
				break;
			case 1:
				fout << "  " << element[1] << "    ";
				break;
			case 2:
				fout << "  " << element[2] << "    ";
				break;

			}
			for (int j = 0; j < 3; j++) {
				fout << setprecision(9) << cell_pattern.p[i][j] << "    ";
			}
			fout << "1     1     1" << endl;
		}

	}
	else if (flag == "type10") {
		for (int i = 0; i < this->num; i++) {
			if (i < 6) {
				fout << "  " << element[0] << "    ";
			}
			else if (i < 8) {
				fout << "  " << element[1] << "    ";
			}
			else
				fout << "  " << element[2] << "    ";
			for (int j = 0; j < 3; j++) {
				fout << setprecision(9) << cell_pattern.p[i][j] << "    ";
			}
			fout << "1     1     1" << endl;
		}
	}

	fout.close();
	return true;
}

vector<int> calculate_n123(double a, double b, double c, int ECUT2) {
	//用于产生n123的数值
	vector<int> res(3, 0);
	res[0] = ceil(pow(2 * ECUT2 / 2, 0.5) / PI * 1.889725989*a);
	res[1] = ceil(pow(2 * ECUT2 / 2, 0.5) / PI * 1.889725989*b);
	res[2] = ceil(pow(2 * ECUT2 / 2, 0.5) / PI * 1.889725989*c);
	//然后控制n1和n2是4的倍数，n3是整10的倍数
	res[0] = (res[0] / 2 + 1) * 2;
	res[1] = (res[1] / 2 + 1) * 2;
	res[2] = (res[2] / 10 + 1) * 10;
	return res;
}


double get_vecror_mo(double*a, int n) {
	double sum = 0;
	for (int i = 0; i < n; i++) {
		sum += a[i] * a[i];
	}
	return pow(sum, 0.5);
}

void fullfill_inter_real_location(double**inter, int kpt, cell &cell_a, string cell_path,
	int full_fill_index, int band, int spin, string flag)
{
	//根据输入的相关信息，获得转化后的坐标x，y，z
	report rep;
	if (rep.read_report(const_cast<char*>((cell_path + "/REPORT").c_str()), 1) == 1)
	{
		//说明没这个文件
		cout << "no report!" << endl;
		cin.get();
		return;
	}
	//分别根据vb和cb填充上坐标
	double test_kpt[3] = { 0.37582,0.6096,0 };
	double ka = ((2 * PI*get_vecror_mo(cell_a.letice[1], 3) / BOHR) / (get_dot_product(cell_a.letice[0], cell_a.letice[1], 3) / BOHR / BOHR));
	double kb = ((2 * PI*get_vecror_mo(cell_a.letice[0], 3) / BOHR) / (get_dot_product(cell_a.letice[0], cell_a.letice[1], 3) / BOHR / BOHR));
	//投影后的x和y,注意先是VB然后是CB
	double X_real = 0, Y_real = 0;
	/*X_real = ka * 1.7321 / 2 * rep.kp_p[kpt][0];
	Y_real = kb * rep.kp_p[kpt][1] + ka * 0.5*rep.kp_p[kpt][0];*/
	X_real = ka * 1.7321 / 2 * test_kpt[0];
	Y_real = kb * test_kpt[1] + ka * 0.5*test_kpt[0];
	double energy = rep.eig_kns[kpt][band][spin] / HARTREE;
	cout << flag << "loca" << rep.kp_p[kpt][0] << "," << rep.kp_p[kpt][1] << "," << rep.kp_p[kpt][2] << endl;
	cout << X_real << "," << Y_real << endl;
	cin.get();
	if (flag == "vb") {

		inter[0][full_fill_index] = X_real;
		inter[1][full_fill_index] = Y_real;
		inter[2][full_fill_index] = energy;
	}
	else if (flag == "cb") {
		inter[3][full_fill_index] = X_real;
		inter[4][full_fill_index] = Y_real;
		inter[5][full_fill_index] = energy;
	}
	else {
		cout << "wrong!flag!" << endl;
		cin.get();
		return;
	}

	return;

}

vector<double> Cal_cell::get_all_magnetic() {
	//获得该结构下所有自旋的情况可能
	int types = 0;
	int metal_count = 0;
	for (int i = 0; i < this->num; i++) {
		if (this->my_classify[i] == 1) {
			metal_count++;
		}
	}
	types = pow(2, metal_count);
	cout << "for count:" << metal_count << ", may have :" << types << " situations!" << endl;
	vector<double>res(1);
	res[0] = types;
	return res;
}

int Cal_cell::get_all_magnetic_for_all(string name, string path) {
	//批量化的产生所有电子的所有情况
	int all_count = 0;
	vector<string>names = get_file_name(name);
	for (auto singel : names) {
		cout << "generate config:" << singel << endl;
		Cal_cell* cell = new Cal_cell(singel,const_cast<char*>((path+singel).c_str()),1);
		all_count+=cell->get_all_magnetic()[0];
	}
	cout << "all situations:"<<all_count<<",all total work done" << endl;
	cin.get();
	return 1;

}
void Cal_cell::get_fromation_energy_peratom() {
	//获得该结构的形成能的方法
	//认为体结构已经绑定好了，直接计算相关能量就行

	double res = -100;
	double towd_averge = this->total_energy / this->num;
	double bulk_everge = this->bulk->total_energy / this->bulk->num;
	this->formation_energy_peratom= towd_averge - bulk_everge;
	return;
}
void Cal_cell::get_formation_energy_persquare() {
	//获得结构的形成能
	//计算方法是基于平面的结合能
	int yu = this->bulk->num%this->num;
	double times = this->bulk->num / this->num;
	if (yu != 0) {
		cout << "cell:" << this->name << ": do not sataisify the chemistry pair!" << endl;
	}

	if (CALCULATION_MODE == 2) {
		//注意这里所谓的剥离能，应该是各自取自己对应的面积		
		//2d一定是ab平面
		double a = pow(letice[0][0] * letice[0][0] + letice[0][1] * letice[0][1] + letice[0][2] * letice[0][2], 0.5);
		double b = pow(letice[1][0] * letice[1][0] + letice[1][1] * letice[1][1] + letice[1][2] * letice[1][2], 0.5);
		this->suqre = a * b;		
		//如果正好是符合化学计量比,2d的能量除以relax后的面积
		/*double twod_everge = this->total_energy / (max(this->suqre,this->org_square));*/
		double twod_everge = this->total_energy / this->suqre;
		//体结构的用relax之前的面积
		//double bulk_everge = this->bulk->total_energy / times / min(this->org_square,this->suqre);
		//cout << "org_square:" << this->org_square << ",after square:" << this->suqre << endl;
		double bulk_everge = this->bulk->total_energy / times / this->suqre;
		this->formation_energy_square = (twod_everge - bulk_everge)/2;
		return;		
	}
	else if (CALCULATION_MODE == 1) {
		//这个时候需要等效求解圆柱体表面积
		//v=(a*b)·c
		double* cross_res = get_cross_product(this->bulk->letice[0], this->bulk->letice[1],3);
		double volume = get_dot_product(cross_res, this->bulk->letice[2], 3);
		cout << "volume:" << volume << endl;
		
		this->ridus = pow((volume / times) / PI / this->h_length, 0.5);
		cout << "ridus:" << this->ridus << endl;
		this->suqre = 2 * PI*this->ridus*this->h_length;
		double energy_delta = this->total_energy - this->bulk->total_energy / times;
		this->formation_energy_square = energy_delta / this->suqre;
		return;

	}
	else if (CALCULATION_MODE == 0) {
		//等效为球体
		double* cross_res = get_cross_product(this->bulk->letice[0], this->bulk->letice[1], 3);
		double volume = get_dot_product(cross_res, this->bulk->letice[2], 3);
		this->ridus = pow((3 * volume / times / 4 / PI), 1 / 3);
		this->suqre = 4 * PI*ridus*ridus;
		double energy_delta = this->total_energy - this->bulk->total_energy / times;
		this->formation_energy_square = energy_delta / this->suqre;
		return;
	}
	else {
		cout << "Now cauculation mode do not work!" << endl;
		cin.get();
		return;
	}
	
}

void Cal_cell::get_org_cell(string path) {
	this->cell_org = new cell(path.c_str(),2);
	double a = pow(cell_org->letice[0][0] * cell_org->letice[0][0] + cell_org->letice[0][1] * cell_org->letice[0][1] + cell_org->letice[0][2] * cell_org->letice[0][2], 0.5);
	double b = pow(cell_org->letice[1][0] * cell_org->letice[1][0] + cell_org->letice[1][1] * cell_org->letice[1][1] + cell_org->letice[1][2] * cell_org->letice[1][2], 0.5);
	this->org_square = a * b;//如果是2维情况下，求解出原来的面积
	this->h_length = get_length(cell_org->letice[0]);//一定是沿着a方向切下来的
	return;
}

//输入是一个report的name，输出是spinup和down的不同带隙情况
vector<double> report::get_half_situation(string report_name) {
	vector<double>res(3, 0.0);
	double spin_up = 0;
	double spin_down = 0;
	string comm = "grep \"diff\" " + report_name + "| tail -1";
	FILE* fin;
	char buf[100];
	char allbuf[1024];
	if ((fin = popen(comm.c_str(), "r")) != NULL)
	{
		while (fgets(buf, sizeof(buf)-1, fin) != NULL)
		{
			strcat(allbuf, buf);
		}
		pclose(fin);
		fin = NULL;
	}
	//获得了输出
	string result = string(allbuf).substr(3);
	int index = result.find_first_of("=");
	string next = result.substr(index+1);
	//cout << result << endl;
	//按照空格分隔字符串
	vector<string>string_list = cut_string(next);
	//cout << next << "," << string_list[0] << ","<<string_list[1]<<endl;
	spin_up = stod(string_list[0].c_str());
	spin_down = stod(string_list[1].c_str());
	//判断取整的精度
	bool session = true;
	int up_number = 0, down_unmber = 0;
	int up, down = 0;
	up_number = ceil(spin_up);
	down_unmber = floor(spin_up);
	if (abs(spin_up - up_number) > Half_band_check && abs(spin_up - down_unmber) > Half_band_check) {
		cout << "spin_up:" << spin_up << ": do not sataixfy the need!" << endl;
		session = false;
	}
	up = (abs(spin_up - up_number) < abs(spin_up - down_unmber)) ? up_number : down_unmber;
	up_number = ceil(spin_down);
	down_unmber = floor(spin_down);
	if (abs(spin_down - up_number) > Half_band_check && abs(spin_down - down_unmber) > Half_band_check) {
		cout << "spin_down:" << spin_down << ": do not sataixfy the need!" << endl;
		session = false;
	}
	down = (abs(spin_down - up_number) < abs(spin_down - down_unmber)) ? up_number : down_unmber;
	session = true;
	if (session == false) {
		//如果是有错误的话，打印两个diff
		res[0] = spin_up;
		res[1] = spin_down;
		res[2] = -1;
		return res;
	}
	//下面开始返回各自的带隙情况	
	//1 返回第一个的需要的两个vector
	double now_fermi = (this->eig[num_kpt * num_electron * this->spin / 2 - 1].energy+ this->eig[num_kpt * num_electron * this->spin / 2 ].energy)/2;
	if (report_name.find("86554") != string::npos) {
		up = 44;
		down = 42;
	}
	else {
		up = 16;
		down = 13;
	}
	vector<double>::iterator max_positon, min_position;
	vector<double>spin1_low = get_allk_energy(0, up - 1);
	vector<double>spin1_high = get_allk_energy(0, up);
	max_positon = max_element(spin1_low.begin(), spin1_low.end());
	min_position = min_element(spin1_high.begin(), spin1_high.end());
	res[0] = *min_position - *max_positon;
	cout << "spin_up gap:" << *min_position - now_fermi << endl;
	//2 返回第2个需要的vector
	vector<double>spin2_low = get_allk_energy(1, down - 1);
	vector<double>spin2_high = get_allk_energy(1, down);
	max_positon = max_element(spin2_low.begin(), spin2_low.end());
	min_position = min_element(spin2_high.begin(), spin2_high.end());
	res[1] = *min_position - *max_positon;
	res[2] = 1;
	cout << "spin_down gap:" << *min_position - now_fermi << endl;
	return res;

}

//获取指定spin和band的所有的本正能量值
vector<double> report::get_allk_energy(int spin, int band) {
	vector<double> res;
	for (int i = 0; i < this->num_kpt; i++) {
		res.push_back(this->eig_kns[i][band][spin]);
	}
	return res;
}

//构造函数读取相关文件
Bandstru::Bandstru(string path,int flag) {
	this->read_report(const_cast<char*>((path + "/REPORT").c_str()),2);
	this->k_num = this->num_kpt;
	this->band_num = this->num_band;
	return;
	ifstream fin;
	if (flag == 1) {
		fin.open(path + "/bandstructure_1.txt", ios::in);
	}
	else if (flag == 2) {
		fin.open(path + "/bandstructure_2.txt", ios::in);
	}
	else {
		cout << "unkonwn the flag band txt!" << endl;
		cin.get();
		return;
	}
	if (!fin.is_open()) {
		cout << "i can not find the file:" << path << endl;
		cin.get();
		return;
	}
	this->k_step.resize(k_num);
	this->band_eig.resize(k_num);
	for (int i = 0; i < k_num; i++) {
		band_eig[i].resize(band_num);
	}
	bool finish = false;
	char buff[200];


	double correct_value = 0;//修正值
	double read_fermi = 0;//读取的OUT.FERMI能量值
	bool button = false;//文件是否存在得标致
	//然后再读取费米能级
	ifstream ifer;
	ifer.open(path + "/OUT.FERMI", ios::in);
	if (!ifer.is_open()||!ifer.good()) {
		//说明没有这个文件，那么直接减去计算得到的fermi能级
		button = false;
	}
	else
	{
		button = true;
		ifer.getline(buff, sizeof(buff) - 1);
		cout << "get:" << buff << endl;
		int index = string(buff).find_first_of("=");
		string latter = string(buff).substr(index + 1);
		index = latter.find_first_of("e");
		string last = latter.substr(0, index - 1);
		read_fermi = stod(last.c_str());
		ifer.close();
	}
	
	//主流程错了，有fermi文件时已经读取了相关结果
	
	if (Fermi_orig == true&&button==true) {
		correct_value = 0;
		this->fermi = read_fermi;
		cout << "fermi:" << this->fermi << endl;
	}
	else if(Fermi_orig==false && button==true){
		correct_value = read_fermi - this->orig_fermi;
		cout << "my cal fermi:" << this->orig_fermi << endl;
		this->fermi = this->orig_fermi;
	}
	else if (Fermi_orig==true&&button==false) {
		cout << "error! No fermi file exists!" << endl;
		cin.get();
	}
	else {
		//那就只剩一种情况，就是用自己的，但是也找不到fermi文件
		correct_value = -this->orig_fermi;
		this->fermi = this->orig_fermi;
	}
	//循环把能量值读取进去
	double k_show = 0;
	double energy_show = 0;
	while (fin.good() && fin.peek() != EOF) {
		for (int i = 0; i < band_num; i++) {
			for (int j = 0; j < k_num; j++) {
				fin >> k_show;
				fin >> energy_show;
				fin.getline(buff, sizeof(buff) - 1);
				/*cout << k_show << ":" << energy_show << endl;
				if (j > 140) {
					cout << "a" << endl;
				}*/
				k_step[j] = k_show;
				if(button==true)
					band_eig[j][i] = energy_show + correct_value;
				else
					band_eig[j][i] = energy_show -this->fermi;
			}
			fin.getline(buff, sizeof(buff) - 1);
		}
		finish = true;
	}
	if (!finish) {
		cout << "read the data unfinished!" <<path<< endl;
		cin.get();
	}
	fin.close();
	
	return;
}
//绘制输出能带绘制的数据
int Bandstru::drawBandstru(int start, int expand) {
	
	int real_start = 0;
	int real_end = 0;
	
	if (start == -1) {
		//start = this->vbm_band-1;
		init_start_band_num();
		start = this->orig_start_band;
		cout << "sure start:" << start << endl;
	}
	//设定范围
	real_start = start - expand;
	real_end = start + expand;
	if (start - expand < 0) {
		cout << "start <0!cut the range automatically" << endl;
		real_start = 0;
		//cin.get();
		//return 0;
		
	}
	if (start + expand > this->band_num) {
		cout << "end >0!cut the range automatically" << endl;
		real_end = band_num;
		/*cin.get();
		return 0;*/
	}
	this->real_start = real_start;
	this->real_end = real_end;

	
	return 1;
}
void Bandstru::output_band(string out_name, vector<Bandstru*> bandlist) {
	//输出band文件
	ofstream fout;
	fout.open(out_name, ios::out);
	fout.setf(ios::fixed, ios::floatfield);  // 设定为 fixed 模式，以小数点表示浮点数
	fout.precision(6);
	//首先打印第一行
	fout << "k_point" << "\t";
	int label = 1;
	for (auto band : bandlist) {
		for (int i = 0; i < (band->real_end - band->real_start); i++)
			fout << label << "\t";
		label++;
	}
	fout << endl;


	//然后打印拼的内容
	for (int step = 0; step <bandlist[0]->k_num; step++) {
		//首先打印表头
		fout << bandlist[0]->k_step[step] << "\t";
		for (int i = 0; i < bandlist.size(); i++) {
			for (int j = bandlist[i]->real_start; j < bandlist[i]->real_end; j++) {
				fout << bandlist[i]->band_eig[step][j] << "\t";
			}
		}
		
		fout << endl;
	}
	fout.close();
	cout << "draw the band :" << out_name << " finished!" << endl;
	return;
}

//绘制修改之后我们需要的dos的数据
int Bandstru::drawDos(string path,string out_name) {
	if (Fermi_orig == false) {
		cout << "use read the fermi :" << this->orig_fermi << endl;
	}
	else {
		cout << "use automitail the fermi :" << this->orig_fermi << endl;
	}
	//1 读取相关的文件
	int lines = 0;
	ifstream fin;
	fin.open(path + "/DOS.totalspin", ios::in);
	if (!fin.is_open()) {
		cout << "i can not find the file!total spin" << endl;
		cin.get();
		return 0;
	}
	fin.close();
	//1.1 说明有这个文件，输出有多少行
	string com = "cat " + path + "/DOS.totalspin" + "| wc -l >line";
	system(com.c_str());
	fin.open("line", ios::in);
	fin >> lines;
	lines -= 1;
	fin.close();
	com = "rm line";
	system(com.c_str());
	//1.2 建立vector的数组,分别是total，up以及down的数值
	int start = 0;
	double temp;
	vector<vector<double>> save(lines, vector<double>(3));
	char buffer[300];
	fin.open(path + "/DOS.spinup", ios::in);
	if (!fin.is_open()) {
		cout << "i can not find the file!spin up" << endl;
		cin.get();
		return 0;
	}
	while (fin.good() && fin.peek() != EOF) {
		fin >> temp;		
		save[start][0] = temp - this->orig_fermi;
		fin >> temp;
		save[start++][1] = abs(temp);
		fin.getline(buffer,299);		
	}	
	fin.close();

	//1.3 读取down的文件
	start = 0;
	fin.open(path + "/DOS.spindown", ios::in);
	if (!fin.is_open()) {
		cout << "i can not find the file!spin up" << endl;
		cin.get();
		return 0;
	}
	while (fin.good() && fin.peek() != EOF) {
		fin >> temp;
		fin >> temp;
		save[start++][2] = abs(temp);
		fin.getline(buffer, 299);
	}
	fin.close();
	cout << "读取产生的dos结果完毕：" << path << endl;

	//2 然后输出相关的结果
	ofstream fout;
	fout.open(path + "/" + out_name, ios::out);
	fout.setf(ios::fixed);
	setprecision(6);
	fout << "energy   up   down" << endl;
	for (int i = 0; i < save.size(); i++) {
		for (int j = 0; j < 3; j++) {
			fout << save[i][j] << " ";
		}
		fout << endl;
	}
	fout.close();
	cout << "产生结果dos完毕！" << path<<endl;
	return 1;
}

//批量化的将一组bandstru或者dos产生出来
int Bandstru::generate_all_bandstu(string name, string path) {
	vector<string>names = get_file_name(name);
	for (auto singel : names) {
		cout << "generate band for config:" << singel << endl;
		vector<Bandstru*>band_room;
		Bandstru* bandstru = new Bandstru(path + singel,1);
		//bandstru->drawDos(path + singel, "DOS_save");
		if (bandstru->drawBandstru(-1, 5) == 0) {
			cout << "draw the band:" << singel << ",encounter error!" << endl;
			cin.get();
			continue;
		}
		band_room.push_back(bandstru);
		cout << "get bandstructure1 finished!" << endl;
		//////需不需要打印第二个
		Bandstru* bandstru2 = new Bandstru(path + singel, 2);
		if (bandstru2->drawBandstru(-1,5) == 0) {
			cout << "draw the band:" << singel << ",encounter error!" << endl;
			cin.get();
			continue;
		}
		band_room.push_back(bandstru2);
		cout << "get bandstructure2 finished!" << endl;
		//然后才开始进行输出
		//输出和生产数据分离
		Bandstru::output_band(path + singel + "/my_band", band_room);
	}
	cout << "all total work done!" << endl;
	cin.get();
	return 1;
}
//初始化开始画的相关数值
void Bandstru::init_start_band_num() {
	//cout << this->orig_start_band << endl;
	double target = 0;
	while (this->orig_start_band == 0) {
		for (int i = 0; i < band_num; i++) {
			for (int j = 0; j <k_num; j++) {
				if (band_eig[j][i] > target) {
					//cout << band_eig[j][i] << endl;
					this->orig_start_band = i;
					return;
				}
			}
		}
		target = target - 0.3;
		//cout << "now target is:" << target << endl;
	}
	if (this->orig_start_band != 0) {
		return;
	}
	else
	{
		cout << "no value is larger than value!" << endl;
		cin.get();
		return;
	}
	
}

//将vasp的in，kpt的格式转变为pwmat的相关格式
void report::vaspKpt2pwmat(string file_name, string out_path,int cutnum) {
	ifstream fin;
	vector<vector<double>> save;
	fin.open(file_name.c_str(),ios::in);
	if (!fin.is_open()) {
		cout << "i can not find the file!" << file_name << endl;
		cin.get();
		return;
	}
	char buffer[100];
	vector<double> temp;
	double data;
	bool flag = false;
	int line = 0;
	while (fin.good() && fin.peek() != EOF) {
		fin.getline(buffer, 99);
		cout << buffer << endl;
		if(line%2==0&&line!=0)
			fin.getline(buffer, 99);
		if (strstr(buffer, "Reciprocal") != NULL||flag==true) {
			//说明下面开始读取真正的数值了
			//每行读取三个数字
			line += 1;
			flag = true;
			for (int i = 0; i < 3; i++) {
				fin >> data;
				cout << data << ",";
				temp.push_back(data);
			}
			save.push_back(temp);
			temp.clear();
			cout << endl;
			
		}
	}
	fin.close();
	cout << save.size() << endl;

	//下面开始输出相关结果
	ofstream fout;
	fout.open("gen.kpt", ios::out);
	fout << "Band" << endl;	
	fout.setf(ios::fixed);
	setprecision(6);
	for (int i = 0; i +1< save.size(); i += 2) {
		fout << cutnum << endl;
		for (int j = 0; j < 3; j++) {
			fout << save[i][j] << "    ";
		}
		fout << endl;
		for (int j = 0; j < 3; j++) {
			fout << save[i+1][j] << "    ";
		}
		fout << endl;
	}
	fout.close();



	//最后执行脚本，产生IN,KPT的结果
	system("split_kp.x < gen.kpt");
	string com = "mv IN.KPT " + out_path;
	system(com.c_str());
	cout << "has finished " + file_name << endl;
	return;
}


void report::allvaspKpt2pwmat(string filename, string inpath, string outpath) {
	vector<string> name = get_file_name(filename);
	for (auto singel : name) {
		report::vaspKpt2pwmat(inpath + singel, outpath + singel);
	}
	cout << "output kpath finished!" << endl;
	return;
}

bool cell::check_contain_supoerO() {
	//判断结构是否含有超氧结构
	int count = 0;
	for (int i = 0; i < this->num; i++) {
		if (this->type[i] == 8) {
			count++;
		}
	}
	if (count < 2)
		return false;
	double distance = 0;
	//下面开始检查任意两个O原子距离的情况了
	for (int i = 0; i < num; i++) {
		if (this->type[i] == 8) {
			for (int j = 0; j < yanshen; j++) {
				for (int k = 0; k < num; k++) {
					if (this->type[k] == 8) {
						distance = dis(this->real_position[CENTER][i], this->real_position[j][k]);
						if (distance > 0.1 && distance < cell::SuperO_rule) {
							return true;
						}
					}
				}
			}
		}
	}
	return false;
}

void Cal_cell::output_inforamtion(ofstream &fout) {
	fout << " "<<this->num << "  " << this->total_energy <<  "  " << this->bulk->num << "  " << this->bulk->total_energy << "  " << this->formation_energy_square << endl;
	//fout << " "<<this->num << "  " << this->total_energy << "  " << this->check_contain_supoerO()<<" "<<this->get_formula_name();
	return;
}

//判断一些指定路径的能带信息
vector<double> report::get_diff(string report_name) {
	//获得diff的相关信息，就三个数值
	vector<double>res(3, -1);		
	string comm = "grep \"diff\" " + report_name + "| tail -1";
	FILE* fin;
	char buf[100];
	char allbuf[1024];
	if ((fin = popen(comm.c_str(), "r")) != NULL)
	{
		while (fgets(buf, sizeof(buf) - 1, fin) != NULL)
		{
			strcat(allbuf, buf);
		}
		pclose(fin);
		fin = NULL;
	}
	//获得了输出
	string result = string(allbuf).substr(3);
	int index = result.find_first_of("=");
	string next = result.substr(index + 1);
	//cout << result << endl;
	//按照空格分隔字符串
	vector<string>string_list = cut_string(next);
	//cout << next << "," << string_list[0] << ","<<string_list[1]<<endl;
	res[0] = stod(string_list[0].c_str());
	res[1] = stod(string_list[1].c_str());
	res[2]= stod(string_list[2].c_str());
	return res;
}
vector<vector<double>> report::get_band_value_onkpt(vector<int>&klist,string file_name) {
	//获取指定k点列表对应的所有的本征能量值
	//考虑到两个方向，所以有两列数据
	//考虑到自旋占据问题，我们扩大一下能带的条数
	int num = klist.size();
	vector<vector<double>> target_eigen(num, vector<double>(4));
	vector<double> diff_info = get_diff(file_name);
	int up_band = round(diff_info[0]);
	int down_band = round(diff_info[1]);
	for (int i = 0; i < this->spin; i++) {
		for (int j = 0; j < num; j++) {
			if (i == 0) {
				target_eigen[j][0] = this->eig_kns[klist[j]][up_band-1][i];
				target_eigen[j][1] = this->eig_kns[klist[j]][up_band][i];
			}
				
			else {
				target_eigen[j][2] = this->eig_kns[klist[j]][down_band-1][i];
				target_eigen[j][3] = this->eig_kns[klist[j]][down_band][i];
			}
				
		}
	}
	/*for (int i = 0; i < num; i++)
		cout << klist[i] << ":" << target_eigen[i][0] << "," << target_eigen[i][1] << "," << target_eigen[i][2] << "," << target_eigen[i][3] << endl;*/
	return target_eigen;

}
vector<int>report::get_all_kptlist(string flag) {
	//获取指定标记下对应k点的序号
	vector<int>result_klist;
	cout << flag << endl;
	for (int i = 0; i < this->num_kpt; i++) {
		int test = this->base_kpt_flag[i];
		if (flag == "a") {
			if (test == 0 || test == 1)
				result_klist.push_back(i);
		}
		else if (flag == "b") {
			if (test == 0 || test == 2)
				result_klist.push_back(i);
		}
		else if (flag == "c") {
			if (test == 0 || test == 3)
				result_klist.push_back(i);
		}
		else {
			cout << "unknown flags,please check!" << endl;
			cin.get();
		}
	}
	
	return result_klist;
}

bool report::if_band_cross(vector<double>&eigen) {
	//判断指定的一系列的本征能量值是否穿过了能带
	/*cout << eigen.size() << endl;
	for (int i = 0; i < eigen.size(); i++)
		cout << eigen[i] << ",";*/
	vector<double>::iterator max_index, min_index;
	max_index = max_element(eigen.begin(), eigen.end());
	min_index = min_element(eigen.begin(), eigen.end());
	return ((*max_index - this->orig_fermi)*(*min_index - this->orig_fermi)) < 0;
}

bool report::read_base_kpt(string file_name)
{
	//读取指定路径下面的分数坐标的k点
	this->base_kpt_flag = vector<int>(this->num_kpt,-1);
	this->base_kpp = vector<vector<double>>(this->num_kpt, vector<double>(3));
	ifstream fin;
	fin.open(file_name.c_str(), ios::in);
	if (!fin.is_open()) {
		cout << "i can not find the file:" << file_name << endl;
		return false;
	}
	char temp[100];
	//先读取前两行
	for(int i=0;i<2;i++)
		fin.getline(temp, 99);
	int start = 0;
	while (fin.good() && fin.peek() != EOF) {
		if (start == this->num_kpt)
			break;
		for(int i=0;i<3;i++)
			fin >> this->base_kpp[start][i];
		this->base_kpt_flag[start] = judge_base_kpt_singel(this->base_kpp[start]);
		fin.getline(temp, 99);
		start++;
	}
	fin.close();
	return true;


}

int report::judge_base_kpt_singel(vector<double>&input) {
	//输入k点三个坐标，输出这个k点所属方向是什么
	double sameZero = 1e-4;
	bool a = abs(input[0]) < sameZero ? true : false;
	bool b = abs(input[1]) < sameZero ? true : false;
	bool c = abs(input[2]) < sameZero ? true : false;
	//cout << input[0] << "," << input[1] << "," << input[2] << endl;
	//cout << a << "," << b << "," << c << endl;
	
	if (a&&b&&c)
		return 0;
	else if (!a && b && c)
		return 1;
	else if (a&&!b &&c)
		return 2;
	else if (a && b&&!c)
		return 3;
	else		
		return 4;
}
//客户函数，调用此函数简单获得带隙穿过情况
//1表示复合情况，2个没穿过，1个穿过
//11表示特殊情况，上下自旋存在不同的情况
//2 表示不是我们注意的情况
int report::check_direction_gap(vector<string>flag, string report_name, string cell_name, string out_kpt_name) {
	//0 建立储存的结构
	vector<vector<bool>> judge_res = vector<vector<bool>>(flag.size(), vector<bool>(2));

	//1读取输入文件OUT.KPT
	if (!read_base_kpt(out_kpt_name)) {
		cout << "read OUT.KPT failed!" << endl;
		//读取失败情况下，进行k点的转化
		//cin.get();
		cout << "proceed k point transform!" << endl;
		report_k2_fractional(cell_name);

	}
	//2 读取不同标记，返回k点的集合
	for (int i = 0; i < flag.size(); i++) {

		vector<bool> if_cross;
		string target = flag[i];
		vector<int>klist = get_all_kptlist(target);
		vector<vector<double>> eigen = get_band_value_onkpt(klist, report_name);
		//然后构造其中的某列对应数据
		vector<double>choose(eigen.size(), -1);

		for (int j = 0; j < eigen[0].size(); j++) {
			for (int index = 0; index < eigen.size(); index++) {
				choose[index]=eigen[index][j];
			}
			if_cross.push_back(if_band_cross(choose));
		}
		judge_res[i][0] = if_cross[0] || if_cross[1];
		judge_res[i][1] = if_cross[1] || if_cross[2];
	}

	//3 判断不同方向的穿过没穿过情况，返回判断值
	/*cout << this->orig_fermi << endl;
	for (int i = 0; i < flag.size(); i++) {
		cout << flag[i] << ":" << judge_res[i][0] << "," << judge_res[i][1] << endl;
	}*/


	int cross_num_up = 0;
	int cross_num_down = 0;
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < flag.size(); j++) {
			if (i == 0 && judge_res[j][i] == true)
				cross_num_up++;
			else if (i == 1 && judge_res[j][i] == true)
				cross_num_down++;
		}
	}
	if (cross_num_up != 1 || cross_num_down != 1)
		return 2;
	for (int i = 0; i < flag.size(); i++) {
		if (judge_res[i][0] != judge_res[i][1])
			return 11;
	}
	return 1;
}

void report::report_k2_fractional(string file_name) {
	//读取晶胞结构文件，将report的真实坐标还原为分数坐标
	cell cell_a(file_name.c_str(), 0);	
	double**d_temp = new double*[3];
	double**rd = new double*[3];
	for (int i = 0; i < 3; i++) {
		d_temp[i] = new double[3];
		rd[i] = new double[3];
	}
	double det = cell_a.letice[0][0] * cell_a.letice[1][1] * cell_a.letice[2][2] + cell_a.letice[0][1] * cell_a.letice[1][2] * cell_a.letice[2][0] + cell_a.letice[0][2] * cell_a.letice[1][0] * cell_a.letice[2][1] - cell_a.letice[0][0] * cell_a.letice[1][2] * cell_a.letice[2][1] - cell_a.letice[0][1] * cell_a.letice[1][0] * cell_a.letice[2][2] - cell_a.letice[0][2] * cell_a.letice[1][1] * cell_a.letice[2][0];
	rd[0][0] = (cell_a.letice[1][1] * cell_a.letice[2][2] - cell_a.letice[1][2] * cell_a.letice[2][1]) / det * 2 * PI*BOHR;
	rd[0][1] = (cell_a.letice[1][2] * cell_a.letice[2][0] - cell_a.letice[1][0] * cell_a.letice[2][2]) / det * 2 * PI*BOHR;
	rd[0][2] = (cell_a.letice[1][0] * cell_a.letice[2][1] - cell_a.letice[2][0] * cell_a.letice[1][1]) / det * 2 * PI*BOHR;

	rd[1][0] = (cell_a.letice[2][1] * cell_a.letice[0][2] - cell_a.letice[2][2] * cell_a.letice[0][1]) / det * 2 * PI*BOHR;
	rd[1][1] = (cell_a.letice[2][2] * cell_a.letice[0][0] - cell_a.letice[2][0] * cell_a.letice[0][2]) / det * 2 * PI*BOHR;
	rd[1][2] = (cell_a.letice[2][0] * cell_a.letice[0][1] - cell_a.letice[2][1] * cell_a.letice[0][0]) / det * 2 * PI*BOHR;

	rd[2][0] = (cell_a.letice[0][1] * cell_a.letice[1][2] - cell_a.letice[1][1] * cell_a.letice[0][2]) / det * 2 * PI*BOHR;
	rd[2][1] = (cell_a.letice[0][2] * cell_a.letice[1][0] - cell_a.letice[1][2] * cell_a.letice[0][0]) / det * 2 * PI*BOHR;
	rd[2][2] = (cell_a.letice[0][0] * cell_a.letice[1][1] - cell_a.letice[1][0] * cell_a.letice[0][1]) / det * 2 * PI*BOHR;
	inv(3, rd, d_temp);

	//下面开始转换k点的坐标了
	//this->kp_fractional = vector<vector<double>>(this->num_kpt,vector<double>(3));
	double location[3];
	for (int i = 0; i < this->num_kpt; i++) {
		for (int j = 0; j < 3; j++) {
			location[j] = kp_p[i][0] * d_temp[0][j] + kp_p[i][1] * d_temp[1][j] + kp_p[i][2] * d_temp[2][j];
		}	
		//cout << kp_p[i][0] << "," << kp_p[i][1] << "," << kp_p[i][2] << endl;
		//cout <<location[0] << "," << location[1] << "," << location[2] << endl;

		this->kp_fractional.push_back(vector<double>{ location[0],location[1],location[2] });
		this->base_kpt_flag[i] = judge_base_kpt_singel(this->kp_fractional[i]);
	}
	for (int i = 0; i < 3; i++) {
		delete[]d_temp[i];
		delete[]rd[i];
	}
	delete[]d_temp;
	delete[]rd;
	return;
}

void report::check_direction_gap_multi(vector<string>flag, string name_file, string path, string kpt_name ) {
	//批量化的产生结构各个方向的统计信息
	ofstream fout;
	fout.open("structure_check", ios::out);
	vector<string>all_name = get_file_name(name_file);
	for (auto singel : all_name) {
		cout << "check structure:" << singel << endl;
		fout << singel << " ";
		report rtemp;
		rtemp.read_report(const_cast<char*>((path+singel+"/REPORT").c_str()),2);
		int res = rtemp.check_direction_gap(flag, path + singel + "/REPORT", path + singel + "/atom.config");
		rtemp.inrich_source_information(singel);
		fout << res << " "<<rtemp.soucrce<<endl;
	}
	fout.close();
	return;
}

void report::inrich_source_information(string config_name){
	//补充数据来源的信息
	string icsd_path = "/share/home/wangz/icsd_2019/ICSD/Expe.inorganic/Expe.inorganic.atmos/";
	string mp_path = "/share/home/database/BDM/files/";
	string number_name = "";
	number_name = config_name.substr(0, config_name.find_first_of("."));
	int check_line = 5;//负责检查前5行的内容
	char buffer[500];
	char temp[100];
	ifstream fin;
	if (config_name.find("cif") == string::npos) {
		//说明这个是mp的结构
		fin.open(mp_path + number_name + "/" + number_name + ".cif", ios::in);		
	}
	else {
		//说明是icsd的数据
		fin.open(icsd_path +  number_name + ".cif", ios::in);
	}
	//正式读取相关文件进行检查
	for (int i = 0; i < check_line; i++) {
		if (!fin.good() || fin.peek() == EOF)
			break;
		fin.getline(temp, 100);
		strcat(buffer, temp);
	}
	if (strstr(buffer, "generated") != NULL) {
		this->soucrce = "Calculation";
	}
	else
		this->soucrce = "Experiment";
	fin.close();
	return;
}


void bothTestAfsituation(string calculationpath,string name) {
	//测试反铁磁的情况
	//输入路径是原始的scf的路径
	if (IF_PULS_MARNETIC&&IF_AFTEST) {
		chdir(calculationpath.c_str());
		string com;
		cout << "prepare for:" << name << " ,af scf test!";
		com = "cp -r " + name + " " + name + "_aftest";
		system(com.c_str());
		puls_magnetic(calculationpath+name+"_aftest","AF");
		chdir((calculationpath + name+"_aftest").c_str());
		system("qsub aaa.pbs");
		cout << "config:" << name << ",af test proceed!" << endl;
		return;
		
	}
}

void Lobister::runCohp() {
	//执行每个结构下面的cohp的计算
	//1:准备lobsiterin文件
	if (this->sucessful == false)
		return;
	ofstream fout;
	fout.open(this->path + "/lobsterin_org", ios::out);	
	if (this->COHP_GET == "Normal") {
		fout << "COHPstartEnergy  -10" << endl;
	}
	else if (this->COHP_GET == "Special") {
		fout << "COHPstartEnergy  " << this->startRange-2 << endl;
	}
	else {
		cout << "unknown the cohp cal falg!" << endl;
		cin.get();
	}
	fout << "COHPendEnergy	  6" << endl;
	fout << "usebasisset pbeVaspFit2015" << endl;
	//这里不使用我们手动统计过的东西
	//直接采用全部统计的方法
	//fout << "includeOrbitals s p d" << endl;
	////然后循环对于每个元素情况列举出来
	for (auto ele : this->element_info) {
		fout << "basisfunctions ";
		fout << ele.first << " " << getBaseFunctions(ele.first,true) << endl;
	}
	//最后添加所有可能的半径的统计信息
	//直接将该目录之前产生的输入信息拿过来
	fout.close();
	//这里进行修改，针对每一组数据我们单独进行分析
	//因此要产生多组的cohp_in来进行lobster的计算
	string temp_file;
	chdir(this->path.c_str());
	for (int i = 0; i < this->ridusElment.size(); i++) {
		temp_file = "lobsterin_" + to_string(i + 1);
		string singel_plus = getLineString(i + 1, this->path + "/" + "cohpin");
		system((string("cp lobsterin_org lobsterin_") + to_string(i+1)).c_str());
		fout.open("lobsterin_" + to_string(i + 1), ios::app);
		fout << singel_plus << endl;
		fout.close();
		//然后改成正式输入文件
		system(("mv " + temp_file + " " + "lobsterin").c_str());
		system("lobster-4.0.0");
		//运行完之后修改一些文件的名字，防止覆盖
		system("ls |grep \"lobster$\" > change_name");
		vector<string> all_name = get_file_name("change_name");
		for (auto temp_name : all_name) {
			system(("mv " + temp_name + " " + temp_name + "_"+to_string(i + 1)).c_str());
		}
		system("rm change_name");
		system(("mv lobsterin lobsterin_" + to_string(i + 1)).c_str());
		system(("mv lobsterout lobsterout_" + to_string(i + 1)).c_str());
	}
	/*string com = "cat " + this->path + "/" + "cohpin >>" + this->path + "/" + "lobsterin";
	system(com.c_str());*/
	
	
	return;	
}

//获得对应偃师文件的轨道信息
string Lobister::getBaseFunctions(string name,bool org_flag) {
	//建立数字和字符的记录
	vector<string>orbatal;//记录的轨道的字符串
	vector<int>orbat_num;//对应的数字
	string com;
	string res;
	//之前的做法
	bool changed = true;
	if (changed) {
		//在这里就是直接读取手写的相关信息
		//作为输入文件
		string all_base_path = "/share/home/wangz/0522_1d_cohp/sample/All_basefunctions";
		com.clear();
		string in="\"";
		com = "grep " +in+name +":\""+ " " + all_base_path + " >pot_test";
		//cout << com << endl;
		system(com.c_str());
		com = "awk -F\":\" \'{print $2}\' pot_test >out_test";
		system(com.c_str());
		ifstream fin;
		fin.open("out_test", ios::in);
		char buff[100];
		fin.getline(buff, 99);
		res = string(buff);
		fin.close();
		system("rm pot_test out_test");

	}
	else
	{
		if (org_flag == false) {
			string potcar_path = "/share/home/wangz/POTCAR/PBE/";
			string com = "grep \"VRHFIN = \" " + potcar_path + name + " >pot_test";
			system(com.c_str());
		}
		else {
			//这时候就读取自己的就行了
			chdir(this->path.c_str());
			string com = "grep \"VRHFIN =\"" + name + " " + "POTCAR" + " >pot_test";
			system(com.c_str());
		}
		//切割出来需要的字符串
		com = "awk -F\": \" '{print $2}' pot_test > out_test";
		system(com.c_str());
		ifstream fin;
		fin.open("out_test", ios::in);
		if (!fin.is_open()) {
			cout << "i can not find the file:out_test" << endl;
			cin.get();
		}
		char temp_orb;
		int temp_num;
		while (fin.good() && fin.peek() != EOF) {
			fin >> temp_orb;
			fin >> temp_num;
			orbatal.push_back("" + temp_orb);
			orbat_num.push_back(temp_num);
		}
		fin.close();

		//反构建字符串
		system("rm pot_test out_test");
		
		for (int i = 0; i < orbat_num.size() - 1; i++) {
			res += to_string(orbat_num[i]) + string(orbatal[i]) + " ";
		}
	}

	return res;
}

void Lobister::getElementList() {
	//获得元素列表的信息
	cell cell_a((this->path+"/atom.config").c_str(),2);
	
	for (int i = 0; i < cell_a.type_num; i++) {
		pair<string,int> temp=make_pair(a[cell_a.type_save[i]], cell_a.type_save[i]);
		this->element_info.push_back(temp);
	}	
	//然后读取元素匹配的信息
	ifstream fin;
	fin.open(this->path + "/cohp_band", ios::in);
	int count = 0;
	fin >> count;
	char buff[100];
	fin.getline(buff, 99);
	while (fin.good() && fin.peek() != EOF) {
		fin.getline(buff,99);
		int mao_index = string(buff).find(":");
		string tt = string(buff).substr(0, mao_index);
		int index = tt.find("-");
		vector<string> temp(2);
		temp[0] = tt.substr(0, index);
		temp[1] = tt.substr(index + 1, tt.size());
		this->ridusElment.push_back(temp);		
	}
	fin.close();
	getElectronNum();
	return;

}
void Lobister::getRidus(){
	//获得每个结构的键长信息统计
	//这个实现已经产生了
	ifstream fin;
	fin.open(this->path + "/cohp_band", ios::in);
	if (!fin.is_open()) {
		cout << "previous band generate failed!" << this->path << endl;
		this->sucessful = false;
		cin.get();
	}
	//现在产生的ridus信息可能有多个了
	//然后读取元素匹配的信息
	try {
		int count_size = 0;
		fin >> count_size;
		char buff[100];
		fin.getline(buff, 99);
		while (fin.good() && fin.peek() != EOF) {
			fin.getline(buff, 99);
			int mao_index = string(buff).find(":");
			string tt = string(buff).substr(mao_index + 1, string(buff).size());
			this->ridus.push_back(atof(tt.c_str()));
		}
		fin.close();
		this->time = count_size;
	}
	catch (...) {
		this->sucessful = false;
		this->comm = "get ridus failed!";
		cin.get();
	}
	return;
}
void Lobister::getCohp() {
	//获得每个结构的cohp的数值
	if (this->sucessful == false)
		return;
	if (Lobister::COHP_GET == "Special") {
		//对应我们自己读取相应cohp的积分结果，并且进行积分		
		chdir(this->path.c_str());		
		
		//2:获取自旋向上和自旋向下对应范围的cohp数值
		//注意，这个文件能级已经减去了费米能级
		//改进，对于多种连键情况怎么办,对于多种连键情况单独进行处理
		bool start_get = false, end_get = false;		
		vector<double>up_start, up_end, down_start, down_end;
		double move_start = this->startRange - this->ferimi;		
		char buff[2048];
		double same_falg = 4e-1;
		int temp_no = 0;
		try {
			//对于每一组成键的数值都是由对应输出文件给出的
			//单独进行处理
			for (int i = 0; i < this->ridusElment.size(); i++) {
				//1:第一步获取cohp结果对应的N数值
				//首先获取No的数目
				string target_file = "COHPCAR.lobster_"+to_string(i+1);
				system(("grep \"No.\" "+target_file+" |wc -l>nn").c_str());
				ifstream fin;
				fin.open("nn", ios::in);
				fin >> temp_no;
				this->interactions.push_back(temp_no);
				fin.close();
				system("rm nn");



				start_get = false;
				end_get = false;
				fin.open(target_file, ios::in);
				if (!fin.is_open()) {
					cout << "in path:" << this->path << ", no file:" << target_file << endl;
					//cin.get();
					this->sucessful = false;
					return;
				}
				for (int j = 0; j < this->interactions[i] + 3; j++)
					fin.getline(buff, 2047);
				//下面开始正式的读取
				while (fin.good() && fin.peek() != EOF) {
					if (start_get&&end_get)
						break;
					fin.getline(buff, 2047);
					vector<string>res = cut_string(string(buff));
					double temp = atof(res[0].c_str());
					if (temp < move_start&&abs(temp - move_start) < same_falg && !start_get) {
						//捕捉到了start的数值
						//然后我们把不同元素对匹配的cohp数值都拿回来
						cout << temp << "," << move_start << endl;
						//int number = getNoForElementPair(this->ridusElment[i]);
						int number = 1;
						up_start.push_back(atof(res[1 + number].c_str()));
						down_start.push_back(atof(res[2 * this->interactions[i] + 3 + number].c_str()));
						start_get = true;
					}
					else if (abs(temp) < 1e-3 && !end_get) {
						//说明此时捕捉到了费米能级对应的数值
						cout << temp << ":close 0" << endl;
						//int number = getNoForElementPair(this->ridusElment[i]);
						int number = 1;
						up_end.push_back(atof(res[1 + number].c_str()));
						down_end.push_back(atof(res[2 * this->interactions[i] + 3 + number].c_str()));
						end_get = true;						
					}

				}
				fin.close();
			}
			if (start_get&&end_get) {
				//获得了对应的数值，开始进行数据整理
				for (int i = 0; i < this->ridusElment.size(); i++) {
					this->up_differ.push_back(up_end[i] - up_start[i]);
					this->down_differ.push_back(down_end[i] - down_start[i]);
				}
			}
			else {
				this->sucessful = false;
				this->comm = "Intergeted failed!";
				return;
			}
			
		}
		catch (...) {
			this->sucessful = false;
			cout << "get cohp failed!" << endl;
			cin.get();
			this->comm = "get cohp failed!";
			return;
		}
	}
	else if (Lobister::COHP_GET == "Normal") {
		char buff[100];
		//下面针对每一种元素配对的方式进行计算相关信息
		for (int i = 0; i < this->ridusElment.size(); i++) {
			vector<double> data;//当前这组元素对的对应cohp数值
			//对每一种元素匹配进行单独分析
			ifstream fin;
			string target_file = "ICOHPLIST.lobster";
			string back_target_file = "ICOOPLIST.lobster";
			bool coop_flag = false;
			fin.open(this->path + "/" + target_file, ios::in);
			if (!fin.is_open()) {
				cout << this->path << " no COHP outfile!" << endl;
				fin.open(this->path + "/" + back_target_file, ios::in);
				if (!fin.is_open()) {
					cout << "no coop file!" << endl;
					this->sucessful = false;
					return;
				}
				else {
					this->comm = "only coop";
					coop_flag = true;
				}

			}
			while (fin.good() && fin.peek() != EOF) {
				fin.getline(buff, 99);
				if (strstr(buff, "spin") != NULL)
					continue;
				if (strstr(buff, this->ridusElment[i][0].c_str()) != NULL && strstr(buff, this->ridusElment[i][1].c_str()) != NULL) {
					//说明是找到了这个结构
					vector<string> cutRes = cut_string(string(buff));
					data.push_back(atof(cutRes[cutRes.size() - 1].c_str()));

				}
			}

			fin.close();
			//然后正式开始书写相关数据
			////还是自己写求平均值的东西吧……
			double temp_sum = 0;
			for (int i = 0; i < data.size() - 1; i++)
				temp_sum += data[i];

			this->cohp.push_back(temp_sum / (data.size()));
			//这里我们试着给出一个极值
			this->min_cohp.push_back(*min_element(data.begin(), data.end()));
			this->max_cohp.push_back(*max_element(data.begin(), data.end()));

		}
	}
	else {
		cout << "unknown the cohp flag!" << endl;
		cin.get();
	}
	//string com;
	//if(coop_flag==false)
	//	com = "awk \'{print $8}\' " + this->path + "/" + target_file + " > temp_lob_res";
	//else
	//	com = "awk \'{print $8}\' " + this->path + "/" + back_target_file + " > temp_lob_res";
	//system(com.c_str());
	//
	//fin.open("temp_lob_res", ios::in);
	//if (!fin.is_open()) {
	//	cout << "no file!temp_lob_res" << endl;
	//	cout << this->path << endl;
	//	this->sucessful = false;
	//	cin.get();
	//}
	//
	//double temp_res;
	//fin.getline(buff, 99);
	//while (fin.good() && fin.peek() != EOF) {
	//	fin >> temp_res;
	//	data.push_back(temp_res);
	//}
	//fin.close();
	
	//Scout << this->cohp << endl;
	return;
}

//静态函数，最外层调用
//批量计算相关结果
void Lobister::runAllProcesses(string file_name,string path,string out_name) {
	//由于存在多个半径值
	//同样的结构会存在多行
	ofstream fout;
	fout.open(out_name, ios::app);
	vector<string> all_name = get_file_name(file_name);
	for (auto name : all_name) {
		cout << "start the :" << name << endl;		
		Lobister lobb(name,path+"/"+name);
		lobb.runProcess();
		for (int i = 0; i < lobb.time; i++) {
			fout << name << " ";
			if (lobb.sucessful == true) {
				//修改输出的数据
				//fout << lobb.ridus[i] << " " << lobb.cohp[i] << " " << lobb.max_cohp[i] << " " << lobb.min_cohp[i]<<" "<<lobb.comm;
				fout << lobb.ridus[i] << " " << lobb.ridusElment[i][0]<<"-"<<lobb.ridusElment[i][1]<<" "<<lobb.up_differ[i] << " " << lobb.down_differ[i] << " " << (lobb.up_differ[i]+lobb.down_differ[i])/2 << " " << lobb.comm;
			}
			else {
				fout << lobb.ridus[i] << " " << lobb.comm;
			}
			fout << endl;
		}
		
	}
	fout.close();
	cout << "all total work done!please check:" << out_name << endl;
	return;
}
void Lobister::getElectronNum() {
	chdir(this->path.c_str());
	string file_name = "/share/home/wangz/0522_1d_cohp/sample/All_outerElectron";
	cell cell_a((this->path + "/atom.config").c_str(), 2);
	int total = 0;
	char buff[100];
	for (int i = 0; i < cell_a.num; i++) {
		ifstream fin;	
		fin.open(file_name, ios::in);
		if (!fin.is_open()) {
			cout << "no file:" << file_name << endl;
			cin.get();
		}
		int temp=0;
		while (fin.good() && fin.peek() != EOF) {
			fin.getline(buff, 99);
			if (strstr(buff, (string(atom_name[cell_a.type[i]])+":").c_str()) != NULL) {
				int index = string(buff).find(":");
				temp = atoi((string(buff).substr(index + 1, string(buff).size())).c_str());
				break;
			}
		}
		total += temp;
		fin.close();
	}
	this->electron_num = total;
	getEferAndStep();
}


void Lobister::getEferAndStep() {
	chdir(this->path.c_str());
	ifstream fin;
	system("grep \"Error\" output.relax |wc -l>temp_error");
	int num = 0;
	fin.open("temp_error", ios::in);
	fin >> num;
	if (num >= 1)
	{
		this->sucessful = false;
		this->comm = "calculation error!";
		return;
	}
	fin.close();
	//如果这个计算已经报错了，就直接跳出
	//检查当前计算是否报错
	


	//获取对应的费米能级和步长信息
	
	fin.open("DOSCAR", ios::in);
	char buff[100];
	//先读取前5行
	for (int i = 0; i < 6; i++)
		fin.getline(buff, 99);
	try {
		vector <string>res = cut_string(string(buff));
		this->step = (atof(res[0].c_str()) - atof(res[1].c_str())) / (atoi(res[2].c_str()) - 1);
		this->ferimi = atof(res[3].c_str());
	}
	catch (...) {
		this->sucessful = false;
		this->comm = "Calcualtion_failed!";
	}	
	fin.close();
	return;
}

void Lobister::getStartRange() {
	//获得能量的初始下限数值
	//只有指定要求，分能级计算情况下才计算该数值
	if (this->sucessful == false)
		return;
	if (this->COHP_GET == "Special") {
		chdir(this->path.c_str());
		ifstream fin;
		fin.open("DOSCAR", ios::in);
		char buff[300];
		double same_flag = 2e-1;//当判断电子数目几乎近似的条件
		for (int i = 0; i < 6; i++)
			fin.getline(buff, 299);
		//正式开始读取
		//注意这里我们需要反向积分
		stack<vector<double>> energy_save;
		vector<double>temp(3, 0);
		bool find_flag = false;
		double now_energy = 0;
		double temp_sum=0;
		double last_energy = -999;
		try {
			while (fin.good() && fin.peek() != EOF) {
				fin >> now_energy;
				if (now_energy > this->ferimi) {
					//说明找到了，注意这里寻找费米能级我们找到第一个大于他的吧
					find_flag = true;
				}
				temp[0] = now_energy;
				fin >> temp[1];
				fin >> temp[2];
				energy_save.push(temp);
				if (find_flag == true)
					break;
				fin.getline(buff, 299);
			}
			fin.close();
			//下面开始反向积分
			while (!energy_save.empty()) {
				vector<double> temp = energy_save.top();
				temp_sum += temp[1] * this->step + temp[2] * this->step;
				//cout << temp_sum << endl;
				if (abs(temp_sum - this->electron_num) < same_flag) {
					cout << temp_sum << "," << this->electron_num << ",range:"<<temp[0]<<endl;
					this->startRange = temp[0];
					break;
				}
				energy_save.pop();
			}
		}
		catch (...) {
			this->sucessful = false;
			this->comm = "back jifen failed!";
			cin.get();
			return;
		}
	}
	return;
}
int Lobister::getNoForElementPair(vector<string>&ele_info) {
	//对于指定的元素连接，得到是NO.第几个的对应数据
	chdir(this->path.c_str());
	int result = -1;
	string file_name = "COHPCAR.lobster";
	ifstream fin;
	fin.open(file_name, ios::in);
	char buff[500];
	for (int i = 0; i < 3; i++)
		fin.getline(buff, 499);
	//然后就是读取的NO对应
	while (fin.good() && fin.peek() != EOF) {
		fin.getline(buff, 499);
		if (strstr(buff, ele_info[0].c_str()) != NULL && strstr(buff, ele_info[1].c_str()) != NULL) {
			int index1,index2;
			index1 = string(buff).find(ele_info[0])+ele_info[0].size();
			index2 = string(buff).find_last_of(ele_info[1])+ele_info[1].size();
			if (buff[index1] - '0' >= 1 && buff[index1] - '0' <= 9 && buff[index2] - '0' >= 1 && buff[index2] - '0' <= 9) {
				//这样才是最终结果
				index1 = string(buff).find_first_of(".");
				index2 = string(buff).find_first_of(":");
				result = atoi(string(buff).substr(index1 + 1, index2 - index1 - 1).c_str());
				break;
			}
		}
	}
	fin.close();
	return result;
}
string Lobister::getLineString(int line, string file_name) {
	//获取指定文件的指定第几行
	ifstream fin;
	fin.open(file_name, ios::in);
	int line_tmep = 0;
	char buff[200];
	while (fin.good() && fin.peek() != EOF) {
		fin.getline(buff, 199);
		line_tmep++;
		if (line_tmep == line)
			break;
	}
	fin.close();
	return string(buff);
}