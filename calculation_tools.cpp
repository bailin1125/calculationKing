
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
	//��Ԫ���Ե�code
	report rep;
	string test_path = "/share2/wangz/0527_got1T/test/calculation/GaAs.config/REPORT";
	double near_location[3];
	int cbm_kpt_xuhao[] = { 0,1,5,36,180,6,30 };
	int vbm_kpt_xuhao[] = { 0,1,5,36,180,6,30 };
	//��ӡһ����Ҫ�Ķ���
	rep.read_report(const_cast<char*>(test_path.c_str()), 1);
	//CBM������
	//����
	cout << "start cbm!" << endl;
	cout << rep.cbm_kpt << endl;
	cout << rep.eig_kns[rep.cbm_kpt][rep.cbm_band - 1][0] / HARTREE << endl;
	for (int x : cbm_kpt_xuhao) {
		cout << "kpt: " << x << endl;
		cout << "energy:" << rep.eig_kns[x][rep.cbm_band - 1][0] / HARTREE << endl;
		cout << "location: " << rep.kp_p[x][0] << "," << rep.kp_p[x][1] << "," << rep.kp_p[x][2] << endl;
	}


	//����
	cout << rep.cbm_kppt[0] << "," << rep.cbm_kppt[1] << "," << rep.cbm_kppt[2] << endl;
	//cout << dis(near_location, rep.cbm_kppt);


	//VBM������
	//����
	cout << "start vbm!" << endl;
	cout << rep.vbm_kpt << endl;
	cout << rep.eig_kns[rep.vbm_kpt][rep.vbm_band - 1][0] / HARTREE << endl;
	for (int x : vbm_kpt_xuhao) {
		cout << "kpt: " << x << endl;
		cout << "energy:" << rep.eig_kns[x][rep.vbm_band - 1][0] / HARTREE << endl;
		cout << "location: " << rep.kp_p[x][0] << "," << rep.kp_p[x][1] << "," << rep.kp_p[x][2] << endl;
	}
	//����
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
		//���������������diff��ֵ
		//tt << setprecision(12) <<cell_a->check_energy(calculation_path + x, "REPORT_scf_test")<<" ";
		//cin.get();
		//tt << cell_a->check_diff(calculation_path + x, "REPORT_scf_test") << " ";
		//tt << setprecision(12) << cell_a->check_energy(calculation_path + x, "REPORT_spin2_test") << " ";
		//tt << setprecision(12) << cell_a->check_diff(calculation_path +x, "REPORT_scf_test")<<" ";
		//tt << setprecision(12) << cell_a->check_energy(calculation_path + x, "REPORT") - cell_a->check_energy(calculation_path + x, "REPORT_spin2_test");
		//���centr����ΧԪ��
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


	//�йز��Է������ݵı���
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
	cout << "��麬��ˮ��2D�ṹ������" << endl;
	return 0;
	//�Ȳ��������ļ�
	string yanshi_name;
	string new_yanshi_path;
	full_fill_yanshi(yanshi_name, new_yanshi_path);
	cout << "fullfille yanshi compleste!" << endl;
	cin.get();
}

//������ճ�������Ҫ�����
void test(element*e ) {
	
	//����DOS���
	Bandstru* bandstru = new Bandstru("/share2/wangz/0329_superO/0419_testAF/calculation/23163.cif.config_nou/pbe_dos", 2);
	bandstru->drawDos("/share2/wangz/0329_superO/0419_testAF/calculation/23163.cif.config_nou/pbe_dos", "DOS_save");
	Bandstru* bandstru1 = new Bandstru("/share2/wangz/0329_superO/0419_testAF/calculation/23163.cif.config_nou/hse_dos", 2);
	bandstru1->drawDos("/share2/wangz/0329_superO/0419_testAF/calculation/23163.cif.config_nou/hse_dos", "DOS_save");
	

	//������������ܴ��ṹ
	Bandstru::generate_all_bandstu("/share2/wangz/0329_superO/0410_dirac/calculation/soc/name", 
		"/share2/wangz/0329_superO/0410_dirac/calculation/soc/");

	//��������IN,KPT
	/*report::allvaspKpt2pwmat("/share2/wangz/0329_superO/0404_plusconfig/kpath/narrow_path/name","/share2/wangz/0329_superO/0404_plusconfig/kpath/narrow_path/","/share2/wangz/0329_superO/0404_plusconfig/kpath/narrow_path/inkpt/");
	cin.get();*/

	//����������ͬ��������ĵ���
	//Cal_cell::get_all_magnetic_for_all("/share2/wangz/1216_half_metal/0102/name", "/share2/wangz/1216_half_metal/0102/config/");

	//�����������нṹ�İ������Ϣ
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
	////��������ĳ���ṹ�ǲ��ǰ����
	/*report rt;
	rt.read_report("/share2/wangz/0329_superO/0331_conductor/calculation/hse_dos/23163.config/REPORT_scf",2);
	cout << rt.orig_fermi << endl;
	vector<double>res = rt.get_half_situation("/share2/wangz/1216_half_metal/0114/hse_band_v2/calculation/432515.config/REPORT_scf");
	cout << res[0] << "," << res[1] << "," << res[2];
	cin.get();*/



	////���ĳ���ṹ���б���
	//cell cell_a("/share2/wangz/1216_half_metal/0114/all_spin_situation/config/432515.config", 0);
	//cell_a.supper_cell(2, 2, 1,"/share2/wangz/1216_half_metal/0114/all_spin_situation/config/");
	//cin.get();

	//ofstream fout;
	//
	//vector<double> vacum = { -1,12,12};//���¶�����ղ�
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
	////ÿ���ṹ��Ӧһ����ά���飬�洢һ��K�����꼰�������
	////����������ص�k��Ͷ�Ӧ������
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
	////��ε������ǣ���ӡ��diff������ܣ��Լ��ṹԪ�ص��ж����
	//ϣ�����ⲿ�ֳ����������Ҫ�ظ��Ĵ���
	


	string twod_name = "/share/home/wangz/high/1122_1d_formation_square";
	string bulk_cal_path = "/share2/wangz/1007_formation_energy/1d_new_bulk/calculation/";//�����ṹ�ļ���·��
	string twod_cal_path_a = "/share2/wangz/1007_formation_energy/1d_new/calculation/";//���2d���ּ���·��
	string twod_cal_path_b = "/share2/wangz/1007_formation_energy/2d/calculation/";//���2d���ּ���·��
	string twod_config_path = "/share/home/wangz/2d_search/ridus_cut/result/result_ionicandval/0103/my_1d/config/";//ԭʼ���������Ľṹ����������һЩ�������


	vector<string> name = get_file_name(twod_name);
	ofstream fout;
	fout.open("formtion_1d", ios::out);
	double diff = 0;//diff��ֵ
	int element_judge = 0;//Ԫ���ж����
	for (auto singel : name) {
		cout << "start the :" << singel << endl;
		try {
			Cal_cell cal(singel, const_cast<char*>((twod_cal_path_a + singel+"-0_1d.config/atom.config").c_str()), 0);
			if (cal.Valid == false)
				continue;
			//���Ȼ��Ԫ���ж����
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
				cout << singel << "what the hell!��ô��û�ҵ�������" << endl;
				cout << "������û��yanshi�����������⣬�Ȳ�����" << endl;
				cal.formation_energy_peratom = -100;

			}

			//Ȼ��ʼ����ṹ
			Cal_cell* bulk = NULL;
			bulk = new Cal_cell(singel + "bulk", const_cast<char*>((bulk_cal_path + singel + ".cif.config" + "/atom.config").c_str()), 0);
			bulk->total_energy = cell::check_energy(bulk_cal_path + singel + ".cif.config", "REPORT");
			cal.bulk = bulk;
			cout << cal.bulk->total_energy << endl;
			//�󶨳�ʼ�ṹ
			cal.get_org_cell(twod_config_path + singel + "-0_1d.config");
			double res = 0;
			if (cal.formation_energy_peratom != -100) {
				cout << "��������~" << endl;
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
	//�й��γ��ܵ���ؽ��
	cout << "all total work done!" << endl;
	cin.get();

}














double report::get_accuracy(string file_name) {
	//��ȡ��·���µľ�ȷ��
	//��ֻд1D�����

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
				//˵����ʼ������
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
		//Ȼ��ʼ��д��ʼ�ž�
		fin << "MAGNETIC" << endl;
		if (flag == "F") {
			for (int i = 0; i < cell_a.num; i++) {
				//���¸ı����ӴžصĲ���
				//������ȫ��Ĭ��Ϊ3			
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
				//���¸ı����ӴžصĲ���
				//������ȫ��Ĭ��Ϊ3			
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

//���ṹ������ղ�
void cell::plus_vacum(vector<double>& vacum, string path,string file_name) {
	//�������ÿ����ղ㣬ÿ����ֵС��1e-2��Ϊ���ı�
	vector<double> now_vacum(3, 0);
	//��һ������Сֵ���ڶ��������ֵ,��¼����ԭ�ӵ����
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
	//���Ȼ�õ�ǰ����ղ��Ƕ���
	//��Ҫ�޸���ղ�Ķ��壬������ͶӰ���Ǹ���ľ���
	for (int i = 0; i < 3; i++) {
		if (i == 0)
			now_vacum[i] = (1 + a_save[0] - a_save[1])*length_axle[i];
		else if (i == 1)
			now_vacum[i] = (1 + b_save[0] - b_save[1])*length_axle[i];
		else if (i == 2)
			now_vacum[i] = (1 + c_save[0] - c_save[1])*length_axle[i];
	}
	//��ӡ��ǰ����ղ��С
	
	for (double x : now_vacum)
		cout <<"now_cacum:"<< x << ",";
	
	//���ȿ������ڵ�letice��ʲô
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
	//Ȼ����ڷ������任��ȥ
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
	//�������ڵ�letice��ʲô
	/*cout << "origin: letice:" << endl;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			cout << this->letice[i][j] << " ";
		}
		cout << endl;
	}*/


	//����������ļ�
	output_config(path,file_name);
	//cin.get();
	return;

}
double cell::get_positive_distance() {
	//����ṹ�е����ھ�����������֮��ľ���
	double res = 1e5;
	int center = (yanshen - 1) / 2;
	double temp_dis = 0;
	//���ȼ����û��������
	bool posi_flag = false;
	int index = -1;//���Ҫ��¼��������Ԫ����ʲô
	for (int i = 0; i < num; i++) {
		if (my_classify[i] == 1) {
			posi_flag = true;
			index = type[i];
			true;
		}			
	}
	if(!posi_flag)
	{
		//���û�������ӣ���ôָ��һ��Ԫ����������
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
	//��report�ļ��У���ȡ������Ϣ
	string comm = "grep 'E_tot(eV)' " + file_name + " |tail -1";

	string check_comm = "ls |grep '" + file_name + "' |wc -l";

	char buf[200];
	char result[1024];
	FILE* fp = NULL;
	int i = 0;

	chdir(path.c_str());
	//system("pwd");
	//������Ҫȷ��������ļ�
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
	//Ȼ���Ǵ��ַ����ж�������
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
	//������Ҫȷ��������ļ�
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
	//Ȼ���Ǵ��ַ����ж������һ������

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
	//��ν������Ԫ�ؾ���A�Ľ���+֮ǰ��Ϊ�ĸ���Ԫ��
	vector<int>main_element(Agroup_matal, Agroup_matal + 22);
	for (int i = 0; i < main_groupnum; i++) {
		main_element.push_back(main_group_element[i]);
	}
	sort(main_element.begin(), main_element.end());
	//�����ж��ǲ��Ǻ�������Ԫ��
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
	//Ȼ�����жϽ��
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
		cout << "�ṹԪ���ж��������飡" << endl;
		cin.get();
		return 4;
	}
}

bool cell::check_contain_transition_metal()
{
	//�жϸýṹ�ǲ��Ǻ��й��ɽ���
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
	//�����ԭ·��������������ص��ļ���Ŀ��·�����Լ�����Ľṹ�ļ�����
	vector<string> config_name = get_file_name(name);
	Wannier_cell::get_source_dependent();//����·������
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
	//�������·��������ؽű�����
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
	//����cb��
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
	//Ȼ������еĶ��ƶ���CB�ļ���
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
	//Ȼ������еĶ��ƶ���CB�ļ���
	system("mv UNK* *.xsf need wannier90.win wannier90.eig wannier90.amn wannier90.mmn wannier90.nnkp prefix.wf* VB/");
	system("rm check_cb check_vb");
	cout << "VB got wannier function successfuly!" << endl;
}

vector<int> Wannier_cell::get_left_and_right(int start, int end, int num, int band, string flag) {
	vector<int> res(4, 0);
	res[0] = start;
	res[3] = end;
	if (flag == "CBM") {
		//�����CBM�Ļ���Ҫ��cbm��ʼ����
		res[1] = band - 1;
		res[2] = band + num;
		if (res[2] >= end) {
			res[2] = -100;
		}
	}
	else if (flag == "VBM") {
		//VBM�Ļ���Ҫ������
		res[1] = band - num;
		if (res[1] <= start) {
			res[1] = -100;
		}
		res[2] = band + 1;
	}
	else
	{
		cout << "wrong flag!please check��" << endl;
	}
	return res;
}
void Wannier_cell::get_wain90(string path, string source_path) {
	//�������һ���ǽṹ����·����һ���������õ��ļ�·��
	chdir(path.c_str());

	//system("sleep 2");
	//���滻��������
	string com;
	ofstream fout;
	//�����ƿ�����ȥ
	com = "cp *.UPF prefix.save/";
	system(com.c_str());


	//���ȹ���cb��,��Ҫд�������ӵ�λ��
	com = "cp " + source_path + this->source_name + "_cb" + " wannier90.win_cb";
	system(com.c_str());
	//����cb��Ҫ����report��ȡ��õ������Ϣ
	fout.open("wannier90.win_cb", ios::app);
	report rp;
	rp.read_report("REPORT", 1);
	double efer = rp.get_singel_number("E_Fermi(eV)=");
	double cbm = rp.get_singel_number("CBM");

	//���ȴ�ӡexcludes_band
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

	//�����������
	if (this->source_name.find("SnBr2") != string::npos) {
		//������SnBr2��һ��ר��
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
	//��д������Ϣ
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
	//��д��������͸������
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



	//�����ǲ���vb��
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
	//��ӡ�ϲ��Ǻ���Ҫ��
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
		//������SnBr2��һ��ר��
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
	//��д������Ϣ
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
	//��д��������͸������
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


	//���������k������
	//int line_num= Wannier_cell::get_start_line();//���ص��Ǵ��ű�ǵ���һ���к�
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
	//���ȫ��������
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
	//�����������ļ���������ǽ���õ�ilmbda������ļ�����

	ofstream fout;
	fout.open(file_name, ios::out);
	if (fout.is_open()) {
		//��map�־û����ļ�
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
	//��ȡһ��Ŀ¼���������report������lambda
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
	//ͨ����ȡ��ص��ļ�,��������Ƕ�ȡ·��
	chdir(path.c_str());

	//�����ǻ�ȡcb�����
	//��Ҫע�⽫��õ�ilmbda��䵽������
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

	//���Ŷ�ȡvb�����
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

	//���ֵ���vector����
	if (this->source_name.find("Ti") != string::npos) {
		//������cb��
		ordered_lambda_cb.push_back(dict2tofind_double(lambda, "CB " + metal_name + "_s"));
		ordered_lambda_cb.push_back(dict2tofind_double(lambda, "CB " + metal_name + "_dz2"));
		ordered_lambda_cb.push_back(dict2tofind_double(lambda, "CB " + metal_name + "_dxz"));
		ordered_lambda_cb.push_back(dict2tofind_double(lambda, "CB " + metal_name + "_dxz"));
		ordered_lambda_cb.push_back(dict2tofind_double(lambda, "CB " + metal_name + "_dx2-y2"));
		ordered_lambda_cb.push_back(dict2tofind_double(lambda, "CB " + metal_name + "_dx2-y2"));

		//Ȼ����VB��
		vector<double >temp_vb;
		temp_vb.push_back(dict2tofind_double(lambda, "VB " + main_group_name + "_s"));
		temp_vb.push_back(dict2tofind_double(lambda, "VB " + main_group_name + "_p"));
		temp_vb.push_back(dict2tofind_double(lambda, "VB " + main_group_name + "_p"));
		temp_vb.push_back(dict2tofind_double(lambda, "VB " + main_group_name + "_p"));
		ordered_lambda_vb = temp_vb;
		ordered_lambda_vb.insert(ordered_lambda_vb.end(), temp_vb.begin(), temp_vb.end());

	}
	else if (this->source_name.find("Mg") != string::npos) {
		//������cb��
		vector<double >temp_cb;
		temp_cb.push_back(dict2tofind_double(lambda, "CB " + main_group_name + "_s"));
		temp_cb.push_back(dict2tofind_double(lambda, "CB " + main_group_name + "_p"));
		temp_cb.push_back(dict2tofind_double(lambda, "CB " + main_group_name + "_p"));
		temp_cb.push_back(dict2tofind_double(lambda, "CB " + main_group_name + "_p"));
		ordered_lambda_cb.push_back(dict2tofind_double(lambda, "CB " + metal_name + "_s"));
		ordered_lambda_cb.insert(ordered_lambda_cb.end(), temp_cb.begin(), temp_cb.end());
		ordered_lambda_cb.insert(ordered_lambda_cb.end(), temp_cb.begin(), temp_cb.end());

		//Ȼ����vb��
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
		//����Ǻ�Zn����һ��
		//������cb��
		vector<double >temp_cb;
		temp_cb.push_back(dict2tofind_double(lambda, "CB " + main_group_name + "_s"));
		temp_cb.push_back(dict2tofind_double(lambda, "CB " + main_group_name + "_p"));
		temp_cb.push_back(dict2tofind_double(lambda, "CB " + main_group_name + "_p"));
		temp_cb.push_back(dict2tofind_double(lambda, "CB " + main_group_name + "_p"));
		ordered_lambda_cb.push_back(dict2tofind_double(lambda, "CB " + metal_name + "_s"));
		ordered_lambda_cb.insert(ordered_lambda_cb.end(), temp_cb.begin(), temp_cb.end());
		ordered_lambda_cb.insert(ordered_lambda_cb.end(), temp_cb.begin(), temp_cb.end());

		//Ȼ����vb��
		double p = dict2tofind_double(lambda, "VB " + main_group_name + "_p");
		for (int i = 0; i < 6; i++) {
			ordered_lambda_vb.push_back(p);
		}

	}
	return;
}
double Wannier_cell::dict2tofind_double(map<string, double> mat, string key)
{
	//�����ص�double��ֵ����map��
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
	//ע���ȴ�ӡvb���ٴ�ӡcb
	report* rp = new report();
	rp->read_report("REPORT", 1);//�������ע��һ��
	fout << " " << 4 * rp->n_x << " " << 4 * rp->n_y << " " << rp->n_z << endl;
	fout << " " << this->ordered_lambda_vb.size() + this->ordered_lambda_cb.size() << " " << "1" << endl;
	fout << " " << "------------------" << endl;

	for (int i = 0; i < 2; i++) {
		//ע��Ҫд����һģһ����
		//���Ȳ���vb��
		for (int i = 0; i < this->ordered_lambda_vb.size(); i++) {
			fout << " " << i + 1 << " " << setprecision(6) << ordered_lambda_vb[i] << endl;
		}
		//Ȼ����cb��
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
	//�������ļ���ĺ���
	if (this->lmbda_valid == true) {
		chdir(now_path.c_str());
		string com;
		get_wannier_param();
		cout << "got wannier param completed!" << endl;
		//��Ҫ׼��in.wannier���ļ���ע���ʽ
		//�����Ȱ�vb���ù���
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
		//Ȼ����CB���ù���
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
		//���濪ʼ��׼�������ļ��Լ��ű�

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
		cout << "got etot for last scf completed��" << endl;
		string default_name = "gpu3";
		//cout << real_name << endl;
		com = real_name;
		int detect_flag = detect_free("gpu2");
		if (detect_flag == 1)
			pbs_got(com, "scf", now_path, "gpu2");//������ɣ���Ҫ����
		else if (detect_flag == 2)
			pbs_got(com, "scf", now_path, "gpu3");//������ɣ���Ҫ����
		else
			pbs_got(com, "scf", now_path, default_name);
		chdir(now_path.c_str());
		system("qsub aaa.pbs");
		return;
	}


}


void Wannier_cell::get_valid_m_new(string now_path, string name) {
	//�µķ�����ȡ��Ч����
	//��Ҫ�ǽ��ж���etot.input�����޸ĺ���������õ�
	//��Ҫ׼��IN.KPT���ļ�
	string com;

	chdir(now_path.c_str());
	//����׼��IN.KPT
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
	int inter_num = 1;//��ÿһ�������ϲ�����ٸ���
	kpt_num += 8 * inter_num * 2;
	ofstream fout;
	fout.open("IN.KPT", ios::out);
	fout << "            " << kpt_num << endl;
	fout << "            2           0" << endl;
	fout.close();
	system("cat temp>>IN.KPT");
	report rp;
	rp.read_report(const_cast<char*>((now_path + "/REPORT").c_str()), 1);
	//ע���������ֵ
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
	//ͬ���������VB��	
	for (int i = 0; i < new_num; i++) {
		for (int j = 0; j < new_num; j++) {
			fout << setprecision(8) << "     " << rp.vbm_kppt[0] + i * inter_x_length << "     " << rp.vbm_kppt[1] + j * inter_y_length << "      " << "0.00000000" << "      " << "0.00000000" << endl;
		}
	}
	fout.close();
	system("rm temp");
	fin.close();

	//Ȼ�������etot.input
	system("cp etot.input  etot_last");
	system("sed -i '7d' etot.input");
	system("echo '    IN.KPT=T' >>etot.input");

	//����ǽ�report����һ��
	system("cp REPORT REPORT_last");
	system("qsub aaa.pbs");
	return;


}
void Wannier_cell::get_valid_m(string now_path, string out_name) {
	//���wkm�µ���Ч����
	//���������ǱȽ�pbe�µ��ҵ���cbm��vbm��λ�ã����Ӧ����wkm��report��Ѱ��
	//��Ҫע�⼸��k��������12 12 1 ��
	//���������ܣ�ע������ĸ�ʽ����Ҫ�����ص�������ֵ����ֵ
	int colom = 20;
	if (this->real_name.find("BaI2") != string::npos || this->real_name.find("CdI2") != string::npos || real_name.find("HgI2") != string::npos || real_name.find("SrI2") != string::npos) {
		colom = 12;
	}
	chdir(now_path.c_str());
	report rp;
	rp.read_report(const_cast<char*>((Wannier_cell::pbe_path + real_name + "/REPORT").c_str()), 1);;
	//�൱�ڻ�ȡ��λ��
	report rp_now;
	rp_now.read_report("REPORT", 2);
	//��Ҫ��õڼ���k��͵ڼ���band��������ӽ���
	//����vb��
	int kk = find_close_kppt(rp.vbm_kppt, rp_now);
	int band = find_close_band(rp.vbm_energy, kk, rp_now);
	int save_band = band;
	string base_path = "/share/home/wangz/high/";
	ofstream fout, fout2;
	fout.open(base_path + out_name, ios::app);
	fout2.open(base_path + out_name + "_data", ios::app);
	fout << real_name << " " << endl;
	//fout << band << " ";
	//�ȼ���x�������Ч����
	double a1 = 0, a2 = 0, a3 = 0, b = 0;
	double a4 = 0, a5 = 0, a6 = 0, a7 = 0;//�ֱ��Ӧ�����������������µ�����
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
	//����������������ĸ�����ֵ




	valid_vb.push_back(a4);
	valid_vb.push_back(a5);
	valid_vb.push_back(a6);
	valid_vb.push_back(a7);
	//Ȼ����CB��
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
	//�ȼ���x�������Ч����
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
	//Ϊ1���ļ���׼�����е������ļ�
	string default_name = "gpu3";
	chdir(dir_path.c_str());
	system(("cp " + base_path + "/atom_super.config " + " atom.config").c_str());//���atom.config
	system(("cp " + base_path + "/*.UPF " + "./").c_str());//��������ļ�
	system(("cp " + base_path + "/REPORT " + "./REPORT_stepa").c_str());//��������ļ�
	report* rp = new report();
	rp->read_report("REPORT_stepa", 1);
	if (flag == "wkm_vb")
		sole_generate_jobs("atom.config", "./", "wkm_vb", rp);//���etot.input
	else if (flag == "wkm_cb")
		sole_generate_jobs("atom.config", "./", "wkm_cb", rp);//���etot.input
	int detect_flag = detect_free("gpu2");
	if (detect_flag == 1)
		pbs_got(real_name + x, "scf", dir_path, "gpu2");//������ɣ���Ҫ����
	else if (detect_flag == 2)
		pbs_got(real_name + x, "scf", dir_path, "gpu3");//������ɣ���Ҫ����
	else
		pbs_got(real_name + x, "scf", dir_path, default_name);
	//���aaa.pbs
	chdir(dir_path.c_str());
	system(("cp " + source_path + "IN.S_WKM_" + x + " " + "IN.S_WKM").c_str());//���in��swkm
	//����ҪҪ���IN,WANNIER�ļ�
	return;
}
void Wannier_cell::run_wkm(string nowpath, string source_path) {
	//ֱ������
	string com;
	vector<string> test = { "0.0","0.5","1.0" };
	chdir(nowpath.c_str());
	this->supper_cell(3, 3, 3, nowpath + "/");//���������ļ�,�������ļ�����atom_super.config
	system("mkdir CB_wkm VB_wkm");
	string temp_despath;

	//�ȴ���VB�Ķ���
	if (this->source_name.find("Si") != string::npos) {
		//�������ص��йز��Ե�
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
		//�������ص��йز��Ե�
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








	//Ȼ����CB�Ķ���	
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





	//�����ǹ����Ķ���,��s�����ͬ
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
	//���������k�������ж���ӽ���k������
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
	//�����������������ָ��k����Ѱ�Ҷ�Ӧ��band��Ŀ
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
	//�������
	string com = "source /share/home/wangz/pwmat.intel.cuda.env";
	system(com.c_str());
	com = "source /share/home/wangz/PWMAT/20200414/pwmat.env";
	system(com.c_str());
	com = "source ~/.bashrc";
	system(com.c_str());
	return;
}
int Wannier_cell::get_start_line(string flag) {
	//��ȡһ���ļ�
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
	//��ȡԪ���б�ȷ����ǰ����Դ��ʲô

	this->source_name = "SnBr2";
	return;

	bool falg_s = false, flag_f = false, flag_mg = false;
	bool flag_si = false, flag_as = false;
	bool flag_snbr = false;
	for (int i = 0; i < this->type_num; i++) {

		//�����ϲ��Ե������ṹ
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
	//Ȼ��ʼ����
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
		//������Ҫȷ��������ļ�
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
		//Ȼ���Ǵ��ַ����ж������һ������

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
	//������һ���ṹ��������жϸýṹ���޴���ˮ����
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
	//��ȡ��Ӧ�������ļ�,ȱ�ٵĴ�path��Ѱ�Ҳ������
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
			//˵��û�У���Ҫ��path�и���
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
	//����һ�����ذ汾
	int i = 0, j = 0;
	string com;
	ifstream fin;
	//��һ����ȡ����
	vector<string> all_name = get_file_name(config_name);


	//�ڶ�����ȡstatus�ļ�
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
	//��󣬸���log_all���в���
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

	//��Ҫע������ֻ��Ҫ��ȡ����е���Ϣ����
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
		//Ȼ����ݲ�ͬ���������flagֵ
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
			//����ע�⣬71��61�Ķ����ԣ���Ҫ���м��
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
				//ÿ�ζ�ȡ���mp��Ҫ����
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
	//��������ĸ߶ԳƵ����������õ��ܴ�ͼ
	chdir(path.c_str());
	//���ȵõ�cbm��vbm�������ߵ�5����k����
	//��Ϊƫ��Ķ���-45�ȣ��������6��������
	//K	0.5052386279	0.8750989735	0.0000000000
	//K'	-0.5052386279	-0.8750989735	-0.0000000000
	//	L	0.7578579418	0.4375494867	0.1306455197
	//	L'	-0.7578579418	-0.4375494867	-0.1306455197
	//	M	0.7578579418	0.4375494867	0.0000000000
	vector<vector<double>> vbm, cbm;
	vector<int> vbm_k, cbm_k;
	report rp;
	rp.read_report("REPORT", 1);
	//�����ص���ţ�����vbm,Ȼ�����vbm
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
	//Ȼ����Ǹ�����Ž�����չ��
	for (int xuhao : vbm_k)
		vbm.push_back(all_k[xuhao]);
	for (int xuhao : cbm_k)
		cbm.push_back(all_k[xuhao]);



	vector<vector<double>> vbm_plus, cbm_plus;
	double rule = 1.1;//��׼��1.0244���ﳤ��һ��
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

	//Ȼ��ʼ����IN.KPT�ļ�
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
	//��������Ĳ����������㣬����һ���µĹ��ߵĵ�
	double x = abs(a[0] - b[0]);
	double y = abs(a[1] - b[1]);
	vector<double> res(3);
	//����Ҫ������
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
	//�򵥵�����ļ�
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
	//��������ļ���
	ofstream fout;
	fout.open(path + file_name, ios::out);
	fout << "    " << this->num << endl;
	fout << " Lattice vector\n";

	int temp = 0;
	fout.setf(ios::fixed);
	//���濪ʼ���
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
	//�Ծ����������������������config�ļ�
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
	//��ȡ·������Ľṹ��Ϊmd��������
	if (CALCULATION_MODE == 1) {
		double realLength = pow(pow(this->letice[0][0], 2) + pow(this->letice[0][1], 2) + pow(this->letice[0][2], 2), 0.5);
		int times = ceil(MD_LENGTH / realLength);
		//��֮ǰ�޸�Դ�ļ�
		chdir(path.c_str());
		system("mv atom.config atom_org.config");
		supper_cell(times,1,1,path);
		//�޸ĺ�Ļ�atom.config
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
	//����Ӧ��aaa.pbs���Ƶ���Ӧ·��
	//���������falg�Ͷ�������ȷ��aaa.pbs
	//�ű�ͬĿ¼Ӧ����׼��relax.pbs,scf.pbs,nonscf.pbs
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

	//Ȼ�����mp123
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
	//�����õ�����Ľű�
	//�����ǽṹ���֣��Լ��ű�Ŀ¼����Ķ�Ӧ�ṹ�ļ�
	//�������ָ��·�����棬������Ӧ�ļ��У����Ҹ���flag������Ӧ�Ľű��ļ�
	//Ŀǰ��Ե���2d�ṹ
	//ѭ���������֣�����һЩ�в���

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
		//�����ļ����Լ�atom.config
		fin >> temp_name;
		cout << temp_name << endl;
		com = "mkdir " + dir_path + temp_name;
		system(com.c_str());
		com.clear();
		com = "cp " + data_path + temp_name + " " + dir_path + temp_name + "/atom.config";
		system(com.c_str());
		//��ʼ��ȡ�ṹ�ļ��������ű�
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


		//���������ļ���ָ��Ŀ¼��
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
		//��ʼдetot.input

		fclose(atom);

		//������ݲ�ͬ����Ҫ��д��ͬ������
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


	//���������ļ���ָ��Ŀ¼��
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
	//�����õ�����Ľű�
	//�����ǽṹ���֣��Լ��ű�Ŀ¼����Ķ�Ӧ�ṹ�ļ�
	//�������ָ��·�����棬������Ӧ�ļ��У����Ҹ���flag������Ӧ�Ľű��ļ�
	//Ŀǰ��Ե���2d�ṹ
	//ѭ���������֣�����һЩ�в���

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
		//�����ļ����Լ�atom.config
		fin >> temp_name;
		cout << temp_name << endl;
		/*com = "mkdir " + dir_path + temp_name;
		system(com.c_str());
		com.clear();*/
		com = "cp " + data_path + temp_name + " " + dir_path + temp_name + "/atom.config";
		system(com.c_str());
		//��ʼ��ȡ�ṹ�ļ��������ű�
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


		//���������ļ���ָ��Ŀ¼��
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
		//��ʼдetot.input
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

		//Ȼ�����mp123
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

		//������ݲ�ͬ����Ҫ��д��ͬ������
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
	//����һ���ļ����������ű������ڹ���
	//����etot.input
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


	//��ʼ��ȡ�ṹ�ļ��������ű�,��ʱ�ṹ�ļ�����Ŀ¼����
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


	//��ʼдetot.input
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


	//Ȼ�����mp123
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
	//�����relax�Ļ��������κ�Լ��������ֱ��1000
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

	//������ݲ�ͬ����Ҫ��д��ͬ������
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
		//�����0D�Ļ��Ͳ��ٿ���ʲô�����������
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
		//�����Ҫ�����ԭʼ��n123
		fout << "    N123 = " << rp->n_x << " " << rp->n_y << " " << rp->n_z << endl;
	}
	else if (flag == "scf")
	{
		if (Soc_on == true) {
			//����soc���������ϵ��������k��ȷ��������������
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
			fout << "    SPIN = 1\n";//�������wkm�Ļ���һ��ֻ����1
		}
		else if (IF_MD == true) {
			//�йؽ���md�ļ���
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
			//���������
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




		//����ͨ�õ�����
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
		//������wkm�ļ������
		fout << "    JOB=WKM" << endl;
		fout << "    CHARGE_DECOMP = T\n";
		fout << "    SPIN = 2\n";
		fout << "    ECUT = 50\n";
		fout << "    ECUT2 = 100\n";
		//�����������
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
		//����dos����
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
	//��������Ľṹ������Ӧ�Ľṹ�ļ�
	//���ݹؼ��ʲ�����Ӧ�����ṹ
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
	//������һ��·��
	//����Ǽ����Ƿ������0��ʾû��reportû������1��ʾ�������� ��2 ��ʾ����������
	//�����ж���û��report�ļ�
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
		//���������REPORT�ļ�����Ҫע���ʱ������������
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

		//���ﻹҪ������relax�����
		if (flag == 1)
		{
			ifstream fin;
			fin.open(path + "/RELAXSTEPS", ios::in);
			if (!fin.is_open())
			{
				cout << "can not find the file :RELAXSTEPS��" << endl;
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

	//�ȳ�ʼ����Ϣ
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
	//��ȡ���ۼ����»����ͽ�����

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
	//Ȼ���ȡ���Ӽ������ֵ����Сֵ
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
	//����ж��ٿ���ڵ㣬���ؿ���ڵ�ĸ���
	//�����Ǽ���ĸ�����
	//����0��ʾû�п��У�1��ʾgpu2���У�2��ʾgpu3����
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
	//ɾ����صĴ���ڵ�
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
		i++;//��ʱ�����ǵ�i-1���ڵ�
		if (strstr(line, "free") != NULL)
		{
			for (j = 0; j < nouse.size(); j++)
			{
				if (nouse[j] == (i - 1))
				{
					cout << i-1 << endl;
					num++;
					//��������һ�����ٰ汾
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
			else if (line[0] == 's')//����s��ʱ��˵���ǵ���
			{
				fscanf(in, "%s", line);
				fscanf(in, "%lf", &ridus[i][0]);
			}
			else//�����negative�Ļ����Ͷ�����ridusֵ
			{
				positive[i] = -1;
				fscanf(in, "%lf", &ridus[i][0]);//��һ���Ǹ���̬��Ӧ�İ뾶
				fscanf(in, "%lf", &ridus[i][1]);//��Ӧ����̬�İ뾶
			}
		}
		fgets(temp, 300, in);
	}
	//int x_xishu = 0;
	//int y_xishu = 0;
	//int z_zishu = 0;

	//���Ӽ�¼type�ļ�¼
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
	//��ʼ����֮ǰ�Ľ�����б궨
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
	//��Ҫ�Լ���д���εĲ���
	/*for (int i = 0; i < num; i++) {
		if (my_classify[i] == 1) {
			this->type[i] = 29;
		}
		else if (my_classify[i] == 2) {
			this->type[i] = 35;
		}
	}*/

	//��ȡԪ���滻���ԣ���ȡ����vector���ֱ�ָ��˳��ĵ�i��ԭ��Ӧ��ΪʲôԪ��
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
	//����ԭ������С��������
	sort(vec.begin(), vec.end(), map_com);
	//cout << vec[0].first << "," << vec[0].second;
	//cout << vec[1].first << "," << vec[1].second;
	//�������Ԫ�ص��м�ͻ���Ԫ��
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
		//��Ҫ����Ƿ��ǽ���Ԫ�ؽ��б궨
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
	//���һ�µ�ǰ�Ļ�ȡԪ��
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
	//����������Ǹ�ָ���İ���·�����ļ�����
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
	//������ʽ��ʼ�����̬��ƽ��ֵ

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
	//����һ����ά�����֣����������ϵ����֣�ֻ�������
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
		//�����ļ����������ϵĿռ�Ⱥ��Ϣ
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
	//�޸ķ����ܼ�����ֵ
	//this->orig_fermi = (eig[num_kpt * num_electron * this->spin / 2 - 1].energy + eig[num_kpt * num_electron * this->spin / 2].energy) / 2;//��ȡ���ݶ��ߵõ��ķ����ܼ�
	this->orig_fermi = eig[num_kpt * num_electron * this->spin / 2 - 1].energy;
	cbm_band = eig[num_kpt * num_electron * this->spin / 2].band + 1;
	vbm_band = eig[num_kpt * num_electron * this->spin / 2 - 1].band + 1;

	vbm_kppt = kp_p[eig[num_kpt * num_electron * this->spin / 2 - 1].kpt];
	cbm_kppt = kp_p[eig[num_kpt * num_electron * this->spin / 2].kpt];
	vbm_energy = eig[num_kpt * num_electron * this->spin / 2 - 1].energy;
	cbm_energy = eig[num_kpt * num_electron *this->spin / 2].energy;

	cbm_kpt = eig[num_kpt*num_electron*this->spin / 2].kpt;//Ӧ����12
	vbm_kpt = eig[num_kpt*num_electron*this->spin / 2 - 1].kpt;//Ӧ����37
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
		//ע��Ҫ��¼��Ӧ�Ĳ�����
		step[0] = to_string(kp_p[eig[num_kpt*num_electron*spin / 2].kpt][0]) + ","+to_string(kp_p[eig[num_kpt*num_electron*spin / 2].kpt][1]) + ","+to_string(kp_p[eig[num_kpt*num_electron*spin / 2].kpt][2]);
		step[1] = to_string(kp_p[eig[num_kpt*num_electron*spin / 2-1].kpt][0]) + "," + to_string(kp_p[eig[num_kpt*num_electron*spin / 2-1].kpt][1]) + "," + to_string(kp_p[eig[num_kpt*num_electron*spin / 2-1].kpt][2]);

		printf("CBM step original is:%lf,%lf,%lf\n", kp_p[eig[num_kpt*num_electron*spin / 2].kpt][0], kp_p[eig[num_kpt*num_electron*spin / 2].kpt][1], kp_p[eig[num_kpt*num_electron*spin / 2].kpt][2]);
		printf("VBM step original is:%lf,%lf,%lf\n", kp_p[eig[num_kpt*num_electron*spin / 2 - 1].kpt][0], kp_p[eig[num_kpt*num_electron*spin / 2 - 1].kpt][1], kp_p[eig[num_kpt*num_electron*spin / 2 - 1].kpt][2]);
	
	
	}
	return gap;

}
double report::get_singel_number(string flag) {
	//������һ����һ��
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
		//��������
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
		//�ҵ���1/52��λ�����ֵ��Ȼ�����1/104��������
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



		//Ȼ����cbm
		kp_temp[0] = rep.kp_p[rep.eig[rep.num_kpt * rep.num_electron * rep.spin / 2 - 1].kpt][0] * d_temp[0][0] + rep.kp_p[rep.eig[rep.num_kpt * rep.num_electron * rep.spin / 2 - 1].kpt][1] * d_temp[1][0] + rep.kp_p[rep.eig[rep.num_kpt * rep.num_electron * rep.spin / 2 - 1].kpt][2] * d_temp[2][0];
		kp_temp[1] = rep.kp_p[rep.eig[rep.num_kpt * rep.num_electron * rep.spin / 2 - 1].kpt][0] * d_temp[0][1] + rep.kp_p[rep.eig[rep.num_kpt * rep.num_electron * rep.spin / 2 - 1].kpt][1] * d_temp[1][1] + rep.kp_p[rep.eig[rep.num_kpt * rep.num_electron * rep.spin / 2 - 1].kpt][2] * d_temp[2][1];
		kp_temp[2] = rep.kp_p[rep.eig[rep.num_kpt * rep.num_electron * rep.spin / 2 - 1].kpt][0] * d_temp[0][2] + rep.kp_p[rep.eig[rep.num_kpt * rep.num_electron * rep.spin / 2 - 1].kpt][1] * d_temp[1][2] + rep.kp_p[rep.eig[rep.num_kpt * rep.num_electron * rep.spin / 2 - 1].kpt][2] * d_temp[2][2];

		//��������
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
	//�����Ǽ�鵱ǰ��vbm��cbm�Ƿ����м��λ���ϣ��Ƿ���1�����Ƿ���0
	//��Ҫ���ֲ���k�㽻�������
	//�������������һ����֮ǰ�Ķ�ȡ��һ�������Ķ�ȡ
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
	//ע�����������˹��Ѿ������˾��ȣ�������Ų��Ե�����������ȥ��
	int vbm = 5, cbm = 5;
	vbm = rep.eig[rep.num_kpt*rep.num_electron*rep.spin / 2].kpt;//Ӧ����12
	cbm = rep.eig[rep.num_kpt*rep.num_electron*rep.spin / 2 - 1].kpt;//Ӧ����37
	//������������Ҫ��ʵ�ǲ��Ƕ��ַ���ȷ��
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
	//���������ǳ���ʹ��map���ֵ��
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
	//���ָ��k���Ӧ�Ĵ�϶��С
	double energy_1 = this->eig_kns[kpt][this->num_electron/2][0];
	double energy_2 = this->eig_kns[kpt][this->num_electron/2-1][0];
	return energy_1 - energy_2;
}
void report::get_valid_m(int calculation_mode) {
	//�����Ч����,ע������1D��2D����3D�����
	int band_vbm = eig[num_kpt*num_electron*spin / 2].band;
	int band_cbm = eig[num_kpt*num_electron*spin / 2 - 1].band;
	int spin_vbm = eig[num_kpt*num_electron*spin / 2].spin;
	int spin_cbm = eig[num_kpt*num_electron*spin / 2 - 1].spin;
	int vbm = 5, cbm = 5;
	vbm = eig[num_kpt*num_electron*spin / 2].kpt;//Ӧ����12
	cbm = eig[num_kpt*num_electron*spin / 2 - 1].kpt;//Ӧ����37
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
		//����Ҫ�������������

		top = vbm - 5;
		botom = vbm + 5;
		left = vbm - 1;
		right = vbm + 1;
		left_top = top - 1;
		right_top = top + 1;
		left_end = botom - 1;
		right_end = botom + 1;
		int vbm_location[9] = { left_end,botom,right_end,left,vbm,right,left_top,top,right_top };

		//����ÿһ��inter����أ����ں������
		for (int i = 0; i < 9; i++) {
			inter[0][i] = kp_p[vbm_location[i]][0];
			inter[1][i] = kp_p[vbm_location[i]][1];
		}

		//�����ǿ��ٻ��x��y�����
		inter[2][0] = eig_kns[left_end][band_vbm][spin_vbm] / HARTREE;
		inter[2][1] = eig_kns[botom][band_vbm][spin_vbm] / HARTREE;
		inter[2][2] = eig_kns[right_end][band_vbm][spin_vbm] / HARTREE;

		inter[2][3] = eig_kns[left][band_vbm][spin_vbm] / HARTREE;
		inter[2][4] = eig_kns[vbm][band_vbm][spin_vbm] / HARTREE;
		inter[2][5] = eig_kns[right][band_vbm][spin_vbm] / HARTREE;

		inter[2][6] = eig_kns[left_top][band_vbm][spin_vbm] / HARTREE;
		inter[2][7] = eig_kns[top][band_vbm][spin_vbm] / HARTREE;
		inter[2][8] = eig_kns[right_top][band_vbm][spin_vbm] / HARTREE;

		//ע���������Ǹ����µ��޸ĺ�ķ�����������ʵ�������������
		/*a1 0.2 0.2 0->ʵ������
			a2 0.2 0.3 0
			a3 0.2 0.1 0

			a1 0.2 0.1 0
			a2 0.2 0.2 0
			a3 0.2 0.3 0
			bΪת��Ϊʵ������ľ���*/
			//ͬʱ����任֮�������λ��

			//Ȼ����x�ľ��룬y�ľ���
		averge_x = (dis(kp_p[vbm], kp_p[left]) + dis(kp_p[vbm], kp_p[right])) / 2;
		averge_y = (dis(kp_p[vbm], kp_p[top]) + dis(kp_p[vbm], kp_p[botom])) / 2;
		//ע��۴����ǿ�Ѩ���������ǵ���
		//�����ǵ��ӵ���Ч��������Ϊ����������y��x
		//������Ϊ��ͬһ��band

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

		//Ȼ��ͬ�������cbm��
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


		//Ȼ�����inter		
		inter[5][0] = eig_kns[left_end][band_cbm][spin_cbm] / HARTREE;
		inter[5][1] = eig_kns[botom][band_cbm][spin_cbm] / HARTREE;
		inter[5][2] = eig_kns[right_end][band_cbm][spin_cbm] / HARTREE;

		inter[5][3] = eig_kns[left][band_cbm][spin_cbm] / HARTREE;
		inter[5][4] = eig_kns[cbm][band_cbm][spin_cbm] / HARTREE;
		inter[5][5] = eig_kns[right][band_cbm][spin_cbm] / HARTREE;

		inter[5][6] = eig_kns[left_top][band_cbm][spin_cbm] / HARTREE;
		inter[5][7] = eig_kns[top][band_cbm][spin_cbm] / HARTREE;
		inter[5][8] = eig_kns[right_top][band_cbm][spin_cbm] / HARTREE;


		//Ȼ����x�ľ��룬y�ľ���
		//ע�������x��y�ľ���ע���˻���ʵ�ռ�ȥ��



		//ͶӰ���x��y

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
		//ע����1d�����ʱ��û�в�ֵ��������������

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
		//�����ά����µĻ�ȡ��Ч����
		//������vbm��
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

		//ͬ��ʼ��cbm��
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


//���صİ汾������·��������ֻ��ӡ������ؽ��,���ظ������ˣ���Ҫע���ˡ�������
void get_last_result(string path, vector<string>&step, vector<vector<double>>&table,int flag)
{
	//�м�������ǲ����ļ�¼

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
			//˵��û����ļ�
			cout << "no report!" << endl;
			return;
		}

	}
	else if (flag == 2)
	{
		//��2�ĵĻ�����Ϊ30��ֻ��ȡ��϶
		if (rep.read_report(const_cast<char*>((path + "/REPORT").c_str()), spin_check) == 1)
		{
			//˵��û����ļ�
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
	vbm = rep.eig[rep.num_kpt*rep.num_electron*rep.spin / 2].kpt;//Ӧ����12
	cbm = rep.eig[rep.num_kpt*rep.num_electron*rep.spin / 2 - 1].kpt;//Ӧ����37
	//cout << vbm << "," << cbm << endl;
	//���濪ʼ�����ؽ��
	if (CALCULATION_MODE == 1) {
		int top = 0, botom = 0, left = 0, right = 0, left_top = 0, right_top = 0, left_end = 0, right_end = 0;
		
		//Ȼ��ֱ����vbm��cbm�������
		int band_vbm = rep.eig[rep.num_kpt*rep.num_electron*rep.spin / 2].band;
		int band_cbm = rep.eig[rep.num_kpt*rep.num_electron*rep.spin / 2 - 1].band;
		int spin_vbm = rep.eig[rep.num_kpt*rep.num_electron*rep.spin / 2].spin;
		int spin_cbm = rep.eig[rep.num_kpt*rep.num_electron*rep.spin / 2 - 1].spin;
		//��ʼ���м�¼table
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
	//����ӡ����
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
			//˵��û����ļ�
			cout << "no report!" << endl;
			cin.get();
			gap = ex = ey = vy = vx = -100;
			return;
		}

	}
	else if (flag == 2)
	{
		//��2�ĵĻ�����Ϊ30��ֻ��ȡ��϶
		if (rep.read_report(const_cast<char*>((path + "/REPORT").c_str()), spin_check) == 1)
		{
			//˵��û����ļ�
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
	vbm = rep.eig[rep.num_kpt*rep.num_electron*rep.spin / 2].kpt;//Ӧ����12
	cbm = rep.eig[rep.num_kpt*rep.num_electron*rep.spin / 2 - 1].kpt;//Ӧ����37
	cout << vbm << "," << cbm << endl;


	//������������Ҫ��ʵ�ǲ��Ƕ��ַ���ȷ��
	//ֻ���ڷ�3D������²Ż����IN.KPT�ļ�,���ҵ�������Ϊtrueʱ��ż��
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
		//���������ǳ���ʹ��map���ֵ��
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

	//���濪ʼ���,����������ǣ�Ҫô���ֵ���ʣ�Ҫô�������õ���
	if ((vbm_ok&& cbm_ok) || IF_ITERATION==false)
	{
		int top = 0, botom = 0, left = 0, right = 0, left_top = 0, right_top = 0, left_end = 0, right_end = 0;

		//Ȼ��ֱ����vbm��cbm�������
		int band_vbm = rep.eig[rep.num_kpt*rep.num_electron*rep.spin / 2].band;
		int band_cbm = rep.eig[rep.num_kpt*rep.num_electron*rep.spin / 2 - 1].band;
		int spin_vbm = rep.eig[rep.num_kpt*rep.num_electron*rep.spin / 2].spin;
		int spin_cbm = rep.eig[rep.num_kpt*rep.num_electron*rep.spin / 2 - 1].spin;

		if (CALCULATION_MODE == 0) {
		
			//0d�Ļ�ʲô�����øɣ�gap�����Ѿ���������

			return;
		}
		else if (CALCULATION_MODE == 2)
		{
			//����Ҫ�������������
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
			//����Щ�����ӳ�䣬��ֹ����߽�
			for (int i = 0; i < 9; i++) {
				if (vbm_location[i] < 0) {
					vbm_location[i] = 0;
				}
				else if (vbm_location[i] > (max_num - 1)) {
					vbm_location[i] = max_num - 1;
				}
			}
			//����ÿһ��inter����أ����ں������
			for (int i = 0; i < 9; i++) {
				inter[0][i] = rep.kp_p[vbm_location[i]][0];
				inter[1][i] = rep.kp_p[vbm_location[i]][1];
			}
			//�����ǿ��ٻ��x��y�����
			inter[2][0] = rep.eig_kns[vbm_location[0]][band_vbm][spin_vbm] / HARTREE;
			inter[2][1] = rep.eig_kns[vbm_location[1]][band_vbm][spin_vbm] / HARTREE;
			inter[2][2] = rep.eig_kns[vbm_location[2]][band_vbm][spin_vbm] / HARTREE;

			inter[2][3] = rep.eig_kns[vbm_location[3]][band_vbm][spin_vbm] / HARTREE;
			inter[2][4] = rep.eig_kns[vbm_location[4]][band_vbm][spin_vbm] / HARTREE;
			inter[2][5] = rep.eig_kns[vbm_location[5]][band_vbm][spin_vbm] / HARTREE;

			inter[2][6] = rep.eig_kns[vbm_location[6]][band_vbm][spin_vbm] / HARTREE;
			inter[2][7] = rep.eig_kns[vbm_location[7]][band_vbm][spin_vbm] / HARTREE;
			inter[2][8] = rep.eig_kns[vbm_location[8]][band_vbm][spin_vbm] / HARTREE;

				//Ȼ����x�ľ��룬y�ľ���
			averge_x = (dis(rep.kp_p[vbm], rep.kp_p[vbm_location[3]]) + dis(rep.kp_p[vbm], rep.kp_p[vbm_location[5]])) / 2;
			averge_y = (dis(rep.kp_p[vbm], rep.kp_p[vbm_location[7]]) + dis(rep.kp_p[vbm], rep.kp_p[vbm_location[1]])) / 2;
			//ע��۴����ǿ�Ѩ���������ǵ���
			//�����ǵ��ӵ���Ч��������Ϊ����������y��x
			//������Ϊ��ͬһ��band

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


			//Ȼ��ͬ�������cbm��
			top = cbm - hang;
			botom = cbm + hang;
			left = cbm - 1;
			right = cbm + 1;
			left_top = top - 1;
			right_top = top + 1;
			left_end = botom - 1;
			right_end = botom + 1;

			int cbm_location[9] = { left_end,botom,right_end,left,cbm,right,left_top,top,right_top };
			//����Щ�����ӳ�䣬��ֹ����߽�
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
			
			//������������Ϣ
			/*fullfill_inter_real_location(inter, left_end, cell_a, path, 0, band_cbm, spin_cbm, "cb");
			fullfill_inter_real_location(inter, botom, cell_a, path, 1, band_cbm, spin_cbm, "cb");
			fullfill_inter_real_location(inter, right_end, cell_a, path, 2, band_cbm, spin_cbm, "cb");
			fullfill_inter_real_location(inter, left, cell_a, path, 3, band_cbm, spin_cbm, "cb");
			fullfill_inter_real_location(inter, cbm, cell_a, path, 4, band_cbm, spin_cbm, "cb");
			fullfill_inter_real_location(inter, right, cell_a, path, 5, band_cbm, spin_cbm, "cb");
			fullfill_inter_real_location(inter, left_top, cell_a, path, 6, band_cbm, spin_cbm, "cb");
			fullfill_inter_real_location(inter, top, cell_a, path, 7, band_cbm, spin_cbm, "cb");
			fullfill_inter_real_location(inter, right_top, cell_a, path, 8, band_cbm, spin_cbm, "cb");*/
			//Ȼ�����inter		
			inter[5][0] = rep.eig_kns[cbm_location[0]][band_cbm][spin_cbm] / HARTREE;
			inter[5][1] = rep.eig_kns[cbm_location[1]][band_cbm][spin_cbm] / HARTREE;
			inter[5][2] = rep.eig_kns[cbm_location[2]][band_cbm][spin_cbm] / HARTREE;

			inter[5][3] = rep.eig_kns[cbm_location[3]][band_cbm][spin_cbm] / HARTREE;
			inter[5][4] = rep.eig_kns[cbm_location[4]][band_cbm][spin_cbm] / HARTREE;
			inter[5][5] = rep.eig_kns[cbm_location[5]][band_cbm][spin_cbm] / HARTREE;

			inter[5][6] = rep.eig_kns[cbm_location[6]][band_cbm][spin_cbm] / HARTREE;
			inter[5][7] = rep.eig_kns[cbm_location[7]][band_cbm][spin_cbm] / HARTREE;
			inter[5][8] = rep.eig_kns[cbm_location[8]][band_cbm][spin_cbm] / HARTREE;





			//Ȼ����x�ľ��룬y�ľ���
			//ע�������x��y�ľ���ע���˻���ʵ�ռ�ȥ��



			//ͶӰ���x��y

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
			//ע����1d�����ʱ��û�в�ֵ��������������

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
			//�����ά����µĻ�ȡ��Ч����
			//������vbm��
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

			//ͬ��ʼ��cbm��
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

	//����ͬ����Ҫ�޸ģ���Ϊ���ǵ��ٽ���ֵ������



}


void get_all_output(string &out_name, string& chazhi_name, string & calculation, string *config, int**flag, int num)
{
	//���������صĽ����������϶��Ч�����ȡ�

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
	//���д洢vbm��cbm����������
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
		fout << config[i] << "\t";//���org_name
		fout << cell_a.get_formula_name() << "\t";//���formula_name
		fout << cell_a.get_positive_distance() << "\t";//���������֮��ľ���
		//���濪ʼ��ӡ��commment���Լ�֮�������

		
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
			
			//�йذ���������
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


	//��ʽ����õ��ľŸ������ֵ�����
	//cout << "now generate the inter_place number:" << endl;
	//fout2.open(chazhi_name, ios::out);
	////���ö���6λ��Ч����
	//fout2.setf(ios::fixed);
	//fout2.precision(10);

	//for (i = 0; i < num; i++)
	//{
	//	if (flag[i][0] == 99)
	//	{
	//		cout << "generate the inter use number:" << config[i] << endl;
	//		//����ÿ���ṹ����6��,ÿ����9������
	//		//������vbm������,�ֱ���x�����꣬y�������z������
	//		//ע�������ص�����
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
	//��������Ǵ�����ԭʼ�ṹ�ļ�λ��,ָ��·����falg��ʶʲô���ͽṹ����Ӧ��ͬ�Ĳ�������
	//�������ָ��·�������Ӧ�ṹ
	char* target = const_cast<char*>(filename.c_str());
	vector<vector<int>> elements;
	Pattern_cell *pattern = new Pattern_cell(flag, target, 0);
	if (flag == "1T") {

		//����Ԫ���б�		
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
			double factor = 1;//�������������ϵ�������������Ӱ뾶���з���

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
		//��������࣬��������,ע������Ԫ��
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
			double factor = 1;//�������������ϵ�������������Ӱ뾶���з���
			//factor = (e[elements[i][1]].nega_ionic_ridus + e[elements[i][0]].posi_ionic_ridus) / now;
			if (factor < 0.8)
				factor = 1;
			pattern->get_pattern(dir_name, elements[i], c_pattern, factor, "type6");

		}


	}
	else if (flag == "type10") {
		//ͬ��Ҳ����Ԫ�ģ��ֳ�����
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
		//Ȼ���ǵڶ���
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
			double factor = 1;//�������������ϵ�������������Ӱ뾶���з���
			//factor = (e[elements[i][1]].nega_ionic_ridus + e[elements[i][0]].posi_ionic_ridus) / now;
			if (factor < 0.8)
				factor = 1;
			pattern->get_pattern(dir_name, elements[i], c_pattern, factor, "type10");

		}

	}
	else if (flag == "SnBr2_three") {
		//��Ӧ������Ԫ�����
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
			double factor = 1;//�������������ϵ�������������Ӱ뾶���з���

			//factor = (e[elements[i][1]].nega_ionic_ridus + e[elements[i][0]].posi_ionic_ridus) / now;
			if (factor < 0.8)
				factor = 1;
			//��һ��Ԫ��������Ԫ�أ��ڶ����ǻ���Ԫ��
			pattern->get_pattern(dir_name, elements[i], *pattern, factor);

		}

	}
	else if (flag == "1_4_A"||flag=="1_4_B"||flag=="2_1"||flag=="3_1"||flag=="SnBr2") {
		//��Ӧ�����1:4�����
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
			double factor = 1;//�������������ϵ�������������Ӱ뾶���з���

			factor = (e[elements[i][1]].nega_ionic_ridus + e[elements[i][0]].posi_ionic_ridus) / now;
			if (factor < 0.8)
				factor = 1;
			//��һ��Ԫ��������Ԫ�أ��ڶ����ǻ���Ԫ��
			pattern->get_pattern(dir_name, elements[i], *pattern, factor);

		}

	}
	else if (flag == "2_1") {
		//����ǵڶ����
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
			double factor = 1;//�������������ϵ�������������Ӱ뾶���з���

			factor = (e[elements[i][1]].nega_ionic_ridus + e[elements[i][0]].posi_ionic_ridus) / now;
			if (factor < 0.8)
				factor = 1;
			//��һ��Ԫ��������Ԫ�أ��ڶ����ǻ���Ԫ��
			pattern->get_pattern(dir_name, elements[i], *pattern, factor);

		}
	}
	else
		return;

}

bool Pattern_cell::get_pattern(string& name, vector<int>element, cell &cell_pattern, double factor, string flag) {
	//��������Ľṹ���нṹ�������½ṹ
	ofstream fout;
	fout.open(name, ios::out);
	if (!fout.is_open()) {
		cout << "can not generate the file!" << endl;
		return false;
	}
	//��ʼ��ʽд��ṹ�ļ�
	//������ʾ��ʽ
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
	//����������޸ģ�ֱ���滻Ԫ�ص������Ϣ
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
		//�������ǹ涨����һ�������ģ��ڶ����Ǹ���
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
	//���ڲ���n123����ֵ
	vector<int> res(3, 0);
	res[0] = ceil(pow(2 * ECUT2 / 2, 0.5) / PI * 1.889725989*a);
	res[1] = ceil(pow(2 * ECUT2 / 2, 0.5) / PI * 1.889725989*b);
	res[2] = ceil(pow(2 * ECUT2 / 2, 0.5) / PI * 1.889725989*c);
	//Ȼ�����n1��n2��4�ı�����n3����10�ı���
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
	//��������������Ϣ�����ת���������x��y��z
	report rep;
	if (rep.read_report(const_cast<char*>((cell_path + "/REPORT").c_str()), 1) == 1)
	{
		//˵��û����ļ�
		cout << "no report!" << endl;
		cin.get();
		return;
	}
	//�ֱ����vb��cb���������
	double test_kpt[3] = { 0.37582,0.6096,0 };
	double ka = ((2 * PI*get_vecror_mo(cell_a.letice[1], 3) / BOHR) / (get_dot_product(cell_a.letice[0], cell_a.letice[1], 3) / BOHR / BOHR));
	double kb = ((2 * PI*get_vecror_mo(cell_a.letice[0], 3) / BOHR) / (get_dot_product(cell_a.letice[0], cell_a.letice[1], 3) / BOHR / BOHR));
	//ͶӰ���x��y,ע������VBȻ����CB
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
	//��øýṹ�������������������
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
	//�������Ĳ������е��ӵ��������
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
	//��øýṹ���γ��ܵķ���
	//��Ϊ��ṹ�Ѿ��󶨺��ˣ�ֱ�Ӽ��������������

	double res = -100;
	double towd_averge = this->total_energy / this->num;
	double bulk_everge = this->bulk->total_energy / this->bulk->num;
	this->formation_energy_peratom= towd_averge - bulk_everge;
	return;
}
void Cal_cell::get_formation_energy_persquare() {
	//��ýṹ���γ���
	//���㷽���ǻ���ƽ��Ľ����
	int yu = this->bulk->num%this->num;
	double times = this->bulk->num / this->num;
	if (yu != 0) {
		cout << "cell:" << this->name << ": do not sataisify the chemistry pair!" << endl;
	}

	if (CALCULATION_MODE == 2) {
		//ע��������ν�İ����ܣ�Ӧ���Ǹ���ȡ�Լ���Ӧ�����		
		//2dһ����abƽ��
		double a = pow(letice[0][0] * letice[0][0] + letice[0][1] * letice[0][1] + letice[0][2] * letice[0][2], 0.5);
		double b = pow(letice[1][0] * letice[1][0] + letice[1][1] * letice[1][1] + letice[1][2] * letice[1][2], 0.5);
		this->suqre = a * b;		
		//��������Ƿ��ϻ�ѧ������,2d����������relax������
		/*double twod_everge = this->total_energy / (max(this->suqre,this->org_square));*/
		double twod_everge = this->total_energy / this->suqre;
		//��ṹ����relax֮ǰ�����
		//double bulk_everge = this->bulk->total_energy / times / min(this->org_square,this->suqre);
		//cout << "org_square:" << this->org_square << ",after square:" << this->suqre << endl;
		double bulk_everge = this->bulk->total_energy / times / this->suqre;
		this->formation_energy_square = (twod_everge - bulk_everge)/2;
		return;		
	}
	else if (CALCULATION_MODE == 1) {
		//���ʱ����Ҫ��Ч���Բ��������
		//v=(a*b)��c
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
		//��ЧΪ����
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
	this->org_square = a * b;//�����2ά����£�����ԭ�������
	this->h_length = get_length(cell_org->letice[0]);//һ��������a������������
	return;
}

//������һ��report��name�������spinup��down�Ĳ�ͬ��϶���
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
	//��������
	string result = string(allbuf).substr(3);
	int index = result.find_first_of("=");
	string next = result.substr(index+1);
	//cout << result << endl;
	//���տո�ָ��ַ���
	vector<string>string_list = cut_string(next);
	//cout << next << "," << string_list[0] << ","<<string_list[1]<<endl;
	spin_up = stod(string_list[0].c_str());
	spin_down = stod(string_list[1].c_str());
	//�ж�ȡ���ľ���
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
		//������д���Ļ�����ӡ����diff
		res[0] = spin_up;
		res[1] = spin_down;
		res[2] = -1;
		return res;
	}
	//���濪ʼ���ظ��ԵĴ�϶���	
	//1 ���ص�һ������Ҫ������vector
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
	//2 ���ص�2����Ҫ��vector
	vector<double>spin2_low = get_allk_energy(1, down - 1);
	vector<double>spin2_high = get_allk_energy(1, down);
	max_positon = max_element(spin2_low.begin(), spin2_low.end());
	min_position = min_element(spin2_high.begin(), spin2_high.end());
	res[1] = *min_position - *max_positon;
	res[2] = 1;
	cout << "spin_down gap:" << *min_position - now_fermi << endl;
	return res;

}

//��ȡָ��spin��band�����еı�������ֵ
vector<double> report::get_allk_energy(int spin, int band) {
	vector<double> res;
	for (int i = 0; i < this->num_kpt; i++) {
		res.push_back(this->eig_kns[i][band][spin]);
	}
	return res;
}

//���캯����ȡ����ļ�
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


	double correct_value = 0;//����ֵ
	double read_fermi = 0;//��ȡ��OUT.FERMI����ֵ
	bool button = false;//�ļ��Ƿ���ڵñ���
	//Ȼ���ٶ�ȡ�����ܼ�
	ifstream ifer;
	ifer.open(path + "/OUT.FERMI", ios::in);
	if (!ifer.is_open()||!ifer.good()) {
		//˵��û������ļ�����ôֱ�Ӽ�ȥ����õ���fermi�ܼ�
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
	
	//�����̴��ˣ���fermi�ļ�ʱ�Ѿ���ȡ����ؽ��
	
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
		//�Ǿ�ֻʣһ��������������Լ��ģ�����Ҳ�Ҳ���fermi�ļ�
		correct_value = -this->orig_fermi;
		this->fermi = this->orig_fermi;
	}
	//ѭ��������ֵ��ȡ��ȥ
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
//��������ܴ����Ƶ�����
int Bandstru::drawBandstru(int start, int expand) {
	
	int real_start = 0;
	int real_end = 0;
	
	if (start == -1) {
		//start = this->vbm_band-1;
		init_start_band_num();
		start = this->orig_start_band;
		cout << "sure start:" << start << endl;
	}
	//�趨��Χ
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
	//���band�ļ�
	ofstream fout;
	fout.open(out_name, ios::out);
	fout.setf(ios::fixed, ios::floatfield);  // �趨Ϊ fixed ģʽ����С�����ʾ������
	fout.precision(6);
	//���ȴ�ӡ��һ��
	fout << "k_point" << "\t";
	int label = 1;
	for (auto band : bandlist) {
		for (int i = 0; i < (band->real_end - band->real_start); i++)
			fout << label << "\t";
		label++;
	}
	fout << endl;


	//Ȼ���ӡƴ������
	for (int step = 0; step <bandlist[0]->k_num; step++) {
		//���ȴ�ӡ��ͷ
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

//�����޸�֮��������Ҫ��dos������
int Bandstru::drawDos(string path,string out_name) {
	if (Fermi_orig == false) {
		cout << "use read the fermi :" << this->orig_fermi << endl;
	}
	else {
		cout << "use automitail the fermi :" << this->orig_fermi << endl;
	}
	//1 ��ȡ��ص��ļ�
	int lines = 0;
	ifstream fin;
	fin.open(path + "/DOS.totalspin", ios::in);
	if (!fin.is_open()) {
		cout << "i can not find the file!total spin" << endl;
		cin.get();
		return 0;
	}
	fin.close();
	//1.1 ˵��������ļ�������ж�����
	string com = "cat " + path + "/DOS.totalspin" + "| wc -l >line";
	system(com.c_str());
	fin.open("line", ios::in);
	fin >> lines;
	lines -= 1;
	fin.close();
	com = "rm line";
	system(com.c_str());
	//1.2 ����vector������,�ֱ���total��up�Լ�down����ֵ
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

	//1.3 ��ȡdown���ļ�
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
	cout << "��ȡ������dos�����ϣ�" << path << endl;

	//2 Ȼ�������صĽ��
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
	cout << "�������dos��ϣ�" << path<<endl;
	return 1;
}

//�������Ľ�һ��bandstru����dos��������
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
		//////�費��Ҫ��ӡ�ڶ���
		Bandstru* bandstru2 = new Bandstru(path + singel, 2);
		if (bandstru2->drawBandstru(-1,5) == 0) {
			cout << "draw the band:" << singel << ",encounter error!" << endl;
			cin.get();
			continue;
		}
		band_room.push_back(bandstru2);
		cout << "get bandstructure2 finished!" << endl;
		//Ȼ��ſ�ʼ�������
		//������������ݷ���
		Bandstru::output_band(path + singel + "/my_band", band_room);
	}
	cout << "all total work done!" << endl;
	cin.get();
	return 1;
}
//��ʼ����ʼ���������ֵ
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

//��vasp��in��kpt�ĸ�ʽת��Ϊpwmat����ظ�ʽ
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
			//˵�����濪ʼ��ȡ��������ֵ��
			//ÿ�ж�ȡ��������
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

	//���濪ʼ�����ؽ��
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



	//���ִ�нű�������IN,KPT�Ľ��
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
	//�жϽṹ�Ƿ��г����ṹ
	int count = 0;
	for (int i = 0; i < this->num; i++) {
		if (this->type[i] == 8) {
			count++;
		}
	}
	if (count < 2)
		return false;
	double distance = 0;
	//���濪ʼ�����������Oԭ�Ӿ���������
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

//�ж�һЩָ��·�����ܴ���Ϣ
vector<double> report::get_diff(string report_name) {
	//���diff�������Ϣ����������ֵ
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
	//��������
	string result = string(allbuf).substr(3);
	int index = result.find_first_of("=");
	string next = result.substr(index + 1);
	//cout << result << endl;
	//���տո�ָ��ַ���
	vector<string>string_list = cut_string(next);
	//cout << next << "," << string_list[0] << ","<<string_list[1]<<endl;
	res[0] = stod(string_list[0].c_str());
	res[1] = stod(string_list[1].c_str());
	res[2]= stod(string_list[2].c_str());
	return res;
}
vector<vector<double>> report::get_band_value_onkpt(vector<int>&klist,string file_name) {
	//��ȡָ��k���б��Ӧ�����еı�������ֵ
	//���ǵ�����������������������
	//���ǵ�����ռ�����⣬��������һ���ܴ�������
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
	//��ȡָ������¶�Ӧk������
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
	//�ж�ָ����һϵ�еı�������ֵ�Ƿ񴩹����ܴ�
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
	//��ȡָ��·������ķ��������k��
	this->base_kpt_flag = vector<int>(this->num_kpt,-1);
	this->base_kpp = vector<vector<double>>(this->num_kpt, vector<double>(3));
	ifstream fin;
	fin.open(file_name.c_str(), ios::in);
	if (!fin.is_open()) {
		cout << "i can not find the file:" << file_name << endl;
		return false;
	}
	char temp[100];
	//�ȶ�ȡǰ����
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
	//����k���������꣬������k������������ʲô
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
//�ͻ����������ô˺����򵥻�ô�϶�������
//1��ʾ���������2��û������1������
//11��ʾ��������������������ڲ�ͬ�����
//2 ��ʾ��������ע������
int report::check_direction_gap(vector<string>flag, string report_name, string cell_name, string out_kpt_name) {
	//0 ��������Ľṹ
	vector<vector<bool>> judge_res = vector<vector<bool>>(flag.size(), vector<bool>(2));

	//1��ȡ�����ļ�OUT.KPT
	if (!read_base_kpt(out_kpt_name)) {
		cout << "read OUT.KPT failed!" << endl;
		//��ȡʧ������£�����k���ת��
		//cin.get();
		cout << "proceed k point transform!" << endl;
		report_k2_fractional(cell_name);

	}
	//2 ��ȡ��ͬ��ǣ�����k��ļ���
	for (int i = 0; i < flag.size(); i++) {

		vector<bool> if_cross;
		string target = flag[i];
		vector<int>klist = get_all_kptlist(target);
		vector<vector<double>> eigen = get_band_value_onkpt(klist, report_name);
		//Ȼ�������е�ĳ�ж�Ӧ����
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

	//3 �жϲ�ͬ����Ĵ���û��������������ж�ֵ
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
	//��ȡ�����ṹ�ļ�����report����ʵ���껹ԭΪ��������
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

	//���濪ʼת��k���������
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
	//�������Ĳ����ṹ���������ͳ����Ϣ
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
	//����������Դ����Ϣ
	string icsd_path = "/share/home/wangz/icsd_2019/ICSD/Expe.inorganic/Expe.inorganic.atmos/";
	string mp_path = "/share/home/database/BDM/files/";
	string number_name = "";
	number_name = config_name.substr(0, config_name.find_first_of("."));
	int check_line = 5;//������ǰ5�е�����
	char buffer[500];
	char temp[100];
	ifstream fin;
	if (config_name.find("cif") == string::npos) {
		//˵�������mp�Ľṹ
		fin.open(mp_path + number_name + "/" + number_name + ".cif", ios::in);		
	}
	else {
		//˵����icsd������
		fin.open(icsd_path +  number_name + ".cif", ios::in);
	}
	//��ʽ��ȡ����ļ����м��
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
	//���Է����ŵ����
	//����·����ԭʼ��scf��·��
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
	//ִ��ÿ���ṹ�����cohp�ļ���
	//1:׼��lobsiterin�ļ�
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
	//���ﲻʹ�������ֶ�ͳ�ƹ��Ķ���
	//ֱ�Ӳ���ȫ��ͳ�Ƶķ���
	//fout << "includeOrbitals s p d" << endl;
	////Ȼ��ѭ������ÿ��Ԫ������оٳ���
	for (auto ele : this->element_info) {
		fout << "basisfunctions ";
		fout << ele.first << " " << getBaseFunctions(ele.first,true) << endl;
	}
	//���������п��ܵİ뾶��ͳ����Ϣ
	//ֱ�ӽ���Ŀ¼֮ǰ������������Ϣ�ù���
	fout.close();
	//��������޸ģ����ÿһ���������ǵ������з���
	//���Ҫ���������cohp_in������lobster�ļ���
	string temp_file;
	chdir(this->path.c_str());
	for (int i = 0; i < this->ridusElment.size(); i++) {
		temp_file = "lobsterin_" + to_string(i + 1);
		string singel_plus = getLineString(i + 1, this->path + "/" + "cohpin");
		system((string("cp lobsterin_org lobsterin_") + to_string(i+1)).c_str());
		fout.open("lobsterin_" + to_string(i + 1), ios::app);
		fout << singel_plus << endl;
		fout.close();
		//Ȼ��ĳ���ʽ�����ļ�
		system(("mv " + temp_file + " " + "lobsterin").c_str());
		system("lobster-4.0.0");
		//������֮���޸�һЩ�ļ������֣���ֹ����
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

//��ö�Ӧ��ʦ�ļ��Ĺ����Ϣ
string Lobister::getBaseFunctions(string name,bool org_flag) {
	//�������ֺ��ַ��ļ�¼
	vector<string>orbatal;//��¼�Ĺ�����ַ���
	vector<int>orbat_num;//��Ӧ������
	string com;
	string res;
	//֮ǰ������
	bool changed = true;
	if (changed) {
		//���������ֱ�Ӷ�ȡ��д�������Ϣ
		//��Ϊ�����ļ�
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
			//��ʱ��Ͷ�ȡ�Լ��ľ�����
			chdir(this->path.c_str());
			string com = "grep \"VRHFIN =\"" + name + " " + "POTCAR" + " >pot_test";
			system(com.c_str());
		}
		//�и������Ҫ���ַ���
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

		//�������ַ���
		system("rm pot_test out_test");
		
		for (int i = 0; i < orbat_num.size() - 1; i++) {
			res += to_string(orbat_num[i]) + string(orbatal[i]) + " ";
		}
	}

	return res;
}

void Lobister::getElementList() {
	//���Ԫ���б����Ϣ
	cell cell_a((this->path+"/atom.config").c_str(),2);
	
	for (int i = 0; i < cell_a.type_num; i++) {
		pair<string,int> temp=make_pair(a[cell_a.type_save[i]], cell_a.type_save[i]);
		this->element_info.push_back(temp);
	}	
	//Ȼ���ȡԪ��ƥ�����Ϣ
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
	//���ÿ���ṹ�ļ�����Ϣͳ��
	//���ʵ���Ѿ�������
	ifstream fin;
	fin.open(this->path + "/cohp_band", ios::in);
	if (!fin.is_open()) {
		cout << "previous band generate failed!" << this->path << endl;
		this->sucessful = false;
		cin.get();
	}
	//���ڲ�����ridus��Ϣ�����ж����
	//Ȼ���ȡԪ��ƥ�����Ϣ
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
	//���ÿ���ṹ��cohp����ֵ
	if (this->sucessful == false)
		return;
	if (Lobister::COHP_GET == "Special") {
		//��Ӧ�����Լ���ȡ��Ӧcohp�Ļ��ֽ�������ҽ��л���		
		chdir(this->path.c_str());		
		
		//2:��ȡ�������Ϻ��������¶�Ӧ��Χ��cohp��ֵ
		//ע�⣬����ļ��ܼ��Ѿ���ȥ�˷����ܼ�
		//�Ľ������ڶ������������ô��,���ڶ�����������������д���
		bool start_get = false, end_get = false;		
		vector<double>up_start, up_end, down_start, down_end;
		double move_start = this->startRange - this->ferimi;		
		char buff[2048];
		double same_falg = 4e-1;
		int temp_no = 0;
		try {
			//����ÿһ��ɼ�����ֵ�����ɶ�Ӧ����ļ�������
			//�������д���
			for (int i = 0; i < this->ridusElment.size(); i++) {
				//1:��һ����ȡcohp�����Ӧ��N��ֵ
				//���Ȼ�ȡNo����Ŀ
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
				//���濪ʼ��ʽ�Ķ�ȡ
				while (fin.good() && fin.peek() != EOF) {
					if (start_get&&end_get)
						break;
					fin.getline(buff, 2047);
					vector<string>res = cut_string(string(buff));
					double temp = atof(res[0].c_str());
					if (temp < move_start&&abs(temp - move_start) < same_falg && !start_get) {
						//��׽����start����ֵ
						//Ȼ�����ǰѲ�ͬԪ�ض�ƥ���cohp��ֵ���û���
						cout << temp << "," << move_start << endl;
						//int number = getNoForElementPair(this->ridusElment[i]);
						int number = 1;
						up_start.push_back(atof(res[1 + number].c_str()));
						down_start.push_back(atof(res[2 * this->interactions[i] + 3 + number].c_str()));
						start_get = true;
					}
					else if (abs(temp) < 1e-3 && !end_get) {
						//˵����ʱ��׽���˷����ܼ���Ӧ����ֵ
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
				//����˶�Ӧ����ֵ����ʼ������������
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
		//�������ÿһ��Ԫ����Եķ�ʽ���м��������Ϣ
		for (int i = 0; i < this->ridusElment.size(); i++) {
			vector<double> data;//��ǰ����Ԫ�ضԵĶ�Ӧcohp��ֵ
			//��ÿһ��Ԫ��ƥ����е�������
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
					//˵�����ҵ�������ṹ
					vector<string> cutRes = cut_string(string(buff));
					data.push_back(atof(cutRes[cutRes.size() - 1].c_str()));

				}
			}

			fin.close();
			//Ȼ����ʽ��ʼ��д�������
			////�����Լ�д��ƽ��ֵ�Ķ����ɡ���
			double temp_sum = 0;
			for (int i = 0; i < data.size() - 1; i++)
				temp_sum += data[i];

			this->cohp.push_back(temp_sum / (data.size()));
			//�����������Ÿ���һ����ֵ
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

//��̬��������������
//����������ؽ��
void Lobister::runAllProcesses(string file_name,string path,string out_name) {
	//���ڴ��ڶ���뾶ֵ
	//ͬ���Ľṹ����ڶ���
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
				//�޸����������
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
	//�����������Ѿ������ˣ���ֱ������
	//��鵱ǰ�����Ƿ񱨴�
	


	//��ȡ��Ӧ�ķ����ܼ��Ͳ�����Ϣ
	
	fin.open("DOSCAR", ios::in);
	char buff[100];
	//�ȶ�ȡǰ5��
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
	//��������ĳ�ʼ������ֵ
	//ֻ��ָ��Ҫ�󣬷��ܼ���������²ż������ֵ
	if (this->sucessful == false)
		return;
	if (this->COHP_GET == "Special") {
		chdir(this->path.c_str());
		ifstream fin;
		fin.open("DOSCAR", ios::in);
		char buff[300];
		double same_flag = 2e-1;//���жϵ�����Ŀ�������Ƶ�����
		for (int i = 0; i < 6; i++)
			fin.getline(buff, 299);
		//��ʽ��ʼ��ȡ
		//ע������������Ҫ�������
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
					//˵���ҵ��ˣ�ע������Ѱ�ҷ����ܼ������ҵ���һ���������İ�
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
			//���濪ʼ�������
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
	//����ָ����Ԫ�����ӣ��õ���NO.�ڼ����Ķ�Ӧ����
	chdir(this->path.c_str());
	int result = -1;
	string file_name = "COHPCAR.lobster";
	ifstream fin;
	fin.open(file_name, ios::in);
	char buff[500];
	for (int i = 0; i < 3; i++)
		fin.getline(buff, 499);
	//Ȼ����Ƕ�ȡ��NO��Ӧ
	while (fin.good() && fin.peek() != EOF) {
		fin.getline(buff, 499);
		if (strstr(buff, ele_info[0].c_str()) != NULL && strstr(buff, ele_info[1].c_str()) != NULL) {
			int index1,index2;
			index1 = string(buff).find(ele_info[0])+ele_info[0].size();
			index2 = string(buff).find_last_of(ele_info[1])+ele_info[1].size();
			if (buff[index1] - '0' >= 1 && buff[index1] - '0' <= 9 && buff[index2] - '0' >= 1 && buff[index2] - '0' <= 9) {
				//�����������ս��
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
	//��ȡָ���ļ���ָ���ڼ���
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