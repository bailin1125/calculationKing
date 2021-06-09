#include "calculation_tools.h"
#include "Mobility.h"
#include "others.h"
#include "config.h"
#include <fstream>
#include <iostream>
#include <math.h>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#define _CRT_SECURE_NO_WARNINGS
#pragma warning(disable : 4996)
using namespace std;
const char aa[120][3] = { " ", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg" };

const int TASK_NUM =1;

int run_gpu2 = 0;
int run_gpu3 = 0;
string Wannier_cell::pbe_path = "/share2/wangz/0527_got1T/valid_mass_0615/";
string cell::cif_path = "/share/home/wangz/icsd_2019/ICSD/Expe.inorganic/Expe.inorganic.atmos/";
const double cell::SuperO_rule = 1.33;
string Lobister::COHP_GET = "Special";//"Special是分能级，Normal是普通模式，直接读取相应的平均值"
//本程序用来进行辅助材料计算，包括任务提交，产生结构,跟踪report的计算精度，最终输出带隙和有效质量等所有pwmat相关的高通量计算任务
//注意我dtmde，32节点交不上去任务
int main(int argv, char*argc[])
{
	int i = 0, j = 0, k = 0;
	int break_num = 0;
	//第一步：产生atom.config,或者直接提供config
	//第二步：对产生的结构进行relax计算
	//第三步：relax计算结束后进行scf计算，hse
	//第四步：nonscf计算
	//第五步：检测report来决定是否进行进一步的计算


	//必要的准备的路径、文件名等变量
	system("pwd");	
	const string base_path = "/share/home/wangz/high/";
	element *e = new element[120];
	string file_ridus = base_path+"ridus";
	string ionic = base_path + "max_ionic";
	read_element(e, file_ridus,ionic);
	//string calculation_path = "/share/home/wangz/king_calculation/double_1t/";//计算结果放在哪里
	//string  data_path = "/share/home/wangz/king_calculation/double_1t/config/";//储存的结构文件位置

	string calculation_path = "/share2/wangz/1007_formation_energy/1d_new_0422/md/md_calculation/";//计算结果放在哪里
	//string calculation_path = "/share2/wangz/0708_newcreate/calculation_moreaccurate_166/";
	//string calculation_path = "/share2/wangz/0527_got1T/valid_mass_0615/";
	//string  data_path = "/share2/wangz/0708_newcreate/config/166_after_relax/";//储存的结构文件位置
	string  data_path = "/share2/wangz/1007_formation_energy/1d_new_0422/md/recal_finish_relax_config/";//储存的结构文件位置
	string new_yanshi_path = "/share/home/wangz/pwmat/yanshi/PSEUDO/PD04.PBE/";

		
	const static string config_name = base_path+"0423_md";//储存的结构文件名字
	
	string log_name = base_path + "log_all";//日志的名
	string status_name = base_path + "status_all";//任务状态记录的文件名
	string out_putname = base_path + "out_all";//输出最终结果的文件名,这个输出的是x和y方向的统计
	string chazhi_out_name = base_path + "cubic_interploc";//输出的用于插值的文件，每个结构输出三行，9个点的x，y，z
	string yanshi_name = base_path + "yanshi";//储存了拥有的偃势的文件名
	ofstream fout,fout2;
	fout.open(log_name, ios::app);//日志文件记录每个结构到了哪一步
	string gpu2 = "gpu2";
	string gpu3 = "gpu3";	
	string com;
	
	

	//首先获得有关自旋问题的结构名字	
	/*string out_name_spin = "sure_spin";
	detect_spin_eualone(calculation_path, config_name, out_name_spin, flag);
	cout << "spin check finished!" << endl;
	cin.get();*/


	//计算迁移率的模块
	//string mobilty_name = "mobil_name";
	///*cout << "calculation:" << get_mobility(0.97, 0.45, 4.52) << endl;
	//cin.get();*/
	//mobility_client(mobilty_name, calculation_path);
	//


	////产生wannier的模块	
	/*string source_path = "/share2/wangz/0527_got1T/calculation/wannier90/";	
	string wan_name = base_path+"0916_name";	
	generate_winf90(source_path, calculation_path, wan_name);
	cout << "generate wan complete!" << endl;
	cin.get();
	return 0;*/
	


	//批量产生结构ls

	/*string Got_path = "/share2/wangz/0916_dos_and_wannier/bulk_config/";
	string pattern_org = "/share2/wangz/0916_dos_and_wannier/bulk_config/SnBr2.config";
	generate_pattern_configs(pattern_org, Got_path,"SnBr2_three", e);
	cout << "generate config complete!" << endl;
	cin.get();*/


	//所有稀奇古怪的测试放在这里
	Lobister::runAllProcesses("/share/home/wangz/high/0525_name","/share/home/wangz/0522_1d_cohp/calculation/","/share/home/wangz/high/0525_out");
	test(e);
	
	//申请空间
	string *config = new string[TASK_NUM];
	//这里可以修改flag值，中间修改流程
	int **flag = new int*[TASK_NUM];
	for (i = 0; i < TASK_NUM; i++)
	{
		flag[i] = new int[8];
	}
	for (i = 0; i < TASK_NUM; i++)
	{
		for (j = 0; j < 8; j++)
		{
			flag[i][j] = 0;
		}
		flag[i][0] = 2;
		flag[i][1] = 1;//使用gpu23	
		flag[i][2] = 1;//上一步的迭代步数默认从1开始
	}
	ifstream fin;
	fin.open(config_name, ios::in);
	if (!fin.is_open())
	{
		cout << "i can not find the file " << config_name << endl;
		cin.get();
		return 0;
	}
	i = 0;
	while (fin.good() && fin.peek() != EOF && i < TASK_NUM)
	{
		fin >> config[i];
		
		//config[i] = data_path + config[i];
		i++;
	}
	fin.close();

	//然后记录目前拥有的元素偃势
	int yanshi[120] = { 0 };
	fin.open(yanshi_name, ios::in);
	if (!fin.is_open())
	{
		cout << "i can not find the file " << yanshi_name << endl;		
		cin.get();
		return 0;
	}
	
	string temp;
	i = 0;
	while (fin.good() && fin.peek() != EOF)
	{
		fin >> temp;
		for (i = 0; i < 120; i++)
		{
			if (strcmp(aa[i], temp.c_str()) == 0)
			{
				yanshi[i] = 1;
				cout << temp << "  " << aa[i] << endl;
				//cin.get();
				break;
			}
			
		}
		if (i == 120)
		{
			cout << "ukonwn element!" << temp << endl;
			cin.get();
			return 0;
		}
	}
	int mp[3];
	int** mp_all = new int*[TASK_NUM];
	for (i = 0; i < TASK_NUM; i++)
	{
		mp_all[i] = new int[3];
	}
	fin.close();
	//储存计算获得的mp123
	//建立日志文件，监控每个结构的计算情况
	//0 未开始计算
	//1 relax进行，2 relax结束
	//3 scf进行，4 scf结束
	//5 nonscf第一步进行 51nonscf第一步结束
	//61上一步还没有提交任务，71上一步提交了任务
	//8表示wkm计算

	//出错和其他标记
	//7表示对应结构的带隙
	//87表示spin检查中不是spin=1的情况
	//20表示nonscf中发现没有带隙
	//30是二分法错误
	//99 表示最终完成，每一位后面+7表示遇到未知错误
	//第二位标记上一步使用哪个队列，1表示gpu2,2表示gpu3
	////第三位标记仅对nonscf有效，数值表示这个结构任务迭代的次数


	
	//cin.get();
	////然后直接输出最终文件了	
	////cout << config[646] << "," << flag[646][0] << endl;
	//generate_last_file(out_putname, flag,TASK_NUM, config_name);
	//cout << "last result got complete!" << endl;
	//cin.get();
	
	////////先建立relax的文件夹以及脚本
	if (generate_jos(config_name, data_path,calculation_path, TASK_NUM,"generate",yanshi,flag,fout) != 0)
	{
		cout << "generate errors!" << endl;
		cin.get();
	}
	
	////首先重新确认缺少堰势的是哪些结构
	//for (i = 0; i < TASK_NUM; i++)
	//{
	//	if (if_no_yanshi((calculation_path + config[i]), yanshi) == 1)
	//	{
	//		flag[i][0] = 7;
	//		fout << config[i] << " :lack yanshi file!" << endl;
	//		cout << config[i] << " :lack yanshi file!" << endl;
	//	}
	//}


	//for (i = 0; i < TASK_NUM; i++)
	//{		
	//	if (i % 10 == 0)
	//	{
	//		cout << "has got mp_all " << i << endl;
	//	}
	//	give_mp123(calculation_path + config[i], mp_all[i]);	
	//	//sole_generate_jobs(calculation_path + config[i]+"/atom.config", calculation_path + config[i],"relax");
	//}



	////////添加功能，程序中断后根据status和log文件读取当时的状态
	//注意这里需要好好修改status file,注意手动修改后一定要将修改后的顺序命名为新的name文件
	/*string backup_name1 =base_path+ "status_all";
	string backup_log = base_path+"log_all_backup";
	if (backup_status(config_name,backup_name1, backup_log, flag, TASK_NUM, mp_all, -1)==0)
	{
		cout << "backup status successful!" << endl;		
	}*/
	
	////跳出来之后先储存一下结果
	/*chdir("/share/home/wangz/projects/calculation_a/bin/x64/Debug");
	com = "mv " + status_name + " " + status_name + "_backup";
	system(com.c_str());
	fout2.open(status_name, ios::out);
	com.clear();
	for (j = 0; j < TASK_NUM; j++)
	{
		fout2 << config[j] << "\t" << flag[j][0] << "\t" << flag[j][1] << "\t" << flag[j][2] << "\t";
		for (int i = 0; i < 3; i++)
		{
			fout2 << mp_all[j][i] << "\t";
		}
		fout2 << endl;
	}
	fout2.close();
	cout << "has save the save!" << endl;
	cin.get();*/
	//return 0;

	//最后再确认下一下当前的剩余数量	
	
	for (int i = 0; i < TASK_NUM; i++)
	{
		if (flag[i][0] == 99 || flag[i][0] % 10 == 7 || flag[i][0] == 20 || flag[i][0] == 30)
		{
			break_num++;
		}
	}	
	/*cout << "now break_num is:" << break_num << ", others " << TASK_NUM - break_num << endl;
	get_all_output(out_putname, chazhi_out_name, calculation_path, config, flag, TASK_NUM);
	cin.get();
	return 0;*/
	
	/*sole_generate_jobs("/share2/wangz/1007_formation_energy/1d_new_0422/calculation/23262_1d.config/atom.config",
		"/share2/wangz/1007_formation_energy/1d_new_0422/calculation/23262_1d.config", "relax");
	
	sole_generate_jobs("/share2/wangz/1007_formation_energy/1d_new_0422/calculation/23262_1d.config/atom.config",
		"/share2/wangz/1007_formation_energy/1d_new_0422/calculation/23262_1d.config", "scf");*/

	
	i = 0;
	break_num = 0;
	while (true)
	{
		for (int i = 0; i < TASK_NUM; i++)
		{
			if (flag[i][0] == 99 || flag[i][0] % 10 == 7 || flag[i][0]==20 ||flag[i][0]==30)
			{
				break_num++;
			}
		}
		cout << "now break_num is:" << break_num << endl;
		if (break_num == TASK_NUM)
		{
			cout << "already done！please make sure!" << endl;
			//cin.get();	
			break;
		}
		else
			break_num = 0;


		for (i = 0; i < TASK_NUM; i++)
		{
			//每隔一定时间，记录所有任务的标志号，上次使用gpu，以及迭代次数，还有当前精度的mp			
			//cin.get();
			
			if ((i + 1)%SAVE_STEP==0)
			{
				chdir(base_path.c_str());				
				com = "mv " + status_name+" "+status_name+"_backup";
				system(com.c_str());
				fout2.open(status_name, ios::out);
				com.clear();
				for (j = 0; j < TASK_NUM; j++)
				{
					fout2 << config[j] << "\t" << flag[j][0] << "\t" << flag[j][1] << "\t" << flag[j][2] << "\t";
					for (int i = 0; i < 3; i++)
					{
						fout2 << mp_all[j][i] << "\t";
					}
					fout2 << endl;
				}
				cout << "has save the save!" << endl;
			}

			fout2.close();
			fout2.clear();
			if (flag[i][0] > 10 && flag[i][0] % 10 == 7)
			{
				//这种情况表示这个任务发生错误
				cout << config[i] << ": encounter error:" << flag[i][0] << endl;				
				continue;
			}
			else if (flag[i][0] == 99) {

			}
			else
				system("sleep 3");
			if (flag[i][0] == 0)//表示还未开始计算
			{
				
				//重新修改逻辑，监测到
				int detect_flag = detect_free(gpu2);
				if (detect_flag > 0)
				{
					string name = calculation_path + config[i] + "/atom.config";
					
					if (sole_generate_jobs(name, calculation_path + config[i], "relax") != 0)
					{
						cout << "generate etot wrong!" << endl;
						flag[i][0] = 17;
						continue;
						//cin.get();
					}
					if (detect_flag == 1)
					{
						pbs_got(config[i], "relax", calculation_path + config[i], gpu2);
						flag[i][1] = 1;
					}
					else if (detect_flag == 2)
					{
						pbs_got(config[i], "relax", calculation_path + config[i], gpu3);
						flag[i][1] = 2;
					}
					else
					{
						cout << "no free points! wait 20 seconds!" << endl;
						system("sleep 20");
					}
					chdir((calculation_path + config[i]).c_str());
					if (IF_STRESS == true)
					{
						plus_stress_mask(calculation_path + config[i] + "/atom.config");
					}
					if (IF_PULS_MARNETIC == true) {
						puls_magnetic(calculation_path + config[i]);
					}
					
					string comp = "qsub aaa.pbs";
					system(comp.c_str());
					flag[i][0] = 1;
					cout << config[i] << ":relax start!" << endl;
					fout << config[i] << ":relax start!" << endl;					
				}				
				else
				{
					cout << "no free points! wait 20 seconds!" << endl;
					system("sleep 20");
				}
			}


			if (flag[i][0] == 1)//relax进行
			{
				cout << config[i] << ": runing relax!" << endl;
				if (if_finish(calculation_path + config[i], 1) == 1)//正常算完
				{
					//还必须relaxsteps符合步数内才行
					if (check_relax(calculation_path + config[i]) == true)
					{
						cout << config[i] << ": relax finished!" << endl;
						fout << config[i] << ": relax finished!" << endl;
						flag[i][0] = 2;						
					}
					else
					{
						cout << config[i] << ": relax error!please check!" << endl;
						fout << config[i] << ": relax error!please check!" << endl;
						flag[i][0] = 17;
					}
				}
				else if (if_finish(calculation_path + config[i], 1) == 2)//算的有问题
				{
					cout << config[i] << ": relax error!please check!" << endl;
					fout << config[i] << ": relax error!please check!" << endl;
					flag[i][0] = 17;
				}
				else//没有算完
				{
					//没有算完那就继续算呗
				}
			}		

			if (flag[i][0] == 2)//relax结束
			{				
				int detect_flag = detect_free(gpu2);
				if (detect_flag > 0)
				{
					//改变先检查有无提交可能，再修改
					string com;				
					
					/*com = "mv " + calculation_path + config[i] + "/atom.config" + " " + calculation_path + config[i] + "/org.config";
					system(com.c_str());
					com = "mv " + calculation_path + config[i] + "/aaa.pbs" + " " + calculation_path + config[i] + "/aaa_ralax.pbs";
					system(com.c_str());
					com.clear();
					com = "mv " + calculation_path + config[i] + "/final.config" + " " + calculation_path + config[i] + "/atom.config";
					system(com.c_str());
					com.clear();
					com = "mv " + calculation_path + config[i] + "/REPORT" + " " + calculation_path + config[i] + "/REPORT_relax";
					system(com.c_str());
					com.clear();
					com = "mv " + calculation_path + config[i] + "/output" + " " + calculation_path + config[i] + "/output_relax";
					system(com.c_str());

					com = "mv " + calculation_path + config[i] + "/etot.input" + " " + calculation_path + config[i] + "/etot_relax";
					system(com.c_str());*/

								
					string name = calculation_path + config[i] + "/atom.config";
					if (sole_generate_jobs(name, calculation_path + config[i], "scf") != 0)
					{
						cout << "generate scf wrong!" << endl;
						flag[i][0] = 17;
						//continue;
						cin.get();
					}
					plusRelaxQdiv(calculation_path + config[i], "");
					bothTestAfsituation(calculation_path, config[i]);//针对AF进行测试
					plus_etot_error(calculation_path + config[i], "etot.input", "REPORT_scf_test");
					puls_magnetic(calculation_path + config[i]);

					//生成aaa.pbs,并且提交任务
					if (detect_flag == 1)
					{
						pbs_got(config[i], "scf", calculation_path + config[i], gpu2);
						flag[i][1] = 1;
					}
					else if (detect_flag == 2)
					{
						pbs_got(config[i], "scf", calculation_path + config[i], gpu3);
						flag[i][1] = 2;
					}		
										
					//此时对应普通scf的计算
					if (IF_MD == false) {				
											
						cout << config[i] << ":scf start!" << endl;
						fout << config[i] << ":scf start!" << endl;
					}
					else
					{	
						//还需要复制一个文件到路径下面
						com = "cp ~/high/IN.MDOPT " + calculation_path + config[i];
						system(com.c_str());
						//还要对结构进行扩胞操作
						cell cell_md((calculation_path+config[i]+"/atom.config").c_str());
						cell_md.supper_cell_md(calculation_path+config[i]+"/");
						cout << config[i] << ":md start!" << endl;
						fout << config[i] << ":md start!" << endl;
					}
					com.clear();
					chdir((calculation_path + config[i]).c_str());
					com = "qsub aaa.pbs";
					system(com.c_str());
					flag[i][0] = 3;
					
				}				
				else
				{
					//说明没有空闲节点
					cout << "no free points! wait 20 seconds!" << endl;
					system("sleep 20");
				}


			}

			////只想着提交上scf就先停

			if (flag[i][0] == 3)//scf进行G
			{
				if (if_finish(calculation_path + config[i]) == 1)//正常算完
				{
					cout << config[i] << ": scf finished!" << endl;
					fout << config[i] << ": scf finished!" << endl;
					flag[i][0] = 4;
					if(IF_ITERATION==false&&IF_DOS==false)
						flag[i][0] = 99;//算到scf结束就行了
				}
				else if (if_finish(calculation_path + config[i]) == 2)//算的有问题
				{
					cout << config[i] << ": scf error!please check!" << endl;
					fout << config[i] << ": scf error!please check!" << endl;
					flag[i][0] = 37;
				}
				else//没有算完
				{
					//没有算完那就继续算呗
				}
			}


			if (flag[i][0] == 4)//scf结束
			{
				int detect_flag = detect_free(gpu2);
				if (detect_flag > 0) {

					string com = "mv " + calculation_path + config[i] + "/REPORT" + " " + calculation_path + config[i] + "/REPORT_scf";
					system(com.c_str());
					com.clear();
					com = "mv " + calculation_path + config[i] + "/output" + " " + calculation_path + config[i] + "/output_scf";
					system(com.c_str());
					com.clear();
					com = "mv " + calculation_path + config[i] + "/OUT.KPT" + " " + calculation_path + config[i] + "/OUT.KPT_scf";
					system(com.c_str());
					com = "mv " + calculation_path + config[i] + "/aaa.pbs" + " " + calculation_path + config[i] + "/scf.pbs";
					system(com.c_str());
					com = "mv " + calculation_path + config[i] + "/etot.input" + " " + calculation_path + config[i] + "/etot_scf.input";
					system(com.c_str());
					com = "mv " + calculation_path + config[i] + "/OUT.SYMM" + " " + calculation_path + config[i] + "/OUT.SYMM_scf";
					system(com.c_str());

					if (IF_DOS == false) 
					{
						//然后是读取report产生IN.KPT
						/*report rep;
						if (rep.read_report(const_cast<char*>((calculation_path + config[i] + "/REPORT_scf").c_str()), spin_check) == 1)
						{
							flag[i][0] = 47;
							fout << config[i] << " report wrong！" << endl;
							cout << config[i] << " report wrong！" << endl;
							continue;
						}
						
						mp_all[i][0] = rep.mp_x;
						mp_all[i][1] = rep.mp_y;
						mp_all[i][2] = rep.mp_z;*/

						//if (make_efm(rep, calculation_path + config[i], mp_all[i][0], mp_all[i][1], mp_all[i][2]) == 1)
						//{
						//	//没有带隙
						//	flag[i][0] = 47;
						//	fout << config[i] << " scf do not have band gap!" << endl;
						//	cout << config[i] << " scf do not have band gap!" << endl;
						//	continue;
						//}

						//然后生成nonscf的脚本
						/*com = "rm " + calculation_path + config[i] + "/etot.input";
						system(com.c_str());*/
						string name = calculation_path + config[i] + "/atom.config";
						if (sole_generate_jobs(name, calculation_path + config[i], "nonscf") != 0)
						{
							cout << "generate nonscf wrong!" << endl;
							flag[i][0] = 47;
							continue;
							//cin.get();
						}

						//生成aaa.pbs,并且提交任务
						if (detect_flag == 1)
						{
							pbs_got(config[i], "nonscf", calculation_path + config[i], gpu2);
							flag[i][1] = 1;
						}
						else if (detect_flag == 2)
						{
							pbs_got(config[i], "nonscf", calculation_path + config[i], gpu3);
							flag[i][1] = 2;
						}
						com.clear();
						chdir((calculation_path + config[i]).c_str());
						
						com = "qsub aaa.pbs";
						system(com.c_str());
						flag[i][0] = 5;
						fout << config[i] << ":1stnonscf start!" << endl;						
					}
					else//对应DOS的情况
					{
						//业务性代码，做好计算dos的准备工作

						com = "mv " + calculation_path + config[i] + "/OUT.WG "+ calculation_path + config[i] + "/IN.WG";
						system(com.c_str());
						com = "mv " + calculation_path + config[i] + "/OUT.WG_2 " + calculation_path + config[i] + "/IN.WG_2";
						system(com.c_str());
						
						
						string name = calculation_path + config[i] + "/atom.config";
						if (sole_generate_jobs(name, calculation_path + config[i], "dos") != 0)
						{
							cout << "generate dos wrong!" << endl;
							flag[i][0] = 47;
							continue;
							//cin.get();
						}
						if (detect_flag == 1)
						{
							pbs_got(config[i], "dos", calculation_path + config[i], gpu2);
							flag[i][1] = 1;
						}
						else if (detect_flag == 2)
						{
							pbs_got(config[i], "dos", calculation_path + config[i], gpu3);
							flag[i][1] = 2;
						}
						com.clear();
						chdir((calculation_path + config[i]).c_str());
						com = "qsub aaa.pbs";
						system(com.c_str());
						flag[i][0] = 5;
						fout << config[i] << "dos start!" << endl;
					}
			    }
				else
				{
					//说明没有空闲节点
					cout << "no free points! wait 20 seconds!" << endl;
					system("sleep 20");					
				}				

			}

			if (flag[i][0] == 5)//nonscf或者dos进行
			{
				//在这种情况下，需要进行多次迭代的nonscf运算来使其计算完成

				if (if_finish(calculation_path + config[i]) == 1)//正常算完
				{
					cout << config[i] << ": nonscf 1st finished!" << endl;
					fout << config[i] << ": nonscf 1st finished!" << endl;
					flag[i][0] = 51;
					flag[i][0] = 99;
				}
				else if (if_finish(calculation_path + config[i]) == 2)//算的有问题
				{
					cout << config[i] << ": nonscf 1st error!please check!" << endl;
					fout << config[i] << ": nonscf 1st error!please check!" << endl;
					flag[i][0] = 57;
				}
				else//没有算完
				{
					//没有算完那就继续算呗
				}
			}



			//重要的循环nonscf迭代

			if (flag[i][0] > 10 && flag[i][0] % 10 == 1)//第一步nonscf结束，开始进行迭代运算
			{
				//终止的条件是要么有错，要么结束达到了需要的精度
				//开始写循环体了
				string com;
				while (1)
				{
					if (flag[i][0] == 61)
					{
						//还没提交呢
						int detect_flag = detect_free(gpu2);
						if (detect_flag > 0)
						{
							if (detect_flag == 1)
							{
								pbs_got(config[i], "nonscf", calculation_path + config[i], gpu2);
								flag[i][1] = 1;
							}
							else if (detect_flag == 2)
							{
								pbs_got(config[i], "nonscf", calculation_path + config[i], gpu3);
								flag[i][1] = 2;
							}
							com.clear();
							chdir((calculation_path + config[i]).c_str());
							com = "qsub aaa.pbs";
							system(com.c_str());
							flag[i][0] = 71;
							fout << config[i] << ":nonscf start!" << flag[i][2] << endl;
							break;
						}						
						else
						{
							//说明没有空闲节点
							cout << "no free points! wait 20 seconds!" << endl;
							system("sleep 20");
							break;
						}

					}
					else if (flag[i][0] == 71 || flag[i][0] == 51)
					{
						//上一步算完了没有？？
						if (if_finish(calculation_path + config[i]) == 1)//正常算完
						{
							cout << config[i] << ": nonscf  finished!" << endl;
							fout << config[i] << ": nonscf  finished!" << endl;
							flag[i][0] = 71;
							int max = mp_all[i][0];
							if (mp_all[i][1] > mp_all[i][0])
							{
								max = mp_all[i][1];
							}
							//这里需要有两种去判断是不是收敛了
							chdir((calculation_path + config[i]).c_str());
							double now_space_rule = report::get_accuracy("REPORT");
							if((RULE_FLAG_NUM==1 && max >= RULE_NUM)|| (RULE_FLAG_NUM==2&&now_space_rule<Reverse_space_rule))							
							{
								if ((check_now_erfen_pk(calculation_path + config[i], 1, flag[i][2]) != 0))
								{
									//达到了收敛标准
									flag[i][0] = 99;
									fout << config[i] << ":succeed!!" << endl;
									cout << config[i] << ": succeed!" << endl;
									break;
								}
								else
								{
									fout << config[i] << ": two_cut method ecounter wrong!please check by hand!" << endl;
									cout << config[i] << ": two_cut method ecounter wrong!please check by hand!" << endl;
									flag[i][0] = 30;
									break;
								}
								
							}

							else
							{
								int cycle_num = flag[i][2];
								//保存上一步的IN.KPT
								com = "mv " + calculation_path + config[i] + "/IN.KPT " + calculation_path + config[i] + "/IN.KPT_" + to_string(cycle_num);
								system(com.c_str());
								//保存这步的report等
								/*if (check_now_erfen_pk(calculation_path + config[i],0,flag[i][2]) == 0)
								{
									fout << config[i] << ": two_cut method ecounter wrong!please check by hand!" << endl;
									cout << config[i] << ": two_cut method ecounter wrong!please check by hand!" << endl;
									flag[i][0] = 30;
									break;
								}*/
								com = "mv " + calculation_path + config[i] + "/REPORT " + calculation_path + config[i] + "/REPORT_" + to_string(cycle_num);
								system(com.c_str());
								com = "mv " + calculation_path + config[i] + "/output " + calculation_path + config[i] + "/output_" + to_string(cycle_num);
								system(com.c_str());
								com = "mv " + calculation_path + config[i] + "/OUT.KPT " + calculation_path + config[i] + "/OUT.KPT_" + to_string(cycle_num);
								system(com.c_str());
								flag[i][2]++;

								//更新新的mp数值，尝试产生kpt
								if (CALCULATION_MODE ==2)
								{
									for (int j = 0; j < 2; j++)
									{
										mp_all[i][j] = 2 * mp_all[i][j];
									}
								}
								else if (CALCULATION_MODE == 1)
								{
									mp_all[i][0] = 2 * mp_all[i][0];
									mp_all[i][1] = 1;
								}
								mp_all[i][2] = 1;
								report temp_rep;
								temp_rep.read_report(const_cast<char*>((calculation_path + config[i] + "/REPORT_" + to_string(cycle_num)).c_str()),spin_check);
								if (make_efm(temp_rep, calculation_path + config[i], mp_all[i][0], mp_all[i][1], mp_all[i][2]) == 1)
								{
									//没有带隙
									flag[i][0] = 20;
									cout << config[i] << " nonscf do not have band gap!" << endl;
									fout << config[i] << " nonscf do not have band gap!" << endl;
									break;
								}

								//然后检查有没有空闲节点
								int detect_flag = detect_free(gpu2);
								if (detect_flag > 0)
								{
									if (detect_flag == 1)
									{
										pbs_got(config[i], "nonscf", calculation_path + config[i], gpu2);
										flag[i][1] = 1;
									}
									else if (detect_flag == 2)
									{
										pbs_got(config[i], "nonscf", calculation_path + config[i], gpu3);
										flag[i][1] = 2;
									}
									com.clear();
									chdir((calculation_path + config[i]).c_str());
									com = "qsub aaa.pbs";
									system(com.c_str());
									flag[i][0] = 71;
									fout << config[i] << ":nonscf start!" << flag[i][2] - 1 << endl;
									cout << config[i] << ":nonscf cycle!:" << flag[i][2] - 1 << endl;
									break;
								}								
								else
								{
									//说明没有空闲节点，标记为这一步还没有提交
									fout << config[i] << ":nonscf wait!" << flag[i][2] - 1 << endl;
									cout << config[i] << ":nonscf cycle!wait!" << flag[i][2] - 1 << endl;
									flag[i][0] = 61;
									break;
								}

							}

						}
						else if (if_finish(calculation_path + config[i]) == 2)//算的有问题
						{
							cout << config[i] << ": nonscf  error!please check!" << endl;
							fout << config[i] << ": nonscf  error!please check!" << endl;
							flag[i][0] = 67;
							break;
						}
						else//没有算完
						{
							//没有算完那就跳出去下一个
							break;
						}

					}
					else
					{
						cout << "wrong!" << endl;
						cin.get();
					}
				}
			}
		}
		
		cout << "sleep 400s ,wait long time!" << endl;
		system("sleep 800");
		/*get_all_output(out_putname, chazhi_out_name, calculation_path, config, flag, TASK_NUM);
		cout << "generate the result file completed!" << endl;*/
	}
	//跳出来之后先储存一下结果
	chdir(base_path.c_str());
	com = "mv " + status_name + " " + status_name + "_backup";
	system(com.c_str());
	fout2.open(status_name, ios::out);
	com.clear();
	for (j = 0; j < TASK_NUM; j++)
	{
		fout2 << config[j] << "\t" << flag[j][0] << "\t" << flag[j][1] << "\t" << flag[j][2] << "\t";
		for (int i = 0; i < 3; i++)
		{
			fout2 << mp_all[j][i] << "\t";
		}
		fout2 << endl;
	}
	cout << "has save the save!" << endl;
	//到这里跳出了大循环
	fout.close();
	fout.clear();
	

	//最后在这里根据相应的结果文件，我们输出有效质量和带隙
	//这里加一个功能，输出化学式,完整的输出需要输出三行，每行九个数字,vbm和cbm都需要进行生成

	//                                                                                                                                                                                          get_all_output(out_putname, chazhi_out_name, calculation_path, config, flag, TASK_NUM);
	get_all_output(out_putname, chazhi_out_name, calculation_path, config, flag, TASK_NUM);
	cout << "generate the result file completed!" << endl;

	delete[]config;
	for (i = 0; i < TASK_NUM; i++)
	{
		delete[]flag[i];
		delete[]mp_all[i];		
	}	
	delete[]mp_all;
	delete[]flag;		
	cout << "all calculation finished!" << endl;
	return 0;
}
