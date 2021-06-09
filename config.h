#pragma once
#ifndef CONFIG
#define CONFIG
//主要是相关的常量

#define BOHR 0.529177
#define PI  3.1415926535897932
#define HARTREE 27.2114
const int Electro_CAL_NUM = 22;//尝试用更多的点计算二次偏导,注意这个必须是2+4的倍数
const int cengshu = 3;								 //��������������������ھ�����չʱ�Ĳ���
const int yanshen = cengshu * cengshu * cengshu;
const int CENTER = (yanshen - 1) / 2;
const int CALCULATION_MODE = 1;//是否针对1D结构进行计算，默认为2D,1为1d，2为2D，3为3D,0为0D


//临时性的计算设置,以及k点的设置
const bool IF_SURE_KPOINT =false;//目前只有wkm为true时以及2Drelax时候才生效，强制一个给定的k点
const bool IF_K_GENERATE_DIMENTIONAL =false;//对于按照倒格矢生成的k点，是否区分维度
const int default_kpoint[3] = {50,1,1};//强制的k点是什么
const double K_rule = 1000;//控制k点生成的密度,默认的k点生成精度
const double K_rule_for_towD =1500;//控制k点生成的精度，对于2D

//计算监控的设置
const int RULE_NUM = 100;//控制迭代的收敛精度
const double Reverse_space_rule = 0.005;//收敛精度，判据是倒空间的真实坐标
const int RULE_FLAG_NUM = 2;//1表示使用分数坐标，2表示使用倒空间坐标
const int SAVE_STEP = 5;//隔多少任务进行状态的重新记录
const bool IGNORE_WAIT_LINE = false;//无限上交还是有空上交，fasle为有空上交
const bool IF_USE_ALL_QUEUE = true;//是否提交任务的时候占用所有队列,也就是占用gpu2和gpu3
const bool IF_QUEUE_CHECK = false;//是否忽视指定的节点，预留出gpu3的部分节点




//主要用于检查的参数，检查spin 情况
const bool IF_HSE = false;//是否使用hse，reelax无效
const int spin_check = 2;//自旋确认,wkm last无效

const bool IF_PULS_MARNETIC = false;//是否在结构文件末尾补充magnetic，一般用于测试，默认会初始化一个铁磁态
const bool IF_VDW_CHECK = false;//是否在计算中开启vdw修正
const bool IF_ETOT_PLUS_EROR = false;//是否在etot文件中添加上error来进行进一步的计算
const bool IF_PLUS_U = false;//是否使用+U的方法,全部计算均会用到
const bool IF_STRESS = true;//是否设置优化晶格，仅relax使用
const bool IF_AFTEST = false;//是否需要对比测试反铁磁的情况。仅在plus_magnetic开启时才有效，会在scf阶段复制一个文件，定义初始反铁磁态
const bool IF_ACCERATION_SCF = false;//是否读取relax的QDIV，放置到结构上，来一定程度加速scf



//计算中其他的特殊设置
const double SPIN_check_rule = 0.1;//检查自旋为1的rule
const int RE_STEP = 500;//relax的步数限制


//有关半金属的设置
const bool IF_output_half =false;//是否输出结构的半金属的信息
const double Half_band_check = 0.3;//对于half金属情况下，确定bspin_up的阈值
const bool Fermi_orig =false;//是否使用原始fermi，即自动产生的out.fermi


//主要用于dos和wkm的计算参数
const bool IF_WKM = false;//是否使用wkm方法，仅在scf中有用
const bool IF_DOS =false;//控制在scf之后进行什么，true表示进行dos计算
const bool IF_HAND_N123 = false;//是否根据需要认为提高设置n123，仅ifwkm为真时有用，注意只有当你觉得有必要打开此项时再打开
const bool IF_MD = true;//是否在relax之后执行md计算，注意会替换scf计算！！！！！！！！！
const double MD_LENGTH = 12;//md的扩胞的标准
const bool Soc_on = false;//!!!!!!!!!!!!!!!重要，这个会替换scf计算！！以及替换偃师文件！是否在scf计算中开启soc

//计算电子结构是否需要迭代
const bool IF_ITERATION =false;
const int K_point_time = 2;//仅在不需要迭代情况下有效，控制scf的精度是常见的多少倍
const double HSE_minaval = 1;//控制hse的scf精度是一般情况下的多少倍

#endif // !CONFIG
