//韦布分布 潜在扇区错误  磁盘周期清洗
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <iostream>
#include <vector>
#include <set>
#include <math.h>
using std::vector;
using std::set;


// Total max disk numbers
#define MAX_DISKS 50000
#define Block 1000000

// Disk event type
#define DISK_FAIL 0  // disk failure
#define DAILY_REPAIR 1 // disk daily repair
#define DISK_REPAIRED 2 // disk failure repaired over 
#define DISK_WARN 3 // disk is predicted to fail
#define DISK_LD 4 // disk latent section failure	//硬盘潜在部分故障 潜在错误 跟 DISK_WARN有什么不一样？
#define DISK_SCRUB 5 // disk scrub					//磁盘清理
#define NODE_FAIL 6 // node failure
#define NODE_REPAIRED 7 // node failure repaired over

// data coding mode
#define Twocopy 0
#define Threecopy 1
#define RShorizontal 2 // 水平放置
#define RSvertical 3

// event and disk structure
typedef struct
{
	int type;
	long time;
	int disk_no;
	int node_no;
	int rack_no;
} event;
event event_heap[MAX_DISKS * 2];
int events;

struct disk_struct
{
	int disk_no;
	int node_no;
	int rack_no;
	// 用于两个Disk结构体的比较
	bool operator == (const struct disk_struct& rhs) {
		return disk_no == rhs.disk_no && node_no == rhs.node_no && rack_no == rhs.rack_no;
	}
};

typedef struct
{
	disk_struct disks[MAX_DISKS];
	int front = 0, end = 0;
}disk_array;

// gsl variables
const gsl_rng_type *T;
gsl_rng *r;

// extern functions 
// typedef int (*array_failed)(int d, int n, int r);
typedef int (*array_failed)(int b, int d, int n, int r);
extern int __checkheap();
extern void heap_insert(long time, int type, int disk_no, int node_no, int rack_no);
extern int heap_search(int type, int disk_no, int node_no, int rack_no);
extern void heap_delete(int i);

inline void testCombine();

// log utility functions
#define TDEBUG
#ifdef TDEBUG
#define _TRACE(d, time, str) printf("%s[%d %d %d] time=%ld\n", str, d.disk_no, d.node_no, d.rack_no, time)
#else
#define _TRACE(d, time, str) {}
#endif // DEBUG
#define TRACE_LD(d, time) _TRACE(d, time, "DISK_LD")
#define TRACE_SCRUB(d, time) _TRACE(d, time, "DISK_SCRUB")
#define TRACE_FAIL(d, time) _TRACE(d, time, "DISK_FAIL")
#define TRACE_WARN(d, time) _TRACE(d, time, "DISK_WARN")
#define TRACE_REPAIR(d, time) _TRACE(d, time, "DISK_REPAIR")
#define TRACE_DLOSS(d, time) _TRACE(d, time, "DATA_LOSS")
#define TRACE_NODE_FAIL(d, time) _TRACE(d, time, "NODE_FAIL")
#define TRACE_NODE_REPAIR(d, time) _TRACE(d, time, "NODE_REPAIR")

// data loss functions
// extern int Threecopy_failed(int d, int n, int r);
int Threecopy_failed(int b, int d, int n, int r);
int is_cross_multiple_rack(int coding, const vector<disk_struct>& failed_d_num,
	const vector<disk_struct>& failed_n_num, disk_struct& disk_target, int rack_size);
int which_situation_is(const vector<disk_struct>& failed_d_num,
	const vector<disk_struct>& failed_n_num, disk_struct& disk_target);
array_failed func[] = { Threecopy_failed, Threecopy_failed, Threecopy_failed,Threecopy_failed,Threecopy_failed };

/*****************************************************
* Global variables which can be used anywhere
****************************************************/
int block, disk, node, rack, coding, rs_k, rs_m;
double Lyita_op,  Lbita_op, Lyita_tia, Lbita_tia, fdr; // failure & prediction
double Lyita_ld,  Lbita_ld,  Lgama_scrub, Lyita_scrub, Lbita_scrub; // ld & scrub
double Lyita_node, Lbita_node, Lgama_node_rest, Lyita_node_rest, Lbita_node_rest; // node failure & repair
double Lgama_rest, Lyita_rest, Lbita_rest; // disk repair
										   // warn_disks的存在感觉没有价值，并不会用来判断data_loss
disk_array failed_disks, warn_disks, ld_disks, failed_nodes;


double percent_avail_bandwidth = 1;						// 可用的网络带宽百分比
double data_per_disk = pow(2, 20);						// 1TB 2^20MB 
#define TOTAL_BANDWIDTH 125 * 60 * 60					// 125MB/s = 1Gb/s 总带宽， 目前系统实际可用的带宽为 percent_avail_bandwidth * TOTAL_BANDWIDTH，最后 * 3600,会换算一下单位（以h为单位）
#define RShorizontal_RepairBandwidth  0.01				// 垂直散布，每次修复消耗的系统带宽的百分比
#define RSvertical_RepairBandwidth    0.01				// 水平散步，每次修复消耗的系统带宽的百分比


long long factorial(long number)
{
	//if (number <= 1)
	//	return 1;
	//else
	//	return number*factorial(number - 1);
	vector<long long> res(number + 1, 1);
	if (number == 1) {
		return 1;
	}
	for (int i = 2; i <= number; i++) {
		res[i] = res[i - 1] * i;
	}
	return res[number];
}

long long combinator(int n, int m)
{
	if (n == 0 || n == m) {
		return 1;
	}
	int temp;
	if (n < m)
	{
		temp = n;
		n = m;
		m = temp;
	}

	long long res = 1;
	for (int i = n; i >= (n - m); i--) {				// 优化阶乘运算
		res *= i;
	}
	long long temp1 = factorial(m);
	return res / temp1;
}



/*****************************************************
* Data loss functions
****************************************************/
template <typename T>
void combine_inner(T &data, int start, int n, int m, int depth, T temp, vector<T> &result)
{
	if (depth == m - 1)
	{
		//最内层循环 将temp加入result
		for (int i = start; i < n - (m - depth - 1); ++i)
		{
			temp[depth] = data[i];
			result.push_back(temp);
		}
	}
	else
		for (int i = start; i < n - (m - depth - 1); ++i)
		{
			temp[depth] = data[i];//每层输出一个元素
			combine_inner(data, i + 1, n, m, depth + 1, temp, result);
		}
}

// 返回从k中选出m个的下标值
//combine(4, -1); combine(4, 0); combine(4, 1); combine(4, 2); combine(4, 3); combine(4, 4); combine(4, 5);
vector<vector<int>> combine(int k, int m)
{
	if (m <= 0 || m > k)
		return{};
	vector<int> data;
	for (int i = 0; i < k; ++i) data.push_back(i);	// 56
	if (k == m) return{ data };
	int depth = 0;
	vector<vector<int>> result;
	vector<int> temp(m, 0);
	combine_inner(data, 0, data.size(), m, depth, temp, result);
	printf("k:%d, m:%d, size:%d\n", k, m, result.size());	// 56 5 3819816 溢出 异常值

	//for (auto it : result) {
	//	for (auto i : it) printf("%d ", i);
	//	printf("\n");
	//}
	return result;
}
// 从故障数组中选择m个故障
vector<vector<disk_struct>> combine_darray(disk_array& darray, int m) {
	vector<vector<disk_struct>> results;
	int k = (darray.end + MAX_DISKS - darray.front) % MAX_DISKS;			// 故障队列的元素的size会很大？
	vector<vector<int>> coms = combine(k, m);
	for (vector<int> vi : coms) {
		vector<disk_struct> res;
		// 例如index为(0,1,3)，那么将drray中第0,1,3个元素拿出来放入res中
		for (int index : vi) {
			res.push_back(darray.disks[(index + darray.front) % MAX_DISKS]);
		}
		results.push_back(res);
	}
	return results;
}

//判断dataloss 调用
//根据传入的编码方案，以及故障的磁盘数目failed_d_num和节点数据failed_n_num，计算数据丢失概率
double func_probability(int coding, const vector<disk_struct>& failed_d_num, 
	const vector<disk_struct>& failed_n_num, disk_struct& disk_target) {
	// 跳过不包含当前disk的磁盘
	bool is_find = false;
	for (auto it : failed_d_num) {
		if (it == disk_target) {
			is_find = true;
			break;
		}
	}
	for (auto it : failed_n_num) {
		if (it == disk_target) {
			is_find = true;
			break;
		}
	}
	// 1. 我们需要跳过不包含当前故障的disk的磁盘
	if (!is_find) {
		return 0.0; 
	}

	// 2. 我们需要跳过不可能引起故障的情形，例如对于RShorizontal，如果不跨越rs_m+1个机架，就不可能引起数据丢失
	if (coding == RShorizontal && !is_cross_multiple_rack(coding, failed_d_num, failed_n_num, disk_target, rs_m)) { // 跨越 m + 1 个机架
		return 0.0;
	}
	else if (coding == RSvertical && !is_cross_multiple_rack(coding, failed_d_num, failed_n_num, disk_target, 1)) { // 跨越 1 个机架
		return 0.0;
	}
	else if (coding == Threecopy && !is_cross_multiple_rack(coding, failed_d_num, failed_n_num, disk_target, 1)) {  // 跨越 1 个机架
		return 0.0;
	}


	
	// 3. 针对不同的编码方案和故障情形，计算相应的概率值  -> && 判断相应错误的 是不是返回真 如果是的话 计算丢失的概率
	// if (coding == RShorizontal && RShorizontal_failed(block, disk, node, rack)) {
	if (coding == RShorizontal) {
		// 计算相应的概率
		// 分两种情况， 如果故障磁盘的数据为0，全部为节点故障 ljf
		if (failed_d_num.size() == 0) {
			//double temp1 = combinator(node - 1, 1);
			//double temp2 = combinator(rack - 1, 1);
			//double temp3 = combinator(node, 2);

			double temp1 = combinator(rack - rs_m - 1, rs_k - 1);			// n 很大，因为rack很大，可能会有溢出的问题
			double temp2 = combinator(rack - 1, rs_k + rs_m - 1);
			double p_node = temp1 * pow(node, rs_k - 1) / (temp2 * pow(node, rs_k + rs_m - 1));
			return 1 - pow((1 - p_node), block * disk);
		}
		
		else{	// 磁盘故障不为0，又是一种情况 ljf
			int x = failed_d_num.size();
			double temp1 = combinator(rack - 1 - rs_m, rs_k - 1);
			double temp2 = combinator(rack - 1, rs_k + rs_m - 1);

			double p_disk = pow(disk, rs_m + 1 - x) * temp1 * pow(node, rs_k - 1) * pow(disk, rs_k - 1) / (temp2 * pow(node, rs_k + rs_m - 1) * pow(disk, rs_k + rs_m - 1));
			return 1 - pow((1 - p_disk), block);
		}
	}
	else if (coding == RSvertical) {	// 怎么判断什么样的情况呢？ 3 1 还是 2 2 还是 2 1 1？
		//int type = RSvertical_failed(block, disk, node, rack);
		int type = which_situation_is(failed_d_num, failed_d_num, disk_target);
		if (type == 0) {
			return 0.0;
		}
		else if (type == 1) { // 1个机架上3个节点发生了故障，另外1个机架1个节点发生了故障
			// 分两种情况， 如果故障磁盘的数据为0,都是节点故障
			if (failed_d_num.size() == 0) {
				double temp1 = combinator(node - 1, 2);
				double temp2 = combinator(rack - 2, 1);
				double temp3 = combinator(node, 3);
				double temp4 = combinator(rack - 1, 2);

				double H_node = temp1 *  temp2 * temp3 / (temp1 *  temp4 * temp3 * temp3);
				return 1 - pow((1 - H_node), block * disk);
			}else{
				int x = failed_d_num.size();
				double temp1 = combinator(node - 1, 2);
				double temp2 = combinator(rack - 2, 1);
				double temp3 = combinator(node, 3);
				double temp4 = combinator(rack - 1, 2);
				double temp5 = combinator(node, 3);

				double H_disk = temp1 * temp2 * temp3 * pow(disk, 9 - x) / (temp1 * temp4 * temp5 * temp5 * pow(disk, 8)) ;
				return 1 - pow((1 - H_disk), block);
			}
		}
		else if (type == 2) { //两个机架上2个节点发生了故障
			// 分两种情况， 如果故障磁盘的数据为0,都是节点故障
			if (failed_d_num.size() == 0) {
				double temp1 = combinator(node - 2, 1);
				double temp2 = combinator(rack - 2, 1);
				double temp3 = combinator(node, 3);
				double temp4 = combinator(node - 1, 2);
				double temp5 = combinator(rack - 1, 2);

				double HH_node = temp1 * temp1 *  temp2 * temp3 / (temp4 * temp5 * temp3 * temp3);
				return 1 - pow((1 - HH_node), block * disk);
			}else{
				int x = failed_d_num.size();
				double temp1 = combinator(node - 2, 1);
				double temp2 = combinator(rack - 2, 1);
				double temp3 = combinator(node, 3);
				double temp4 = combinator(node - 1, 2);
				double temp5 = combinator(rack - 1, 2);
				double temp6 = combinator(node, 3);

				double HH_disk = temp1 * temp1 * temp2 * temp3 * pow(disk, 9 - x) / (temp4 * temp5 * temp6 * temp6 * pow(disk, 8)) ;
				return 1 - pow((1 - HH_disk), block);
			}
		}
		else if (type == 3) { // 1个机架的2个节点发生了故障，另外两个机架的1个节点发生故障
			// 分两种情况， 如果故障磁盘的数据为0,都是节点故障
			if (failed_d_num.size() == 0) {
				double temp1 = combinator(node - 1, 2);
				double temp2 = combinator(rack - 1, 2);
				double temp3 = combinator(node, 3);
				double temp4 = combinator(node - 2, 1);

				double HHH_node = temp1 * temp1 *  temp4 / (temp1 * temp2 * temp3 * temp3);
				return 1 - pow((1 - HHH_node), block * disk);
			}else{
				int x = failed_d_num.size();
				double temp1 = combinator(node - 1, 2);
				double temp2 = combinator(rack - 1, 2);
				double temp3 = combinator(node, 3);
				double temp4 = combinator(node - 2, 1);

				double HHH_disk = temp1 * temp1 * temp4 * pow(disk, 9 - x) / (temp1 * temp2 * temp3 * temp3 * pow(disk, 8)) ;
				return 1 - pow((1 - HHH_disk), block);
			}
		}
	}
	else if (coding == Threecopy) {
		// 分两种情况， 如果故障磁盘的数据为0,都是节点故障
		if (failed_d_num.size() == 0) {
			double temp1 = combinator(node - 1, 2);
			double temp2 = combinator(rack - 1, 1);
			double temp3 = combinator(node, 2);

			double copy_node = 1 / (temp1 * temp2 * node  + (rack - 1) * temp3);
			return 1 - pow((1 - copy_node), block * disk);
		}else{
			int x = failed_d_num.size();
			double temp1 = combinator(node - 1, 1);
			double temp2 = combinator(rack - 1, 1);
			double temp3 = combinator(node, 2);

			double copy_disk = pow(disk, 3 - x) / ((temp1 * temp2 * node +  temp2 * temp3) * pow(disk, 2)) ;
			return 1 - pow((1 - copy_disk), block);
		}
		
	}
	else {
		printf("ERROR coding!");
		exit(-1);
	}
	
	return 0.0;
}
// 精度 根据概率判断是否发生数据丢失？
bool data_loss_with_prob(double prob) {
	if ((rand() % 100) < prob * 100) {
		return true;
	}
	return false;
	// return true;
}


// 
// 从node故障和disk故障数组中，每次选择m个故障，然后判定故障情形，并计算这种故障导致的概率 -> 调用的时候，替代 fun
disk_struct failed[MAX_DISKS];
bool is_data_loss_with_m_fails(disk_struct disk, int event_type) {
	// 其实有一种不需要使用跳过disk的方法，也就是说，在求combine的时候修改diskarray的end位置，进行--操作
	// 然后再把当前的disk加入进去。由于时间关系，先简单的寻找m+1个failed，然后跳过不包含当前故障的disk。

	int m = rs_m + 1;
	// x个磁盘故障，x >= 0, m-x个节点故障
	for (int x = 0; x <= m; ++x) {
		int n_nums = m - x;
		vector<vector<disk_struct>> disks_vec = combine_darray(failed_disks, x);
		vector<vector<disk_struct>> nodes_vec = combine_darray(failed_nodes, n_nums);

		// 选择x个磁盘和m-x个节点故障，并判断是否发生数据丢失
		for (auto i : disks_vec) {
			for (auto j : nodes_vec) {
				// 如果选择的故障磁盘+节点个数总和不等于m，那么不可能引起数据丢失
				if (i.size() + j.size() != m) continue;

				// ljf 增加计算概率函数 这个函数是个累加操作？
				double prob = func_probability(coding, i, j, disk);
				if (data_loss_with_prob(prob)) return true;
			}
		}
	}

	return false;
}


/*****************************************************
* Basic event functions for disk and node
****************************************************/
inline void events_for_new_disk(disk_struct d)
{
	int j = d.disk_no, k = d.node_no, i = d.rack_no;
	long time_op, time_ld, time_scrub, time_tia;
	// disk failure and prediction
	time_op = (long)gsl_ran_weibull(r, Lyita_op, Lbita_op);
	heap_insert(time_op, DISK_FAIL, j, k, i);
	if ((rand() % 10000) < fdr * 100)
	{
		time_tia = 1 + (long)gsl_ran_weibull(r, Lyita_tia, Lbita_tia);
		heap_insert((time_op - time_tia) >= 0 ? (time_op - time_tia) : 0, DISK_WARN, j, k, i);
	}
	// ld and scrub will be ommited in the insert
	time_ld = (long)gsl_ran_weibull(r, Lyita_ld, Lbita_ld);
	heap_insert(time_ld, DISK_LD, j, k, i);
	time_scrub = (long)(Lgama_scrub + gsl_ran_weibull(r, Lyita_scrub, Lbita_scrub));
	heap_insert(time_scrub, DISK_SCRUB, j, k, i);
}
inline void event_failure_for_node(disk_struct d)
{
	long time_op = (long)gsl_ran_weibull(r, Lyita_node, Lbita_node);
	heap_insert(time_op, NODE_FAIL, -1, d.node_no, d.rack_no);
}
inline void event_repair_for_node(disk_struct d, long time)
{
	// long time_mttr = (long)(Lgama_node_rest + gsl_ran_weibull(r, Lyita_node_rest, Lbita_node_rest));
	long time_mttr = 0;
	//if (coding == RShorizontal && percent_avail_bandwidth > 0) {
	//	time_mttr = disk * data_per_disk *  rs_k / (percent_avail_bandwidth * TOTAL_BANDWIDTH);			// RS水平散布，生成修复事件，总共需要跨rs_k个机架，需要修复disk个，每次传输per_disk，传输的数据总量为disk * data_per_disk *  rs_k
	//	percent_avail_bandwidth -= RShorizontal_RepairBandwidth;			
	//}
	//else if (coding == RSvertical && percent_avail_bandwidth > 0){
	//	time_mttr = disk * data_per_disk *  2 / (percent_avail_bandwidth * TOTAL_BANDWIDTH);	// RS垂直散布，假定了数据可以先在机架内部做计算，RS（6,3）数据散布在两个机架上，只需要传输2个磁盘的数据量就可以了
	//	percent_avail_bandwidth -= RSvertical_RepairBandwidth;	
	//}
	//else if(percent_avail_bandwidth <= 0){
	//	printf("The current network bandwidth is not enough to fix the fault and exit... \n");
	//	return;
	//}



	if (coding == RShorizontal) {
		time_mttr = disk * data_per_disk *  rs_k / (percent_avail_bandwidth * TOTAL_BANDWIDTH);			// RS水平散布，生成修复事件，总共需要跨rs_k个机架，需要修复disk个，每次传输per_disk，传输的数据总量为disk * data_per_disk *  rs_k	
	}
	else if (coding == RSvertical){
		time_mttr = disk * data_per_disk *  2 / (percent_avail_bandwidth * TOTAL_BANDWIDTH);	// RS垂直散布，假定了数据可以先在机架内部做计算，RS（6,3）数据散布在两个机架上，只需要传输2个磁盘的数据量就可以了
	}

	heap_insert(time + time_mttr, NODE_REPAIRED, -1, d.node_no, d.rack_no);
}
inline void event_scrub_for_disk(disk_struct d, long time) 
{
	long time_scrub = (long)(time + Lgama_scrub + gsl_ran_weibull(r, Lyita_scrub, Lbita_scrub));
	heap_insert(time_scrub, DISK_SCRUB, d.disk_no, d.node_no, d.rack_no);
}
inline void event_ld_for_disk(disk_struct d, long time)
{
	long time_ld = time + (long)gsl_ran_weibull(r, Lyita_ld, Lbita_ld);
	heap_insert(time_ld, DISK_LD, d.disk_no, d.node_no, d.rack_no);
}
inline void event_repair_for_disk(disk_struct d, long time)
{
	// long time_mttr = (long)(Lgama_rest + gsl_ran_weibull(r, Lyita_rest, Lbita_rest));
	long time_mttr = 0;
	//if (coding == RShorizontal && percent_avail_bandwidth > 0) {
	//	time_mttr = data_per_disk *  rs_k / (percent_avail_bandwidth * TOTAL_BANDWIDTH);			// RS水平散布，生成修复事件，总共需要跨rs_k个机架，需要修复disk个，每次传输per_disk，传输的数据总量为disk * data_per_disk *  rs_k
	//	percent_avail_bandwidth -= RShorizontal_RepairBandwidth;			
	//}
	//else if (coding == RSvertical && percent_avail_bandwidth > 0){
	//	time_mttr = data_per_disk * 2 / (percent_avail_bandwidth * TOTAL_BANDWIDTH);				// RS垂直散布，有多少个rack,就传多少个磁盘
	//	percent_avail_bandwidth -= RSvertical_RepairBandwidth;	
	//}
	//else if(percent_avail_bandwidth <= 0){
	//	printf("The current network bandwidth is not enough to fix the fault and exit... \n");
	//	return;
	//}

	if (coding == RShorizontal) {
		time_mttr = data_per_disk *  rs_k / (percent_avail_bandwidth * TOTAL_BANDWIDTH);			// RS水平散布，生成修复事件，总共需要跨rs_k个机架，需要修复disk个，每次传输per_disk，传输的数据总量为disk * data_per_disk *  rs_k	
	}
	else if (coding == RSvertical){
		time_mttr = data_per_disk * 2 / (percent_avail_bandwidth * TOTAL_BANDWIDTH);				// RS垂直散布，有多少个rack,就传多少个磁盘
	}


	heap_insert(time + time_mttr, DISK_REPAIRED, d.disk_no, d.node_no, d.rack_no);
} 
// C++引用
inline void add_disk_to_array(disk_array& darray, disk_struct d)
{
	darray.disks[darray.end] = d;
	darray.end++;
	darray.end %= MAX_DISKS;
}
inline void remove_disk_from_array(disk_array& darray, disk_struct d)
{
	int j = 0;
	for (int i = darray.front; i != darray.end; i++, i %= MAX_DISKS)
	{
		if (d == darray.disks[i])
		{
			++j;
		}
		else
		{
			darray.disks[(i + MAX_DISKS - j) % MAX_DISKS] = darray.disks[i];
		}
	}
	darray.end = (darray.end + MAX_DISKS - j) % MAX_DISKS;
}
// true for truly deleting
inline bool remove_event_from_heap(disk_struct d, int ev_type)
{
	int number = heap_search(ev_type, d.disk_no, d.node_no, d.rack_no);
	if (number > 0)
	{
		heap_delete(number);
		return true;
	}
	return false;
}

void initialize()
{
	events = 0;
	event_heap[0].time = 0;

	printf("the failed_disk.size is %d \n" , (failed_disks.end + MAX_DISKS - failed_disks.front ) % MAX_DISKS);
	printf("the failed_node.size is %d \n" , (failed_nodes.front + MAX_DISKS - failed_nodes.end ) % MAX_DISKS);
	if ((failed_nodes.end + MAX_DISKS - failed_nodes.front) % MAX_DISKS == 49999) {
		printf("the failed_node has error ... \n");
	}

	failed_disks.front = failed_disks.end = 0;
	warn_disks.front = warn_disks.end = 0;
	ld_disks.front = ld_disks.end = 0;
	failed_nodes.front = failed_nodes.end = 0;

	for (int i = 0; i < rack; i++)
		for (int k = 0; k < node; k++)
		{
			for (int j = 0; j < disk; j++)
			{
				events_for_new_disk({ j, k, i });
			}
			event_failure_for_node({-1, k, i});
		}
}
int test_mttdl;
double test_mttdl_cycle;

int sim_num_loss(double t)
{
	int type;
	int num_loss = 0; //数据丢失数量
	//   int test_mttdl_num = 1;
	long time;

	heap_delete(1);
	while (event_heap[0].time <= t)
	{
		type = event_heap[0].type;
		time = event_heap[0].time;
		disk_struct d ={event_heap[0].disk_no, event_heap[0].node_no, event_heap[0].rack_no};

		if (type == DISK_LD)
		{
			//记录当前潜在错误，并启动下次潜在错误生成
			TRACE_LD(d, time);
			add_disk_to_array(ld_disks, d);
			event_ld_for_disk(d, time);

			//if (func[coding](disk, node, rack))
			//{
			//	num_loss = num_loss + 1;
			//}
			if (is_data_loss_with_m_fails(d, type)) {
				num_loss = num_loss + 1;
				TRACE_DLOSS(d, time);
			}
		}
		else if (type == DISK_SCRUB)
		{
			//清除该磁盘的所有扇区错误
			TRACE_SCRUB(d, time);
			remove_disk_from_array(ld_disks, d);
			//磁盘清洗完成，清除所有已发生的潜在错误。启动下一次清洗。
			event_scrub_for_disk(d, time);
		}
		else if (type == NODE_FAIL)
		{
			// 清除该节点上的所有磁盘的相关事件，并加入节点故障标记和修复事件
			TRACE_NODE_FAIL(d, time);

			// 1. 遍历该node下的所有磁盘，逐个从array和heap中清除
			for (int i = 0; i < disk; ++i ) {
				disk_struct to_delete = {i, d.node_no, d.rack_no};
				remove_disk_from_array(warn_disks, to_delete);
				remove_disk_from_array(ld_disks, to_delete);
				remove_disk_from_array(failed_disks, to_delete);

				remove_event_from_heap(to_delete, DISK_LD);
				remove_event_from_heap(to_delete, DISK_SCRUB);
				remove_event_from_heap(to_delete, DISK_WARN);
				remove_event_from_heap(to_delete, DISK_REPAIRED);
			}

			// 2. 加入node故障标记，并启动node的下一次修复
			add_disk_to_array(failed_nodes, d);

			// 先 生成修复事件，插入到事件堆中，然后判断此时系统状态会不会丢失数据，即系统在这个盘、节点丢失了数据之后，仍会出现故障修复事件（即修复晚了）
			event_repair_for_node(d, time);

			//判断是否dataloss
			//if (func[coding](disk, node, rack))
			//{
			//	num_loss = num_loss + 1;
			//}
			if (is_data_loss_with_m_fails(d, type)) {
				num_loss = num_loss + 1;
				TRACE_DLOSS(d, time);
				printf("the num_loss is : %d in the past %ld hours...\n", num_loss, time);
			}
		}
		else if (type == NODE_REPAIRED)
		{
			// 删除该node下的全部事件
			TRACE_NODE_REPAIR(d, time);
			// 1. 遍历该node下的所有磁盘，逐个从array和heap中清除
			for (int i = 0; i < disk; ++i ) {
				disk_struct to_delete = {i, d.node_no, d.rack_no};
				remove_disk_from_array(warn_disks, to_delete);
				remove_disk_from_array(ld_disks, to_delete);
				remove_disk_from_array(failed_disks, to_delete);

				remove_event_from_heap(to_delete, DISK_LD);
				remove_event_from_heap(to_delete, DISK_SCRUB);
				remove_event_from_heap(to_delete, DISK_WARN);
			}
			// 2. 删除node故障事件
			remove_disk_from_array(failed_nodes, d);

			// 重新生成该node下的全部事件
			for (int i = 0; i < disk; ++i ) {
				disk_struct to_new = {i, d.node_no, d.rack_no};
				events_for_new_disk(to_new);
			}
			event_failure_for_node(d);

			// 记得归还带宽
			//if (coding == RShorizontal) {
			//	percent_avail_bandwidth += RShorizontal_RepairBandwidth;			
			//}
			//else if (coding == RSvertical){
			//	percent_avail_bandwidth += RSvertical_RepairBandwidth;	
			//}
		}
		else if (type == DISK_WARN)
		{ 
			//发生预警，生成修复事件
			TRACE_WARN(d, time);
			add_disk_to_array(warn_disks, d);
			event_repair_for_disk(d, time);
		}
		else if (type == DISK_REPAIRED)
		{
			// 磁盘修复完成，删除磁盘上的所有事件
			TRACE_REPAIR(d, time);
			remove_disk_from_array(failed_disks, d);
			remove_disk_from_array(warn_disks, d);
			remove_disk_from_array(ld_disks, d);

			// 删除堆上关于这个磁盘的所有事件 ?
			remove_event_from_heap(d, DISK_LD);
			remove_event_from_heap(d, DISK_SCRUB);
			remove_event_from_heap(d, DISK_FAIL);
			remove_event_from_heap(d, DISK_WARN);

			// 重新生成新的磁盘上的所有事件
			events_for_new_disk(d);

			// 记得归还带宽
			//if (coding == RShorizontal) {
			//	percent_avail_bandwidth += RShorizontal_RepairBandwidth;			
			//}
			//else if (coding == RSvertical){
			//	percent_avail_bandwidth += RSvertical_RepairBandwidth;	
			//}
		}
		else if (type == DISK_FAIL)
		{ 
			// 记录故障，并从array和heap中移除相关对象
			TRACE_FAIL(d, time);
			add_disk_to_array(failed_disks, d);
			// 因为当前的盘已经故障，所有的LD,SCRUB,WARN就没有存在的意义
			remove_disk_from_array(ld_disks, d);
			remove_disk_from_array(warn_disks, d);

			remove_event_from_heap(d, DISK_LD);
			remove_event_from_heap(d, DISK_SCRUB);
			remove_event_from_heap(d, DISK_WARN);

			// 启动当前盘的修复（这个与dataloss无关，不管什么情况下，都应该准备修复此磁盘）
			event_repair_for_disk(d, time);

			//判断是否dataloss -> 应该换成 is_data_loss
			//if (func[coding](disk, node, rack))
			//{
			//	num_loss = num_loss + 1;
			//}
			if (is_data_loss_with_m_fails(d, type)) {
				num_loss = num_loss + 1;
				TRACE_DLOSS(d, time);
				printf("the num_loss is : %d in the past %ld hours...\n", num_loss, time);
			}
		}
		else
		{
			printf("BUG\n");
			exit(-1);
		}
		heap_delete(1);
	}

	return num_loss;
}

int main(int argc, char *argv[])
{
	//combine(4, -1); combine(4, 0); combine(4, 1); combine(4, 2); combine(4, 3); combine(4, 4); combine(4, 5);
	//combine(56, 5);
	//testCombine();
	/*int result = 0;
	result = combinator(4, -1);
	printf("result is : %d\n", result);
	result = combinator(4, 0); 
	printf("result is : %d\n", result);
	result = combinator(4, 1); 
	printf("result is : %d\n", result);
	result = combinator(4, 2); 
	printf("result is : %d\n", result);
	result = combinator(4, 3); 
	printf("result is : %d\n", result);
	result = combinator(4, 4); 
	printf("result is : %d\n", result);
	result = combinator(4, 5);
	printf("result is : %d\n", result);*/



	
	double result[10000];
	int i, reps1, num_dloss, reps2, add_reps1 = 0, add_loss = 0;
	double pro_loss, t;

	if (argc < 19)
	{
		printf("Usage: mttdlsim coding (Twocopy/Threecopy/RS) disks nodes racks  Lyita_op Lbita_op Lgama_rest Lyita_rest Lbita_rest Lyita_ld, Lbita_ld, Lgama_scrub, Lyita_scrub, Lbita_scrub Lyita_tia Lbita_tia FDR t reps rs_k rs_m\n");
		return -1;
	}

	if (!strcmp(argv[1], "Twocopy"))
		coding = Twocopy;
	else if (!strcmp(argv[1], "Threecopy"))
		coding = Threecopy;
	else if (!strcmp(argv[1], "RShorizontal"))
		coding = RShorizontal;
	else if (!strcmp(argv[1], "RSvertical"))
		coding = RSvertical;
	else
	{
		printf("coding BUG\n");
		return -1;
	}
	//./sim_predict Threecopy 4 15 60 461386 1.1 6 12 2 9259 1 6 168 3 354 1 80 43800 10
	//  0.000897

	// new version
	//./sim_predict Threecopy 4 15 60 461386 1.1 6 12 2 9259 1 6 168 3 354 1 461386 1.1  6 12 2  80 43800 10
	//./sim_predict RS 4 15 60 461386 1.1 6 12 2 9259 1 6 168 3 354 1 461386 1.1  6 12 2  80 43800 10 6 3 
	//	  RShorizontal 4 15 60 461386 1.1 6 12 2 9259 1 6 168 3 354 1 461386 1.1  6 12 2  80 8760   1 6 3 100000 

	//	  RShorizontal 4 15 60 10000 1.5 6 12 2 9259  1 6 168 3 354 1 10000 1.1  6 12 2  0 4380  1 6 3 100000 
	//  80 预测准确率 -》 改为0 ？
	// 8760 改为 87600 ？


	//  0.000897										// 上一篇论文 可靠性公式 结果 看一下逻辑对不对，体现一种情况，跟原来差不多。。
	//block = Block;
	disk = atoi(argv[2]); // 4
	node = atoi(argv[3]); // 15
	rack = atoi(argv[4]); // 60							//
	Lyita_op = atof(argv[5]);	//461386				//MTTF? 小时h 越小越容易坏。。。。
	Lbita_op = atof(argv[6]);	//1.1
	Lgama_rest = atof(argv[7]);	//6
	Lyita_rest = atof(argv[8]);	//12
	Lbita_rest = atof(argv[9]);	//2
	Lyita_ld = atof(argv[10]);	//9259
	Lbita_ld = atof(argv[11]);	//1
	Lgama_scrub = atof(argv[12]);	//6
	Lyita_scrub = atof(argv[13]);	//168
	Lbita_scrub = atof(argv[14]);	//3
	Lyita_tia = atof(argv[15]); // 354 14.75天
	Lbita_tia = atof(argv[16]); // 1
	Lyita_node = atof(argv[17]); // node 461386 
	Lbita_node = atof(argv[18]); // 1.1
	Lgama_node_rest = atof(argv[19]) * disk; // node repair 6
	Lyita_node_rest = atof(argv[20]) * disk; //12
	Lbita_node_rest = atof(argv[21]); //2
	fdr = atof(argv[22]); // 80									// 改为0 。。。。如果为100%，存储系统根本不会发生数据丢失，都能预警出来。。？？
	t = atof(argv[23]); // 43800 5年 小时 1年 8760
	reps2 = atoi(argv[24]); // 10
	rs_k = atoi(argv[25]); // rs_k 6
	rs_m = atoi(argv[26]);		// rs_m 3
	block = atoi(argv[27]);		// 100000	// 新增
	fprintf(stderr, "block:%d disk:%d node:%d rack:%d fdr:%.0f coding:%s(%d,%d) t:%.1fyear reps:%d\n", 
		block, disk, node, rack, fdr, argv[1], rs_k, rs_m, t/24/365, reps2);
	srand((unsigned)time(NULL));

	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);

	gsl_rng_set(r, (unsigned long)time(NULL));

	for (i = 0; i < reps2; i++)
	{
		num_dloss = 0;
		reps1 = 0;
		while (num_dloss < 5)
		{
			initialize();
			num_dloss = num_dloss + sim_num_loss(t);
			reps1++; //出现数据丢失的循环次数
			if( (reps1+1)%10 == 0){
				printf("...repeat %d...\n",reps1);
				printf("...num_dloss %d...\n",num_dloss);
				//fflush(stdout);
			} 
		}
		add_loss = add_loss + num_dloss;
		add_reps1 = add_reps1 + reps1;
		result[i] = (double)num_dloss / (double)reps1;
		printf("...pro_dloss %f\n", result[i]);
	}

	pro_loss = gsl_stats_mean(result, 1, reps2);
	printf("\naverage = %f\n", pro_loss);
	printf("...add_loss %d....add_reps1 %d...pro_loss %lf \n", add_loss, add_reps1, 1.0 * add_loss / add_reps1);
	printf("\nvariance = %f\n", gsl_stats_variance_with_fixed_mean(result, 1, reps2, pro_loss));
	//    printf("\n mttdl = %f\n", test_mttdl/test_mttdl_cycle);

	//printf("\num_dataloss = %d\n", num_dloss);
	//printf("\resp = %d\n", reps1);
	gsl_rng_free(r);
	
	return 0;
}


/*****************************************************************
*            The Heap structure functions
*  heap_insert: 插入event
*  heap_search: 搜索event
*  heap_delete: 删除给定索引的event
****************************************************************/
inline void heap_insert(long time, int type, int disk_no, int node_no, int rack_no)
{
	int i, p;
	int j;
	if (type == DISK_LD || type == DISK_SCRUB)
	{
		return;
	}

	for (j = 1; j <= events; j++)
		if (event_heap[j].type == type && event_heap[j].disk_no == disk_no && event_heap[j].node_no == node_no && event_heap[j].rack_no == rack_no)
			return;

	events++;
	i = events;
	p = i / 2;
	while (p > 0 && event_heap[p].time > time)
	{
		event_heap[i] = event_heap[p];
		i = p;
		p = i / 2;
	}
	event_heap[i].time = time;
	event_heap[i].type = type;
	event_heap[i].disk_no = disk_no;
	event_heap[i].node_no = node_no;
	event_heap[i].rack_no = rack_no;
}

inline int heap_search(int type, int disk_no, int node_no, int rack_no)
{
	int i;
	for (i = 1; i <= events; i++)
	{
		if (event_heap[i].type == type && event_heap[i].disk_no == disk_no && event_heap[i].node_no == node_no && event_heap[i].rack_no == rack_no)
		{
			return i;
		}
	}
	return -1;
}

inline void __MinHeapIfy(int i)
{
	int min = i;
	int left = 2 * i;
	int right = left + 1;

	if (left <= events && event_heap[i].time > event_heap[left].time)
	{
		min = left;
	}
	else
	{
		min = i;
	}
	if (right <= events && event_heap[min].time > event_heap[right].time)
	{
		min = right;
	}
	if (i != min)
	{
		event tmp = event_heap[i];
		event_heap[i] = event_heap[min];
		event_heap[min] = tmp;
		__MinHeapIfy(min);
	}
}

//x节点的time减小到更小的值
inline void __HeapDecreaseKey(int x)
{
	event tmp;
	while (x > 1 && event_heap[x >> 1].time > event_heap[x].time)
	{
		tmp = event_heap[x];
		event_heap[x] = event_heap[x >> 1];
		event_heap[x >> 1] = tmp;
		x >>= 1;
	}
}

inline int __checkheap()
{
	int i, left, right;
	int end = events >> 1;
	for (i = 1; i <= end; ++i)
	{
		if (event_heap[i].time < 0)
			return 0;
		left = 2 * i;
		right = left + 1;
		if (left <= events && event_heap[i].time > event_heap[left].time)
			return 0;
		if (right <= events && event_heap[i].time > event_heap[right].time)
			return 0;
	}
	return 1;
}

//删除堆中指定元素i
inline void heap_delete(int i)
{

	event_heap[0] = event_heap[i];
	event tmp = event_heap[events];

	if (event_heap[i].time == tmp.time)
	{
		event_heap[i] = tmp;
		events--;
	}
	else if (event_heap[i].time < tmp.time) //i为根节点的分支的最小堆性质可能遭到破坏，要进行调整
	{
		event_heap[i] = tmp;
		events--;
		__MinHeapIfy(i);
	}
	else //i节点的值减小到更小的tmp
	{
		event_heap[i] = tmp;
		events--;
		__HeapDecreaseKey(i);
	}
}


// check whether data loss happened
// lp师兄写的
int Threecopy_failed(int b, int d, int n, int r)
{
	int i, j, k;
	int rack_no1, rack_no2, rack_no3;
	int node_no1, node_no2, node_no3;

	if (((failed_disks.end + MAX_DISKS - failed_disks.front) % MAX_DISKS) <= 2) {
		return 0;
	}

	for (i = failed_disks.front; i != (failed_disks.end + MAX_DISKS - 2) % MAX_DISKS; i++, i %= MAX_DISKS)
	{
		rack_no1 = failed_disks.disks[i].rack_no;
		node_no1 = failed_disks.disks[i].node_no;
		for (j = (i + 1) % MAX_DISKS; j != (failed_disks.end + MAX_DISKS - 1) % MAX_DISKS; j++, j %= MAX_DISKS)
		{
			rack_no2 = failed_disks.disks[j].rack_no;
			node_no2 = failed_disks.disks[j].node_no;
			for (k = (j + 1) % MAX_DISKS; k != failed_disks.end; k++, k %= MAX_DISKS)
			{
				rack_no3 = failed_disks.disks[k].rack_no;
				node_no3 = failed_disks.disks[k].node_no;
				if ((rack_no1 == rack_no2 && node_no1 != node_no2 && rack_no1 != rack_no3) 
					|| (rack_no1 == rack_no3 && node_no1 != node_no3 && rack_no1 != rack_no2) 
					|| (rack_no2 == rack_no3 && node_no2 != node_no3 && rack_no2 != rack_no1))
					return 1;
			}
		}
	}
	return 0;
}


//// check whether data loss happened Rs（）
//int RShorizontal_failed(int b, int d, int n, int r)	//RS水平
//{
//	int i;
//	int rack_no;
//	int temp_MAXDISKS[MAX_DISKS] = { 0 };
//	// 这里检查的是 故障磁盘 跨越的机架个数 会不会有问题？
//	if (((failed_disks.end + MAX_DISKS - failed_disks.front) % MAX_DISKS) <= rs_m) {
//		return 0;
//	}
//	int count_error = 0;
//	for (i = failed_disks.front; i != (failed_disks.end + MAX_DISKS) % MAX_DISKS; i++, i %= MAX_DISKS)
//	{
//		rack_no = failed_disks.disks[i].rack_no;
//		temp_MAXDISKS[rack_no]++;
//		if (temp_MAXDISKS[rack_no] > rs_m) {
//			return 1;
//		}
//	}
//	return 0;
//}
//
//
//// check whether data loss happened Rs（）
//int RSvertical_failed(int b, int d, int n, int r)	//RS 垂直
//{
//	int i, j, k, m;
//	int rack_no1, rack_no2, rack_no3, rack_no4;
//	int node_no1, node_no2, node_no3, node_no4;
//	int temp_MAXDISKS[MAX_DISKS] = { 0 };
//	if (((failed_disks.end + MAX_DISKS - failed_disks.front) % MAX_DISKS) <= rs_m) {
//		return 0;
//	}
//
//	//RS 
//	for (i = failed_disks.front; i != (failed_disks.end + MAX_DISKS - 3) % MAX_DISKS; i++, i %= MAX_DISKS)
//	{
//		rack_no1 = failed_disks.disks[i].rack_no;
//		node_no1 = failed_disks.disks[i].node_no;
//		for (j = (i + 1) % MAX_DISKS; j != (failed_disks.end + MAX_DISKS - 2) % MAX_DISKS; j++, j %= MAX_DISKS)
//		{
//			rack_no2 = failed_disks.disks[j].rack_no;
//			node_no2 = failed_disks.disks[j].node_no;
//			for (k = (j + 1) % MAX_DISKS; k != (failed_disks.end + MAX_DISKS - 1) % MAX_DISKS; k++, k %= MAX_DISKS)
//			{
//				rack_no3 = failed_disks.disks[k].rack_no;
//				node_no3 = failed_disks.disks[k].node_no;
//				for (m = (k + 1) % MAX_DISKS; m != failed_disks.end; m++, m %= MAX_DISKS) {
//					rack_no4 = failed_disks.disks[m].rack_no;
//					node_no4 = failed_disks.disks[m].node_no;
//					//这个是不是应该用failed_node数组，可是看师兄判断三副本时用的是 failed_disk数组？这个是不是有问题？
//					//共分为3种情况 1个机架上3个节点发生了故障，另外1个机架1个节点发生了故障；两个机架上2个节点发生了故障；1个机架的2个节点发生了故障，另外两个机架的1个节点发生故障
//					if ((rack_no1 == rack_no2 && rack_no1 == rack_no3  &&  node_no1 != node_no2 && node_no1 != node_no3 && node_no3 != node_no2 && rack_no1 != rack_no4) ||		// 1==2==3 4
//						(rack_no1 == rack_no2 && rack_no1 == rack_no4  &&  node_no1 != node_no2 && node_no1 != node_no4 && node_no4 != node_no2 && rack_no1 != rack_no3) ||		// 1==2==4 3
//						(rack_no1 == rack_no3 && rack_no1 == rack_no4  &&  node_no1 != node_no3 && node_no1 != node_no4 && node_no3 != node_no4 && rack_no1 != rack_no2) ||		// 1==3==4 2
//						(rack_no2 == rack_no3 && rack_no3 == rack_no4  &&  node_no2 != node_no3 && node_no3 != node_no4 && node_no2 != node_no4 && rack_no2 != rack_no1))		// 2==3==4 1
//						return 1;
//
//					// 两个机架上各自2个节点发生了故障
//					if ((rack_no1 == rack_no2 && node_no1 != node_no2 && rack_no3 == rack_no4 && node_no3 != node_no4) ||	// 1==2 3==4
//						(rack_no1 == rack_no3 && node_no1 != node_no3 && rack_no2 == rack_no4 && node_no2 != node_no4) ||	// 1==3 2==4
//						(rack_no1 == rack_no4 && node_no1 != node_no4 && rack_no3 == rack_no2 && node_no3 != node_no2))		// 1==4 2==3
//						return 2;
//
//					// 一个机架的2个节点放生了故障，另外两个机架上的1个节点发生了故障
//					if ((rack_no1 == rack_no2 && node_no1 != node_no2 && rack_no1 != rack_no3 && rack_no1 != rack_no4 && rack_no3 != rack_no4) ||	// 1==2 3 4
//						(rack_no1 == rack_no3 && node_no1 != node_no3 && rack_no1 != rack_no2 && rack_no1 != rack_no4 && rack_no2 != rack_no4) ||	// 1==3 2 4
//						(rack_no1 == rack_no4 && node_no1 != node_no4 && rack_no1 != rack_no2 && rack_no1 != rack_no3 && rack_no2 != rack_no3) ||	// 1==4 2 3
//						(rack_no2 == rack_no3 && node_no2 != node_no3 && rack_no2 != rack_no1 && rack_no2 != rack_no4 && rack_no1 != rack_no4) ||	// 2==3 1 4
//						(rack_no2 == rack_no4 && node_no2 != node_no4 && rack_no2 != rack_no1 && rack_no2 != rack_no3 && rack_no1 != rack_no3) ||	// 2==4 1 3
//						(rack_no3 == rack_no4 && node_no3 != node_no4 && rack_no3 != rack_no1 && rack_no3 != rack_no2 && rack_no1 != rack_no2)		// 3==4 1 2
//						)
//						return 3;
//					return 0;
//				}
//			}
//		}
//	}
//
//	return 0;
//}

int is_cross_multiple_rack(int coding, const vector<disk_struct>& failed_d_num, 
	const vector<disk_struct>& failed_n_num, disk_struct& disk_target, int rack_size)
{
	set<int> s_rack;
	s_rack.insert(disk_target.rack_no);
	if (failed_d_num.size() > 0) {
		for (auto it : failed_d_num) {
			s_rack.insert(it.rack_no);
			if (s_rack.size() > rack_size) {
				return 1;
			}
		}
	}
	if (failed_n_num.size() > 0) {
		for (auto it : failed_n_num) {
			s_rack.insert(it.rack_no);
			if (s_rack.size() > rack_size) {
				return 1;
			}
		}
	}

	return 0;
}


int which_situation_is(const vector<disk_struct>& failed_d_num,
	const vector<disk_struct>& failed_n_num, disk_struct& disk_target)
{
	int node_size = failed_n_num.size();
	int rack_no1, rack_no2, rack_no3, rack_no4;
	int node_no1, node_no2, node_no3, node_no4;
	for (int i = 0; i < node_size - 3; i++) {
		rack_no1 = failed_n_num[i].rack_no;
		node_no1 = failed_n_num[i].node_no;
		for (int j = i + 1; j < node_size - 2; j++) {
			rack_no2 = failed_n_num[j].rack_no;
			node_no2 = failed_n_num[j].node_no;
			for (int k = j + 1; k < node_size - 1; k++) {
				rack_no3 = failed_n_num[k].rack_no;
				node_no3 = failed_n_num[k].node_no;
				for (int m = k + 1; m < node_size ; m++) {
					rack_no4 = failed_n_num[m].rack_no;
					node_no4 = failed_n_num[m].node_no;
					//共分为3种情况 1个机架上3个节点发生了故障，另外1个机架1个节点发生了故障；两个机架上2个节点发生了故障；1个机架的2个节点发生了故障，另外两个机架的1个节点发生故障
					if ((rack_no1 == rack_no2 && rack_no1 == rack_no3  &&  node_no1 != node_no2 && node_no1 != node_no3 && node_no3 != node_no2 && rack_no1 != rack_no4) ||		// 1==2==3 4
						(rack_no1 == rack_no2 && rack_no1 == rack_no4  &&  node_no1 != node_no2 && node_no1 != node_no4 && node_no4 != node_no2 && rack_no1 != rack_no3) ||		// 1==2==4 3
						(rack_no1 == rack_no3 && rack_no1 == rack_no4  &&  node_no1 != node_no3 && node_no1 != node_no4 && node_no3 != node_no4 && rack_no1 != rack_no2) ||		// 1==3==4 2
						(rack_no2 == rack_no3 && rack_no3 == rack_no4  &&  node_no2 != node_no3 && node_no3 != node_no4 && node_no2 != node_no4 && rack_no2 != rack_no1))		// 2==3==4 1
						return 1;

					// 两个机架上各自2个节点发生了故障
					if ((rack_no1 == rack_no2 && node_no1 != node_no2 && rack_no3 == rack_no4 && node_no3 != node_no4) ||	// 1==2 3==4
						(rack_no1 == rack_no3 && node_no1 != node_no3 && rack_no2 == rack_no4 && node_no2 != node_no4) ||	// 1==3 2==4
						(rack_no1 == rack_no4 && node_no1 != node_no4 && rack_no3 == rack_no2 && node_no3 != node_no2))		// 1==4 2==3
						return 2;

					// 一个机架的2个节点放生了故障，另外两个机架上的1个节点发生了故障
					if ((rack_no1 == rack_no2 && node_no1 != node_no2 && rack_no1 != rack_no3 && rack_no1 != rack_no4 && rack_no3 != rack_no4) ||	// 1==2 3 4
						(rack_no1 == rack_no3 && node_no1 != node_no3 && rack_no1 != rack_no2 && rack_no1 != rack_no4 && rack_no2 != rack_no4) ||	// 1==3 2 4
						(rack_no1 == rack_no4 && node_no1 != node_no4 && rack_no1 != rack_no2 && rack_no1 != rack_no3 && rack_no2 != rack_no3) ||	// 1==4 2 3
						(rack_no2 == rack_no3 && node_no2 != node_no3 && rack_no2 != rack_no1 && rack_no2 != rack_no4 && rack_no1 != rack_no4) ||	// 2==3 1 4
						(rack_no2 == rack_no4 && node_no2 != node_no4 && rack_no2 != rack_no1 && rack_no2 != rack_no3 && rack_no1 != rack_no3) ||	// 2==4 1 3
						(rack_no3 == rack_no4 && node_no3 != node_no4 && rack_no3 != rack_no1 && rack_no3 != rack_no2 && rack_no1 != rack_no2)		// 3==4 1 2
						)
						return 3;
				}
			}
		}
	}

	return 0;
}


inline void testCombine() {
	disk_struct d ={0, 0, 0};
	add_disk_to_array(failed_disks, d);
	add_disk_to_array(failed_disks, d);
	vector<vector<disk_struct>> disks_vec = combine_darray(failed_disks, 1);
	for (auto it : disks_vec) {
		for (auto i : it) printf("%d ", i);
		printf("\n");
	}
}
