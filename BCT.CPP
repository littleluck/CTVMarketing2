#include <iostream>
#include <cmath>
#include <algorithm>

#include <stdlib.h>
#include <time.h> 

//读取文件所需要
#include <fstream>
#include <cassert>
#include <string>

#include "Graph.h"

//BCT算法实现影响最大化
using namespace std;

const double C = 2*(exp(1) - 2) ;
const double coefficient = 8*C*(1 - 1/(2*exp(1)))*(1 - 1/(2*exp(1))) ;

	
/**
 *  log(C(n,k)) = log(n!) - log((n-k)!) - log(k!)
 * @return C(n,k)组合数的对数
 */
 double logCombination(int n, int k){
	cout<<"logCombination start"<<endl ;
	if(n < k){
		cout<<"n要大于等于k"<<endl ;
		exit(1) ;
	}
	double sum = 0 ;
	for(int i=1;i<=k; i++){
		sum = sum + log(n-k+i) - log(i) ;
	}
	return sum ;
}


/**
 * |S|=k, 且c(S)<=budget，返回最大的k
 * @param cost
 * @param budget
 * @return
 */
int getKMAX(double cost[], int length, double budget){
	cout<<"getKMAX start"<<endl ;
	sort(cost, cost+length) ;
	double sum = 0 ;
	int i ;
	for(i=0; i<length; i++){
		sum += cost[i] ;
		if(sum > budget) break ;
	}
	return i ;
}

/**
 * Algorithm 2 BSA - Benefit Sampling Algorithm for LT model
 * @param graph
 * @return 采样的节点集合
 */
void benefitSampling(vector<int> &seed, Graph graph){
	//cout<<"benefitSampling start"<<endl ;
	seed.clear() ;
	int n = graph.n ;
	double sumBenefits = 0 ;
	//Vertex vertex[] = graph.vertex ;
	//double adj[][] = graph.adj ;
	int i ;
	for(i=0; i<n; i++){	//得到benefit总和
		sumBenefits += graph.vertex[i].benefit ;
	}
	
	//注意：此句代码不能写在这里，因为每次进入此方法调用srand()的时间间隔极为短暂，
	//所以会导致每次生成的随机数p都一样
	//srand(time(NULL));

	double p = rand()/double(RAND_MAX) ;	//随机生成一个0~1之间的小数
	//cout<<"p: "<<p<<endl ;
	double sum1 = 0 ;
	int u = n-1 ;	//将u初始化为结点n-1，因为当p=1时，下面的if语句不能执行，不知道为啥
	for(i=0; i<n; i++){
		sum1 += graph.vertex[i].benefit/sumBenefits ;
		if(p <= sum1){	//选中了一个节点
			u = i ;
			break ;
		}
	}
	//cout<<"p="<<p<<" sum1="<<sum1<<" u="<<u<<endl ;
	
	//只要节点u不在seed中就循环
	vector<int>::iterator it;
	it = find(seed.begin(), seed.end(), u) ;
	while(it == seed.end()){	//在seed中没找到u
		seed.push_back(u) ;	//添加结点u到seed中
		
		//Attempt to select an edge (v, u) using live-edge model
		double p2 = rand()/double(RAND_MAX) ;
		int v = -1 ;
		double sum2 = 0 ;
		int size = graph.adj[u].size() ;
		for(i=0; i<size; i++){	//这里要遍历所有指向u的结点
			//if(graph.adj[u][i] != 0){	//存在边(i,u)
				sum2 += graph.adj[u][i] ;
				if(p2 <= sum2){	//选中了一条边
					v = graph.tag[u][i] ;
					break ;
				}
			//}
		}
		if(v != -1){	//选中边(v,u)
			u = v ;
		}

		it = find(seed.begin(), seed.end(), u) ;	//因为u可能已更新，所以要重新定位
	}
	//cout<<"benefitSampling end"<<endl ;
}

int degHS(vector<vector<int> > hypergraph, vector<int> s) ;

/**
 * Algorithm 3 Weighted-Max-Coverage Algorithm
 * @param n
 * @param hypergraph
 * @param B
 * @return 临时种子集合
 */
void weightedMaxCoverage(int n, vector<int> &seed, Graph graph, vector<vector<int> > hypergraph, double B)
{	
	cout<<"weightedMaxCoverage start"<<endl ;
	seed.clear() ;	//将seed集合置空
	bool flag = true ;
	int i ;
	int currentPosition = 0 ;
	double sumCost = 0 ;	//保存seed集合中所有节点的cost总和
	vector<int>::iterator it;
	while(flag){
		flag = false ;	
		
		for(i=currentPosition; i<n; i++){	//判断是否还有节点满足条件，若有则将flag设置为true，此处每次都从i=currentPosition开始遍历
			it = find(seed.begin(), seed.end(), i) ;
			if((it == seed.end()) && (graph.vertex[i].cost<=(B-sumCost))){
				flag = true ;
				currentPosition = i ;
				break ;
			}
		}
		if(flag){
			double max = -1 ;
			int maxNumber = 0 ;
			int deg = degHS(hypergraph, seed) ;
			for(i=currentPosition; i<n; i++){	//找到具有最大(degH(S∪{v})-degH(S))/c(v)的节点
				it = find(seed.begin(), seed.end(), i) ;	//此处仍要保证结点i不在seed中，不然同一个节点可能会被加入到seed中多次
				if((it == seed.end()) && (graph.vertex[i].cost<=(B-sumCost))){
					seed.push_back(i) ;
					int deg2 = degHS(hypergraph, seed) ;
					double temp = (deg2 - deg)/graph.vertex[i].cost ;
					if(max < temp){
						max = temp ;
						maxNumber = i ;
					}
					seed.pop_back() ;	//记得将结点i删除
				}
			}//end for
			seed.push_back(maxNumber) ;	//将具有最大(degH(S∪{v})-degH(S))/c(v)的节点加入到seed中
			sumCost += graph.vertex[maxNumber].cost ;	//要更新sumCost
			cout<<"在weightedMaxCoverage中每次选择的结点u："<<maxNumber<<endl ;
		}
	}//end while
	
	vector<int> v ;
	int u = -1 ;	//保存最大degH(v)的结点编号
	int maxDeg = -1 ;
	for(i=0; i<n; i++){
		if(graph.vertex[i].cost <= B){
			v.push_back(i) ;
			int temp = degHS(hypergraph, v) ;
			if(maxDeg < temp){
				maxDeg = temp ;
				u = i ;
			}
			v.pop_back() ;
		}
	}
	if(degHS(hypergraph, seed) < maxDeg){
		seed.clear() ;
		seed.push_back(u) ;
	}
	cout<<"weightedMaxCoverage end"<<endl ;
}
	
/**
 * 求degH(S)，应该是求集合s与hypergraph中的vector有交集的个数
 * @param hypergraph
 * @param s
 * @return	
 */
int degHS(vector<vector<int> > hypergraph, vector<int> s){
	int intersections = 0 ;	//保存交集个数
	int i, j ;
	int m = hypergraph.size() ;
	int n ;
	vector<int>::iterator it ;
	for(i=0; i<m; i++)
	{
		n = hypergraph[i].size() ;
		for(j=0; j<n; j++)
		{
			it = find(s.begin(), s.end(), hypergraph[i][j]) ;
			if(it != s.end())	//hypergraph[i][j]在s中
			{
				intersections++ ;
				break ;
			}
		}
	}
	
	return intersections ;
}


void mainBCT(vector<int> &s, Graph graph, double budget, double epi, double q){
	cout<<"mainBCT start"<<endl ;
	if(budget<=0 || (epi<=0 || epi>=1)||(q<=0 || q>=1)) exit(1) ;	//程序异常退出
	
	int n = graph.n ;
	double *cost = new double[n] ;
	for(int i=0; i<n; i++){
		cost[i] = graph.vertex[i].cost ;
	}
	int kMax = getKMAX(cost, n, budget) ;

	delete[] cost ;	//释放cost内存
	
	double L = coefficient*(-log(q) + logCombination(n, kMax) + 2.0/n)/(epi*epi) ;	//for uniform cost
//	double L = coefficient*(-Math.log(q) + kMax*Math.log(n) + 2.0/n)/(epi*epi) ;	//otherwise
	
	double AL = (1 + exp(1)*epi/(2*exp(1) - 1))*L ;
	double Nt = AL ;
	vector<vector<int> > hypergraph ;
	
	do{
		for(int j=1; j<=Nt-hypergraph.size(); j++){
			vector<int> set ;
			benefitSampling(set, graph) ;
			hypergraph.push_back(set) ;
		}
		Nt = 2*Nt ;
		weightedMaxCoverage(n, s, graph, hypergraph, budget) ;
	}while(degHS(hypergraph, s) < AL) ;
	cout<<"mainBCT end"<<endl ;
}

/**
分割字符串，s为主串，c为分隔符
*/
void splitString(const string& s, vector<string>& v, const string& c)
{
	//cout<<"splitString start"<<endl ;
	v.clear() ;
	string::size_type pos1, pos2;
	pos2 = s.find(c);
	pos1 = 0;
	while(string::npos != pos2)
	{
		v.push_back(s.substr(pos1, pos2-pos1));
	
		pos1 = pos2 + c.size();
		pos2 = s.find(c, pos1);
	 }
	if(pos1 != s.length())
		v.push_back(s.substr(pos1));
}


int main()
{
	//记得先把graph中的vertex和adj都初始化为0
	//节点下标要从0开始
	struct Graph graph ;

	string file = "./graphfile.txt" ;
	ifstream infile ;
	infile.open(file.data()) ;	//将文件流对象与文件连接起来 
	assert(infile.is_open()) ;	//若失败，则输出错误信息，并终止程序

	string str ;
	getline(infile, str) ;	//第一行为结点总数，先读取出来
	cout<<str<<endl ;
	vector<string> result ;	//保存分割后的字符串
	//splitString(str, result, "	") ;	//这里以Tab键为分隔符
	graph.n = atoi(str.c_str()) ;	//只要第一个

	//下面初始化adj、tag和vertex
	int n = graph.n ;
	graph.adj.resize(n) ;
	graph.tag.resize(n) ;
	int i ;
	//for(i=0; i<n; i++)
	//{
		//graph.adj[i].resize(n) ;
	//}
	cout<<"----------"<<endl ;

	graph.vertex.resize(n) ;
	for(i=0; i<n; i++)
	{
		graph.vertex[i].benefit = 1 ;
		graph.vertex[i].cost = 1 ;
	}

	int a, b ;
	double weight ;
	while(getline(infile, str))
	{
		//cout<<"----------"<<endl ;
		splitString(str, result, " ") ;	//这里以Tab键为分隔符
		a = atoi(result[0].c_str()) ;
		b = atoi(result[1].c_str()) ;
		weight = atof(result[2].c_str()) ;
		graph.adj[b].push_back(weight) ;	//反向存储边，graph.adj[b]存储的都是指向b的结点
		graph.tag[b].push_back(a) ;	//将实际的首节点编号保存
		//graph.adj[a][b] = weight ;
	}
	infile.close() ;


	clock_t start, finish ;
	double duration ;
	start = clock() ;

	srand(time(NULL));	//根据当前时间初始化随机种子，一定要在main()方法中初始化
	vector<int> s ;
	mainBCT(s, graph, 3, 0.1, 1.0/n) ;

	finish = clock() ;
	duration = (double)(finish - start)/CLOCKS_PER_SEC ;
	
	cout<<"用时"<<duration<<"s"<<endl ;
	int size = s.size() ;
	cout<<"种子集合大小为："<<size<<endl ;
	cout<<"种子集合为：" ;
	for(i=0; i<size; i++)
	{
		cout<<s[i]<<" " ;
	}
	cout<<endl ;


	//cout<<"read end!"<<endl ;
	return 0 ;
}