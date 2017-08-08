#include<vector>
#include "Vertex.h"

using namespace std;


struct Graph
{
	int n ;
	vector<Vertex> vertex ;	//顶点数组
	vector<vector<double> > adj ;	//反向邻接表存储图，先初始化第一维，第二维用push_back
	vector<vector<int> > tag ;	//因为边<a,b> push_back到adj之后，边的首结点编号a与vector编号不对应
	//所以用tag来存储实际编号a

	//vector<double *> adj ;	//邻接矩阵存储图，这样也可以表示二维数组

};