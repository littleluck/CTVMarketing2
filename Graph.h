#include<vector>
#include "Vertex.h"

using namespace std;


struct Graph
{
	int n ;
	vector<Vertex> vertex ;	//��������
	vector<vector<double> > adj ;	//�����ڽӱ�洢ͼ���ȳ�ʼ����һά���ڶ�ά��push_back
	vector<vector<int> > tag ;	//��Ϊ��<a,b> push_back��adj֮�󣬱ߵ��׽����a��vector��Ų���Ӧ
	//������tag���洢ʵ�ʱ��a

	//vector<double *> adj ;	//�ڽӾ���洢ͼ������Ҳ���Ա�ʾ��ά����

};