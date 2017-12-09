#include<iostream>
#include<vector>
#include <omp.h> 
#include<time.h>
#include<cstdlib>
#include<list>
#include<iterator>
#include<omp.h>
#include<cstdio>
#define rep(i,a,b) for(i=a;i<b;++i)
#define o(a) puts(a)
using namespace std;
const int cache_size=64/sizeof(int);
//vertex:vector
//
//set randval for each of the vertex two arrays one contains existing vertices id,another store random value;
//How could 
//for each edge:two ends
void graph_construction(int n,vector<pair<int,int> > &edges,double p)
{
	int i,j;
	rep(i,0,n)
		rep(j,0,i)
		{
			if(rand()*1.0/RAND_MAX<p)
			{
				edges.push_back(make_pair(i,j));
			}
		}
}
int main()
{
	cout<<cache_size<<endl;
	FILE* ftime=fopen("time4.txt","w");
	const int n=30000;
	double poss=0.2;
	vector<pair<int,int> > ori_edges;
	vector<int> independent_set;
	graph_construction(n,ori_edges,poss);
	int round;
	int round_num=20;
	rep(round,0,round_num){
		vector<int> existing_vertices(n,0);
		int randval[n]={0};
		int label[n]={0};
		vector<pair<int,int> > edges(ori_edges);
		int i;
		rep(i,0,n)
			existing_vertices[i]=i;
		vector<double> delta_time(5,0.0);
		while(!existing_vertices.empty())
		{
			//set randval
			double time=omp_get_wtime();
			int i;
			#pragma omp parallel for private(i) schedule(static,cache_size*2)
			rep(i,0,existing_vertices.size())
			{
				randval[existing_vertices[i]]=rand();
			}
			delta_time[0]+=omp_get_wtime()-time;
			//cout<<"rand val assigned"<<endl;
			//traverse the edge list,-1 remaining,1 represent neighbor 0 represent choose
			#pragma omp parallel for private(i) schedule(static,cache_size)
			rep(i,0,edges.size())
			{
				int start=edges[i].first;
				int end=edges[i].second;
				if(label[start]==-1 && label[end]==-1)
					continue;
				if(randval[start]<=randval[end])
				{
					label[end]=-1;
				}
				else
					label[start]=-1;
			}
			delta_time[1]+=omp_get_wtime()-time;
			//mark the neighbors
			#pragma omp parallel for private(i)
			rep(i,0,edges.size())
			{
				int start=edges[i].first;
				int end=edges[i].second;
				//if(label[start]==0 && label[end]==0)
					//cout<<"error"<<endl;
				if(label[start]==0)
				{
					label[end]=1;
				}
				else if(label[end]==0)
					label[start]=1;			
			}
			delta_time[2]+=omp_get_wtime()-time;
			//rep(i,0,n)
			//	cout<<label[i]<<endl;
			//getchar();
			//add to independent set and make new existing_vertices and delete edges
			int j=0;
			int* new_existing_vertices=new int[existing_vertices.size()];
			#pragma omp parallel for private(i)
			rep(i,0,existing_vertices.size())
			{
				if(label[existing_vertices[i]]==0)
				{
					#pragma omp critical
					independent_set.push_back(existing_vertices[i]);
				}
				if(label[existing_vertices[i]]<0 && i!=j)
				{
					#pragma omp critical
					new_existing_vertices[j++]=existing_vertices[i];
				}
			}
			cout<<j<<endl;
			//getchar();
			#pragma omp parallel for private(i)
			rep(i,0,j)
				existing_vertices[i]=new_existing_vertices[i];
			existing_vertices.resize(j);
			delete [] new_existing_vertices;
			delta_time[3]+=omp_get_wtime()-time;
			j=0;
			std::vector<pair<int,int> > new_edges;
			#pragma omp parallel for private(i)
			rep(i,0,edges.size())
			{
				if(label[edges[i].first]<0 && label[edges[i].second]<0)
				{
					#pragma omp critical
					new_edges.push_back(edges[i]);				
				}
			}
			j=new_edges.size();
			#pragma omp parallel for private(i)
			rep(i,0,j)
				edges[i]=new_edges[i];
			edges.resize(j);
			#pragma omp parallel for private(i)
			rep(i,0,existing_vertices.size())
			{
				#pragma omp critical
				label[existing_vertices[i]]=0;
			}
			delta_time[4]+=omp_get_wtime()-time;
		}
			fprintf(ftime,"%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",n,delta_time[0],delta_time[1],delta_time[2],delta_time[3],delta_time[4]);
			o("independent_set");
			rep(i,0,independent_set.size())
					cout<<independent_set[i]<<endl;
			cout<<independent_set.size()<<endl;
			independent_set.clear();
	}
	//set randval
}
