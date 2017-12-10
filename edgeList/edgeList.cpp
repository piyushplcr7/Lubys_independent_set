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
using std::vector;
using std::cout;
using std::pair;
using std::endl;
using std::make_pair;
const int n=10000;
const int cache_size=64/sizeof(int);
//vertex:vector
//
//set randval for each of the vertex two arrays one contains existing vertices id,another store random value;
//How could 
//for each edge:two ends
void ER_graph_construction(int n,vector<pair<int,int> > &edges,double p)
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
void rmat_graph_construction(char* filename,vector<pair<int,int> > &edges)
{
	FILE* fgraph=fopen(filename,"r");
	int start,end;
	while(fscanf(fgraph,"%d%d",&start,&end)!=EOF)
	{
		//cout<<start<<end<<endl;
		//getchar();
		if(start>=n || end>=n)
			exit(-1);
		edges.push_back(make_pair(start,end));
	}
	printf("construction end\n");
	fclose(fgraph);
}
int main(int argc,char* argv[])
{
	cout<<cache_size<<endl;
	cout<<argv[1]<<endl;
	FILE* ftime=fopen(argv[1],"w");
	double poss=0.2;
	vector<pair<int,int> > ori_edges;
	vector<int> independent_set;
	char filename[]="out.txt";
	ER_graph_construction(n,ori_edges,poss);
	//rmat_graph_construction(filename,ori_edges);
	int round;
	int round_num=20;
	srand(time(NULL));
	vector<double> total_delta_time(5,0.0);
	rep(round,0,round_num){
		vector<int> existing_vertices(n,0);
		int randval[n][cache_size]={0};
		int label[n][cache_size]={0};
		vector<pair<int,int> > edges(ori_edges);
		int i;
		rep(i,0,n)
			existing_vertices[i]=i;
		vector<double> delta_time(6,0.0);
		while(!existing_vertices.empty())
		{
			//set randval
			unsigned seed;
			int i;			
			double curr_time=omp_get_wtime();
			#pragma omp parallel private(seed)
			{
				seed=time(NULL)+19*omp_get_thread_num();
				#pragma omp for private(i)
				rep(i,0,existing_vertices.size())
				{
					randval[existing_vertices[i]][0]=rand_r(&seed);
				}
			}
			delta_time[0]+=omp_get_wtime()-curr_time;
			//cout<<"rand val assigned"<<endl;
			//traverse the edge list,-1 remaining,1 represent neighbor 0 represent choose
			curr_time=omp_get_wtime();
			#pragma omp parallel
			{
				//int chunk=edges.size()/omp_get_num_threads();
				//int end=(omp_get_thread_num()==omp_get_num_threads()-1)?edges.size():chunk*(omp_get_thread_num()+1);
				#pragma omp for
				rep(i,0,edges.size())
				{
					int start=edges[i].first;
					int end=edges[i].second;
					if(label[start][0]==-1 && label[end][0]==-1)
						continue;
					if(randval[start][0]<=randval[end][0])
					{
						label[end][0]=-1;
					}
					else
						label[start][0]=-1;
				}
			}
			delta_time[1]+=omp_get_wtime()-curr_time;
			//mark the neighbors
			curr_time=omp_get_wtime();
			#pragma omp parallel for private(i)
			rep(i,0,edges.size())
			{
				int start=edges[i].first;
				int end=edges[i].second;
				//if(label[start]==0 && label[end]==0)
					//cout<<"error"<<endl;
				if(label[start][0]==0)
				{
					label[end][0]=1;
				}
				else if(label[end][0]==0)
					label[start][0]=1;			
			}
			delta_time[2]+=omp_get_wtime()-curr_time;
			//add to independent set and make new existing_vertices and delete edges

			//delete vertices
			curr_time=omp_get_wtime();
			int j=0;
			//int* new_existing_vertices=new int[existing_vertices.size()];
			vector<int> new_existing_vertices;
			/*#pragma omp parallel for private(i)
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
			cout<<j<<endl;*/
			#pragma omp parallel
			{
				vector<int> private_vertices;
				vector<int> private_chosen;
				int chunk=existing_vertices.size()/omp_get_num_threads();
				int i;
				int end=(omp_get_thread_num()==omp_get_num_threads()-1)?existing_vertices.size():chunk*(omp_get_thread_num()+1);
				rep(i,chunk*omp_get_thread_num(),end)
				{
					if(label[existing_vertices[i]][0]==0)
					{
						private_chosen.push_back(existing_vertices[i]);
					}
					if(label[existing_vertices[i]][0]<0)
					{
						private_vertices.push_back(existing_vertices[i]);
					}
				}
				#pragma omp critical
				{
					new_existing_vertices.insert(new_existing_vertices.end(),private_vertices.begin(),private_vertices.end());
					independent_set.insert(independent_set.end(),private_chosen.begin(),private_chosen.end());
				}
			} 
			j=new_existing_vertices.size();
			cout<<j<<endl;
			//getchar();
			#pragma omp parallel for private(i)
			rep(i,0,j)
				existing_vertices[i]=new_existing_vertices[i];
			existing_vertices.resize(j);
			//delete [] new_existing_vertices;
			delta_time[3]+=omp_get_wtime()-curr_time;


			//delete edges
			j=0;
			std::vector<pair<int,int> > new_edges;
			curr_time=omp_get_wtime();
			#pragma omp parallel
			{
				vector<pair<int,int> > private_edges;
				int chunk=edges.size()/omp_get_num_threads();
				int i;
				int end=(omp_get_thread_num()==omp_get_num_threads()-1)?edges.size():chunk*(omp_get_thread_num()+1);
				//#pragma omp for private(i) nowait 
				rep(i,chunk*omp_get_thread_num(),end)
				{
					if(label[edges[i].first][0]<0 && label[edges[i].second][0]<0)
					{
						//#pragma omp critical
						private_edges.push_back(edges[i]);				
					}
				}
				#pragma omp critical
					new_edges.insert(new_edges.end(),private_edges.begin(),private_edges.end());
			}
			delta_time[4]+=omp_get_wtime()-curr_time;
			curr_time=omp_get_wtime();
			j=new_edges.size();
			//getchar();
			#pragma omp parallel for private(i)
			rep(i,0,j)
				edges[i]=new_edges[i];
			edges.resize(j);
			#pragma omp parallel for private(i) schedule(static,cache_size)
			rep(i,0,existing_vertices.size())
			{
				label[existing_vertices[i]][0]=0;
			}
			delta_time[5]+=omp_get_wtime()-curr_time;
		}
			fprintf(ftime,"%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",n,delta_time[0],delta_time[1],delta_time[2],delta_time[3],delta_time[4],delta_time[5]);
			o("independent_set");

			//rep(i,0,independent_set.size())
			//		cout<<independent_set[i]<<endl;
			cout<<independent_set.size()<<endl;
			//getchar();
			independent_set.clear();
	}
	fclose(ftime);
	//set randval
}
