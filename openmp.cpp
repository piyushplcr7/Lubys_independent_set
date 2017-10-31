#include<iostream>
#include<vector>
#include <omp.h> 
#include<time.h>
#include<cstdlib>
#define rep(i,a,b,c) for(int i=a;i<b;i+=c)
using std::vector;
struct Graph {
	public:
	int n_vertices;
	//Set of Vertices; contains 1 if vertex is in the graph and 0 if not in the graph
	vector<int> vertices;
	//std::array<int,n_vertices> vertices;
	//set of Edges
	//std::array<std:array<int,n_vertices>,n_vertices> edges;
	vector<vector<int> > edges;
	//constructor
	Graph(int n):n_vertices(n),vertices(n,0),edges(n,vector<int>(n,0))
	{
			//create the matrix
		rep(i,0,n,1)
			vertices[i]=i;
	}
	Graph(int n,double prob):n_vertices(n),vertices(n,0),edges(n,vector<int>(n,0))//random graph G(n_vertices,prob) where prob is the probability that an edge exists
	{
		rep(i,0,n,1)
			vertices[i]=i;
		srand(time(NULL));
		rep(i,0,n,1)
		{
			rep(j,0,n,1)
				edges[i][j]=(rand()*1.0/RAND_MAX<prob);
			edges[i][i]=0;
		}
	}
		//create the matrix
		/*for (int i=0;i<n_vertices;++i)
			{
			edges[i]=new int[n_vertices];
			vertices[i]=1;
			}
		};*/
	//destructor
	/*~Graph(){
		delete [] vertices;
		for (int i = 0; i < n_vertices; ++i)
    		delete [] edges[i];
		delete [] edges;
		};*/
};

int main(int argc, char** argv)
{
       //Test Graph
	/*Graph G(8);
	for(int i=0;i<8;i++)
              for(int j=0;j<8;j++)
                     G.edges[i][j]=0;

	G.edges[1][4]=1; G.edges[4][1]=1;
	G.edges[1][2]=1; G.edges[2][1]=1;
	G.edges[1][0]=1; G.edges[0][1]=1;
	G.edges[4][5]=1; G.edges[5][4]=1;
	G.edges[4][7]=1; G.edges[7][4]=1;
	G.edges[0][7]=1; G.edges[7][0]=1;
	G.edges[5][2]=1; G.edges[2][5]=1;
	G.edges[3][0]=1; G.edges[0][3]=1;
	G.edges[6][3]=1; G.edges[3][6]=1;
	G.edges[3][2]=1; G.edges[2][3]=1;
	G.edges[5][6]=1; G.edges[6][5]=1;
	G.edges[6][7]=1; G.edges[7][6]=1;*/
	int n_vertices=100;
	double prob=0.5;
	Graph G(n_vertices,prob);
	vector<int> dv(n_vertices,0);
    double time;
    
    time = -omp_get_wtime();

	#pragma omp parallel for
	       rep(i,0,n_vertices,1)
              {
                     rep(j,0,n_vertices,1)
                            dv[i] = dv[i]+(G.edges[i][j]);
                     int this_thread = omp_get_thread_num();
                     int num_threads = omp_get_num_threads();
                     std::cout<<"thread index is:"<<this_thread<<std::endl;
              }



	rep(i,0,n_vertices,1)
		std::cout<<"dv for index "<<i<<"="<<dv[i]<<std::endl;

	time += omp_get_wtime();

	std::cout<<"it took time:"<<time<<std::endl; 

       return 0;
}
