//implementation of the greedy algorithm for finding a maximal independent set in a given graph:(V,E) where V is the set of vertices and E is the set of edges
#include<iostream>
#include<array>
#include<cstdlib>
#include<fstream>
#include<chrono>
#include<vector>
using namespace std;

using Vlist_t = std::vector<int>;
using Emat_t = std::vector<Vlist_t>;

void mis(Vlist_t &Vertices,Emat_t &Edges,vector<int>& independent_set,const int& N)
	{
	//int r1=rand()%N;
	//int r2=rand()%N;
	int j;
	for(int i=0;i<N;++i) //iterate over all the vertices of the graph
		{
		//i=(it+r1)%N;
		if (Vertices[i]!=0) //if the current vertex is still in the remaining set
			{
        independent_set.push_back(i);
			//greedily remove all its neighbours if they are still in the remaining set
			for(int jt=0 ; jt<Edges[i].size() ; ++jt) //iterate over all the neighbours of the current vertex
				{
				//j=(jt+r2)%Edges[i].size();
        j = Edges[i][jt];
				if(Vertices[j]!=0)
					Vertices[j]=0;
				}
			}
		}
	}

void Construct_graph(Emat_t &Edges, const int& N,char* filename) {
  cout << "Reading graph from: " << filename << endl;
  FILE* fp=fopen(filename,"r");
  int start,end,n_edges = 0;
  while(fscanf(fp,"%d%d",&start,&end)!=EOF)
		{
			Edges[start][end]=1;
			Edges[end][start]=1;
			n_edges++;
		}
    cout<<n_edges<<endl;
		fclose(fp);
}

void Construct_graph(Emat_t &Edges, const int& N, const double& prob) {
  cout << "ER graph construction " << endl;
  int n_edges = 0;
  for(int i = 0; i<N ; ++i)
  {
    for(int j = i+1 ; j<N ; ++j)
    {
      if(rand()*1.0/RAND_MAX<prob){
        ++n_edges;
        Edges[i].push_back(j);
        Edges[j].push_back(i);
      }
    }
  }
}

void print_graph(const Emat_t &Edges, const int& N) {
	for (int i = 0 ; i<N ; ++i) {
		std::cout << "Vertex: " << i <<" |";
		for (int j = 0 ;j<Edges[i].size() ; ++j) {
			std::cout << " -> " << Edges[i][j];
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}


int main() {
  srand(time(NULL));
  int N = 10; //no. of vertices
  double prob = 0.3;
  Vlist_t Vertices(N,1);
  Emat_t Edges(N);
  Construct_graph(Edges,N,prob);
  //print_graph(Edges,N);
	//char filename[]="output.txt";
	//function to find the maximal independent set mis(G) {G is a graph object} which returns a set of vertices which constitute an independent set
	auto start = std::chrono::system_clock::now();
	vector<int> independent_set;
//  cout << "test 1" << endl;
 	mis(Vertices,Edges,independent_set,N);
  //cout << "test 2 " << endl;
	auto end = std::chrono::system_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
//  std::cout << "Independent set vertices: " << std::endl;
//  for (int i =0 ; i<independent_set.size() ; ++i)
  //  std::cout << independent_set[i] << " ";
  //std::cout<<std::endl;
  std::ofstream output;
  output.open("greedy_data.txt",'w');
	output<<"Independent set size: " << independent_set.size()<< " Calculated in: " << "		"<< elapsed.count() << " us"<< std::endl;
	return 0;
}
