//implementation of the greedy algorithm for finding a maximal independent set in a given graph:(V,E) where V is the set of vertices and E is the set of edges
#include<iostream>
#include<array>
#include<cstdlib>
#include<chrono>


struct Graph {
	public:
	//no. of vertices
	int n_vertices;
	int n_edges;
	//Set of Vertices; contains 1 if vertex is in the graph and 0 if not in the graph
	//vector<int> vertices;
	int* vertices;
	//set of Edges
	//vector<vector<int> > edges;
	int** edges;
	//constructor
	Graph(int n,double prob):n_vertices(n)//random graph G(n_vertices,prob) where prob is the probability that an edge exists
		{
		n_edges=0;
		vertices=new int[n_vertices];
		edges=new int*[n_vertices];
		for (int i=0;i<n;++i)
			{			
			edges[i]=new int[n];
			vertices[i]=1;			
			}
		//Try to seed the random number generator with a constant (makes debugging easier as we can obtain the same graph)
		srand(7);
		for (int i=0;i<n;++i)
			{
			for (int j=i+1;j<n;++j)
				{
				
				edges[i][j]=(rand()*1.0/RAND_MAX<prob);
				edges[j][i]=edges[i][j];
				
				n_edges+=edges[i][j];
				}
			edges[i][i]=0;
			
			}
		
		}
	//destructor
	~Graph(){
		delete [] vertices;
		for (int i = 0; i < n_vertices; ++i)
    		delete [] edges[i];
		delete [] edges;		
		}
};

void mis(Graph& G) 
	{
	int r1=rand()%G.n_vertices;
	int r2=rand()%G.n_vertices;
	int i,j;
	for(int it=0;it<G.n_vertices;++it) //iterate over all the vertices of the graph
		{
		i=(it+r1)%G.n_vertices;
		if (G.vertices[i]!=0) //if the current vertex is still in the remaining set
			{ 
			//greedily remove all its neighbours if they are still in the remaining set
			for(int jt=0;jt<G.n_vertices;++jt) //iterate over all the neighbours of the current vertex
				{
				j=(jt+r2)%G.n_vertices;
				if(G.edges[i][j]==1 && G.vertices[j]!=0)
					G.vertices[j]=0;
				}
			}
		}
	}


int main() {
	srand(time(0));
	//initialization or construction of the target graph
	
	/*for (int i=0; i<8;++i)
		G.vertices[i]=1;
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

	
	//function to find the maximal independent set mis(G) {G is a graph object} which returns a set of vertices which constitute an independent set
	for (int n=10;n<1000000;n=n+50)
	{
		Graph G(n,0.5);	
		auto start = std::chrono::system_clock::now();
	 	mis(G);
		auto end = std::chrono::system_clock::now();
		auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
		std::cout <<G.n_vertices<<"	"<<G.n_edges<<"		"<< elapsed.count() << std::endl;
	}

	return 0;
}
