//implementation of the greedy algorithm for finding a maximal independent set in a given graph:(V,E) where V is the set of vertices and E is the set of edges
#include<iostream>
#include<array>
#include<cstdlib>
#include<ctime>

/* IGNORE THIS PART RIGHT NOW, ITS FOR IMPLEMENTING THE ADJACENCY LIST DATA STRUCTURE

//structure for the vertex nodes of the adjacency list
struct node_v{
	public:
	node* down;
	node* up;
	node* right;
	const unsigned int vertex_id;

	public:
	//constructor to initialize the node with vertex_id with v and the pointers with NULL value
	node_v(unsigned int v):vertex_id(v),down(NULL),up(NULL),right(NULL) {};
	//Default destructor
	~node_v{};
};

//structure for edge nodes of the adjacency list
struct node_e{
	public:
	node* left;
	node* right;
	const unsigned int vertex_id;

	public:
	//constructor to initialize the node with vertex_id with v and the pointers with NULL value
	node_e(unsigned int v):vertex_id(v),left(NULL),right(NULL) {};
	//Default destructor
	~node_e{};
};

void rm_v(node_v& v,node_v* head) {
	//remove a vertex 

}


//Template struct for graph where Edges and Vertices can be represented as an adjacency list and matrix or simple 2 dimensional matrix and array
//template <typename T>
*/

struct Graph {
	public:
	int n_vertices;
	//Set of Vertices; contains 1 if vertex is in the graph and 0 if not in the graph
	int* vertices=new int[n_vertices];;	
	//std::array<int,n_vertices> vertices;
	//set of Edges
	//std::array<std:array<int,n_vertices>,n_vertices> edges;
	int** edges=new int*[n_vertices];;	
	//constructor
	Graph(int n):n_vertices(n){
		for (int i=0;i<n_vertices;++i)
			edges[i]=new int[n_vertices];
						};
	//destructor
	~Graph(){
		delete [] vertices;
		for (int i = 0; i < n_vertices; ++i)
    		delete [] edges[i];
		delete [] edges;
		};
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
	Graph G(8);
	for (int i=0; i<8;++i)
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
	G.edges[6][7]=1; G.edges[7][6]=1;

	
	//function to find the maximal independent set mis(G) {G is a graph object} which returns a set of vertices which constitute an independent set
 	mis(G);
	//verification of the correctness of algorithm
	for (int j=0; j<8;++j)
	{
	if (G.vertices[j]!=0)
		std::cout<<j<<std::endl;	
	}	
	return 0;
}
