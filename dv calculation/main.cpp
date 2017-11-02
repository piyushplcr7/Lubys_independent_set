#include<iostream>
#include<mpi.h>
//#include<vector>
#include<time.h>
#include<cstdlib>
//using std::vector;
//main process, initializes V no. of processes which individually calculate the degree of one vertex and return it to the main process

//Using the structure Graph (with array for vertices and matrix for edges)

struct Graph {
	public:
	//no. of vertices
	int n_vertices;
	//Set of Vertices; contains 1 if vertex is in the graph and 0 if not in the graph
	//vector<int> vertices;
	int* vertices;
	//set of Edges
	//vector<vector<int> > edges;
	int** edges;
	//constructor
	Graph(int n,double prob):n_vertices(n)//random graph G(n_vertices,prob) where prob is the probability that an edge exists
		{
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
			for (int j=i;j<n;++j)
				{
				edges[i][j]=(rand()*1.0/RAND_MAX<prob);
				edges[j][i]=edges[i][j];
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



//main function
int main(int argc, char** argv) {
	MPI_Init(&argc,&argv);
	double prob=.5;
	MPI_Status status;
	MPI_Comm intercomm;
	double time;
	int NP,temp;
	for (int n_v=10;n_v<200;n_v=n_v+20)
		{
		//Test Graph
		//n_v: no. of vertices for the random graph
		 //probability for an edge to exist for the random graph
		Graph G(n_v,prob); //random graph generator
		
		NP=n_v; //NP is the no. of Vthreads to be spawned, equal to the no. of vertices
		int errcodes[NP];
		int *dv=new int[NP];
		
		 //take start time
		//MPI_Comm_rank(MPI_COMM_WORLD,&rank);
		//MPI_Comm_size(MPI_COMM_WORLD,&size);
		//dv array to be populated by the spawned processes
		
		//std::cout<<"Creating "<< NP << " Processes for calculation of dv"<<std::endl;
		//spawning Vthreads
		//time = -MPI_Wtime();
		MPI_Comm_spawn("Vthread",MPI_ARGV_NULL,NP,MPI_INFO_NULL,0,MPI_COMM_WORLD,&intercomm,errcodes);
		//sending the data to vthreads (rows of size NP of the edge matrix) for calculation of dv
		time = -MPI_Wtime();
		for (int i=0;i<NP;++i)
			MPI_Send(*(G.edges+i),NP,MPI_INT,i,0,intercomm);
		time += MPI_Wtime();
		//receiving the dv data from the Vthreads and storing it in dv array
		for (int i=0;i<NP;++i)
			{
			MPI_Recv(&temp,1,MPI_INT,MPI_ANY_SOURCE,2,intercomm,&status);
			dv[status.MPI_SOURCE]=temp;
			}
		//for (int i=0;i<100000;++i) {}
		//MPI_Comm_size(intercomm,&size1);
		//for (int i=0;i<NP;++i)
		//	std::cout<<"dv for index "<<i<<"="<<dv[i]<<std::endl;
	       //time += MPI_Wtime(); //take end time
	       std::cout<<time<<"	"<< n_v<<std::endl;
		delete[] dv;
		}
	MPI_Finalize();
	return 0;
}
