#include<iostream>
#include<mpi.h>
//main process, initializes V no. of processes which individually calculate the degree of one vertex and return it to the main process

//Using the structure Graph (with array for vertices and matrix for edges)

struct Graph {
	public:
	int n_vertices;
	//Set of Vertices; contains 1 if vertex is in the graph and 0 if not in the graph
	int* vertices;
	//std::array<int,n_vertices> vertices;
	//set of Edges
	//std::array<std:array<int,n_vertices>,n_vertices> edges;
	int** edges;
	//constructor
	Graph(int n):n_vertices(n)
		{
		vertices=new int[n_vertices];
		edges=new int*[n_vertices];
		//create the matrix
		for (int i=0;i<n_vertices;++i)
			{
			edges[i]=new int[n_vertices];
			vertices[i]=1;
			}
		};
	//destructor
	~Graph(){
		delete [] vertices;
		for (int i = 0; i < n_vertices; ++i)
    		delete [] edges[i];
		delete [] edges;
		};
};



//main function
int main(int argc, char** argv) {
	//Test Graph
	Graph G(8);
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

	//int rank,size;
	double time;
	int NP=8,temp;
	int errcodes[NP];
	MPI_Status status;
	MPI_Comm intercomm;
	MPI_Init(&argc,&argv);
	time = -MPI_Wtime(); //take start time
	//MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	//MPI_Comm_size(MPI_COMM_WORLD,&size);
	//dv array to be populated by the spawned processes
	int *dv=new int[NP];
	std::cout<<"Creating "<< NP << " Processes for calculation of dv"<<std::endl;
	//spawning Vthreads
	MPI_Comm_spawn("Vthread",MPI_ARGV_NULL,NP,MPI_INFO_NULL,0,MPI_COMM_WORLD,&intercomm,errcodes);
	//sending the data to vthreads (rows of size NP of the edge matrix) for calculation of dv
	for (int i=0;i<NP;++i)
		MPI_Send(*(G.edges+i),NP,MPI_INT,i,0,intercomm);
	//receiving the dv data from the Vthreads and storing it in dv array
	for (int i=0;i<NP;++i)
		{
		MPI_Recv(&temp,1,MPI_INT,MPI_ANY_SOURCE,2,intercomm,&status);
		dv[status.MPI_SOURCE]=temp;
		}
	//for (int i=0;i<100000;++i) {}
	//MPI_Comm_size(intercomm,&size1);
	for (int i=0;i<NP;++i)
		std::cout<<"dv for index "<<i<<"="<<dv[i]<<std::endl;
       time += MPI_Wtime(); //take end time
       std::cout<<"this took time:"<<time<<std::endl;
	MPI_Finalize();
	return 0;
}
