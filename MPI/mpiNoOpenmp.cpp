//how to parallel for mpi?
//First we need to compute the 
#include<iostream>
#include<vector>
#include<time.h>
#include<cstdlib>
#include<list>
#include<iterator>
#include<omp.h>
#include<mpi.h>
#include<cstdio>
#define rep(i,a,b) for(int i=a;i<b;++i)
#define o(a) puts(a)
using std::vector;
using std::cout;
using std::pair;
using std::endl;
using std::make_pair;
const int n=30000;
const int cache_size=64/sizeof(int);
const int root=0;
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
	printf("construction begin\n");
	FILE* fgraph=fopen(filename,"r");
	int start,end;
	while(fscanf(fgraph,"%d%d",&start,&end)!=EOF)
	{
		if(start>=n || end>=n)
			exit(-1);
		edges.push_back(make_pair(start,end));
	}
	printf("construction end\n");
	fclose(fgraph);
}
void set_randval(int* randval,const int &begin,const int &end,unsigned int &seed,const int &rank)
{
	//#pragma omp parallel private(seed)
	{
		rep(i,begin,end)
		{
			randval[i]=rand_r(&seed);
		}
	}	
}
void choose_independent_vertice(const vector<pair<int,int> > &edges,int* randval,int* label)
{
	{
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
	}
}

void choose_neighbor_vertices(const vector<pair<int,int> > &edges,int* label,int n)
{
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
}
void reconstruct_vertices(vector<int> &independent_set,vector<int> &existing_vertices,int *label,int begin,int end)
{
		int j=0;
		rep(i,begin,end)
		{
			if(label[existing_vertices[i]]==0)
			{
				independent_set.push_back(existing_vertices[i]);
			}
			else if(label[existing_vertices[i]]<0)
			{
				existing_vertices[j++]=existing_vertices[i];
			}
		}
		existing_vertices.resize(j);
}
void reconstruct_vertices(vector<int>& new_existing_vertices,vector<int> &independent_set,const vector<int> &existing_vertices,int *label,int begin,int end)
{
		rep(i,begin,end)
		{
			if(label[existing_vertices[i]]==0)
			{
				independent_set.push_back(existing_vertices[i]);
			}
			else if(label[existing_vertices[i]]<0)
			{
				new_existing_vertices.push_back(existing_vertices[i]);
			}
		}
}
void reconstruct_edges(vector<pair<int,int> > &edges,int *label)
{
		int j=0;
		rep(i,0,edges.size())
		{
			if(label[edges[i].first]<0 && label[edges[i].second]<0)
			{
				edges[j++]=edges[i];				
			}
		}
		edges.resize(j);

}
void recover_label(int *label,vector<int> &existing_vertices,int begin,int end)
{
	rep(i,begin,end)
	{
		label[existing_vertices[i]]=0;
	}
}
int main(int argc,char* argv[])
{
	MPI_Init(&argc,&argv);
	int rank,size;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	FILE* ftime=NULL;
	if(rank==root)
		ftime=fopen(argv[1],"w");
	vector<pair<int,int> > ori_edges;
	vector<int> independent_set;

	//double poss=0.2;
	//ER_graph_construction(n,ori_edges,poss);
	char filename[]="out.txt";
	rmat_graph_construction(filename,ori_edges);
	srand(time(NULL));
	unsigned seed=time(NULL)+19+17*rank;
	int round_num=20;
	vector<double> total_delta_time(5,0.0);
	rep(round,0,round_num){
		MPI_Barrier(MPI_COMM_WORLD);
		vector<int> existing_vertices(n,0);
		int randval[n]={0};
		int label[n]={0};
		int end_edges;
		int edges_chunk_size=ori_edges.size()/size;
		int begin_edges=rank*edges_chunk_size;
		if(rank==size-1)
			end_edges=ori_edges.size();
		else
			end_edges=(rank+1)*edges_chunk_size;
		vector<pair<int,int> > edges(ori_edges.begin()+begin_edges,ori_edges.begin()+end_edges);
		int i;
		rep(i,0,n)
			existing_vertices[i]=i;
		int vertice_chunk_size=existing_vertices.size()/size;
		int begin_vertices=vertice_chunk_size*rank;
		int end_vertices;
		if(rank==size-1)
			end_vertices=existing_vertices.size();
		else
			end_vertices=(rank+1)*vertice_chunk_size;	
		vector<double> delta_time(6,0.0);
		while(!existing_vertices.empty())
		{
			//set randval
			int i;			
			MPI_Barrier(MPI_COMM_WORLD);
			double curr_time=MPI_Wtime();		
			set_randval(randval,begin_vertices,end_vertices,seed,rank);
			//send randval
			//int this_vertice_chunk_size=end_vertices-begin_vertices;
			int last_vertice_chunk_size=existing_vertices.size()-(size-1)*vertice_chunk_size;
			int *recvcounts=new int[size];
			rep(i,0,size-1)
			{
				recvcounts[i]=vertice_chunk_size;
			}
			recvcounts[size-1]=last_vertice_chunk_size;
			int* displs=new int[size];
			displs[0]=0;
			for(int i=1;i<size;++i)
				displs[i]=displs[i-1]+recvcounts[i-1];

			//If the number of vertices divide number of processes

			MPI_Allgatherv(MPI_IN_PLACE,recvcounts[rank],MPI_INT,&randval[0],recvcounts,displs,MPI_INT,MPI_COMM_WORLD);
			delta_time[0]+=MPI_Wtime()-curr_time;
			//traverse the edge list,-1 remaining,1 represent neighbor 0 represent choose
			curr_time=MPI_Wtime();
			choose_independent_vertice(edges,randval,label);
			MPI_Allreduce(MPI_IN_PLACE,&label[0],n,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
			delta_time[1]+=MPI_Wtime()-curr_time;

			curr_time=MPI_Wtime();
			//mark the neighbors
			choose_neighbor_vertices(edges,label,n);
			//synchronize the labels			
			MPI_Allreduce(MPI_IN_PLACE,&label[0],n,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
			delta_time[2]+=MPI_Wtime()-curr_time;
			//add to independent set and make new existing_vertices and delete edges  
			//delete vertices
			curr_time=MPI_Wtime();
			vector<int> new_existing_vertices;
			reconstruct_vertices(new_existing_vertices,independent_set,existing_vertices,label,begin_vertices,end_vertices);
			int local_new_size=new_existing_vertices.size();
			MPI_Allgather(&local_new_size,1,MPI_INT,recvcounts,1,MPI_INT,MPI_COMM_WORLD);

			displs[0]=0;
			for(int i=1;i<size;++i)
			{
				displs[i]=displs[i-1]+recvcounts[i-1];
			}
			int new_vertice_size=displs[size-1]+recvcounts[size-1];
			cout<<new_vertice_size<<endl;
			existing_vertices.resize(new_vertice_size);
			MPI_Allgatherv(&new_existing_vertices[0],local_new_size,MPI_INT,&existing_vertices[0],
				recvcounts,displs,MPI_INT,MPI_COMM_WORLD);
			delete [] recvcounts;
			delete [] displs;

			//synchronize the edges insert data
			//MPI_Barrier(MPI_COMM_WORLD);
			delta_time[3]+=MPI_Wtime()-curr_time;
			int j=0;
			//std::vector<pair<int,int> > new_edges;
			curr_time=MPI_Wtime();
			reconstruct_edges(edges,label);
			//MPI_Barrier(MPI_COMM_WORLD);			
			delta_time[4]+=MPI_Wtime()-curr_time;
			curr_time=MPI_Wtime();
			//update the vertice chunk size
			vertice_chunk_size=existing_vertices.size()/size;
			begin_vertices=vertice_chunk_size*rank;
			if(rank==size-1)
				end_vertices=existing_vertices.size();
			else
				end_vertices=(rank+1)*vertice_chunk_size;	
			recover_label(label,existing_vertices,begin_vertices,end_vertices);
			MPI_Allreduce(MPI_IN_PLACE,&label[0],n,MPI_INT,MPI_LAND,MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);
			delta_time[5]+=MPI_Wtime()-curr_time;
		}
		int *recvcounts=new int[size];
		int independent_set_size=independent_set.size();
		//cout<<independent_set_size<<endl;
		MPI_Allgather(&independent_set_size,1,MPI_INT,&recvcounts[0],1,MPI_INT,MPI_COMM_WORLD);
		//getchar();
		int* displs=new int[size];
		displs[0]=0;
		for(int i=1;i<size;++i)
			displs[i]=displs[i-1]+recvcounts[i-1];	
		int total_independent_set_size=(displs[size-1]+recvcounts[size-1]);
		if(rank==root)
			independent_set.resize(total_independent_set_size);
		int* total_independent_set=NULL;
		if(rank==root)
			total_independent_set=new int[total_independent_set_size];
		MPI_Gatherv(&independent_set[0],independent_set_size,MPI_INT,&independent_set[0],
			recvcounts,displs,MPI_INT,root,MPI_COMM_WORLD);	
		//for(int i=0;i<independent_set.size();++i)
		//	cout<<independent_set[i]<<endl;
		delete [] recvcounts;
		delete [] displs;		
		delete [] total_independent_set;	
		if(rank==root){
			fprintf(ftime,"%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",n,delta_time[0],delta_time[1],delta_time[2],delta_time[3],delta_time[4],delta_time[5],delta_time[1]+delta_time[2]+delta_time[3]+delta_time[4]+delta_time[5]);
			o("independent_set");
			cout<<total_independent_set_size<<endl;
		}
		//delete [] total_independent_set;
		independent_set.clear();
	}
	if(rank==root)
		fclose(ftime);
	MPI_Finalize();
	//set randval
}
