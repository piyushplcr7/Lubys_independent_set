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
#define rep(i,a,b) for(i=a;i<b;++i)
#define o(a) puts(a)
using std::vector;
using std::cout;
using std::pair;
using std::endl;
using std::make_pair;
const int n=10000;
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
void set_randval(int* randval,int begin,int end,int seed)
{
	#pragma omp parallel private(seed)
	{
		seed=time(NULL)+19*omp_get_thread_num()+17*rank;
		#pragma omp for private(i)
		rep(i,begin_vertices,end_vertices)
		{
			randval[i]=rand_r(&seed);
		}
	}	
}
void choose_independent_vertice(const vector<pair<int,int> > &edges,int* label,int n)
{
	#pragma omp parallel
	{
		//int chunk=edges.size()/omp_get_num_threads();
		//int end=(omp_get_thread_num()==omp_get_num_threads()-1)?edges.size():chunk*(omp_get_thread_num()+1);
		#pragma omp for
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
void recvcounts(int* count)
{
	
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
void reconstruct_vertices(vector<int> &new_existing_vertices,vector<int> &independent_set,const vector<int> &existing_vertices,int chunk_size_this,int rank)
{

	#pragma omp parallel
	{
		vector<int> private_vertices;
		vector<int> private_chosen;
		int chunk=chunk_size_this/omp_get_num_threads();
		int i;
		int end=(omp_get_thread_num()==omp_get_num_threads()-1)?chunk_size_this:chunk*(omp_get_thread_num()+1);
		//beginProcess=chunk_size*rank
		//beginThread=chunk*omp_get_thread_num()
		rep(i,chunk_size*rank+chunk*omp_get_thread_num(),chunk_size*rank+end)
		{
			if(label[existing_vertices[i]]==0)
			{
				private_chosen.push_back(existing_vertices[i]);
			}
			if(label[existing_vertices[i]]<0)
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
}
void reconstruct_edges(vector<pair<int,int> > &new_edges,const vector<pair<int,int> > &edges,int *label)
{
	#pragma omp parallel
	{
		vector<pair<int,int> > private_edges;
		int chunk=edges.size()/omp_get_num_threads();
		int i;
		int end=(omp_get_thread_num()==omp_get_num_threads()-1)?edges.size():chunk*(omp_get_thread_num()+1);
		//#pragma omp for private(i) nowait 
		rep(i,chunk*omp_get_thread_num(),end)
		{
			if(label[edges[i].first]<0 && label[edges[i].second]<0)
			{
				//#pragma omp critical
				private_edges.push_back(edges[i]);				
			}
		}
		#pragma omp critical
			new_edges.insert(new_edges.end(),private_edges.begin(),private_edges.end());
	}
}
void recover_label(int *label,int begin,int end)
{
	#pragma omp parallel for private(i) schedule(static,cache_size)
	rep(i,begin,end)
	{
		label[existing_vertices[i]]=0;
	}
}
int main(int argc,char* argv[])
{
	cout<<cache_size<<endl;
	cout<<argv[1]<<endl;
	FILE* ftime=fopen(argv[1],"w");
	MPI_Init(&argc,&argv);
	int rank,size;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	//double poss=0.2;
	vector<pair<int,int> > ori_edges;
	vector<int> independent_set;
	char filename[]="out.txt";
	//ER_graph_construction(n,ori_edges,poss);
	rmat_graph_construction(filename,ori_edges);
	int round;
	int round_num=20;
	srand(time(NULL));
	vector<double> total_delta_time(5,0.0);

	rep(round,0,round_num){
		vector<int> existing_vertices(n,0);
		int randval[n]={0};
		int label[n]={0};
		int edges_chunk_size;
		int end_edges;
		int edges_chunk_size=ori_edges.size()/size;
		int begin_edges=rank*chunk_size;
		if(rank==size-1)
			end_edges=ori_edges.size();
		else
			end_edges=(rank+1)*edges_chunk_size;
		vector<pair<int,int> > edges(ori_edges.begin()+begin_edges,ori_edges.begin()+end_edges);
		int i;
		rep(i,0,n)
			existing_vertices[i]=i;
		vector<double> delta_time(6,0.0);
		int vertice_chunk_size=existing_vertices.size()/size;
		int begin_vertices=vertice_chunk_size*rank;
		int end_vertices;
		if(rank==size-1)
			end_vertices=existing_vertices.size();
		else
			end_vertices=(rank+1)*vertice_chunk_size;	
		while(!existing_vertices.empty())
		{
			//set randval
			unsigned seed;
			int i;			
			double curr_time=omp_get_wtime();
			int this_vertice_chunk_size=end_vertices-begin_vertices;		
			set_randval(randval,begin_vertices,end_vertices);
			//If the number of vertices divide number of processes
			MPI_Allgather(randval+begin_vertices,this_vertice_chunk_size,MPI_INT,randval,vertice_chunk_size,MPI_INT,MPI_COMM_WORLD);
				//MPI_Recv(randval_send,existing_vertices.size(),MPI_INT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			delta_time[0]+=omp_get_wtime()-curr_time;
			//cout<<"rand val assigned"<<endl;
			//traverse the edge list,-1 remaining,1 represent neighbor 0 represent choose
			curr_time=omp_get_wtime();
			choose_independent_vertice(edges,label,n);
			MPI_Allreduce(label,label,n,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
			delta_time[1]+=omp_get_wtime()-curr_time;
			//mark the neighbors
			curr_time=omp_get_wtime();
			//synchronize the labels
			MPI_Allreduce(label,label,n,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
			delta_time[2]+=omp_get_wtime()-curr_time;
			//add to independent set and make new existing_vertices and delete edges  
			//delete vertices
			curr_time=omp_get_wtime();
			int j=0;
			vector<int> new_existing_vertices;
			reconstruct_vertices(new_existing_vertices,independent_set,existing_vertices,vertice_chunk_size,rank);
			j=new_existing_vertices.size();
			cout<<j<<endl;

			int* recvcounts=malloc(size*sizeof(int));

			MPI_Allgather(&j,1,MPI_INT,recvcounts,1,MPI_INT,MPI_COMM_WORLD);
			//synchronize the new existing vertices insert data
			//should be replaced by MPI_Allgatherv
			int* displs=malloc(size*sizeof(int));
			displs[0]=0;
			for(int i=1;i<size;++i)
				displs[i]=displs[i-1]+recvcounts[i-1];
			MPI_Allgatherv(&new_existing_vertices[0],existing_vertices.size(),MPI_INT,&existing_vertices[0],
				recvcounts,displs,MPI_INT,MPI_COMM_WORLD);
			delete [] recvcounts;
			delete [] displs;
			//synchronize the edges insert data
			existing_vertices.resize(j);
			delta_time[3]+=omp_get_wtime()-curr_time;
			j=0;
			std::vector<pair<int,int> > new_edges;
			curr_time=omp_get_wtime();
			reconstruct_edges(new_edges,edges,label);
			j=new_edges.size();
			#pragma omp parallel for private(i)
			rep(i,0,j)
				edges[i]=new_edges[i];
			edges.resize(j);
			delta_time[4]+=omp_get_wtime()-curr_time;
			curr_time=omp_get_wtime();
			vertice_chunk_size=existing_vertices.size()/size;
			begin_vertices=vertice_chunk_size*rank;
			if(rank==size-1)
				end_vertices=existing_vertices.size();
			else
				end_vertices=(rank+1)*vertice_chunk_size;	
			recover_label(label,begin_vertices,end_vertices);
			MPI_Allreduce(label,label,n,MPI_INT,MPI_LAND,MPI_COMM_WORLD);
			delta_time[5]+=omp_get_wtime()-curr_time;
		}
		fprintf(ftime,"%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",n,delta_time[0],delta_time[1],delta_time[2],delta_time[3],delta_time[4],delta_time[5]);
		o("independent_set");
		cout<<independent_set.size()<<endl;
		independent_set.clear();
	}
	MPI_Finalize();
	fclose(ftime);
	//set randval
}
