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
const int n=10;
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
		//cout<<start<<end<<endl;
		//getchar();
		if(start>=n || end>=n)
			exit(-1);
		edges.push_back(make_pair(start,end));
	}
	printf("construction end\n");
	fclose(fgraph);
}
void set_randval(int* randval,const int &begin,const int &end,unsigned int seed,const int &rank)
{
	//#pragma omp parallel private(seed)
	{
		seed=time(NULL)+19+17*rank;
		rep(i,begin,end)
		{
			randval[i]=rand_r(&seed);
		}
	}	
}
void choose_independent_vertice(const vector<pair<int,int> > &edges,int* randval,int* label)
{
	//#pragma omp parallel
	{
		//int chunk=edges.size()/omp_get_num_threads();
		//int end=(omp_get_thread_num()==omp_get_num_threads()-1)?edges.size():chunk*(omp_get_thread_num()+1);
		//#pragma omp for
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
void reconstruct_vertices(vector<int> &new_existing_vertices,vector<int> &independent_set,const vector<int> &existing_vertices,int *label,int begin,int end)
{

	/*#pragma omp parallel
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
		//#pragma omp critical
		{
			new_existing_vertices.insert(new_existing_vertices.end(),private_vertices.begin(),private_vertices.end());
			independent_set.insert(independent_set.end(),private_chosen.begin(),private_chosen.end());
		}
	} */
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
		//#pragma omp critical	
}
void reconstruct_edges(vector<pair<int,int> > &new_edges,const vector<pair<int,int> > &edges,int *label)
{
	//#pragma omp parallel
	{
		vector<pair<int,int> > private_edges;

		//#pragma omp for private(i) nowait 
		rep(i,0,edges.size())
		{
			if(label[edges[i].first]<0 && label[edges[i].second]<0)
			{
				//#pragma omp critical
				private_edges.push_back(edges[i]);				
			}
		}
		//#pragma omp critical
		new_edges.insert(new_edges.end(),private_edges.begin(),private_edges.end());
	}
}
void recover_label(int *label,vector<int> &existing_vertices,int begin,int end)
{
	//#pragma omp parallel for private(i) schedule(static,cache_size)
	rep(i,begin,end)
	{
		label[existing_vertices[i]]=0;
	}
}
int main(int argc,char* argv[])
{
	//cout<<cache_size<<endl;
	MPI_Init(&argc,&argv);
	int rank,size;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	FILE* ftime=NULL;
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
	MPI_Barrier(MPI_COMM_WORLD);
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
			if(rank==0)
				cout<<begin_vertices<<"\t"<<end_vertices<<endl;
			unsigned seed;
			int i;			
			double curr_time=MPI_Wtime();		
			set_randval(randval,begin_vertices,end_vertices,seed,rank);
			//send randval
			int this_vertice_chunk_size=end_vertices-begin_vertices;
			int *recvcounts=new int[size];
			rep(i,0,size-1)
			{
				recvcounts[i]=vertice_chunk_size;
			}
			recvcounts[size-1]=this_vertice_chunk_size;
			int* displs=new int[size];
			displs[0]=0;
			for(int i=1;i<size;++i)
				displs[i]=displs[i-1]+recvcounts[i-1];

			//If the number of vertices divide number of processes
			int* randval_send=new int[recvcounts[rank]];
			rep(i,0,recvcounts[rank])
			{
				randval_send[i]=randval[begin_vertices+i];
			}
			if(rank==0)
			cout<<"all gather"<<endl;
			rep(i,0,size){
			cout<<i<<"\t"<<recvcounts[rank]<<"\t"<<displs[rank]<<endl;
			MPI_Gatherv(&randval_send[0],recvcounts[rank],MPI_INT,&randval[0],recvcounts,displs,MPI_INT,i,MPI_COMM_WORLD);
			}
			delete [] randval_send;

			//getchar();
				//MPI_Recv(randval_send,existing_vertices.size(),MPI_INT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			delta_time[0]+=MPI_Wtime()-curr_time;
			//cout<<"rand val assigned"<<endl;
			//traverse the edge list,-1 remaining,1 represent neighbor 0 represent choose
			curr_time=MPI_Wtime();
			choose_independent_vertice(edges,randval,label);
			MPI_Barrier(MPI_COMM_WORLD);
			if(rank==0)
			cout<<"end"<<endl;
			MPI_Allreduce(MPI_IN_PLACE,&label[0],n,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
			//if(rank==1)
			//	rep(i,0,n)
			//		cout<<"label"<<label[i]<<endl;
			if(rank==0)
			cout<<"all reduce end"<<endl;
			MPI_Barrier(MPI_COMM_WORLD);
			delta_time[1]+=MPI_Wtime()-curr_time;
			//mark the neighbors
			curr_time=MPI_Wtime();
			//synchronize the labels

			choose_neighbor_vertices(edges,label,n);
			MPI_Allreduce(MPI_IN_PLACE,&label[0],n,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
			//cout<<"labels"<<endl;
			if(rank==0)
			rep(i,0,n)
				cout<<"label"<<label[i]<<endl;
			//cout<<"labels"<<endl;
			MPI_Barrier(MPI_COMM_WORLD);
			delta_time[2]+=MPI_Wtime()-curr_time;
			//add to independent set and make new existing_vertices and delete edges  
			//delete vertices
			curr_time=MPI_Wtime();
			int j=0;
			vector<int> new_existing_vertices;
			MPI_Barrier(MPI_COMM_WORLD);
			reconstruct_vertices(new_existing_vertices,independent_set,existing_vertices,label,begin_vertices,end_vertices);
			j=new_existing_vertices.size();
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Allgather(&j,1,MPI_INT,recvcounts,1,MPI_INT,MPI_COMM_WORLD);
			//cout<<"new vertice size rank"<<rank<<"\t"<<j<<endl;
			//synchronize the new existing vertices insert data
			//should be replaced by MPI_Allgatherv

			displs[0]=0;
			for(int i=1;i<size;++i)
			{
				displs[i]=displs[i-1]+recvcounts[i-1];
				//cout<<displs[i]<<endl;
			}
			int new_vertice_size=displs[size-1]+recvcounts[size-1];

			cout<<"new size"<<new_vertice_size<<endl;
			existing_vertices.resize(new_vertice_size);
			if(rank==0){
				cout<<"vertice size,rank"<<recvcounts[0]<<"\t"<<recvcounts[1]<<endl;
				cout<<new_existing_vertices.size()<<existing_vertices.size()<<endl;
				cout<<displs[0]<<displs[1]<<endl;
				cout<<j<<endl;
				//cout<<&new_existing_vertices[1]-&new_existing_vertices[0]<<endl;
			}
			int* local_existing_vertices=new int[j];
			rep(i,0,j)
				local_existing_vertices[i]=new_existing_vertices[i];
			MPI_Barrier(MPI_COMM_WORLD);
			rep(i,0,size)
			{
			cout<<"rank"<<i<<endl;
			MPI_Gatherv(&local_existing_vertices[0],j,MPI_INT,&existing_vertices[0],
				recvcounts,displs,MPI_INT,i,MPI_COMM_WORLD);
			}
			delete [] local_existing_vertices;
			if(rank==0)
			cout<<"all gather end"<<endl;
			/*if(rank==1){
				for(int i=0;i<existing_vertices.size();++i)
					cout<<existing_vertices[i]<<endl;
				getchar();
			}*/
			delete [] recvcounts;
			delete [] displs;

			//synchronize the edges insert data
			MPI_Barrier(MPI_COMM_WORLD);
			delta_time[3]+=MPI_Wtime()-curr_time;
			j=0;
			std::vector<pair<int,int> > new_edges;
			curr_time=MPI_Wtime();
			MPI_Barrier(MPI_COMM_WORLD);
			reconstruct_edges(new_edges,edges,label);
			j=new_edges.size();
			//#pragma omp parallel for private(i)
			rep(i,0,j)
				edges[i]=new_edges[i];
			edges.resize(j);
			delta_time[4]+=MPI_Wtime()-curr_time;
			curr_time=MPI_Wtime();
			vertice_chunk_size=existing_vertices.size()/size;
			begin_vertices=vertice_chunk_size*rank;
			if(rank==size-1)
				end_vertices=existing_vertices.size();
			else
				end_vertices=(rank+1)*vertice_chunk_size;	
			recover_label(label,existing_vertices,begin_vertices,end_vertices);
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Allreduce(MPI_IN_PLACE,&label[0],n,MPI_INT,MPI_LAND,MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);
			delta_time[5]+=MPI_Wtime()-curr_time;
		}
		int *recvcounts=new int[size];
		int independent_set_size=independent_set.size();
		getchar();
		//cout<<independent_set_size<<endl;
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Allgather(&independent_set_size,1,MPI_INT,&recvcounts[0],1,MPI_INT,MPI_COMM_WORLD);
		cout<<recvcounts[0]<<recvcounts[1]<<endl;
		//getchar();
		int* displs=new int[size];
		displs[0]=0;
		for(int i=1;i<size;++i)
			displs[i]=displs[i-1]+recvcounts[i-1];	
		int* total_independent_set=new int[sizeof(int)*(displs[i-1]+recvcounts[i-1])];
		rep(i,0,independent_set.size())
		{
			cout<<"indepe"<<independent_set[i]<<endl;
		}
		getchar();
		cout<<displs[size-1]+recvcounts[size-1]<<endl;
		cout<<"gather start"<<endl;
		if(rank==root)
			independent_set.resize(displs[size-1]+recvcounts[size-1]);
		MPI_Gatherv(&independent_set[0],independent_set_size,MPI_INT,&total_independent_set[0],
			recvcounts,displs,MPI_INT,root,MPI_COMM_WORLD);	
		//for(int i=0;i<independent_set.size();++i)
		//	cout<<independent_set[i]<<endl;
		delete [] recvcounts;
		delete [] displs;		
		delete [] total_independent_set;	
		if(rank==root){
			FILE* ftime=fopen(argv[1],"a");
			fprintf(ftime,"%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",n,delta_time[0],delta_time[1],delta_time[2],delta_time[3],delta_time[4],delta_time[5]);
			fclose(ftime);
			o("independent_set");
			cout<<independent_set.size()<<endl;
		}
		independent_set.clear();
	}
	MPI_Finalize();
	//set randval
}
