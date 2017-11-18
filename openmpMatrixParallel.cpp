#include<iostream>
#include<vector>
#include <omp.h> 
#include<time.h>
#include<cstdlib>
#include<list>
#include<iterator>
#include<omp.h>
#include<cstdio>
#define rep(i,a,b,c) for(i=a;i<b;i+=c)

//typedef list<int>::iterator Vertex
using namespace std;

	class Edge;
	class Vertex{
	public:
		/*
		class Edge{
		public:
		list<Vertex>::iterator end;
		list<Edge>::iterator same_edge;
		Edge(list<Vertex>::iterator y)
		{
			end=y;
		}
		Edge()
		{
			end=NULL;
			same_edge=NULL
		}
		};*/     
			vector<int> edges;
			int num_edges;//number of existing edges
			int randval;// random value
			int id;//The id for this vertex
			int label;//whether possible for deletion in this round
			bool deleted;//whether the vertex has been deleted
			Vertex(int i,int n)
			{
				deleted=false;
				edges=vector<int>(n,0);
				id=i;
				randval=0;
				label=1;
				num_edges=0;
			}
		 int get_id()
		 {
		 	return id;
		 }
	};

class Graph {
public:

	//int n_vertices;
	//Set of Vertices; contains 1 if vertex is in the graph and 0 if not in the graph
	vector<Vertex> vertices;
	int num_vertices;//number of existing vertex
	int N;//number of vertex
	Graph(int n)
	{
		N=n;
		num_vertices=n;
		int i;
		rep(i,0,n,1){
			vertices.push_back(Vertex(i,n));
			//(vertices+i)->edges=list<Edge>();
		}
	//create the matrix
	}
	Graph(int n,double prob)//random graph G(n_vertices,prob) where prob is the probability that an edge exists
	{
		N=n;
		num_vertices=n;
		int i,source,dest;
		rep(i,0,n,1)
		vertices.push_back(Vertex(i,n));
		rep(source,0,N,1)
		{
			rep(dest,0,source,1){
				if(rand()*1.0/RAND_MAX<prob){
					(vertices[source].edges)[dest]=1;
					(vertices[dest].edges)[source]=1;
					(vertices[source]).num_edges++;
					(vertices[dest]).num_edges++;				
				}
			}
		}
	}
	void showDegree()
	{
		for(vector<Vertex>::iterator iter=vertices.begin();iter!=vertices.end();++iter)
		{
			cout<<"vertex"<<iter->id<<"degree"<<iter->num_edges<<endl;		
		}
	}
	void set_randval()
	{	
		int i;
		#pragma omp parallel for private(i) schedule(dynamic)	
		rep(i,0,vertices.size(),1){
			//cout<<i<<"thread_num"<<omp_get_thread_num()<<endl;
			if(vertices[i].deleted==false){
				vertices[i].randval=rand();
			}
		}
	}
	bool isEmpty()
	{
		return (num_vertices==0);
	}
	void choose(vector<int> &chosen_set)
	{
		int choose_number=0;
		int neighbor_edge;
		int i;
		#pragma omp parallel for private(i,neighbor_edge) schedule(dynamic)
		rep(i,0,vertices.size(),1)
		{
			//cout<<"thread"<<omp_get_thread_num()<<endl;
			//cout<<"max thread"<<omp_get_max_threads()<<endl;
			vector<Vertex>::iterator iter=vertices.begin()+i;
			//cout<<omp_get_thread_num()<<endl;
			if(iter->deleted==true)
				continue;
			if(iter->label==1){
				rep(neighbor_edge,0,N,1){
					if(iter->edges[neighbor_edge]==0 || vertices[neighbor_edge].deleted==true)//edge not exist
					{
						continue;
					}
					if(iter->randval<=(vertices[neighbor_edge]).randval)
					{
						(vertices[neighbor_edge]).label=0;
					}
					else
					{
						iter->label=0;
						break;
					}
				}
				if(iter->label==1){
					#pragma omp atomic
					++choose_number;
					#pragma omp critical
						chosen_set.push_back(iter->id);
				}
			}	
		}
		//chosen_set.resize(choose_number);	
	}
	void deleteNodeNeighbor(vector<Vertex>::iterator pvertex,int signal)
	{
		if(pvertex->deleted==true)
			return;
		#pragma omp critical
		{
		if(pvertex->deleted==true)
			signal=1;
		//num_vertices-=(1-pvertex->deleted);
		else
			pvertex->deleted=true;
		}
		if(signal==1)
			return;
		#pragma omp atomic
		num_vertices--;
		int neighbor_edge;
		#pragma omp parallel for default(shared) private(neighbor_edge) schedule(dynamic)
		for(neighbor_edge=0;neighbor_edge<(pvertex->edges).size();neighbor_edge++)
		{
			if(pvertex->edges[neighbor_edge]==1)
			{
				int signal=0;
				deleteNode(neighbor_edge,signal);
			}
		}
		//cout<<pvertex->deleted<<endl;
		pvertex->num_edges=0;
	}
	void deleteNode(int nvertex,int signal)
	{
		//cout<<"delete Node"<<pvertex->id<<endl;

		if((vertices[nvertex]).deleted==true)
		{
			return;
		}
		#pragma omp critical
		{
		if((vertices[nvertex]).deleted==true)
			signal=1;
		//num_vertices-=(1-pvertex->deleted);
		else
			(vertices[nvertex]).deleted=true;
		}
		if(signal==1)
			return;
		#pragma omp atomic
		num_vertices--;
		int i;
		#pragma omp parallel for default(shared) private(i) schedule(static)
		rep(i,0,vertices.size(),1){
				#pragma omp atomic
				(vertices.begin()+i)->num_edges-=(vertices.begin()+i)->edges[nvertex];
				(vertices.begin()+i)->edges[nvertex]=0;
			}
			(vertices[nvertex]).num_edges=0;

	}
	void deleteSet(vector<int>::iterator begin,int size)
	{
		int i;
		#pragma omp parallel for private(i) schedule(dynamic)
		rep(i,0,size,1)
		{
		    int signal=0;
			deleteNodeNeighbor(vertices.begin()+*(begin+i),signal);
		}
		#pragma omp parallel for private(i)
		rep(i,0,vertices.size(),1)
		{
			if((vertices.begin()+i)->deleted==false)
				(vertices.begin()+i)->label=1;	
		}
	}
};

int main(int argc, char** argv)
{
	FILE* fp=fopen("output3.txt","w");
	FILE* ftime=fopen("time3.txt","w");
	srand(time(NULL));
	int num_round=20;
	for(int n_vertices=1000;n_vertices<=2000;n_vertices+=200){
		double delta=0;
		int avg_iteration=0;
		double avg_iter_float=0;
		for(int round=0;round<num_round;++round){
				double prob=0.5;
				Graph G(n_vertices,prob);
				//G.showDegree();
				//cout<<"Graph initialized"<<endl;
				vector<int> chosen_set;
				vector<int> independent_set;
			    //cout<<"contruction end"<<endl;
			    
			    double time = omp_get_wtime();
				//cout<<"begin"<<endl;
		//#pragma omp parallel for
				while(!G.isEmpty()){
					++avg_iteration;
					//cout<<"round"<<++i<<"num"<<G.num_vertices<<endl;
					G.set_randval();
					//cout<<"Rand val generated"<<endl;
					G.choose(chosen_set);
					//cout<<"The independent set nodes chosen"<<endl;
					int oldsize=independent_set.size();
					independent_set.resize(oldsize+chosen_set.size());
					int i;
					#pragma omp parallel for
					rep(i,0,chosen_set.size(),1){
						independent_set[oldsize+i]=*(chosen_set.begin()+i);
					}
					G.deleteSet(chosen_set.begin(),chosen_set.size());
					//cout<<"num_vertices"<<G.num_vertices<<endl;
					//cout<<" Nodes and neighbors deleted"<<endl;
					chosen_set.clear();
					//cout<<" Vertex size"<<G.vertices.size()<<endl;
	          	}
	          	//gettimeofday(&end, NULL);
				delta += omp_get_wtime()-time;
	    //cout<<"independent set found"<<endl;
	    /*if(round%10==0){
			rep(i,0,independent_set.size(),1)
				std::cout<<"node"<<i<<"in MIS:\t"<<independent_set[i]<<std::endl;
		      }*/
	  	}
	    avg_iter_float=avg_iteration*1.0/num_round;
	    fprintf(ftime,"%d\t%.2fs\n",n_vertices,delta);
	    fprintf(fp,"%d\t%lf\n",n_vertices,avg_iter_float);
	}
	//cout<<independent_set.size()<<endl;
	fclose(fp);

	//time += omp_get_wtime();
	//std::cout<<"it took time:"<<time<<std::endl; 
    return 0;
}
