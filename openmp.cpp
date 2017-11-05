#include<iostream>
#include<vector>
//#include <omp.h>
#include<time.h>
#include<cstdlib>
#include<list>
#include<iterator>
#include<chrono>
#define rep(i,a,b,c) for(int i=a;i<b;i+=c)

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
			list<Edge> edges;
			int randval;
			int id;
			Vertex(int i)
			{
				id=i;
				randval=0;
			}

			//~Vertex()
			//{
			//	delete [] edges;
			//}
		 int get_id()
		 {
		 	return id;
		 }

		/*Vertex(list<Vertex*>::iterator h)
		{
			head=h;
			//num_edges=0;
		}
		//void setlist(list<Vertex*>::iterator begin)
		//{
		//	head=begin;
		//}*/
	};
	class Edge{
		public:
		list<Vertex>::iterator end;
		list<Edge>::iterator same_edge;
		Edge(list<Vertex>::iterator y)
		{
			end=y;
		}
	};

class Graph {
public:

	int n_edges;
	//Set of Vertices; contains 1 if vertex is in the graph and 0 if not in the graph
	list<Vertex> vertices;
	//std::array<int,n_vertices> vertices;
	//set of Edges
	//std::array<std:array<int,n_vertices>,n_vertices> edges;
	//list<list<Vertex*> > edges;
	//constructor
	Graph(int n)
	{
		n_edges=0;
		rep(i,0,n,1){
			vertices.push_back(Vertex(i));
			//(vertices+i)->edges=list<Edge>();
		}
	//create the matrix
	}
	Graph(int n,double prob)//random graph G(n_vertices,prob) where prob is the probability that an edge exists
	{
		n_edges=0;
		//needs to be optimized
		srand(time(NULL));
		list<Vertex>::iterator source;
		list<Vertex>::iterator dest;
		rep(i,0,n,1)
		vertices.push_back(Vertex(i));
		for(source=vertices.begin();source!=vertices.end();++source)
		{
			for(dest=vertices.begin();dest!=source;++dest){
				//edges[i][j]=(rand()*1.0/RAND_MAX<prob);
				if(rand()*1.0/RAND_MAX<prob){
					++n_edges;
					(source->edges).push_back(Edge(dest));
					(dest->edges).push_back(Edge(source));
					(source->edges).back().same_edge=--dest->edges.end();
					(dest->edges).back().same_edge=--source->edges.end();
					//cout<<"add edge"<<((source->edges).back().same_edge->end->id)<<"\t"<<(dest->edges).back().same_edge->end->id<<endl;
				}
			}
			source->randval=rand();
		}
	}
	void showDegree()
	{
<<<<<<< HEAD
		//for(list<Vertex>::iterator iter=vertices.begin();iter!=vertices.end();++iter)
			//cout<<"vertex"<<iter->id<<"degree"<<iter->edges.size()<<endl;		
=======
		for(list<Vertex>::iterator iter=vertices.begin();iter!=vertices.end();++iter)
			cout<<"vertex"<<iter->id<<"degree"<<iter->edges.size()<<endl;
>>>>>>> ce1611e117f41f8764e7058aa5562dd09f345c77
	}
	void set_randval()
	{
		//srand(time(NULL));
		for(list<Vertex>::iterator iter=vertices.begin();iter!=vertices.end();++iter)
			iter->randval=rand();
	}
	bool isEmpty()
	{
		return vertices.empty();
	}
	void choose(vector<list<Vertex>::iterator> &chosen_set)
	{
		for(list<Vertex>::iterator iter=vertices.begin();iter!=vertices.end();++iter)
		{
		//	cout<<"1/d:"<<1.0/(2*iter->edges.size())<<"randval"<<iter->randval<<endl;
			if(iter->randval>0 && (iter->randval*1.0/RAND_MAX)<(1.0/(2*iter->edges.size())))
			{
				for(list<Edge>::iterator eiter=iter->edges.begin();eiter!=iter->edges.end();++eiter)
					eiter->end->randval=-1;
				chosen_set.push_back(iter);
			}
		}
	}
	void deleteNodeNeighbor(list<Vertex>::iterator pvertex)
	{
		//cout<<"deleteNodeneighbor"<<pvertex->id<<endl;
		//for(list<Edge>::iterator neighbor_edge=(pvertex->edges).begin();neighbor_edge!=(pvertex->edges).end();neighbor_edge++)
		//	cout<<neighbor_edge->end->id<<endl;
		for(list<Edge>::iterator neighbor_edge=(pvertex->edges).begin();neighbor_edge!=(pvertex->edges).end();neighbor_edge++)
			deleteNode(neighbor_edge->end,pvertex);
		vertices.erase(pvertex);
	}
	void deleteNode(list<Vertex>::iterator pvertex,list<Vertex>::iterator not_delete)
	{
		//cout<<"delete Node"<<pvertex->id<<endl;

		for(list<Edge>::iterator neighbor_edge=(pvertex->edges).begin();neighbor_edge!=(pvertex->edges).end();neighbor_edge++)
		{
			if(neighbor_edge->end!=not_delete){
			neighbor_edge->end->edges.erase(neighbor_edge->same_edge);
			//cout<<neighbor_edge->end->id<<"\t"<<neighbor_edge->end->edges.size()<<endl;
		}
		}
		//cout<<"delete end"<<endl;
		vertices.erase(pvertex);
	}
	void deleteSet(vector<list<Vertex>::iterator>::iterator begin,vector<list<Vertex>::iterator>::iterator end)
	{
		for(vector<list<Vertex>::iterator>::iterator viter=begin;viter!=end;++viter)
		{
			deleteNodeNeighbor(*viter);
			//cout<<(*viter)->id<<endl;
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
	//~Graph(){
	//	delete [] vertices;
	//	};
};

int main(int argc, char** argv)
{
<<<<<<< HEAD
	//int n_vertices=200;
=======
	int n_vertices=10000;
>>>>>>> ce1611e117f41f8764e7058aa5562dd09f345c77
	double prob=0.5;
	for (int n_vertices=10;n_vertices<1000000;n_vertices=n_vertices+50)
	{
	Graph G(n_vertices,prob);
	//G.showDegree();
	vector<list<Vertex>::iterator > chosen_set;
	vector<int> independent_set;
    //double time;
    //cout<<"contruction end"<<endl;
    //int i=0;
    //time = -omp_get_wtime();
	auto start = std::chrono::system_clock::now();
	//#pragma omp parallel for
	clock_t t_start = clock();
			while(!G.isEmpty()){
				//cout<<"round"<<++i<<endl;
				G.set_randval();
				//cout<<"Rand val generated"<<endl;
				G.choose(chosen_set);
				//cout<<"The independent set nodes chosen"<<endl;
				//rep(i,0,chosen_set.size(),1)
				//	cout<<chosen_set[i]->id<<endl;
				for(vector<list<Vertex>::iterator >::iterator iter=chosen_set.begin();iter!=chosen_set.end();++iter)
					independent_set.push_back((*iter)->id);
				G.deleteSet(chosen_set.begin(),chosen_set.end());
				//cout<<" Nodes and neighbors deleted"<<endl;
				chosen_set.clear();
				//cout<<" Vertex size"<<G.vertices.size()<<endl;
          	}
	auto end = std::chrono::system_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
	std::cout <<n_vertices<<"	"<<G.n_edges<<"		"<< elapsed.count() << std::endl;
	independent_set.clear();
	}
    //cout<<"independent set found"<<endl;
	//rep(i,0,independent_set.size(),1)
		//std::cout<<"node"<<i<<"in MIS:\t"<<independent_set[i]<<std::endl;
	//cout<<chosen_set.size()<<endl;
	//time += omp_get_wtime();
<<<<<<< HEAD
	//std::cout<<"it took time:"<<time<<std::endl; 
=======
	clock_t t_end = clock();
	double exe_time = (double)(t_end - t_start) / CLOCKS_PER_SEC;
	std::cout<<"it took time:"<<exe_time<<std::endl;
>>>>>>> ce1611e117f41f8764e7058aa5562dd09f345c77
    return 0;
}
