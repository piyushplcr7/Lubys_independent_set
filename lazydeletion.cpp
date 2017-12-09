#include<iostream>
#include<vector>
#include <omp.h>
#include<atomic>
#include<time.h>
#include<cstdlib>
#include<list>
#include<vector>
#include<iterator>
#include<chrono>
#include<algorithm>
#define rep(i,a,b,c) for(int i=a;i<b;i+=c)
#define Nthreads 4
using namespace std;


	class Edge;

	class Vertex{
	public:
			std::vector<std::vector<Edge>>::iterator edges;
			int randval;
			//atomic<int> randval;
			int id;
			int label;
			bool status; //true if the vertex is in the graph, false if deleted
			Vertex(int i,std::vector<std::vector<Edge>>::iterator assigned_edge_list):edges(assigned_edge_list)
			{
				id=i;
				randval=0;
				label = 1;
				status = true;
			}
		 int get_id()
		 {
		 	return id;
		 }
	};
  using E_iterator = std::vector<Edge>::iterator;
  using V_iterator = std::list<Vertex>::iterator;

	class Edge{
		public:
		V_iterator end;
		E_iterator same_edge;
		bool status;
		Edge(V_iterator y):end(y),status(true) {}
		Edge(){}
	};

class Graph {
public:
	friend void print_graph(Graph&);
	int n_edges;
	int n_;
	vector<V_iterator> current_vertices; //vector containing the pointers to vertices present in the graph
	list<Vertex> vertices;
  vector<vector<Edge>> EDGES;
	//create the matrix
	Graph(int n,double prob)//random graph G(n_vertices,prob) where prob is the probability that an edge exists
	{
		n_=n;
		n_edges=0;
		//vertices.reserve(n+2);
		current_vertices.reserve(n);
    EDGES.reserve(n);
    //EDGES.reserve(n);
		//srand(time(NULL));
		V_iterator source;
		V_iterator dest;
    vector<vector<Edge>>::iterator it = EDGES.begin();
    //setting the iterators
		for (int i = 0 ; i<n ; ++i)
    {
      vertices.push_back(Vertex(i,it+i));
			V_iterator End = vertices.end(); --End;
			current_vertices.push_back(End);
      (it+i)->reserve(n);
    }
		//#pragma omp parallel for private(source,dest)
		for(source=vertices.begin();source!=vertices.end();++source)
		{
			for(dest=vertices.begin();dest!=source;++dest){
				if(rand()*1.0/RAND_MAX<prob){
					++n_edges;
					source->edges->push_back(Edge(dest));
					dest->edges->push_back(Edge(source));
					source->edges->back().same_edge=--dest->edges->end();
					dest->edges->back().same_edge=--source->edges->end();
				}
			}
			source->randval=rand();
		}
	}

	Graph (const Graph& G) {
		std::cout << "copy constructor called!" << std::endl;
		int n = G.vertices.size();
		n_edges = G.n_edges;
		vector<V_iterator> current_vertices(G.current_vertices); //vector containing the pointers to vertices present in the graph
		list<Vertex> vertices(G.vertices);
		std::list<Vertex>::const_iterator temp = G.vertices.begin();
	  vector<vector<Edge>> EDGES(G.EDGES);
		vector<vector<Edge>>::iterator EDGES_iterator = EDGES.begin();
		vector<V_iterator>::iterator wth = current_vertices.begin();
		for (V_iterator it = vertices.begin() ; it!= vertices.end() ; ++it)
		{
			it->id = temp->id;
			std::cout << it->id << " | " ;
			*wth = it;
			it->edges = EDGES_iterator;
			++EDGES_iterator;
			++wth;
			++temp;
		}

	}

	void showDegree()
	{
		for(vector<V_iterator>::iterator iter=current_vertices.begin();iter!=current_vertices.end();++iter)
			cout<<"vertex"<<(*iter)->id<<"degree"<<(*iter)->edges->size()<<endl;
	}

	void set_randval()
	{
		#pragma omp parallel for schedule(static) num_threads(Nthreads)
		for(int i = 0 ; i<current_vertices.size();++i)
		{
			vector<V_iterator>::iterator iter = current_vertices.begin();
			iter=iter+i;
			(*iter)->randval=rand();
		}
	}

	bool isEmpty()
	{
		return current_vertices.empty();
	}

	void choose(vector<V_iterator> &chosen_set) //check this
	{

		#pragma omp parallel for schedule(static) num_threads(Nthreads)
		for(int i=0 ; i<current_vertices.size() ; ++i)
		{
			vector<V_iterator>::iterator iter=current_vertices.begin();
			iter = iter+i;
			if( (*iter)->label>0 )
			{
				for(E_iterator eiter= (*iter)->edges->begin();eiter!=(*iter)->edges->end();++eiter)
				{
					//#pragma omp critical
					if(eiter->end->randval>(*iter)->randval)
						eiter->end->label=0;
					else if(eiter->end->randval<(*iter)->randval)
					  (*iter)->label=0;
						else{
							if((*iter)->id>eiter->end->id)
							(*iter)->label=0;
							else
							eiter->end->label=0;
						}
						if((*iter)->label==0)
						 break;
				}
				#pragma omp critical
				{
					if((*iter)->label==1 && (*iter)->status == true)
					chosen_set.push_back(*iter);
				}
			}
		}

	}

	void deleteNodeNeighbor(V_iterator pvertex)
	{
		for( E_iterator neighbor_edge=pvertex->edges->begin() ; neighbor_edge!=pvertex->edges->end() ; ++neighbor_edge)
		{
			deleteNode(neighbor_edge->end,pvertex);
		}
		pvertex->status = false;
		//std::cout << "Remaining graph after deletion of: " << pvertex->id << " and its neighbors \n";
		//print_graph(*this);
	}

	void deleteNode(V_iterator pvertex,V_iterator not_delete)
	{
		//std::cout << " deleting the node " << pvertex->id << " Neighbor of node " <<not_delete->id << std::endl;
		int temp = pvertex->id;
		const E_iterator begining = pvertex->edges->begin();
		//#pragma omp parallel
		//{
			#pragma omp parallel for schedule(static) num_threads(Nthreads)
			for( int i = 0 ; i<pvertex->edges->size() ; ++i)
			{
				E_iterator to_be_deleted = begining + i;
					if(to_be_deleted->end->id!=not_delete->id)
					{
						//to_be_deleted->status = false; no need to delete this edge as the node for this is going to be deleted
						to_be_deleted->same_edge->status = false;
					}
			} //end of for loop
		pvertex->status = false;
		//update_current_vertices();
		//std::cout << " delete node ended \n" << std::endl;
		//std::cout << "Remaining graph after deletion of node: " << temp <<" \n";
		//print_graph(*this);
	}

	void deleteSet(vector<V_iterator>::iterator begin,vector<V_iterator>::iterator end)
	{
		for(vector<V_iterator>::iterator viter=begin;viter!=end;++viter)
		{
			deleteNodeNeighbor(*viter);
		}
		update_current_vertices();
		update_edge_list();
	}

void update_current_vertices() {
	vector<V_iterator> new_current_vertices; new_current_vertices.reserve(n_);
	#pragma omp parallel for schedule(static) num_threads(Nthreads)
	for (int i = 0 ; i<current_vertices.size() ; ++i )        	//for ( vector<V_iterator>::iterator it = current_vertices.begin() ; it!=current_vertices.end() ;  )
	{
		if ((current_vertices[i])->status == true && (current_vertices[i])->edges->size()>0 ) //still in the graph
		{
			#pragma omp critical
			new_current_vertices.push_back(current_vertices[i]); //add to new_current_vertices vector
		}
	}
	unsigned int newsize = new_current_vertices.size(); //std::cout << "new current vertices size: " << newsize << std::endl;
	#pragma omp parallel for schedule(static) num_threads(Nthreads)
	for (int i = 0 ; i<newsize ; ++i )
	{
		current_vertices[i] = new_current_vertices[i];
		current_vertices[i]->label = 1;
	}
	current_vertices.resize(newsize);
}

void update_edge_list() {
	int size = current_vertices.size();

	for (int i = 0 ; i<size ; ++i) //going over the current vertices
	{
		vector<Edge> temp; temp.reserve(n_);
		unsigned int esize = current_vertices[i]->edges->size();
		#pragma omp parallel for schedule(static) num_threads(Nthreads)
		for (int j = 0; j<esize ; ++j)
		{
			if(current_vertices[i]->edges->operator[](j).status == true)
			{
				#pragma omp critical
				temp.push_back(current_vertices[i]->edges->operator[](j));
			}
		}
		unsigned int newsize = temp.size();
		#pragma omp parallel for schedule(static) num_threads(Nthreads)
		for (int j = 0; j<newsize ; ++j)
		{
			current_vertices[i]->edges->operator[](j) = temp[j];
			current_vertices[i]->edges->operator[](j).same_edge->same_edge = current_vertices[i]->edges->begin() + j; //fixing the same edge
		}
		current_vertices[i]->edges->resize(newsize);
	}
}

};

void print_graph(Graph& G) {
	for (vector<V_iterator>::iterator it = G.current_vertices.begin() ; it!=G.current_vertices.end() ; ++it) {
		std::cout << "Vertex: " << (*it)->id <<" |";
		for (E_iterator ned = (*it)->edges->begin() ; ned!=(*it)->edges->end() ; ++ned) {
			std::cout << " -> " << ned->end->id;
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

void verify(vector<int>& mis, Graph& G) {
	bool final_decision = true; //assume it is MIS
	int n = G.vertices.size();
	vector<bool> MIS(n,false);
	for (vector<int>::iterator it=mis.begin() ; it!=mis.end() ; ++it)
		MIS[*it] = true;
	//iterating over the graph to verify independent or not
	for ( V_iterator it = G.vertices.begin() ; it!=G.vertices.end() ; ++it) {
		int count = 0;
		for (E_iterator ned = (it)->edges->begin() ; ned!=(it)->edges->end() ; ++ned) {
			if (MIS[it->id] == true && MIS[ned->end->id] == true ) //if the vertex whose edge list is being traversed is in the MIS its neighbor is in the MIS(false)
				final_decision = false;
			else if (MIS[it->id] == false && MIS[ned->end->id] == false)
				++count;
		}
		if (count == (it)->edges->size() && count > 0)
			final_decision = false;
	}
	if (final_decision)
		std::cout << "		MIS" << std::endl;
	else
		std::cout << "		Verified! IT IS NOT MIS" << std::endl;

}

int main(int argc, char** argv)
{
	//for (int trials = 0 ; trials < 1000 ; ++trials) {
	long long int seed = time(NULL);
	srand(seed);
	//std::cout << "SEED VALUE USED = " << seed << std::endl;

		//std::cout << "trial no: " << trials << std::endl;
	double prob=0.5;
	int n_vertices=700,i=0; //for n_vertices 36 the problem begins (same found mis and not mis; no seg fault) (17 for seg fault)
	//for (int n_vertices=10;n_vertices<1000000;n_vertices=n_vertices+50)
	//{
	Graph G(n_vertices,prob);
	//Graph Gcopy(G);
	//std::cout << "Graph Initialized! with seed: " << seed << std::endl;
	//print_graph(G);
	srand(seed);
	Graph Gcopy(n_vertices,prob);
	//std::cout << "Graph copy Initialized! \n" << std::endl;
	//print_graph(Gcopy);
	//V_iterator pvertex=G.vertices.begin(); V_iterator not_delete = pvertex;
	//++not_delete; //vertex 1
	//G.deleteNode(pvertex,not_delete);
	//vector<Edge>::iterator t = pvertex->edges.begin();
	//++t; //edge 1->2
	//std::cout << "Special edge delete called for  "<< pvertex->id << " -> " << t->end->id << "\n";
	//G.special_delete_edge(t);
	//G.showDegree();
	vector<V_iterator > chosen_set;
	vector<int> independent_set;
	auto start = std::chrono::system_clock::now();
	//#pragma omp parallel for
	clock_t t_start = clock();
			while(!G.isEmpty()){
				//cout<<"round"<<++i<<endl;
				G.set_randval();
			//	cout<<"Rand val generated"<<endl;
				G.choose(chosen_set);
				//cout<<"The independent set nodes chosen, size = "<<chosen_set.size()<<endl;
				//rep(i,0,chosen_set.size(),1)
					//cout<<chosen_set[i]->id<<endl;
				for(vector<V_iterator >::iterator iter=chosen_set.begin() ; iter!=chosen_set.end() ; ++iter)
				{
					independent_set.push_back((*iter)->id);
					//std::cout << "TEST 1" << std::endl;
				}
				//std::cout << "TEST 2" << std::endl;
				G.deleteSet(chosen_set.begin(),chosen_set.end());
				//cout<<" Nodes and neighbors deleted"<<endl;
				chosen_set.clear();
				//cout<<" Vertex size"<<G.vertices.size()<<endl;
          }
					//independent_set.push_back(G.current_vertices[0]->id);
	auto end = std::chrono::system_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
	std::cout <<n_vertices<<"	"<<G.n_edges<<"		"<< elapsed.count() << std::endl;
	//independent_set.push_back(0);
	//std::cout << "Independent set: " << std::endl;
	//for (vector<int>::iterator it = independent_set.begin() ; it != independent_set.end() ; ++it)
		//std::cout << *it << " " ;
	//std::cout << std::endl;

	//std::cout << "Graph final! \n" << std::endl;
	//print_graph(G);
	//std::cout << "Graph copy final! \n" << std::endl;
	//print_graph_copy(Gcopy);
	verify(independent_set,Gcopy);
	//independent_set.clear();
//}

    return 0;
}
