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

using namespace std;


	class Edge;

	class Vertex{
	public:
			Edge* edges; //pointer to the edge list for the vertex
			int randval;
			unsigned int size_edges;
			//atomic<int> randval;
			int id; //vertex id
			int label; //label to indicate if the vertex is valid for selection in a certain round
			bool status; //true if the vertex is in the graph, false if deleted
			Vertex(int i,Edge* assigned_edge_list):id(i)//,edges(assigned_edge_list)
			{
				//id=i;
				edges = assigned_edge_list;
				randval=0;
				label = 1;
				status = true;
				size_edges = 0;
			}
		 int get_id()
		 {
		 	return id;
		 }
	};

	using V_iterator = std::list<Vertex>::iterator;

	class Edge{
		public:
		V_iterator end;
		Edge* same_edge;
		int check;
		Edge(V_iterator y):end(y) {check = 10;}
		Edge(){check = 10;}
	};

class Graph {
public:
	friend void print_graph(Graph&);
	int n_edges;
	int n_;
	vector<V_iterator> current_vertices; //vector containing the pointers to vertices present in the graph
	list<Vertex> vertices; //a linked list for the vertices
  Edge** EDGES; //pointer to a pointer for Edge (for creating the 2d structure)
	//create the matrix
	~Graph() {
		for (int i = 0 ; i<n_ ; ++i)
    {
			delete[] *(EDGES+i);
		}
		delete[] EDGES;
	}
	Graph(int n,double prob)//random graph G(n_vertices,prob) where prob is the probability that an edge exists
	{
		n_edges=0;
		n_ = n;
		//vertices.reserve(n+2);
		current_vertices.reserve(n+2);
  	EDGES = new Edge*[n]; //for n edge lists
		V_iterator source;
		V_iterator dest;
    //setting the iterators and pointers
		for (int i = 0 ; i<n ; ++i)
    {
			EDGES[i] = new Edge[n]; //reserving space for the edge list
      vertices.push_back( Vertex(i,EDGES[i]) ); //initializing the vertex
			V_iterator End = vertices.end(); --End; //getting the iterator for initialized vertex
			current_vertices.push_back(End); //adding the pointer for the initialized vertex into the vector current vertices
			//std::cout << "Edge list check: " << EDGES[i]->check << std::endl;
    }
		//std::cout << "ctr1" << std::endl;
		//#pragma omp parallel for private(source,dest)
		for(source=vertices.begin();source!=vertices.end();++source)
		{
			for(dest=vertices.begin();dest!=source;++dest){
				if(rand()*1.0/RAND_MAX<prob){
					++n_edges;
					//source->edges->push_back(Edge(dest));
					Edge temp(dest);
					//std::cout << "ctr2" << std::endl;
					unsigned int source_index = source->size_edges;
					//std::cout << "source index: " << source_index << std::endl;
					//std::cout << "Edge list check: " << source->edges->check << std::endl;
					*(source->edges + source_index) = temp;
					//std::cout << "ctr3" << std::endl;
					//dest->edges->push_back(Edge(source));
					Edge temp1(source);
					dest->edges[dest->size_edges] = temp1;
					source->edges[source->size_edges].same_edge = (dest->edges + dest->size_edges) ;
					dest->edges[dest->size_edges].same_edge = (source->edges + source->size_edges);
					++dest->size_edges ; ++source->size_edges;
				}
			}
			source->randval=rand();
		}
	}

	/*Graph (const Graph& G) {
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
*/
	void showDegree()
	{
		for(vector<V_iterator>::iterator iter=current_vertices.begin();iter!=current_vertices.end();++iter)
			cout<<"vertex"<<(*iter)->id<<"degree"<<(*iter)->size_edges<<endl;
	}

	void set_randval()
	{
		#pragma omp parallel for
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

		#pragma omp parallel for
		for(int i=0 ; i<current_vertices.size() ; ++i)
		{
			vector<V_iterator>::iterator iter=current_vertices.begin();
			iter = iter+i;
			if( (*iter)->label>0 )
			{
				//for(E_iterator eiter= (*iter)->edges->begin();eiter!=(*iter)->edges->end();++eiter)
				for(int j= 0; j < (*iter)->size_edges ; ++j)
				{
					Edge* eiter = (*iter)->edges + j ;
					//#pragma omp critical
					if(eiter->end->randval > (*iter)->randval)
						eiter->end->label=0;
					else if(eiter->end->randval < (*iter)->randval)
					  (*iter)->label=0;
						else{
							if((*iter)->id > eiter->end->id)
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
		//for( E_iterator neighbor_edge=pvertex->edges->begin() ; neighbor_edge!=pvertex->edges->end() ; ++neighbor_edge)
		for( unsigned int j = 0 ; j < pvertex->size_edges ; ++j)
		{
			Edge* neighbor_edge = pvertex->edges + j;
			//vector<Vertex>::iterator temp = neighbor_edge->end;
			deleteNode(neighbor_edge->end,pvertex/*vertices.end()*/);
			//if (pvertex==vertices.end())
				//pvertex = temp;
		}
		//if (vertices.size()>1)
		special_delete_node(pvertex);
		//update_current_vertices();
		//else
			//vertices.clear();
		//std::cout << "Remaining graph after deletion of: " << pvertex->id << " and its neighbors \n";
		//print_graph(*this);
	}

	void deleteNode(V_iterator pvertex,V_iterator not_delete)
	{
		//std::cout << " deleting the node " << pvertex->id << " Neighbor of node " <<not_delete->id << std::endl;
		Edge* begining = pvertex->edges;
			#pragma omp parallel for //num_threads(4)
			for( int i = 0 ; i<pvertex->size_edges ; ++i)
			{
				Edge* to_be_deleted = begining + i;
				//neighbor_edge = neighbor_edge + i ;
				//if (to_be_deleted->end->status == true && to_be_deleted->end->size_edges > 0 )
				{
					if(to_be_deleted->end->id!=not_delete->id)
					{
						//std::cout << "Deleting the edge called for: " << temp << " -> " << to_be_deleted->end->id << "\n";
						Edge* End = to_be_deleted->end->edges + to_be_deleted->end->size_edges-1; //end of the edge list containing same edge, used for swap
						Edge* Betw = to_be_deleted->same_edge;
						//std::cout << "To be swapped: ("<<Betw->end->id << "," << End->end->id << ")" << "\n";
						if (End!=Betw) {
							//implement a manual swap
							Betw->end = End->end;
							Edge** temp = &(Betw->same_edge); Edge* assignment = End->same_edge;

							//Betw->same_edge = End->same_edge;	// this line has to be critical or atomic!!!!!!!!!!!!
							#pragma omp atomic
							*temp = assignment;
							//std::cout << "Pointers swapped to: ("<<Betw->end->id << "," << End->end->id << ")" << "\n";
				      Betw->same_edge->same_edge = Betw;
							//End->same_edge->same_edge = End;
							//print_graph(*this);
							//#pragma omp critical
							//to_be_deleted->end->size_edges = to_be_deleted->end->size_edges - 1;
							--to_be_deleted->end->size_edges;
							//to_be_deleted->end->edges->pop_back();
						}
						else
						{
							//#pragma omp critical
							to_be_deleted->end->size_edges = to_be_deleted->end->size_edges - 1;
							//to_be_deleted->end->edges->pop_back(); //no need to fix the memory locations
						}
						//std::cout << "Remaining graph: " << std::endl;
						//print_graph(*this);
					//}

					}
					//else
						//break;
				}

				//#pragma omp barrier
				//std::cout << "I'm inside parallel region" << std::endl;
			}
			//#pragma omp barrier
		//}

		//std::cout << "I'm outside parallel region" << std::endl;
		special_delete_node(pvertex);
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

	}

	/*void special_delete_edge( E_iterator to_be_deleted ) {
		//std::cout << "deleting the edge : "<< to_be_deleted->end->id << " -> " << to_be_deleted->same_edge->end->id << std::endl;
		//if (to_be_deleted->end->status == true && to_be_deleted->end->edges->size() > 0)
		//{
			E_iterator End = to_be_deleted->end->edges->end(); //end of the edge list containing same edge, used for swap
			End--;
			E_iterator Betw = to_be_deleted->same_edge;
			//std::cout << "To be swapped: ("<<Betw->end->id << "," << End->end->id << ")" << "\n";
			if (End!=Betw) {
				//implement a manual swap
				Betw->end = End->end;
				Betw->same_edge = End->same_edge;
				//std::cout << "Pointers swapped to: ("<<Betw->end->id << "," << End->end->id << ")" << "\n";
	      Betw->same_edge->same_edge = Betw;
				//print_graph(*this);
				to_be_deleted->end->edges->pop_back();
			}
			else
				to_be_deleted->end->edges->pop_back(); //no need to fix the memory locations
			//std::cout << "Remaining graph: " << std::endl;
			//print_graph(*this);
		//}
	}*/

	void special_delete_node(V_iterator to_be_deleted ) {
		/*if (vertices.size() > 1) {
		std::cout << "Special delte node called for node: "<<to_be_deleted->id << std::endl;
		V_iterator End = vertices.end();
		--End;
		if (to_be_deleted != End)
		{
			std::swap(*to_be_deleted,*End);
			vertices.pop_back();
			for (vector<Edge>::iterator iter = to_be_deleted->edges->begin() ; iter!=to_be_deleted->edges->end() ; ++iter)
			{
				iter->same_edge->same_edge = iter;
				iter->same_edge->end = to_be_deleted;
			}
		}
		else
			vertices.pop_back();
		std::cout << "Special delte node ended \n" << std::endl;
	}*/
	to_be_deleted->status = false;
}

void update_current_vertices() {
	for ( vector<V_iterator>::iterator it = current_vertices.begin() ; it!=current_vertices.end() ;  )
	{
		if ((*it)->status == false ) //not in the graph anymore
		{
			vector<V_iterator>::iterator End = current_vertices.end();
			--End;
			std::swap(*it,*End);
			current_vertices.pop_back();
		}
		else
		{
			(*it)->label = 1;
			++it;
		}

	}
}
};

void print_graph(Graph& G) {
	for (vector<V_iterator>::iterator it = G.current_vertices.begin() ; it!=G.current_vertices.end() ; ++it) {
		std::cout << "Vertex: " << (*it)->id <<" |";
		for (unsigned int i = 0; i < (*it)->size_edges ; ++i) {
			Edge* ned = (*it)->edges + i;
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
		for (unsigned int i = 0; i < it->size_edges ; ++i) {
			Edge* ned = it->edges + i;
			if (MIS[it->id] == true && MIS[ned->end->id] == true ) //if the vertex whose edge list is being traversed is in the MIS its neighbor is in the MIS(false)
				final_decision = false;
			else if (MIS[it->id] == false && MIS[ned->end->id] == false)
				++count;
		}
		if (count == it->size_edges && count > 0)
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
	long long int seed = 42;//time(NULL);
	srand(seed);
	//std::cout << "SEED VALUE USED = " << seed << std::endl;

		//std::cout << "trial no: " << trials << std::endl;
	double prob=0.5;
	int n_vertices=700,i=0; //for n_vertices 36 the problem begins (same found mis and not mis; no seg fault) (17 for seg fault)
	//for (int n_vertices=10;n_vertices<1000000;n_vertices=n_vertices+50)
	//{
	//std::cout << "Test1" << std::endl;
	Graph G(n_vertices,prob);
	//std::cout << "Test2" << std::endl;
	//Graph Gcopy(G);
	//std::cout << "Graph Initialized! with seed: " << seed << std::endl;
	//print_graph(G);
	//std::cout << "Test3" << std::endl;
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
