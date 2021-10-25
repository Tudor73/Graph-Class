#include <iostream>
#include <vector>
#include <queue>
#include <fstream>
using namespace std;

class Graph {
    int nrNodes, nrEdges;
    bool oriented;
    vector<vector<int>> edges;
public:
    Graph(int nrNodes, bool oriented = true) ;
    void addEdge(int x, int y); // sets edge between x and y
    void BFS(int startingNode, ostream& g);
    int connectedComponents();  // returns the number of conencted components from the graph
    void DFS(int node, vector<int>& visited);
    ~Graph();
};
Graph::Graph(int nrNodes, bool oriented) {
    this->nrNodes = nrNodes;
    this->oriented = oriented;
    edges.resize(nrNodes + 1);
}
void Graph::addEdge(int x, int y) {
    this->edges[x].push_back(y);
    if(!oriented) {
        this->edges[y].push_back(x);
    }
}
void Graph::BFS(int startingNode, ostream& g) {
    vector<int> cost(this->nrNodes+1, -1);
    vector<int> visited(this->nrNodes+1, 0);
    queue<int> queue;

    cost[startingNode] = 0;
    visited[startingNode] = 1;
    queue.push(startingNode);

    while(!queue.empty()){
        int currentNode = queue.front();
        for(auto node : edges[currentNode]){
            if(!visited[node]){
                visited[node] = 1;
                cost[node] = cost[currentNode] + 1;
                queue.push(node);
            }
        }
        queue.pop();
    }
   for(int i = 1; i<= nrNodes; i++ ){
       g << cost[i] << " ";
   }
}
void Graph::DFS(int currentNode, vector<int> &visited) {
    for(auto node : this->edges[currentNode]) {
        if(!visited[node]) {
            visited[node] = 1;
            DFS(node, visited);
        }
    }
}
int Graph::connectedComponents() {
    int connectedComp = 0;
    vector<int>visited(this->nrNodes + 1, 0);
    for(int i = 1; i<= nrNodes; i++) {
        if(!visited[i]) {
            visited[i] = 1;
            connectedComp++;
            DFS(i, visited);
        }
    }
    return connectedComp;
}
Graph::~Graph(){
    edges.clear();
}
int main()
{
    int nrNodes, nrEdges;
    ifstream f("dfs.in");
    ofstream g("dfs.out");
    f >> nrNodes >> nrEdges;
    Graph G(nrNodes, false); // false - unorientedGraph, true- orientedGraph
    for(int i = 1; i<= nrEdges; i++) {
        int x, y;
        f >> x >> y;
        G.addEdge(x,y);
    }
    g << G.connectedComponents();

}

