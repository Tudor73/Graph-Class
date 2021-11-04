#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <fstream>
#include <algorithm>
using namespace std;

class Graph {
    int nrNodes, nrEdges;
    bool oriented;
    vector<vector<int>> edges;
    vector<vector<int>> transposedGraph;
    stack<int> topoSort;
public:
    Graph(int nrNodes, bool oriented = true) ;
    void addEdge(int x, int y); // sets edge between x and y
    void BFS(int startingNode, ostream& g);
    int connectedComponents();  // returns the number of conencted components from the graph
    void DFS(int node, vector<int>& visited);
    void topologicalSort(int x, vector<int> &visited);
    void stronglyConnectedComponents(ofstream & g);   // outputs the strongly connected components
    ~Graph();
private:
    void DFS2(int currentNode, vector<int> &component, int currentComponent, vector<vector<int>> & solution); // used in stronglyConnComponents

};
Graph::Graph(int nrNodes, bool oriented) {
    this->nrNodes = nrNodes;
    this->oriented = oriented;
    edges.resize(nrNodes + 1);
    this->transposedGraph.resize(this->nrNodes+1);
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
void Graph::topologicalSort(int currentNode, vector<int> &visited) {
    visited[currentNode] = 1;
    for(auto node: this->edges[currentNode]){
        if(!visited[node])
            topologicalSort(node, visited);
    }
    this->topoSort.push(currentNode);
}
void Graph::DFS2(int currentNode, vector<int> &visited, int currentComponent, vector<vector<int>> & solution) { // dfs used for strongly conn comp
    visited[currentNode] = 1;
    solution[currentComponent].push_back(currentNode);
    for(auto i: this->transposedGraph[currentNode]){
        if(visited[i] == 0){
            DFS2(i,visited,currentComponent,solution);
        }
    }
}
void Graph::stronglyConnectedComponents(ofstream &g) {
    vector<vector<int>>solution(nrNodes+1);
    this->transposedGraph.resize(this->nrNodes+1);
    for(int i = 1; i<= this->nrNodes; i++){             // compute transposed graph
        for(auto j : edges[i]){
            transposedGraph[j].push_back(i);
        }
    }
    vector<int>visited(this->nrNodes+1, 0);
    for(int i = 1; i<= nrNodes; i++){
        if(!visited[i])
            topologicalSort(i, visited); // sort in topologcal order
    }

    vector<int>visited2(this->nrNodes+1, 0);

    int currentComponent = 0;
    while(!topoSort.empty()) {
        int i = topoSort.top();
        if(visited2[i] == 0) {
            currentComponent++;
            DFS2(i, visited2, currentComponent,solution);
        }
        topoSort.pop();
    }
    g << currentComponent << "\n";
    for(int i = 1; i<=currentComponent; i++){
        for(auto j : solution[i]){
            g << j << " ";
        }
        g << "\n";
    }
}
Graph::~Graph(){
    edges.clear();
}
int main()
{
    int nrNodes, nrEdges;
    ifstream f("ctc.in");
    ofstream g("ctc.out");
    f >> nrNodes >> nrEdges;
    Graph G(nrNodes, true); // false - unorientedGraph, true- orientedGraph
    for(int i = 1; i<= nrEdges; i++) {
        int x, y;
        f >> x >> y;
        G.addEdge(x, y);
    }
    G.stronglyConnectedComponents(g);
    return 0;


}

