#ifndef MAIN_CPP_GRAPH_H
#define MAIN_CPP_GRAPH_H
#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <fstream>
#include <algorithm>
#include <tuple>

using namespace std;

class Graph {
protected:
    int nrNodes, nrEdges;
    bool oriented;
    vector<vector<int>> edges;
public:
    Graph( const int nrNodes, bool oriented = true) ;
    void addEdge(const int x,const int y); // sets edge between x and y
    int getNumberOfNodes();
    vector<int> BFS(int startingNode);
    vector<int> DFS(int node);
    int connectedComponents();  // returns the number of conencted components from the graph
    int getGraphDiameter();
    vector<int> getTopologicalSortedGraph();
    vector<vector<int>> stronglyConnectedComponents();   // outputs the strongly connected components
    vector<vector<int>> biconnectedComponents(); // returns the array of biconnected components
    vector<pair<int,int>> getCriticalConnections(); // return array of critical edge
    bool havelHakim(vector<int> grades); //  returns true or false if the grades array form a valid graph
    ~Graph();
private:
    void dfsAux(int node, vector<int>& visited);
    void topologicalSort(int x, vector<int> &visited, vector<int>& topoSort);
    void DFS_strongly_conn_comp(int currentNode, vector<int> &component, int currentComponent, vector<vector<int>> & solution, vector<vector<int>>& transposedGraph); // used in stronglyConnComponents
    void DFS_crit(int node, int predecesor, int level, vector<int>& lvl,vector<int>& low, vector<pair<int,int>>& output ); // used in getCriticalConnections
    void DFS_bcc(int node, int parent, int level, vector<int> &lvl, vector<int> &low, stack<int>& nodeStack,
                 vector<vector<int>> &output);

};

#endif //MAIN_CPP_GRAPH_H
