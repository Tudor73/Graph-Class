#ifndef MAIN_CPP_WEIGHTEDGRAPH_H
#define MAIN_CPP_WEIGHTEDGRAPH_H
#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <fstream>
#include <algorithm>
#include <tuple>
#include "Graph.h"

using namespace std;

struct edge {
    int node, cost;
    edge(int node , int cost ) {
        this->node = node;
        this->cost = cost;
    }
};

class WeightedGraph : public Graph
{
private:
    vector<vector<edge>> weigthedEdges;
public:
    WeightedGraph(const int nrNodes, bool oriented = true);
    void addEdge(const int x, const int y,const int cost);
    vector<int> Dijkstra (int node);
    vector<int> BellmanFord(int node);
    vector<int> Prim(int node, int& totalCost);
    void disjoint(vector<tuple<int,int,int>>, vector<int>& roots, ostream& g);

    vector<vector<edge>> FloydWarshall();
private:
    int findRoot(int node, vector<int> &roots);
};

#endif //MAIN_CPP_WEIGHTEDGRAPH_H
