#include "WeightedGraph.h"

const int NMAX = numeric_limits<int>::max();

WeightedGraph::WeightedGraph(const int nrNodes, bool oriented): Graph(nrNodes, oriented) {
    this->weigthedEdges.resize(nrNodes+1);
}

void WeightedGraph::addEdge(const int x,const int y,const int cost) {
    this->weigthedEdges[x].push_back(edge(y,cost));
    if(!oriented) {
        this->weigthedEdges[y].push_back(edge(x,cost));
    }
}

vector<int> WeightedGraph::Dijkstra(int node) {
    vector<int> minDistance(nrNodes+1,NMAX);
    vector<bool> visited(nrNodes+1, false);
    priority_queue<pair<int,int>, vector<pair<int,int>>, greater<pair<int,int>>> pq;
    minDistance[node] = 0;
    pq.push(make_pair(0,node));
    while(!pq.empty()) {
        int current_node = pq.top().second;
        pq.pop();
        if(!visited[current_node]) {
            visited[current_node] = true;
            for(auto i: this->weigthedEdges[current_node]){
                if(minDistance[current_node] + i.cost < minDistance[i.node]){
                    minDistance[i.node] = minDistance[current_node] + i.cost;
                    pq.push(make_pair(minDistance[i.node], i.node));
                }
            }
        }
    }

    return minDistance;
}
vector<int> WeightedGraph::BellmanFord(int node) {
    vector<int> minDistances(nrNodes+1, NMAX);
    vector<int> visited(nrNodes+1, 0);
    vector<bool> inQueue(nrNodes + 1, false);
    queue<int> q;
    q.push(node);
    inQueue[node] = true;
    minDistances[node] = 0;
    while(!q.empty()){
        int current_node = q.front();
        visited[current_node]++;
        q.pop();
        inQueue[current_node] = false;
        if(visited[current_node] >= nrNodes)
            return {};
        for(auto i :this->weigthedEdges[current_node])
        {
            if(minDistances[current_node] + i.cost < minDistances[i.node]){
                minDistances[i.node] = minDistances[current_node] + i.cost;

                if(!inQueue[i.node]){
                    q.push(i.node);
                    inQueue[i.node] = true;
                }
            }
        }
    }
    return minDistances;
}

vector<int> WeightedGraph::Prim(int node, int &totalCost) {
    vector<int> minDistance(nrNodes+1, NMAX);
    vector<int> parent(nrNodes+1, 0);
    vector<bool> selected(nrNodes+1, false);
    priority_queue<pair<int, int>,vector<pair<int, int>> ,greater<>>pq;
    minDistance[node] = 0;
    pq.push(make_pair(0,1));
    while(!pq.empty())
    {
        int currentNode = pq.top().second;
        int currentCost = pq.top().first;
        pq.pop();
        if(!selected[currentNode])
        {
            selected[currentNode] = true;
            totalCost= totalCost + currentCost;
            for(auto i : this->weigthedEdges[currentNode]){
                if(!selected[i.node])
                {
                    if(i.cost < minDistance[i.node])
                    {
                        minDistance[i.node] = i.cost;
                        parent[i.node] = currentNode;
                        pq.push(make_pair(minDistance[i.node], i.node));
                    }
                }
            }
        }
    }
    return parent;
}

int WeightedGraph::findRoot(int node, vector<int>& roots) {
    if(roots[node] != node)
        return roots[node] = findRoot(roots[node], roots);
    return roots[node] ;
}

void WeightedGraph::disjoint(vector<tuple<int,int,int>> v, vector<int>& roots, ostream& g) {
    for(auto i : v) {
        int code = get<0>(i);
        int firstNode = get<1>(i);
        int secondNode = get<2>(i);
        if(code == 1) {
            int firstNodeRoot = findRoot(firstNode, roots);
            int secondNodeRoot = findRoot(secondNode, roots);
            roots[firstNodeRoot] = secondNodeRoot;  // append the root of the first tree to the root of the cost tree
        }
        else if (code == 2)
        {
            if(findRoot(firstNode, roots) == findRoot(secondNode, roots)){
                g << "DA\n";
            }
            else g << "NU\n";
        }
    }
}

vector<vector<edge>> WeightedGraph::FloydWarshall() {
    vector<vector<edge>> distances = this->weigthedEdges;
    for(int k = 0; k < this->nrNodes; k++)
        for(int i = 0; i< this->nrNodes; i++)
            for(int j = 0; j< this->nrNodes; j++) {
                if( i != j && distances[k][j].cost && distances[i][k].cost )
                    if(distances[i][j].cost > distances[i][k].cost + distances[k][j].cost || distances[i][j].cost == 0){
                        distances[i][j].cost = distances[i][k].cost + distances[k][j].cost;
                    }
            }
    return distances;
}
