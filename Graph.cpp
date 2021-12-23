#include "Graph.h"

Graph::Graph(const int nrNodes, bool oriented) {
    this->nrNodes = nrNodes;
    this->oriented = oriented;
    edges.resize(nrNodes + 1);
}
void Graph::addEdge(const int x, const int y) {
    this->edges[x].push_back(y);
    if(!oriented) {
        this->edges[y].push_back(x);
    }
}
int Graph::getNumberOfNodes() {
    return this->nrNodes;
}
vector<int> Graph::BFS(int startingNode) {
    vector<int> distances(this->nrNodes + 1, -1);
    vector<int> visited(this->nrNodes+1, 0);
    queue<int> queue;

    distances[startingNode] = 0;
    visited[startingNode] = 1;
    queue.push(startingNode);

    while(!queue.empty()){
        int currentNode = queue.front();
        for(auto node : edges[currentNode]){
            if(!visited[node]){
                visited[node] = 1;
                distances[node] = distances[currentNode] + 1;
                queue.push(node);
            }
        }
        queue.pop();
    }
    return distances;
}
vector<int> Graph::DFS(int node) {
    vector<int> visited(nrNodes+1, 0);
    dfsAux(node, visited);
    return visited;
}
void Graph::dfsAux(int currentNode, vector<int> &visited) {
    for(auto node : this->edges[currentNode]) {
        if(!visited[node]) {
            visited[node] = 1;
            dfsAux(node, visited);
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
            dfsAux(i,visited);
        }
    }
    return connectedComp;
}

int Graph::getGraphDiameter() {
    vector<int> distance = BFS(1); // bfs from endNode 1
    int maxDistance = -1;
    int endNode;
    for(int i = 1; i < distance.size(); i++) { // the maximum distance to a endNode
        if(distance[i] > maxDistance)
        {
            maxDistance = distance[i];
            endNode = i;
        }
    }
    maxDistance = -1;
    distance = BFS(endNode); // bfs from the endNode with the maximum distance
    for(int i = 1; i< distance.size(); i++) {
        if(distance[i] > maxDistance)
        {
            maxDistance = distance[i];
        }
    }
    return maxDistance;
}
void Graph::topologicalSort(int currentNode, vector<int> &visited, vector<int>& topoSort) {
    visited[currentNode] = 1;
    for(auto node: this->edges[currentNode]){
        if(!visited[node])
            topologicalSort(node, visited, topoSort);
    }
    topoSort.push_back(currentNode);
}
vector<int> Graph::getTopologicalSortedGraph() {
    vector<int> topoSort;
    vector<int>visited(this->nrNodes+1, 0);
    for(int i = 1; i<= nrNodes; i++){
        if(!visited[i])
            topologicalSort(i, visited, topoSort); // sort in topologcal order
    }
    return topoSort;
}
void Graph::DFS_strongly_conn_comp(int currentNode, vector<int> &visited, int currentComponent, vector<vector<int>> & solution, vector<vector<int>>& transposedGraph) { // dfs used for strongly conn comp
    visited[currentNode] = 1;
    solution[currentComponent].push_back(currentNode);
    for(auto i: transposedGraph[currentNode]){
        if(visited[i] == 0){
            DFS_strongly_conn_comp(i,visited,currentComponent,solution, transposedGraph);
        }
    }
}
vector<vector<int>> Graph::stronglyConnectedComponents(){
    vector<vector<int>>solution;
    vector<vector<int>> transposedGraph;
    transposedGraph.resize(this->nrNodes+1);
    for(int i = 1; i<= this->nrNodes; i++){             // compute transposed graph
        for(auto j : edges[i]){
            transposedGraph[j].push_back(i);
        }
    }

    vector<int> v = getTopologicalSortedGraph();
    stack<int> topoSort;

    for(auto i : v)
    {
        topoSort.push(i);
    }
    vector<int>visited2(this->nrNodes+1, 0);

    solution.push_back(vector<int>());
    int currentComponent = 0;
    while(!topoSort.empty()) {
        int i = topoSort.top();
        topoSort.pop();
        if(visited2[i] == 0) {
            currentComponent++;
            solution.push_back(vector<int>());
            DFS_strongly_conn_comp(i, visited2, currentComponent,solution, transposedGraph);
        }
    }
    return solution;
}
bool Graph::havelHakim(vector<int> grades) {
    sort(grades.begin(), grades.end(), greater<>());
    if(grades[0] > grades.size()-1)
        return false;
    while (true){
        int gr = grades[0];
        grades.erase(grades.begin());
        for(int& i: grades){
            gr--;
            i--;
            if(i < 0) return false;
            if(gr == 0) break;
        }
        sort(grades.begin(), grades.end(), greater<>());
        if(grades.empty() || grades[0] == 0){
            return true;
        }
    }
}

vector<pair<int,int>> Graph::getCriticalConnections() {
    vector<int> lvl(this->nrNodes+1, 0);
    vector<int> low(this->nrNodes+ 1, 1);
    vector<pair<int,int>> output(this->nrNodes+1);
    DFS_crit(1,-1,1,lvl, low,output);
    return output;

}
void Graph::DFS_crit(int node, int parent, int level, vector<int> &lvl, vector<int> &low,
                     vector<pair<int, int>> &output) {
    lvl[node] = level;
    low[node] = level;
    for(auto i : this->edges[node]){
        if(lvl[i] == 0) {
            DFS_crit(i,node, level+1, lvl, low, output);
            low[node] = min(low[node], low[i]);
        }
        else if(lvl[i] != 0 && i != parent)
            low[node] = min(low[node], lvl[i]);
    }
    if(low[node] == lvl[node] && node != 0){
        pair<int,int> edge;
        edge.first = node;
        edge.second = parent;
        output.push_back(edge);
    }

}
vector<vector<int>> Graph::biconnectedComponents() {
    vector<int> lvl(this->nrNodes+1, 0);
    vector<int> low(this->nrNodes+ 1, 0);
    vector<vector<int>> output;
    stack<int> nodeStack;
    int level = 1, parent = 0;
    for(int i = 1; i<= this->nrNodes; i++) {
        if(lvl[i] == 0) {
            DFS_bcc(i,parent, level,lvl,low,nodeStack,output);
        }
    }
    return output;
}
void Graph::DFS_bcc(int node, int parent, int level, vector<int> &lvl, vector<int> &low,stack<int>& nodeStack,
                    vector<vector<int>> &output) {
    lvl[node] = level;
    low[node] = level;
    nodeStack.push(node);
    for(auto i : this->edges[node]){
        if(lvl[i] == 0 ){
            DFS_bcc(i,node,level+1,lvl, low,nodeStack, output);
            low[node] =  min(low[node], low[i]);
            if(low[i] >= lvl[node]) {
                vector<int> biconnectedComponent;
                biconnectedComponent.push_back(node);
                while(nodeStack.top() != i ){
                    biconnectedComponent.push_back(nodeStack.top());
                    nodeStack.pop();
                }
                biconnectedComponent.push_back(i);
                nodeStack.pop();
                output.push_back(biconnectedComponent);
            }
        }
        else if(i != parent)
            low[node] = min(low[node], lvl[i]);
    }

}
Graph::~Graph(){
    edges.clear();
}