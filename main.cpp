#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <fstream>
#include <algorithm>
#include <tuple>
using namespace std;

const int NMAX = numeric_limits<int>::max();

class Graph {
protected:
    int nrNodes, nrEdges;
    bool oriented;
    vector<vector<int>> edges;
public:
    Graph(int nrNodes, bool oriented = true) ;
    void addEdge(int x, int y); // sets edge between x and y
    int getNumberOfNodes();
    void BFS(int startingNode, ostream& g);
    int connectedComponents();  // returns the number of conencted components from the graph
    void DFS(int node, vector<int>& visited);
    vector<int> getTopologicalSortedGraph();
    vector<vector<int>> stronglyConnectedComponents();   // outputs the strongly connected components
    bool havelHakim(vector<int> grades); //  returns true or false if the grades array form a valid graph
    vector<vector<int>> biconnectedComponents(); // returns the array of biconnected components
    vector<pair<int,int>> getCriticalConnections(); // return array of critical edge
    ~Graph();
private:
    void topologicalSort(int x, vector<int> &visited, vector<int>& topoSort);
    void DFS_strongly_conn_comp(int currentNode, vector<int> &component, int currentComponent, vector<vector<int>> & solution, vector<vector<int>>& transposedGraph); // used in stronglyConnComponents
    void DFS_crit(int node, int predecesor, int level, vector<int>& lvl,vector<int>& low, vector<pair<int,int>>& output ); // used in getCriticalConnections
    void DFS_bcc(int node, int parent, int level, vector<int> &lvl, vector<int> &low, stack<int>& nodeStack,
                 vector<vector<int>> &output);

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
int Graph::getNumberOfNodes() {
    return this->nrNodes;
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
//        cout << i << " ";
        topoSort.push(i);
    }
    cout << endl;
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
class WeightedGraph : public Graph
{
private:
    vector<vector<pair<int,int>>> weigthedEdges;
public:
    WeightedGraph(int nrNodes, bool oriented = true);
    void addEdge(int x, int y, int cost);
    vector<int> Dijkstra (int node);
    vector<int> BellmanFord(int node);
    vector<int> Prim(int node, int& totalCost);
    void disjoint(vector<tuple<int,int,int>>, vector<int>& roots, ostream& g);
private:
    int findRoot(int node, vector<int> &roots);
};

WeightedGraph::WeightedGraph(int nrNodes, bool oriented): Graph(nrNodes, oriented) {
    this->weigthedEdges.resize(nrNodes+1);
}

void WeightedGraph::addEdge(int x, int y, int cost) {
    this->weigthedEdges[x].push_back(make_pair(y,cost));
    if(!oriented) {
        this->weigthedEdges[y].push_back(make_pair(x,cost));
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
                if(minDistance[current_node] + i.second < minDistance[i.first]){ // i.second is the cost and i.first the node
                    minDistance[i.first] = minDistance[current_node] + i.second;
                    pq.push(make_pair(minDistance[i.first], i.first));
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
            if(minDistances[current_node] + i.second < minDistances[i.first]){ // i.second is the cost and i.first the node
                minDistances[i.first] = minDistances[current_node] + i.second;

                if(!inQueue[i.first]){
                    q.push(i.first);
                    inQueue[i.first] = true;
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
    pq.push(make_pair(0,1)); // first is the cost and second the node
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
                if(!selected[i.first])
                {
                    if(i.second < minDistance[i.first]) // first is the node and second is the cost in the edges array
                    {
                        minDistance[i.first] = i.second;
                        parent[i.first] = currentNode;
                        pq.push(make_pair(minDistance[i.first], i.first));
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
            roots[firstNodeRoot] = secondNodeRoot;  // append the root of the first tree to the root of the second tree
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
int main()
{
    int nrNodes, nrEdges;
    ifstream f("ctc.in");
    ofstream g("ctc.out");
    vector<tuple<int,int,int>> v;
    f >> nrNodes >> nrEdges;
    int x ,y ,z;
    Graph G(nrNodes);
    for(int i = 0; i < nrEdges; i++)
    {
        f >> x >> y;
        G.addEdge(x,y);
    }
    vector<vector<int>> solution = G.stronglyConnectedComponents();
    cout << solution.size() -1 ;
    for(auto i : solution)
    {
        for(auto j : i) {
            cout << j << " ";
        }
        cout << endl;

    }
}

