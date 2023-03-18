

Algorithm: Merge Sort

Code:

#include <iostream>
using namespace std;

void merge(int arr[], int l, int m, int r) {
    int n1 = m - l + 1;
    int n2 = r - m;

    int L[n1], R[n2];

    for (int i = 0; i < n1; i++)
        L[i] = arr[l + i];
    for (int j = 0; j < n2; j++)
        R[j] = arr[m + 1 + j];

    int i = 0;
    int j = 0;
    int k = l;
    while (i < n1 && j < n2) {
        if (L[i] <= R[j]) {
            arr[k] = L[i];
            i++;
        }
        else {
            arr[k] = R[j];
            j++;
        }
        k++;
    }

    while (i < n1) {
        arr[k] = L[i];
        i++;
        k++;
    }

    while (j < n2) {
        arr[k] = R[j];
        j++;
        k++;
    }
}

void mergeSort(int arr[], int l, int r) {
    if (l >= r)
        return;

    int m = l + (r - l) / 2;
    mergeSort(arr, l, m);
    mergeSort(arr, m + 1, r);
    merge(arr, l, m, r);
}

void printArray(int arr[], int size) {
    for (int i = 0; i < size; i++)
        cout << arr[i] << " ";
    cout << endl;
}

int main() {
    int arr[] = { 12, 11, 13, 5, 6, 7 };
    int arr_size = sizeof(arr) / sizeof(arr[0]);

    cout << "Given array is \n";
    printArray(arr, arr_size);

    mergeSort(arr, 0, arr_size - 1);

    cout << "\nSorted array is \n";
    printArray(arr, arr_size);
    return 0;
}

Algorithm: Dijkstra's Shortest Path Algorithm

#include<bits/stdc++.h>
using namespace std;

#define INF 0x3f3f3f3f

typedef pair<int,int> iPair;

class Graph
{
    int V;
    list< pair<int,int> > *adj;

    public:
        Graph(int V);
        void addEdge(int u,int v,int w);
        void shortestPath(int s);
};

Graph::Graph(int V)
{
    this->V=V;
    adj=new list<iPair> [V];
}

void Graph::addEdge(int u,int v,int w)
{
    adj[u].push_back(make_pair(v,w));
    adj[v].push_back(make_pair(u,w));
}

void Graph::shortestPath(int src)
{
    priority_queue< iPair,vector<iPair>,greater<iPair> > pq;
    vector<int> dist(V,INF);

    pq.push(make_pair(0,src));
    dist[src]=0;

    while(!pq.empty())
    {
        int u=pq.top().second;
        pq.pop();

        list<pair<int,int > >::iterator i;
        for(i=adj[u].begin();i!=adj[u].end();++i)
        {
            int v=(*i).first;
            int weight=(*i).second;

            if(dist[v]>dist[u]+weight)
            {
                dist[v]=dist[u]+weight;
                pq.push(make_pair(dist[v],v));
            }

        }
    }
    cout<<"Vertex   Distance from Source"<<endl;
    for(int i=0;i<V;i++)
        cout<<i<<"\t\t"<<dist[i]<<endl;
}

int main()
{

    Graph g(9);

    g.addEdge(0, 1, 4);
    g.addEdge(0, 7, 8);
    g.addEdge(1, 2, 8);
    g.addEdge(1, 7, 11);
    g.addEdge(2, 3, 7);
    g.addEdge(2, 8, 2);
    g.addEdge(2, 5, 4);
    g.addEdge(3, 4, 9);
    g.addEdge(3, 5, 14);
    g.addEdge(4, 5, 10);
    g.addEdge(5, 6, 2);
    g.addEdge(6, 7, 1);
    g.addEdge(6, 8, 6);
    g.addEdge(7, 8, 7);

    g.shortestPath(0);

    return 0;
}

Algorithm: Binary Search

int binarySearch(vector<int> arr, int l, int r, int x) {
   if (r >= l) {
      int mid = l + (r - l)/2;
      if (arr[mid] == x)
         return mid;
      if (arr[mid] > x)
         return binarySearch(arr, l, mid-1, x);
      return binarySearch(arr, mid+1, r, x);
   }
   return -1;
}

Name: Merge Sort Algorithm

Code:

#include <iostream>
using namespace std;

void merge(int arr[], int l, int m, int r) {
    int i, j, k;
    int n1 = m - l + 1;
    int n2 = r - m;

    int L[n1], R[n2];
    
    for (i = 0; i < n1; i++)
        L[i] = arr[l + i];
    for (j = 0; j < n2; j++)
        R[j] = arr[m + 1+ j];

    i = 0; 
    j = 0; 
    k = l; 
    while (i < n1 && j < n2) {
        if (L[i] <= R[j]) {
            arr[k] = L[i];
            i++;
        }
        else {
            arr[k] = R[j];
            j++;
        }
        k++;
    }

    while (i < n1) {
        arr[k] = L[i];
        i++;
        k++;
    }

    while (j < n2) {
        arr[k] = R[j];
        j++;
        k++;
    }
}

void mergeSort(int arr[], int l, int r) {
    if (l < r) {
        int m = l+(r-l)/2;

        mergeSort(arr, l, m);
        mergeSort(arr, m+1, r);

        merge(arr, l, m, r);
    }
}

int main() 
{
    int arr[] = { 12, 11, 13, 5, 6, 7 };
    int arr_size = sizeof(arr)/sizeof(arr[0]);

    cout << "Given array is \n";
    for (int i = 0; i < arr_size; i++)
        cout << arr[i] << " ";

    mergeSort(arr, 0, arr_size - 1);

    cout << "\nSorted array is \n";
    for (int i = 0; i < arr_size; i++)
        cout << arr[i] << " ";

    return 0;
}




One interesting algorithm that a third year university student could implement to sharpen their skills is the QuickSort algorithm. 

Here is an example code snippet in C++:

#include <iostream>
#include <vector>

using namespace std;

void quickSort(vector<int>& arr, int low, int high) {
  if (low < high) {
    int pivot = arr[(low+high)/2];
    int i = low - 1;
    int j = high + 1;
    while (true) {
      do {
        i++;
      } while (arr[i] < pivot);
      do {
        j--;
      } while (arr[j] > pivot);
      if (i >= j) {
        break;
      }
      swap(arr[i], arr[j]);
    }
    quickSort(arr, low, j);
    quickSort(arr, j+1, high);
  }
}

int main() {
  vector<int> arr = {3, 1, 5, 2, 7, 4};
  quickSort(arr, 0, arr.size()-1);
  for (int i = 0; i < arr.size(); i++) {
    cout << arr[i] << " ";
  }
  cout << endl;
  return 0;
} 

This code takes an input vector of integers, and uses the QuickSort algorithm to sort the elements in non-decreasing order. The algorithm works by selecting a pivot element (in this case, we choose the middle element of the array), and partitioning the elements around it such that elements less than the pivot come before it, and elements greater than the pivot come after it. This process is repeated recursively on the two subarrays created by the partition until they are sorted.

Algorithm: Floyd-Warshall algorithm 

Code: 

#include<bits/stdc++.h>
#define ll long long 
#define N 510
using namespace std;
const ll inf=1ll<<60;
ll n,m,dis[N][N],ans;
int main(){
    scanf("%lld%lld",&n,&m);
    for(int i=1;i<=n;i++)
        for(int j=1;j<=n;j++)
            dis[i][j]=i==j?0:inf;
    for(int i=1,u,v,w;i<=m;i++){
        scanf("%d%d%d",&u,&v,&w);
        dis[u][v]=min(dis[u][v],1ll*w);
    }
    for(int k=1;k<=n;k++)
        for(int i=1;i<=n;i++)
            for(int j=1;j<=n;j++)
                dis[i][j]=min(dis[i][j],dis[i][k]+dis[k][j]);
    for(int i=1;i<=n;i++)
        for(int j=1;j<=n;j++)
            ans+=dis[i][j]<inf?dis[i][j]:0;
    printf("%lld\n",ans);
    return 0;
}



Quick Sort: 

#include<bits/stdc++.h>
using namespace std;

int partition(int arr[], int low, int high)
{
    int pivot = arr[high];    
    int i = (low - 1);  
    for (int j = low; j <= high- 1; j++)
    {
        if (arr[j] <= pivot)
        {
            i++;    
            swap(arr[i], arr[j]);
        }
    }
    swap(arr[i + 1], arr[high]);
    return (i + 1);
}

void quickSort(int arr[], int low, int high)
{
    if (low < high)
    {
        int pi = partition(arr, low, high);
        quickSort(arr, low, pi - 1);
        quickSort(arr, pi + 1, high);
    }
}

int main()
{
    int arr[] = {10, 7, 8, 9, 1, 5};  
    int n = sizeof(arr)/sizeof(arr[0]);
    quickSort(arr, 0, n-1);
    cout<<"Sorted array: \n";
    for (int i = 0; i < n; i++)
        cout<<arr[i]<<" ";
    return 0;
}

Algorithm: Prim's Minimum Spanning Tree

void primMST(int graph[V][V]) 
{ 
    int parent[V]; 
    int key[V]; 
    bool mstSet[V]; 

    for (int i = 0; i < V; i++) 
    { 
        key[i] = INT_MAX; 
        mstSet[i] = false; 
    } 
  
    key[0] = 0;     
    parent[0] = -1;  

    for (int count = 0; count < V - 1; count++) 
    { 
        
        int u = minKey(key, mstSet); 

        mstSet[u] = true; 

        for (int v = 0; v < V; v++) 

            if (graph[u][v] && mstSet[v] == false && 
                graph[u][v] < key[v]) 
            { 
                parent[v] = u; 
                key[v] = graph[u][v]; 
            } 
    } 
} 


Algorithm: QuickSort

void quickSort(int arr[], int low, int high) {
    if (low < high) {
        int pivot = arr[high];
        int i = low - 1;
        for (int j = low; j < high; j++) {
            if (arr[j] < pivot) {
                i++;
                int temp = arr[i];
                arr[i] = arr[j];
                arr[j] = temp;
            }
        }
        int temp = arr[i + 1];
        arr[i + 1] = arr[high];
        arr[high] = temp;
        int pi = i + 1;
        quickSort(arr, low, pi - 1);
        quickSort(arr, pi + 1, high);
    }
}