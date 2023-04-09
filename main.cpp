

Algorithm: Dijkstra's Algorithm 

#include<bits/stdc++.h>
using namespace std;

#define pii pair<int,int>
vector<pii> adj[100005];
int dist[100005];

void Dijkstra(int src){
    priority_queue<pii, vector<pii>, greater<pii>> pq;
    pq.push({0,src});

    while(!pq.empty()){
        int d=pq.top().first;
        int u=pq.top().second;
        pq.pop();

        if(dist[u]!=-1) continue;
        dist[u]=d;

        for(auto v:adj[u]){
            if(dist[v.first]!=-1)continue;
            pq.push({v.second+d,v.first});
        }
    }
}

int main(){
    memset(dist,-1,sizeof dist);
    //Add edges and weights to the graph
    int src=1;
    Dijkstra(src);
    //Output shortest distance from vertex src to all other vertices
    return 0;
}

Quick Sort:

```
#include<iostream>
using namespace std;

int Partition(int A[], int start, int end){
    int pivot=A[end];
    int partitionIndex=start;
    for(int i=start;i<end;i++){
        if(A[i]<=pivot){
            swap(A[i],A[partitionIndex]);
            partitionIndex++;
        }
    }
    swap(A[partitionIndex],A[end]);
    return partitionIndex;
}

void QuickSort(int A[], int start, int end){
    if(start<end){
        int partitionIndex=Partition(A,start,end);
        QuickSort(A,start,partitionIndex-1);
        QuickSort(A,partitionIndex+1,end);
    }
}

int main(){
    int A[]={7,6,5,4,3,2,1,0};
    int n=sizeof(A)/sizeof(A[0]);
    QuickSort(A,0,n-1);
    for(int i=0;i<n;i++)
        cout<<A[i]<<" ";
    return 0;
}
```

Algorithm: Dijkstra's Algorithm 

#include<bits/stdc++.h> 
using namespace std; 

typedef pair<int,int> pii; 

const int INF=1e9; 
const int MAXN=1005; 

int n,m,s; 
int dist[MAXN]={0}; 
bool vis[MAXN]={false}; 
int G[MAXN][MAXN]; 

void Dijkstra(){ 
    fill(dist,dist+n+1,INF); 
    dist[s]=0; 
    priority_queue<pii,vector<pii>,greater<pii>> pq; 
    pq.push(make_pair(0,s)); 
    while(!pq.empty()){ 
        int u=pq.top().second; 
        pq.pop(); 
        if(vis[u]){ 
            continue; 
        } 
        vis[u]=true; 
        for(int v=1;v<=n;v++){ 
            if(!vis[v]&&G[u][v]>0&&dist[u]+G[u][v]<dist[v]){ 
                dist[v]=dist[u]+G[u][v]; 
                pq.push(make_pair(dist[v],v)); 
            } 
        } 
    } 
} 

int main(){ 
    scanf("%d%d%d",&n,&m,&s); 
    for(int i=0;i<=n;i++){ 
        for(int j=0;j<=n;j++){ 
            G[i][j]=INF; 
        } 
    } 
    for(int i=1;i<=m;i++){ 
        int u,v,w; 
        scanf("%d%d%d",&u,&v,&w); 
        G[u][v]=w; 
    } 
    Dijkstra(); 
    for(int i=1;i<=n;i++){ 
        printf("%d ",dist[i]); 
    } 
    return 0; 
} 

// Sample Input: 
// 3 3 1 
// 1 2 2 
// 2 3 3 
// 3 1 4 

// Sample Output: 
// 0 2 5

Algorithm: Merge Sort 

Code: 

#include <iostream> 
#include <vector> 

void mergeSort(std::vector<int>& arr, int l, int r); 
void merge(std::vector<int>& arr, int l, int m, int r); 

int main() { 
    std::vector<int> arr {8, 5, 3, 1, 9, 6, 0, 7, 4, 2}; 

    mergeSort(arr, 0, arr.size() - 1); 

    for (int num : arr) { 
        std::cout << num << " "; 
    } 
    std::cout << std::endl; 

    return 0; 
} 

void mergeSort(std::vector<int>& arr, int l, int r) { 
    if (l < r) { 
        int m = (l + r) / 2; 

        mergeSort(arr, l, m); 
        mergeSort(arr, m+1, r); 

        merge(arr, l, m, r); 
    } 
} 

void merge(std::vector<int>& arr, int l, int m, int r) { 
    int i, j, k; 
    int n1 = m - l + 1; 
    int n2 = r - m; 

    std::vector<int> L(n1); 
    std::vector<int> R(n2); 

    for (i = 0; i < n1; i++) { 
        L[i] = arr[l + i]; 
    } 
    for (j = 0; j < n2; j++) { 
        R[j] = arr[m + 1 + j]; 
    } 

    i = 0; 
    j = 0; 
    k = l; 

    while (i < n1 && j < n2) { 
        if (L[i] <= R[j]) { 
            arr[k] = L[i]; 
            i++; 
        } else { 
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

One interesting algorithm that a third year university student could implement to sharpen their skills is the Knapsack Problem. 

Here is an implementation in C++:

#include <iostream> 
using namespace std; 

int max(int a, int b) { return (a > b) ? a : b; }  

int knapSack(int W, int wt[], int val[], int n)  
{  
    if (n == 0 || W == 0)  
        return 0;  

    if (wt[n-1] > W)  
        return knapSack(W, wt, val, n-1);  
  
    else return max( val[n-1] + knapSack(W-wt[n-1], wt, val, n-1), knapSack(W, wt, val, n-1) );  
}  
  
int main()  
{  
    int val[] = {60, 100, 120};  
    int wt[] = {10, 20, 30};  
    int W = 50;  
    int n = sizeof(val)/sizeof(val[0]);  
    cout<< knapSack(W, wt, val, n);  
    return 0;  
}  

Output: 220

Algorithm: Dijkstra's Algorithm

Code:

#include<bits/stdc++.h>
using namespace std;
#define pii pair<int,int>
#define pb push_back
const int N=1e5+5;
int n,m,x,y,c,dis[N];
bool vis[N];
vector<pii> g[N];
priority_queue< pii,vector<pii>,greater<pii> > q;
void dijkstra(int s)
{
    for(int i=0;i<=n;i++) dis[i]=1e9+5; //initializes the distances to infinity
    dis[s]=0;
    q.push({0,s});
    while(!q.empty())
    {
        x=q.top().second;
        q.pop();
        if(vis[x]) continue;
        vis[x]=true;
        for(int i=0;i<g[x].size();i++)
        {
            y=g[x][i].first;
            c=g[x][i].second;
            if(dis[y]>dis[x]+c)
            {
                dis[y]=dis[x]+c;
                q.push({dis[y],y});
            }
        }
    }
}
int main()
{
    cin>>n>>m;
    for(int i=1;i<=m;i++)
    {
        cin>>x>>y>>c;
        g[x].pb({y,c});
        g[y].pb({x,c}); //comment these two lines for directed graph
    }
    dijkstra(1);
    for(int i=1;i<=n;i++) cout<<dis[i]<<" ";
    return 0;
}

Here's an example of a code snippet implementing the Bubble Sort algorithm: 

```
#include <iostream>
using namespace std;

void bubbleSort(int arr[], int n) {
    for (int i = 0; i < n-1; i++) {
        for (int j = 0; j < n-i-1; j++) {
            if (arr[j] > arr[j+1]) {
                swap(arr[j], arr[j+1]);
            }
        }
    }
}

int main() {
    int arr[] = {64, 25, 12, 22, 11};
    int n = sizeof(arr)/sizeof(arr[0]);
    bubbleSort(arr, n);
    cout << "Sorted array: ";
    for (int i = 0; i < n; i++) {
        cout << arr[i] << " ";
    }
    return 0;
}
``` 

This code implements the Bubble Sort algorithm to sort an array of integers in ascending order. The main function declares an array of integers and calls the bubbleSort function to sort the array. The bubbleSort function implements the nested for-loops to compare adjacent elements in the array and swap them if the left element is greater than the right element. The sorted array is outputted at the end of the program.


Algorithm: Kosaraju's algorithm for finding strongly connected components in a directed graph.

Code:

#include<bits/stdc++.h>
using namespace std;

void dfs1(vector<int> graph[], int u, vector<bool> &visited, stack<int> &st){
    visited[u] = true;
    for(int v : graph[u]){
        if(!visited[v]){
            dfs1(graph, v, visited, st);
        }
    }
    st.push(u);
}

void dfs2(vector<int> reverse_graph[], int u, vector<bool> &visited){
    visited[u] = true;
    cout<<u<<" ";
    for(int v : reverse_graph[u]){
        if(!visited[v]){
            dfs2(reverse_graph, v, visited);
        }
    }
}

void kosarajuSCC(vector<int> graph[], vector<int> reverse_graph[], int n){
    vector<bool> visited(n, false);
    stack<int> st;

    for(int i=0;i<n;i++){
        if(!visited[i]){
            dfs1(graph, i, visited, st);
        }
    }

    visited = vector<bool>(n, false);
    while(!st.empty()){
        int u = st.top();
        st.pop();
        if(!visited[u]){
            dfs2(reverse_graph, u, visited);
            cout<<"\n";
        }
    }
}

int main(){
    int n, m;
    cin>>n>>m;

    vector<int> graph[n], reverse_graph[n];

    for(int i=0;i<m;i++){
        int u, v;
        cin>>u>>v;
        graph[u].push_back(v);
        reverse_graph[v].push_back(u);
    }

    kosarajuSCC(graph, reverse_graph, n);

    return 0;
}

One interesting algorithm that a third year university student could implement is the Merge Sort algorithm. Here's an example C++ code snippet:

```
#include <iostream>
using namespace std;

void merge(int arr[], int left, int middle, int right) {
    int n1 = middle - left + 1;
    int n2 = right - middle;

    int L[n1], R[n2];

    for (int i = 0; i < n1; i++) {
        L[i] = arr[left + i];
    }

    for (int j = 0; j < n2; j++) {
        R[j] = arr[middle + 1 + j];
    }

    int i = 0, j = 0, k = left;

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

void mergeSort(int arr[], int left, int right) {
    if (left < right) {
        int middle = left + (right - left) / 2;

        mergeSort(arr, left, middle);
        mergeSort(arr, middle + 1, right);

        merge(arr, left, middle, right);
    }
}

int main() {
    int arr[] = {12, 11, 13, 5, 6, 7};
    int arr_size = sizeof(arr)/sizeof(arr[0]);

    cout << "Original Array: ";
    for (int i = 0; i < arr_size; i++) {
        cout << arr[i] << " ";
    }
    cout << endl;

    mergeSort(arr, 0, arr_size - 1);

    cout << "Sorted Array: ";
    for (int i = 0; i < arr_size; i++) {
        cout << arr[i] << " ";
    }
    cout << endl;

    return 0;
}
```Algorithm: Sieve of Eratosthenes

Code:

#include <iostream>
#include <vector>
using namespace std;

void sieve(int n) {
    vector<bool> primes(n+1, true);

    for(int p=2; p*p<=n; p++) {
        if(primes[p] == true) {
            for(int i=p*p; i<=n; i+=p)
                primes[i] = false;
        }
    }

    for(int p=2; p<=n; p++) {
        if(primes[p])
            cout << p << " ";
    }
}

int main() {
    int n = 30;
    cout << "Prime numbers smaller than or equal to " << n << ":\n";
    sieve(n);
    return 0;
}Algorithm: Merge Sort

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
}Algorithm: Merge Sort

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

    int i = 0, j = 0, k = l;

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
        int m = l + (r - l) / 2;

        mergeSort(arr, l, m);
        mergeSort(arr, m + 1, r);

        merge(arr, l, m, r);
    }
}

int main() {
    int arr[] = { 64, 25, 12, 22, 11 };
    int n = sizeof(arr) / sizeof(arr[0]);

    mergeSort(arr, 0, n - 1);

    for (int i = 0; i < n; i++)
        cout << arr[i] << " ";

    return 0;
}Algorithm: Dijkstra's Algorithm

#include <iostream>
#include <vector>
#include <queue>
using namespace std;

#define INF 0x3f3f3f3f 

typedef pair<int, int> pii;

vector<pair<int, int>> adj[100001];
int dist[100001];

void dijkstra(int s) {
    memset(dist, INF, sizeof(dist));
    dist[s] = 0;
    priority_queue<pii, vector<pii>, greater<pii>> pq;
    pq.push(make_pair(0, s));

    while(!pq.empty()) {
        int u = pq.top().second;
        pq.pop();
        for(auto v : adj[u]) {
            int w = v.first, weight = v.second;
            if(dist[u] + weight < dist[w]) {
                dist[w] = dist[u] + weight;
                pq.push(make_pair(dist[w], w));
            }
        }
    }
}

int main() {
    int n, m, s;
    cin >> n >> m >> s;

    for(int i=0; i<m; i++) {
        int u, v, w;
        cin >> u >> v >> w;
        adj[u].push_back(make_pair(v, w));
        adj[v].push_back(make_pair(u, w)); // If graph is undirected
    }

    dijkstra(s);

    for(int i=1; i<=n; i++) {
        cout << "Distance from " << s << " to " << i << " is " << dist[i] << endl;
    }

    return 0;
}Algorithm: QuickSort

#include<iostream> 
using namespace std; 

void swap(int* a, int* b) 
{ 
    int t = *a; 
    *a = *b; 
    *b = t; 
} 

int partition (int arr[], int low, int high) 
{ 
    int pivot = arr[high]; 
    int i = (low - 1); 
  
    for (int j = low; j <= high- 1; j++) 
    { 
        if (arr[j] <= pivot) 
        { 
            i++;  
            swap(&arr[i], &arr[j]); 
        } 
    } 
    swap(&arr[i + 1], &arr[high]); 
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

void printArray(int arr[], int size) 
{ 
    for (int i=0; i < size; i++) 
        cout<<arr[i]<<" "; 
    cout<<"\n"; 
} 

int main() 
{ 
    int arr[] = {10, 7, 8, 9, 1, 5}; 
    int n = sizeof(arr)/sizeof(arr[0]); 
    quickSort(arr, 0, n-1); 
    cout<<"Sorted array: \n"; 
    printArray(arr, n); 
    return 0; 
}Mergesort:

```
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
        R[j] = arr[m + 1 + j];

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
        int m = l + (r - l) / 2;
        mergeSort(arr, l, m);
        mergeSort(arr, m + 1, r);

        merge(arr, l, m, r);
    }
}

int main() {
    int arr[] = { 12, 11, 13, 5, 6, 7 };
    int arr_size = sizeof(arr) / sizeof(arr[0]);

    cout << "Given array is \n";
    for (int i = 0; i < arr_size; i++)
        cout << arr[i] << " ";

    mergeSort(arr, 0, arr_size - 1);

    cout << "\nSorted array is \n";
    for (int i = 0; i < arr_size; i++)
        cout << arr[i] << " ";
    return 0;
}
```Algorithm: Merge Sort

Code:

```
#include <iostream> 
using namespace std; 

void merge(int arr[], int l, int m, int r) 
{ 
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

void mergeSort(int arr[], int l, int r) 
{ 
    if (l < r) { 
        int m = l + (r - l) / 2; 

        mergeSort(arr, l, m); 
        mergeSort(arr, m + 1, r); 

        merge(arr, l, m, r); 
    } 
} 

void printArray(int arr[], int size) 
{ 
    for (int i = 0; i < size; i++) 
        cout << arr[i] << " "; 
    cout << endl; 
} 

int main() 
{ 
    int arr[] = { 12, 11, 13, 5, 6, 7 }; 
    int arr_size = sizeof(arr) / sizeof(arr[0]); 

    cout << "Given array is \n"; 
    printArray(arr, arr_size); 

    mergeSort(arr, 0, arr_size - 1); 

    cout << "\nSorted array is \n"; 
    printArray(arr, arr_size); 
    return 0; 
} 
```Algorithm: Quick Sort

// C++ implementation of QuickSort
#include <iostream>
using namespace std;

// A utility function to swap two elements
void swap(int* a, int* b)
{
    int t = *a;
    *a = *b;
    *b = t;
}

/* This function takes last element as pivot, places
   the pivot element at its correct position in sorted
    array, and places all smaller (smaller than pivot)
   to left of pivot and all greater elements to right
   of pivot */
int partition (int arr[], int low, int high)
{
    int pivot = arr[high];    // pivot
    int i = (low - 1);  // Index of smaller element
 
    for (int j = low; j <= high- 1; j++)
    {
        // If current element is smaller than the pivot
        if (arr[j] < pivot)
        {
            i++;    // increment index of smaller element
            swap(&arr[i], &arr[j]);
        }
    }
    swap(&arr[i + 1], &arr[high]);
    return (i + 1);
}
 
/* The main function that implements QuickSort
 arr[] --> Array to be sorted,
  low  --> Starting index,
  high  --> Ending index */
void quickSort(int arr[], int low, int high)
{
    if (low < high)
    {
        /* pi is partitioning index, arr[p] is now
           at right place */
        int pi = partition(arr, low, high);
 
        // Separately sort elements before
        // partition and after partition
        quickSort(arr, low, pi - 1);
        quickSort(arr, pi + 1, high);
    }
}
 
/* Function to print an array */
void printArray(int arr[], int size)
{
    int i;
    for (i=0; i < size; i++)
        cout << arr[i] << " ";
    cout << endl;
}
 
// Driver program to test above functions
int main()
{
    int arr[] = {10, 7, 8, 9, 1, 5};
    int n = sizeof(arr)/sizeof(arr[0]);
    quickSort(arr, 0, n-1);
    cout << "Sorted array: \n";
    printArray(arr, n);
    return 0;
}Algorithm: Sieve of Eratosthenes

Code:

#include <iostream>
#include <bitset>

using namespace std;

int main()
{
    const int N = 1000000;
    bitset<N+1> prime;
    prime.set();
    prime[0] = prime[1] = false;

    for(int i=2;i*i<=N;++i)
        if(prime[i])
            for(int j=i*i;j<=N;j+=i)
                prime[j] = false;

    for(int i=1;i<=N;++i)
        if(prime[i])
            cout << i << " ";

    return 0;
}

Algorithm: Insertion Sort

#include <iostream

// Knapsack problem 
#include <bits/stdc

Implementation of Selection Sort Algorithm:
int main() 


Knapsack Problem Algorithm: 

#include<i

//Quick Sort Algorithm
#include <iostream>


Knuth-Morris-Pratt Algorithm:

#

Knapsack Algorithm:

int knapsack(int 