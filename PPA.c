// Achilleas Ntalagiorgos Univesity of Ioannina
// a C Program to demonstrate adjacency list representation of graphs
// find heuristics as minimum distancies from the goal node  with reversed binary heap djkstra algorithm
// and then find pareto optimal set of best paths with BOA*

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <limits.h>
#include <math.h>

#define INF 999999999
#define LeftCHILD(x) 2 * x + 1
#define RightCHILD(x) 2 * x + 2
#define parent(x) (x - 1) / 2

//-------------------STRUCTURES----------------------------------------

// GRAPH STRUCTURES
// A structure to represent an adjacency list node
typedef struct AdjListVertex
{
    int dest;
    int cost1;
    int cost2;
    struct AdjListVertex *next;

} AdjListVertex;

// A structure to represent an adjacency list
typedef struct AdjList
{
    struct AdjListVertex *head;
} AdjList;

// A structure to represent a graph. A graph is an array of adjacency lists. Size of array will be V (number of vertices in graph)
typedef struct Graph
{
    int V;
    struct AdjList *array;

} Graph;

// HEAP  STRUCTURES
// Structure to represent a min heap node
typedef struct MinHeapNode
{
    int v;
    int dist;
} MinHeapNode;

// Structure to represent a min heap
typedef struct MinHeap
{
    // Number of heap nodes present currently
    int size;

    // Capacity of min heap
    int capacity;
    // This is needed for decreaseKey()
    int *pos;

    struct MinHeapNode **array;
} MinHeap;

//BOA STRUCTURES
//node structure
typedef struct BOA_node
{

    int vertex;
    int g_cost[2];
    int f_cost[2];
    int parent;

} BOA_node;

//minheap structure
typedef struct PPA_minHeap
{
    int size;
    struct PPA_node **elem;

} PPA_minHeap;

//PPA STRUCTURES
typedef struct PPA_node
{
    struct BOA_node *PPA_array[2];

} PPA_node;

//----------------------------------------------------------------------

//---------------- functions initiallizations----------------------
struct AdjListVertex *newAdjListNode(int dest, int cost1, int cost2);
struct Graph *createGraph(int V);
void addEdge(struct Graph *graph, int src, int dest, int cost1, int cost2);
void printGraph(struct Graph *graph);
void findHeuristic(struct Graph *lg);

struct MinHeapNode *newMinHeapNode(int v, int dist);
struct MinHeap *createMinHeap(int capacity);
void swapMinHeapNode(struct MinHeapNode **a, struct MinHeapNode **b);
void minHeapify(struct MinHeap *minHeap, int idx);
int isEmpty(struct MinHeap *minHeap);
struct MinHeapNode *extractMin(struct MinHeap *minHeap);
void decreaseKey(struct MinHeap *minHeap, int v, int dist);
int isInMinHeap(struct MinHeap *minHeap, int v);

void dijkstra(struct Graph *graph, int src, int heuristic, int dist[]);
void printArr(int dist[], int n);

void PPA(struct Graph *graph, int start, int goal, int dist1[], int dist2[]);
struct PPA_node *merge(PPA_node *nd1, PPA_node *nd2);
int isBounded(float e1, float e2, PPA_node *ppa_nd);
int isDominated_PPA(PPA_node *nd, int goal_node, int *g2min);
void insert_PPA(PPA_node *ppa_nd, PPA_minHeap *hp, float e1, float e2);
void merge_to_solutions_PPA(PPA_node *goal, PPA_minHeap *sols, float e1, float e2);
struct PPA_node *create_PPA_path(BOA_node *tl, BOA_node *br);
struct BOA_node *create_BOA_node(int vertex, int g[], int f[], int parent);
struct PPA_minHeap *initMinHeap();
void insertNode(struct PPA_minHeap *hp, struct PPA_node *path);
void swap(PPA_node *n1, PPA_node *n2);
void deleteNode(PPA_minHeap *hp, int i);
void heapify(PPA_minHeap *hp, int i);

void readfiles(char *file1, char *file2);

//---------------------------------------------------------------------

//-----------------------CREATE GRAPH FUNCTIONS-------------------------

// A utility function to create a new adjacency list node
struct AdjListVertex *newAdjListNode(int dest, int cost1, int cost2)
{
    struct AdjListVertex *newNode = (struct AdjListVertex *)malloc(sizeof(struct AdjListVertex));
    newNode->dest = dest;
    newNode->cost1 = cost1;
    newNode->cost2 = cost2;
    newNode->next = NULL;
    return newNode;
}

// A utility function that creates a graph of V vertices
struct Graph *createGraph(int V)
{
    struct Graph *graph = (struct Graph *)malloc(sizeof(struct Graph));
    graph->V = V;

    // Create an array of adjacency lists. Size of
    // array will be V
    graph->array = (struct AdjList *)malloc(V * sizeof(struct AdjList));

    // Initialize each adjacency list as empty by
    // making head as NULL
    int i;
    for (i = 0; i < V; ++i)
        graph->array[i].head = NULL;

    return graph;
}

// TODO :: DELETE -1 when i check papper grpah

// Adds an edge to a directed graph
void addEdge(struct Graph *graph, int src, int dest, int cost1, int cost2)
{
    // Add an edge from src to dest. A new node is added to the adjacency list of src. The node is added at the beginning
    struct AdjListVertex *newNode = newAdjListNode(dest - 1, cost1, cost2); // TODO :-1 here
    newNode->next = graph->array[src - 1].head;                             // TODO :-1 here
    graph->array[src - 1].head = newNode;                                   // TODO :-1 here
}

// TODO :: DELETE +1 when i check papper grpah

// A utility function to print the adjacency list
// representation of graph
void printGraph(struct Graph *graph)
{
    FILE *fp;

    fp = fopen("c:\\temp\\1.txt", "w");
    int v;
    for (v = 0; v < graph->V; ++v) //
    {
        struct AdjListVertex *pCrawl = graph->array[v].head;
        fprintf(fp, "\n Adjacency list of vertex %d\n  ", v + 1); // TODO :+1 here
        while (pCrawl)
        {
            fprintf(fp, " \tgoes to -> %d with cost1: %d and cost2: %d", pCrawl->dest + 1, pCrawl->cost1, pCrawl->cost2); // TODO :+1 here
            pCrawl = pCrawl->next;
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}

//-----------------------------------------------------------------------

//-----------------BINARY HEAP FUNCTIONS----------------------------------

// A utility function to create a new Min Heap Node
struct MinHeapNode *newMinHeapNode(int v, int dist)
{
    struct MinHeapNode *minHeapNode = (struct MinHeapNode *)malloc(sizeof(struct MinHeapNode));
    minHeapNode->v = v;
    minHeapNode->dist = dist;
    return minHeapNode;
}

// A utility function to create a Min Heap
struct MinHeap *createMinHeap(int capacity)
{
    struct MinHeap *minHeap = (struct MinHeap *)malloc(sizeof(struct MinHeap));
    minHeap->pos = (int *)malloc(capacity * sizeof(int));
    minHeap->size = 0;
    minHeap->capacity = capacity;
    minHeap->array = (struct MinHeapNode **)malloc(capacity * sizeof(struct MinHeapNode *));
    return minHeap;
}

// A utility function to swap two nodes of min heap.
// Needed for min heapify
void swapMinHeapNode(struct MinHeapNode **a, struct MinHeapNode **b)
{
    struct MinHeapNode *t = *a;
    *a = *b;
    *b = t;
}

// A standard function to heapify at given idx
// This function also updates position of nodes when they are swapped.
// Position is needed for decreaseKey()
void minHeapify(struct MinHeap *minHeap, int idx)
{
    int smallest, left, right;
    smallest = idx;
    left = 2 * idx + 1;
    right = 2 * idx + 2;

    if (left < minHeap->size && minHeap->array[left]->dist < minHeap->array[smallest]->dist)
        smallest = left;

    if (right < minHeap->size && minHeap->array[right]->dist < minHeap->array[smallest]->dist)
        smallest = right;

    if (smallest != idx)
    {
        // The nodes to be swapped in min heap
        struct MinHeapNode *smallestNode = minHeap->array[smallest];
        struct MinHeapNode *idxNode = minHeap->array[idx];

        // Swap positions
        minHeap->pos[smallestNode->v] = idx;
        minHeap->pos[idxNode->v] = smallest;

        // Swap nodes
        swapMinHeapNode(&minHeap->array[smallest], &minHeap->array[idx]);
        minHeapify(minHeap, smallest);
    }
}

// A utility function to check if the given minHeap is empty or not
int isEmpty(struct MinHeap *minHeap)
{
    return minHeap->size == 0;
}

// Standard function to extract minimum node from heap
struct MinHeapNode *extractMin(struct MinHeap *minHeap)
{
    if (isEmpty(minHeap))
        return NULL;

    // Store the root node
    struct MinHeapNode *root = minHeap->array[0];

    // Replace root node with last node
    struct MinHeapNode *lastNode = minHeap->array[minHeap->size - 1];
    minHeap->array[0] = lastNode;

    // Update position of last node
    minHeap->pos[root->v] = minHeap->size - 1;
    minHeap->pos[lastNode->v] = 0;

    // Reduce heap size and heapify root
    --minHeap->size;
    minHeapify(minHeap, 0);

    return root;
}

// Function to decrease/update dist value of a given vertex v.
// This function uses pos[] of min heap to get the current index of node in min heap
void decreaseKey(struct MinHeap *minHeap, int v, int dist)
{
    // Get the index of v in  heap array
    int i = minHeap->pos[v];

    // Get the node and update its dist value
    minHeap->array[i]->dist = dist;

    // Travel up while the complete
    // tree is not hepified.
    // This is a O(Logn) loop
    while (i && minHeap->array[i]->dist <
                    minHeap->array[(i - 1) / 2]->dist)
    {
        // Swap this node with its parent
        minHeap->pos[minHeap->array[i]->v] = (i - 1) / 2;
        minHeap->pos[minHeap->array[(i - 1) / 2]->v] = i;
        swapMinHeapNode(&minHeap->array[i], &minHeap->array[(i - 1) / 2]);
        // move to parent index
        i = (i - 1) / 2;
    }
}

// A utility function to check if a given vertex 'v' is in min heap or not
int isInMinHeap(struct MinHeap *minHeap, int v)
{
    if (minHeap->pos[v] < minHeap->size)
        return 1;
    return 0;
}

//------------------------------------------------------------------

//---------------------DJIKSTRA FUNCTION----------------------------
// The main function that calulates
// distances of shortest paths from src to all
// vertices. It is a O(ELogV) function
//---------------------DJIKSTRA FUNCTION----------------------------
// The main function that calulates
// distances of shortest paths from src to all
// vertices. It is a O(ELogV) function
void dijkstra(struct Graph *graph, int src, int heuristic, int dist[])
{
    // Get the number of vertices in graph
    int V = graph->V;
    // dist values used to pick
    // minimum weight edge in cut

    // minHeap represents set E
    struct MinHeap *minHeap = createMinHeap(V);

    // Initialize min heap with all vertices. dist value of all vertices
    for (int v = 0; v < V; ++v)
    {
        dist[v] = INT_MAX;
        minHeap->array[v] = newMinHeapNode(v, dist[v]);
        minHeap->pos[v] = v;
    }

    // Make dist value of src vertex as 0 so that it is extracted first
    minHeap->array[src]->v = src;
    minHeap->array[src]->dist = dist[src];
    dist[src] = 0;
    decreaseKey(minHeap, src, dist[src]);

    // Initially size of min heap is equal to V
    minHeap->size = V;

    // In the followin loop, minheap contains all node whose shortest distance is not yet finalized.
    if (heuristic == 1)
    {
        while (!isEmpty(minHeap))
        {
            // Extract the vertex with minimum distance value
            struct MinHeapNode *minHeapNode = extractMin(minHeap);

            // Store the extracted vertex number
            int u = minHeapNode->v;

            // Traverse through all adjacent vertices of u (the extracted vertex) and update their distance values
            struct AdjListVertex *pCrawl = graph->array[u].head;
            while (pCrawl != NULL)
            {
                int v = pCrawl->dest;

                // If shortest distance to v is not finalized yet, and distance to v through u is less than its previously calculated distance
                if (isInMinHeap(minHeap, v) && dist[u] != INT_MAX && pCrawl->cost1 + dist[u] < dist[v])
                {
                    dist[v] = dist[u] + pCrawl->cost1;

                    // update distance value in min heap also
                    decreaseKey(minHeap, v, dist[v]);
                }
                pCrawl = pCrawl->next;
            }
            free(minHeapNode);
        }
    }
    else
    {
        while (!isEmpty(minHeap))
        {
            // Extract the vertex with minimum distance value
            struct MinHeapNode *minHeapNode = extractMin(minHeap);

            // Store the extracted vertex number
            int u = minHeapNode->v;

            // Traverse through all adjacent vertices of u (the extracted vertex) and update their distance values
            struct AdjListVertex *pCrawl = graph->array[u].head;
            while (pCrawl != NULL)
            {
                int v = pCrawl->dest;

                // If shortest distance to v is not finalized yet, and distance to v through u is less than its previously calculated distance
                if (isInMinHeap(minHeap, v) && dist[u] != INT_MAX && pCrawl->cost2 + dist[u] < dist[v])
                {
                    dist[v] = dist[u] + pCrawl->cost2;

                    // update distance value in min heap also
                    decreaseKey(minHeap, v, dist[v]);
                }
                pCrawl = pCrawl->next;
            }
            free(minHeapNode);
        }
    }
    free(minHeap->array);
    free(minHeap->pos);
    free(minHeap);

    // print the calculated shortest distances
    //printArr(dist, V);
    //printf("\n\n");
}

// A utility function used to print the solution
void printArr(int dist[], int n)
{
    printf("Vertex   Distance from Source\n");
    for (int i = 0; i < n; ++i)
        printf("%d \t\t %d\n", i, dist[i]);
}

//------------------------------------------------------------------

//-------------------READ FILE FUNCTION------------------------------
//read file call the functions to create graph
//and then call djikstra
void readfiles(char *file1, char *file2)
{

    FILE *fptr1; //pointer for reading the file
    FILE *fptr2;

    char c1[1000];   // array to save the line of file
    char *p1;        // word token
    char *array1[4]; // array to save the intrested lines

    char c2[1000];   // array to save the line of file
    char *p2;        // word token
    char *array2[4]; // array to save the intrested lines

    char *saveptr1, *saveptr2;
    int i = 0;
    int j = 0;

    int start_vertex;
    int goal_vertex;
    int V;
    struct Graph *graph, *rev_graph;

    // open the file to read and check for error
    if ((fptr1 = fopen(file1, "r")) == NULL) //"USA-road-t.NY.gr"
    {
        printf("Error! opening file1");

        // Program exits if the file pointer returns NULL.
        exit(1);
    }

    if ((fptr2 = fopen(file2, "r")) == NULL) //"USA-road-t.NY.gr"
    {
        printf("Error! opening file2");

        // Program exits if the file pointer returns NULL.
        exit(1);
    }

    // read text line by line
    while (((c1[0] = fgetc(fptr1)) != EOF) && (c2[0] = fgetc(fptr2) != EOF))
    {
        //reads text until newline is encountered
        fscanf(fptr1, "%[^\n]", c1);
        fscanf(fptr2, "%[^\n]", c2);

        if ((c1[0] == 'p') && (c2[0] == 'p'))
        {

            //split the string into words
            p1 = strtok_r(c1, " ", &saveptr1);
            p2 = strtok_r(c2, " ", &saveptr2);

            while ((p1 != NULL) && (p2 != NULL))
            {

                array1[i++] = p1;
                array2[j++] = p2;

                p1 = strtok_r(NULL, " ", &saveptr2);
                p2 = strtok_r(NULL, " ", &saveptr1);
            }
            if (strcmp(array1[2], array2[2]) == 0)
                printf("the number of nodes in the graph is: %s\n\n", array2[2]);
            V = atoi(array2[2]);
            graph = createGraph(V);
            rev_graph = createGraph(V);
            break;
        }
    }

    int *dist1 = (int *)malloc(V * sizeof(int));
    int *dist2 = (int *)malloc(V * sizeof(int));

    while (((c1[0] = fgetc(fptr1)) != EOF) && (c2[0] = fgetc(fptr2) != EOF))
    {
        fscanf(fptr1, "%[^\n]", c1);
        fscanf(fptr2, "%[^\n]", c2);

        //check if the data is intrested
        if ((c1[0] == 'a') && (c2[0] == 'a'))
        {

            //reset the array
            i = 0;
            j = 0;

            //split the string into words
            p1 = strtok_r(c1, " ", &saveptr1);
            p2 = strtok_r(c2, " ", &saveptr2);

            while ((p1 != NULL) && (p2 != NULL))
            {

                array1[i++] = p1;
                array2[j++] = p2;

                p1 = strtok_r(NULL, " ", &saveptr1);
                p2 = strtok_r(NULL, " ", &saveptr2);
            }

            if ((strcmp(array1[1], array2[1]) == 0) && (strcmp(array1[2], array2[2]) == 0))
            {

                //create graph
                addEdge(graph, atoi(array1[1]), atoi(array1[2]), atoi(array1[3]), atoi(array2[3]));
                addEdge(rev_graph, atoi(array1[2]), atoi(array1[1]), atoi(array1[3]), atoi(array2[3]));
            }
        }
    }

    printGraph(graph);

    printf("please input starting vertex:\n");
    scanf("%d", &start_vertex);
    while (start_vertex < 1 || start_vertex > V)
    {
        printf("The start vert  ex must have value [1...%d]\n", V);
        printf("Please choose a new start vertex\n");
        scanf("%d", &start_vertex);
    }
    start_vertex--;

    printf("please input goal vertex:\n");
    scanf("%d", &goal_vertex);
    while (goal_vertex < 1 || goal_vertex > V)
    {
        printf("The goal vertex must have value [1...%d]\n", V);
        printf("Please choose a new goal vertex\n");
        scanf("%d", &goal_vertex);
    }
    goal_vertex--;

    //printGraph(graph);

    dijkstra(rev_graph, goal_vertex, 1, dist1);
    dijkstra(rev_graph, goal_vertex, 2, dist2);

    //printArr(dist1, V);
    //printArr(dist2, V);

    PPA(graph, start_vertex, goal_vertex, dist1, dist2);

    //rewind pointers to the start of the file's and close file's
    free(dist1);
    free(dist2);

    //TODO FREE ADJANCY LIST -> NECCESSARY

    free(graph->array);
    free(graph);
    free(rev_graph->array);
    free(rev_graph);

    rewind(fptr1);
    fclose(fptr1);
    rewind(fptr2);
    fclose(fptr2);
}
//--------------------------------------------------------------------

//--------------------- PPA* FUNCTIONS --------------------------------
// BOA MINHEAP FUNCTIONS
struct PPA_minHeap *initMinHeap()
{

    PPA_minHeap *hp = (struct PPA_minHeap *)malloc(sizeof(struct PPA_minHeap));
    hp->size = 0;
    return hp;
}

void insertNode(struct PPA_minHeap *hp, struct PPA_node *path)
{
    // allocating space

    if (hp->size)
    {
        hp->elem = realloc(hp->elem, (hp->size + 1) * sizeof(PPA_node));
    }
    else
    {
        hp->elem = malloc(sizeof(PPA_node));
    }

    int digit_num = 0;
    //Calculate total digits of f2 value
    digit_num = ((long long)(path->PPA_array[1]->f_cost[1]) == 0) ? 1 : (log10((long long)(path->PPA_array[1]->f_cost[1])) + 1);
    // Positioning the node at the right position in the min heap
    int i = (hp->size)++;

    while (i && (path->PPA_array[0]->f_cost[0] * pow((double)10, (double)digit_num) + path->PPA_array[1]->f_cost[1]) < (hp->elem[parent(i)]->PPA_array[0]->f_cost[0] * pow((double)10, (double)digit_num) + hp->elem[parent(i)]->PPA_array[1]->f_cost[1]))
    {
        hp->elem[i] = hp->elem[parent(i)];
        i = parent(i);
    }

    hp->elem[i] = path;
}

void deleteNode(PPA_minHeap *hp, int i)
{
    if (hp->size)
    {
        //printf("Deleting node %d\n\n", hp->elem[0]->vertex+1) ;
        hp->elem[i] = hp->elem[--(hp->size)];
        hp->elem = realloc(hp->elem, hp->size * sizeof(PPA_node));
        heapify(hp, i);
    }
    else
    {
        //printf("\nMin Heap is empty!\n") ;
        free(hp->elem);
    }
}

void heapify(PPA_minHeap *hp, int i)
{

    int digit_num = 0;
    //Calculate total digits
    if (LeftCHILD(i) < hp->size)
        digit_num = ((long long)hp->elem[LeftCHILD(i)]->PPA_array[1]->f_cost[1] == 0) ? 1 : (log10((long long)(hp->elem[LeftCHILD(i)]->PPA_array[1]->f_cost[1])) + 1);

    //&& hp->elem[LeftCHILD(i)].g_cost[1] < hp->elem[i]
    int smallest = (LeftCHILD(i) < hp->size && (hp->elem[LeftCHILD(i)]->PPA_array[0]->f_cost[0] * pow((double)10, (double)digit_num) + hp->elem[LeftCHILD(i)]->PPA_array[1]->f_cost[1]) < (hp->elem[i]->PPA_array[0]->f_cost[0] * pow((double)10, (double)digit_num) + hp->elem[i]->PPA_array[1]->f_cost[1])) ? LeftCHILD(i) : i;

    if (RightCHILD(i) < hp->size)
        digit_num = ((long long)hp->elem[RightCHILD(i)]->PPA_array[1]->f_cost[1] == 0) ? 1 : (log10((long long)(hp->elem[RightCHILD(i)]->PPA_array[1]->f_cost[1])) + 1);

    if ((RightCHILD(i) < hp->size) && ((hp->elem[RightCHILD(i)]->PPA_array[0]->f_cost[0] * pow((double)10, (double)digit_num) + hp->elem[RightCHILD(i)]->PPA_array[1]->f_cost[1]) < (hp->elem[smallest]->PPA_array[0]->f_cost[0] * pow((double)10, (double)digit_num) + hp->elem[smallest]->PPA_array[1]->f_cost[1])))
    {
        smallest = RightCHILD(i);
    }
    if (smallest != i)
    {
        swap(hp->elem[i], hp->elem[smallest]);
        heapify(hp, smallest);
    }
}

void swap(PPA_node *n1, PPA_node *n2)
{
    PPA_node temp = *n1;
    *n1 = *n2;
    *n2 = temp;
}

// a function to create BOA* node
struct BOA_node *create_BOA_node(int vertex, int g[], int f[], int parent)
{
    BOA_node *BOANode = (struct BOA_node *)malloc(sizeof(struct BOA_node));

    BOANode->vertex = vertex;
    BOANode->g_cost[0] = g[0];
    BOANode->g_cost[1] = g[1];
    BOANode->f_cost[0] = f[0];
    BOANode->f_cost[1] = f[1];
    BOANode->parent = parent;

    return BOANode;
}

// a function with input 2 BOA nodes and creates a PPA path
struct PPA_node *create_PPA_path(BOA_node *tl, BOA_node *br)
{
    PPA_node *PPApath = (struct PPA_node *)malloc(sizeof(struct PPA_node));

    PPApath->PPA_array[0] = tl;
    PPApath->PPA_array[1] = br;

    return PPApath;
}

//a function that shows if the PPApath is bounden in e1,e2
int isBounded(float e1, float e2, PPA_node *ppa_merged)
{
    //TODO DELETE
    //printf("ppa merged: tl: (%d, (%d, %d), (%d, %d), %d), br: (%d, (%d, %d), (%d, %d), %d)\n\n", ppa_merged->PPA_array[0]->vertex, ppa_merged->PPA_array[0]->g_cost[0], ppa_merged->PPA_array[0]->g_cost[1], ppa_merged->PPA_array[0]->f_cost[0], ppa_merged->PPA_array[0]->f_cost[1], ppa_merged->PPA_array[0]->parent, ppa_merged->PPA_array[1]->vertex, ppa_merged->PPA_array[1]->g_cost[0], ppa_merged->PPA_array[1]->g_cost[1], ppa_merged->PPA_array[1]->f_cost[0], ppa_merged->PPA_array[1]->f_cost[1], ppa_merged->PPA_array[1]->parent);

    //printf("cost1_br: %d \t cost1_tl %d\n", ppa_merged->PPA_array[1]->g_cost[0], ppa_merged->PPA_array[0]->g_cost[0]);
    //printf("e1: %d\n", (ppa_merged->PPA_array[1]->g_cost[0] - ppa_merged->PPA_array[0]->g_cost[0]) / ppa_merged->PPA_array[0]->g_cost[0]);

    if ((e1 >= ((double)(ppa_merged->PPA_array[1]->f_cost[0] - ppa_merged->PPA_array[0]->f_cost[0]) / (double)ppa_merged->PPA_array[0]->f_cost[0])) && (e2 >= ((double)(ppa_merged->PPA_array[0]->f_cost[1] - ppa_merged->PPA_array[1]->f_cost[1]) / (double)ppa_merged->PPA_array[1]->g_cost[1])))
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

// a function thaty merge 2 PPA paths
struct PPA_node *merge(PPA_node *nd1, PPA_node *nd2)
{

    PPA_node *new_nd = create_PPA_path(nd1->PPA_array[0], nd1->PPA_array[1]);

    if (nd1->PPA_array[0]->f_cost[0] <= nd2->PPA_array[0]->f_cost[0])
    {

        new_nd->PPA_array[0] = nd1->PPA_array[0];
    }
    else
    {

        new_nd->PPA_array[0] = nd2->PPA_array[0];
    }

    if (nd1->PPA_array[1]->f_cost[1] <= nd2->PPA_array[1]->f_cost[1])
    {

        new_nd->PPA_array[1] = nd1->PPA_array[1];
    }
    else
    {
        new_nd->PPA_array[1] = nd2->PPA_array[1];
    }

    return new_nd;
}

// main PPA* algorithm's functions
int isDominated_PPA(PPA_node *nd, int goal_node, int *g2min)
{

    if (nd->PPA_array[1]->f_cost[1] >= g2min[goal_node])
        return 1;

    if (nd->PPA_array[1]->g_cost[1] >= g2min[nd->PPA_array[1]->vertex])
    {
        return 1;
    }
    return 0;
}

void insert_PPA(PPA_node *ppa_nd, PPA_minHeap *hp, float e1, float e2)
{
    for (int i = 0; i < hp->size; i++)
    {
        PPA_node *PPvmerged = merge(hp->elem[i], ppa_nd);
        //TODO DELETE
        printf("PPvmerged : tl: (%d, (%d, %d), (%d, %d), %d), br: (%d, (%d, %d), (%d, %d), %d)\n\n", PPvmerged->PPA_array[0]->vertex + 1, PPvmerged->PPA_array[0]->g_cost[0], PPvmerged->PPA_array[0]->g_cost[1], PPvmerged->PPA_array[0]->f_cost[0], PPvmerged->PPA_array[0]->f_cost[1], PPvmerged->PPA_array[0]->parent + 1, PPvmerged->PPA_array[1]->vertex + 1, PPvmerged->PPA_array[1]->g_cost[0], PPvmerged->PPA_array[1]->g_cost[1], PPvmerged->PPA_array[1]->f_cost[0], PPvmerged->PPA_array[1]->f_cost[1], PPvmerged->PPA_array[1]->parent + 1);
        printf("the bound is %f, %f \n\n", ((double)(PPvmerged->PPA_array[1]->g_cost[0] - PPvmerged->PPA_array[0]->g_cost[0]) / (double)PPvmerged->PPA_array[0]->g_cost[0]), ((double)(PPvmerged->PPA_array[0]->g_cost[1] - PPvmerged->PPA_array[1]->g_cost[1]) / (double)PPvmerged->PPA_array[1]->g_cost[1]));

        if (isBounded(e1, e2, PPvmerged))
        {
            printf("The node -> PPvmerged : tl: (%d, (%d, %d), (%d, %d), %d), br: (%d, (%d, %d), (%d, %d), %d) is bounded \n", PPvmerged->PPA_array[0]->vertex + 1, PPvmerged->PPA_array[0]->g_cost[0], PPvmerged->PPA_array[0]->g_cost[1], PPvmerged->PPA_array[0]->f_cost[0], PPvmerged->PPA_array[0]->f_cost[1], PPvmerged->PPA_array[0]->parent + 1, PPvmerged->PPA_array[1]->vertex + 1, PPvmerged->PPA_array[1]->g_cost[0], PPvmerged->PPA_array[1]->g_cost[1], PPvmerged->PPA_array[1]->f_cost[0], PPvmerged->PPA_array[1]->f_cost[1], PPvmerged->PPA_array[1]->parent + 1);
            printf("the bound is e1>= %f and e2>= %f\n\n", (double)(PPvmerged->PPA_array[1]->g_cost[0] - PPvmerged->PPA_array[0]->g_cost[0]) / (double)PPvmerged->PPA_array[0]->g_cost[0], (double)(PPvmerged->PPA_array[0]->g_cost[1] - PPvmerged->PPA_array[1]->g_cost[1]) / (double)PPvmerged->PPA_array[1]->g_cost[1]);
            //free(ppa_nd);
            deleteNode(hp, i);
            insertNode(hp, PPvmerged);
            return;
        }
        free(PPvmerged);
    }

    //TODO DELETE
    //printf("ppa nd : tl: (%d, (%d, %d), (%d, %d), %d), br: (%d, (%d, %d), (%d, %d), %d)\n\n", ppa_nd->PPA_array[0]->vertex, ppa_nd->PPA_array[0]->g_cost[0], ppa_nd->PPA_array[0]->g_cost[1], ppa_nd->PPA_array[0]->f_cost[0], ppa_nd->PPA_array[0]->f_cost[1], ppa_nd->PPA_array[0]->parent, ppa_nd->PPA_array[1]->vertex, ppa_nd->PPA_array[1]->g_cost[0], ppa_nd->PPA_array[1]->g_cost[1], ppa_nd->PPA_array[1]->f_cost[0], ppa_nd->PPA_array[1]->f_cost[1], ppa_nd->PPA_array[1]->parent);
    insertNode(hp, ppa_nd);
    //TODO DELETE
    /*printf("META TO EXTEND\n");
    for (int i = 0; i < hp->size; i++)
    {
        printf("Openlist elem[%d]: tl: (%d, (%d, %d), (%d, %d), %d), br: (%d, (%d, %d), (%d, %d), %d)\n\n", i, hp->elem[i]->PPA_array[0]->vertex, hp->elem[i]->PPA_array[0]->g_cost[0], hp->elem[i]->PPA_array[0]->g_cost[1], hp->elem[i]->PPA_array[0]->f_cost[0], hp->elem[i]->PPA_array[0]->f_cost[1], hp->elem[i]->PPA_array[0]->parent, hp->elem[i]->PPA_array[1]->vertex, hp->elem[i]->PPA_array[1]->g_cost[0], hp->elem[i]->PPA_array[1]->g_cost[1], hp->elem[i]->PPA_array[1]->f_cost[0], hp->elem[i]->PPA_array[1]->f_cost[1], hp->elem[i]->PPA_array[1]->parent);
    }
    */
    return;
}

void merge_to_solutions_PPA(PPA_node *goal, PPA_minHeap *sols, float e1, float e2)
{
    for (int i = 0; i < sols->size; i++)
    {
        PPA_node *PPvgoalmerged = merge(sols->elem[i], goal);

        if (isBounded(e1, e2, PPvgoalmerged))
        {
            //printf("IM INNNNNNNN");
            free(goal);
            deleteNode(sols, i);
            insertNode(sols, PPvgoalmerged);
            return;
        }
        free(PPvgoalmerged);
    }

    insertNode(sols, goal);
}

void PPA(struct Graph *graph, int start, int goal, int dist1[], int dist2[])
{
    // number of graph verteces
    int V = graph->V;

    //initialize g_values for starting node
    int g_start[2] = {0, 0};

    //initialize f_values for starting node
    int f_start[2];
    f_start[0] = dist1[start];
    f_start[1] = dist2[start];

    // e1 e2 for dominance check
    double e1 = 0;
    double e2 = 0;

    // those variables used to create BOA nodes

    //parent
    int prnt;
    //g_cost
    int g_cost[2];
    //f_cost = g_cost + h
    int f_cost[2];

    PPA_minHeap *solutions_pp = initMinHeap();

    // initialize g2_min array
    int *g2_min = (int *)malloc(V * sizeof(int));
    for (int i = 0; i < V; i++)
        g2_min[i] = INF;

    // create starting BOA_ndoe
    BOA_node *node = create_BOA_node(start, g_start, f_start, start);
    PPA_node *path = create_PPA_path(node, node);

    //TODO DELETE
    //printf("tl: (%d, (%d, %d), (%d, %d), %d), br: (%d, (%d, %d), (%d, %d), %d)\n\n ", path->PPA_array[0]->vertex, path->PPA_array[0]->g_cost[0], path->PPA_array[0]->g_cost[1], path->PPA_array[0]->f_cost[0], path->PPA_array[0]->f_cost[1], path->PPA_array[0]->parent, path->PPA_array[1]->vertex, path->PPA_array[1]->g_cost[0], path->PPA_array[1]->g_cost[1], path->PPA_array[1]->f_cost[0], path->PPA_array[1]->f_cost[1], path->PPA_array[1]->parent);

    // create minheap OpenList
    PPA_minHeap *OpenList = initMinHeap();
    // insert starting node in OpenList

    insertNode(OpenList, path);
    int iter = 0;

    while (OpenList->size > 0)
    {
        //TODO DELETE 3 LINES
        printf("iter: %d\n", iter);
        iter++;
        printf("The expanded node is: tl: (%d, (%d, %d), (%d, %d), %d), br: (%d, (%d, %d), (%d, %d), %d)\n\n ", path->PPA_array[0]->vertex + 1, path->PPA_array[0]->g_cost[0], path->PPA_array[0]->g_cost[1], path->PPA_array[0]->f_cost[0], path->PPA_array[0]->f_cost[1], path->PPA_array[0]->parent + 1, path->PPA_array[1]->vertex + 1, path->PPA_array[1]->g_cost[0], path->PPA_array[1]->g_cost[1], path->PPA_array[1]->f_cost[0], path->PPA_array[1]->f_cost[1], path->PPA_array[1]->parent + 1);

        // delete the node with minimum f value
        deleteNode(OpenList, 0);

        if (isDominated_PPA(path, goal, g2_min))
        {
            //TODO DELETE
            printf("The path (tl: (%d, (%d, %d), (%d, %d), %d), br: (%d, (%d, %d), (%d, %d), %d) is dominated \n\n ", path->PPA_array[0]->vertex, path->PPA_array[0]->g_cost[0], path->PPA_array[0]->g_cost[1], path->PPA_array[0]->f_cost[0], path->PPA_array[0]->f_cost[1], path->PPA_array[0]->parent, path->PPA_array[1]->vertex, path->PPA_array[1]->g_cost[0], path->PPA_array[1]->g_cost[1], path->PPA_array[1]->f_cost[0], path->PPA_array[1]->f_cost[1], path->PPA_array[1]->parent);

            if (OpenList->size > 0)
            {
                free(path->PPA_array[0]);
                free(path->PPA_array[1]);
                free(path);
                path = OpenList->elem[0];
            }

            continue;
        }
        // else update vetrex's g2_min
        g2_min[path->PPA_array[0]->vertex] = path->PPA_array[0]->g_cost[1];
        g2_min[path->PPA_array[1]->vertex] = path->PPA_array[1]->g_cost[1];

        // if node's vertex is goal vertex add it to solutions_pp
        if (path->PPA_array[0]->vertex == goal || path->PPA_array[1]->vertex == goal)
        {
            //TODO DELETE
            printf("The path(tl: (%d, (%d, %d), (%d, %d), %d), br: (%d, (%d, %d), (%d, %d), %d) is a solution \n\n ", path->PPA_array[0]->vertex + 1, path->PPA_array[0]->g_cost[0], path->PPA_array[0]->g_cost[1], path->PPA_array[0]->f_cost[0], path->PPA_array[0]->f_cost[1], path->PPA_array[0]->parent + 1, path->PPA_array[1]->vertex + 1, path->PPA_array[1]->g_cost[0], path->PPA_array[1]->g_cost[1], path->PPA_array[1]->f_cost[0], path->PPA_array[1]->f_cost[1], path->PPA_array[1]->parent + 1);
            merge_to_solutions_PPA(path, solutions_pp, e1, e2);

            if (OpenList->size > 0)
                path = OpenList->elem[0];

            continue;
        }

        // pointer that shows in the beggining of the current's note adjacency list
        struct AdjListVertex *pCrawl0 = graph->array[path->PPA_array[0]->vertex].head;
        struct AdjListVertex *pCrawl1 = graph->array[path->PPA_array[1]->vertex].head;
        // for every vertex in node's vertex's adjancy list

        while (pCrawl0 != NULL && pCrawl1 != NULL)
        {

            //if the neighbor nodes are not dominates add it to the open list
            // crete new BOA nodes
            g_cost[0] = pCrawl0->cost1 + path->PPA_array[0]->g_cost[0];
            g_cost[1] = pCrawl0->cost2 + path->PPA_array[0]->g_cost[1];
            f_cost[0] = g_cost[0] + dist1[pCrawl0->dest];
            f_cost[1] = g_cost[1] + dist2[pCrawl0->dest];
            prnt = path->PPA_array[0]->vertex;

            BOA_node *tl_node = create_BOA_node(pCrawl0->dest, g_cost, f_cost, prnt);
            // printf("tl: (%d, (%d, %d), (%d, %d), %d)\n\n ", tl_node->vertex + 1, tl_node->g_cost[0], tl_node->g_cost[1], tl_node->f_cost[0], tl_node->f_cost[1], tl_node->parent);

            g_cost[0] = pCrawl1->cost1 + path->PPA_array[1]->g_cost[0];
            g_cost[1] = pCrawl1->cost2 + path->PPA_array[1]->g_cost[1];
            f_cost[0] = g_cost[0] + dist1[pCrawl1->dest];
            f_cost[1] = g_cost[1] + dist2[pCrawl1->dest];
            prnt = path->PPA_array[1]->vertex;

            BOA_node *br_node = create_BOA_node(pCrawl0->dest, g_cost, f_cost, prnt);
            //printf("br: (%d, (%d, %d), (%d, %d), %d)\n\n ", br_node->vertex + 1, br_node->g_cost[0], br_node->g_cost[1], br_node->f_cost[0], br_node->f_cost[1], br_node->parent);

            PPA_node *newpath = create_PPA_path(tl_node, br_node);

            // dominance check

            if (isDominated_PPA(newpath, goal, g2_min))
            {
                free(tl_node);
                free(br_node);
                free(newpath);

                pCrawl0 = pCrawl0->next;
                pCrawl1 = pCrawl1->next;

                continue;
            }
            //TODO DELETE
            printf("the new path (tl: (%d, (%d, %d), (%d, %d), %d), br: (%d, (%d, %d), (%d, %d), %d) \n\n ", newpath->PPA_array[0]->vertex + 1, newpath->PPA_array[0]->g_cost[0], newpath->PPA_array[0]->g_cost[1], newpath->PPA_array[0]->f_cost[0], newpath->PPA_array[0]->f_cost[1], newpath->PPA_array[0]->parent + 1, newpath->PPA_array[1]->vertex + 1, newpath->PPA_array[1]->g_cost[0], newpath->PPA_array[1]->g_cost[1], newpath->PPA_array[1]->f_cost[0], newpath->PPA_array[1]->f_cost[1], newpath->PPA_array[1]->parent + 1);

            //add them to the OpenList
            insert_PPA(newpath, OpenList, e1, e2);
            //TODO DELETE
            for (int i = 0; i < OpenList->size; i++)
            {
                printf("AFTER INS tl: (%d, (%d, %d), (%d, %d), %d), br: (%d, (%d, %d), (%d, %d), %d)\n\n ", OpenList->elem[i]->PPA_array[0]->vertex + 1, OpenList->elem[i]->PPA_array[0]->g_cost[0], OpenList->elem[i]->PPA_array[0]->g_cost[1], OpenList->elem[i]->PPA_array[0]->f_cost[0], OpenList->elem[i]->PPA_array[0]->f_cost[1], OpenList->elem[i]->PPA_array[0]->parent + 1, OpenList->elem[i]->PPA_array[1]->vertex + 1, OpenList->elem[i]->PPA_array[1]->g_cost[0], OpenList->elem[i]->PPA_array[1]->g_cost[1], OpenList->elem[i]->PPA_array[1]->f_cost[0], OpenList->elem[i]->PPA_array[1]->f_cost[1], OpenList->elem[i]->PPA_array[1]->parent + 1);
            }

            //next neighboors
            pCrawl0 = pCrawl0->next;
            pCrawl1 = pCrawl1->next;

            //printf("next hop->%d\n", pCrawl0->dest);
        }

        if (OpenList->size > 0)
        {
            //free(path->PPA_array[0]);
            //free(path->PPA_array[1]);
            free(path);
            path = OpenList->elem[0];
        }
    }

    free(g2_min);
    free(OpenList);

    BOA_node **solutions = (BOA_node **)malloc(solutions_pp->size * sizeof(BOA_node *));

    for (int i = 0; i < solutions_pp->size; i++)
    {
        solutions[i] = solutions_pp->elem[i]->PPA_array[0];
        printf("NODE:( %d, (%d,%d), (%d,%d), %d) is a solution\n", solutions[i]->vertex + 1, solutions[i]->g_cost[0], solutions[i]->g_cost[1], solutions[i]->f_cost[0], solutions[i]->f_cost[1], solutions[i]->parent + 1);
    }
    printf("The number of solutions is: %d\n\n", solutions_pp->size);
    free(solutions);
    for (int i = 0; i < solutions_pp->size; i++)
    {
        free(solutions_pp->elem[i]->PPA_array[0]);
        free(solutions_pp->elem[i]->PPA_array[1]);
        free(solutions_pp->elem[i]);
    }

    free(solutions_pp->elem);

    free(solutions_pp);
}

//--------------------------------------------------------------------

// main

//TODO: FIX THAT WITH READFILES
// C program to illustrate
// command line arguments
#include <stdio.h>

int main(int argc, char *argv[])
{

    char *file1 = argv[1];
    char *file2 = argv[2];

    clock_t begin = clock(); //to count time

    readfiles(file1, file2);

    clock_t end = clock(); //to count time

    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC; //find time in sec's

    printf("The time for my PPA* is: %f\n\n", time_spent); //print time

    return 0;
}
