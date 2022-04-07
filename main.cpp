/*
graph.h

Visible structs and functions for graph construction and manipulation.

Skeleton written by Grady Fitzpatrick for COMP20007 Assignment 1 2021
*/

/* Definition of a graph. */
struct graph;

struct solution;

/* A particular solution to a graph problem. */
#ifndef SOLUTION_STRUCT
#define SOLUTION_STRUCT
struct solution {
    int heartsLost;
};
#endif

/* Which part the program should find a solution for. */
#ifndef PART_ENUM
#define PART_ENUM
enum problemPart {
    PART_A = 0,
    PART_B = 1,
    PART_C = 2
};
#endif

/* Creates an undirected graph with the given numVertices and no edges and
returns a pointer to it. NumEdges is the number of expected edges. */
struct graph *newGraph(int numVertices);

/* Adds an edge to the given graph. */
void addEdge(struct graph *g, int start, int end, int cost);

/* Find the number of hearts which will be left if Lonk takes the optimal path.
  In all parts the graph will be the graph formed by regular rooms in the
  dungeon.
  numRooms is the number of rooms in the dungeon.
  Lonk starts in the startingRoom and the boss is in the bossRoom.
  For PART_B, the numShortcuts denotes the number of shortcuts,
    shortcutStarts represents the room on one side of each shortcut which
    Lonk's key could be used to open, and shortuctEnds represents the other side.
  For PART_A and PART_C, numShortcuts will be 0, shortcutStarts and shortcutEnds
    will be NULL.
  For PART_C, numHeartRooms denotes the number of heart rooms, heartRooms is the
    list of rooms which contain hearts.
  For PART_A and PART_B, numHeartRooms will be 0 and heartRooms will be NULL.
 */
struct solution *graphSolve(struct graph *g, enum problemPart part,
                            int numRooms, int startingRoom, int bossRoom, int numShortcuts,
                            int *shortcutStarts, int *shortcutEnds, int numHeartRooms, int *heartRooms);

/* Returns a new graph which is a deep copy of the given graph (which must be
  freed with freeGraph when no longer used). */
struct graph *duplicateGraph(struct graph *g);

/* Frees all memory used by graph. */
void freeGraph(struct graph *g);

/* Frees all data used by solution. */
void freeSolution(struct solution *solution);
/*
list.h

Visible structs and functions for linked lists.

Skeleton written by Grady Fitzpatrick for COMP20007 Assignment 1 2022
*/
/* The linked list. */
struct list;

/* Get a new empty list. */
struct list *newlist(void *item);

/* Add an item to the head of the list. Returns the new list. */
struct list *prependList(struct list *list, void *item);

/* Gets the first item from the list. */
void *peekHead(struct list *list);

/* Takes the first item from the list, updating the list pointer and returns
  the item stored. */
void *deleteHead(struct list **list);

/* Free all list items. */
void freeList(struct list *list);
/*
pq.h

Visible structs and functions for priority queues.

Skeleton written by Grady Fitzpatrick for COMP20007 Assignment 1 2022
*/
/* The priority queue. */
struct pq;

/* Get a new empty priority queue. */
struct pq *newPQ();

/* Add an item to the priority queue - cast pointer to (void *). */
void enqueue(struct pq *pq, void *item, int priority);

/* Take the smallest item from the priority queue - cast pointer back to
  original type. */
void *deletemin(struct pq *pq);

/* Returns 1 if empty, 0 otherwise. */
int empty(struct pq *pq);

/* Remove all items from priority queue (doesn't free) and free the queue. */
void freePQ(struct pq *pq);
/*
utils.h

Visible structs and functions for helper functions to do with reading and
writing.

Skeleton written by Grady Fitzpatrick for COMP20007 Assignment 1 2022
*/
/* Because we use FILE in this file, we should include stdio.h here. */
#include <stdio.h>
/* Because we use struct graph in this file, we should include graph.h here. */
/* The problem specified. */
struct graphProblem;

/* Reads the data from the given file pointer and returns a pointer to this
information. */
struct graphProblem *readProblem(FILE *file, enum problemPart part);

/* Finds a solution for a given problem. */
struct solution *findSolution(struct graphProblem *problem,
                              enum problemPart part);

/* Frees all data used by problem. */
void freeProblem(struct graphProblem *problem);
/*
graph.c

Set of vertices and edges implementation.

Implementations for helper functions for graph construction and manipulation.

Skeleton written by Grady Fitzpatrick for COMP20007 Assignment 1 2022
*/

#include <assert.h>
#include <limits.h>
#include <stdlib.h>

#define INITIALEDGES 32
/*
list.c

Implementations for helper functions for linked list construction and
manipulation.

Skeleton written by Grady Fitzpatrick for COMP20007 Assignment 1 2022
*/

#include <assert.h>
#include <stdlib.h>
/*
utils.c

Implementations for helper functions to do with reading and writing.

Skeleton written by Grady Fitzpatrick for COMP20007 Assignment 1 2022
*/
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

/* Default cost for edges. */
#define DEFAULTCOST 1

struct graphProblem {
    int numRooms;
    int numConnections;
    int startRoom;
    int bossRoom;
    int numShortcuts;
    int *shortcutStart;
    int *shortcutEnd;
    struct graph *graph;
    int numHeartRooms;
    int *heartRooms;
};

struct graphProblem *readProblem(FILE *file, enum problemPart part) {
    int i;
    int startRoom;
    int endRoom;
    /* Allocate space for problem specification */
    struct graphProblem *problem = (struct graphProblem *)
        malloc(sizeof(struct graphProblem));
    assert(problem);

    /* First line of input is the number of rooms. */
    assert(fscanf(file, "%d", &(problem->numRooms)) == 1);
    /* Next line comprises number of connections between rooms. */
    assert(fscanf(file, "%d", &(problem->numConnections)) == 1);
    /* Next line comprises the start room. */
    assert(fscanf(file, "%d", &(problem->startRoom)) == 1);
    /* Next line comprises the boss room. */
    assert(fscanf(file, "%d", &(problem->bossRoom)) == 1);

    /* Build graph number of rooms. */
    problem->graph = newGraph(problem->numRooms);
    /* Add all edges to graph. */
    for (i = 0; i < problem->numConnections; i++) {
        assert(fscanf(file, "%d %d", &startRoom, &endRoom) == 2);
        addEdge(problem->graph, startRoom, endRoom, DEFAULTCOST);
    }

    /* Handle PART_B */
    if (part == PART_B) {
        /* Read number of shortcuts */
        assert(fscanf(file, "%d", &(problem->numShortcuts)) == 1);
        problem->shortcutStart = (int *)malloc(sizeof(int) * problem->numShortcuts);
        assert(problem->shortcutStart || problem->numShortcuts == 0);
        problem->shortcutEnd = (int *)malloc(sizeof(int) * problem->numShortcuts);
        assert(problem->shortcutEnd || problem->numShortcuts == 0);
        for (i = 0; i < problem->numShortcuts; i++) {
            /* Read each shortcut connection. */
            assert(fscanf(file, "%d %d", &startRoom, &endRoom) == 2);
            (problem->shortcutStart)[i] = startRoom;
            (problem->shortcutEnd)[i] = endRoom;
        }
    } else {
        problem->shortcutStart = NULL;
        problem->shortcutEnd = NULL;
        problem->numShortcuts = 0;
    }

    /* Handle PART_C */
    if (part == PART_C) {
        /* Read number of heart rooms */
        assert(fscanf(file, "%d", &(problem->numHeartRooms)) == 1);
        problem->heartRooms = (int *)malloc(sizeof(int) * problem->numHeartRooms);
        assert(problem->heartRooms || problem->numHeartRooms == 0);
        for (i = 0; i < problem->numHeartRooms; i++) {
            /* Read each heart room. */
            assert(fscanf(file, "%d", &(problem->heartRooms[i])) == 1);
        }
    } else {
        problem->numHeartRooms = 0;
        problem->heartRooms = NULL;
    }

    return problem;
}

struct solution *findSolution(struct graphProblem *problem,
                              enum problemPart part) {
    return graphSolve(problem->graph, part, problem->numRooms,
                      problem->startRoom, problem->bossRoom, problem->numShortcuts,
                      problem->shortcutStart, problem->shortcutEnd, problem->numHeartRooms,
                      problem->heartRooms);
}

void freeProblem(struct graphProblem *problem) {
    /* No need to free if no data allocated. */
    if (!problem) {
        return;
    }
    freeGraph(problem->graph);
    if (problem->shortcutStart) {
        free(problem->shortcutStart);
    }
    if (problem->shortcutEnd) {
        free(problem->shortcutEnd);
    }
    if (problem->heartRooms) {
        free(problem->heartRooms);
    }

    free(problem);
}

void freeSolution(struct solution *solution) {
    /* No need to free if no data allocated. */
    if (!solution) {
        return;
    }
    free(solution);
}

/*
pq.c

Unsorted Array Implementation

Implementations for helper functions for priority queue construction and
manipulation.

Skeleton written by Grady Fitzpatrick for COMP20007 Assignment 1 2022
*/
#include <assert.h>
#include <stdlib.h>

#define INITIALITEMS 32

struct pq {
    int count;
    int allocated;
    void **queue;
    int *priorities;
};

struct pq *newPQ() {
    struct pq *pq = (struct pq *)malloc(sizeof(struct pq));
    assert(pq);
    pq->count = 0;
    pq->allocated = 0;
    pq->queue = NULL;
    pq->priorities = NULL;
    return pq;
}

void enqueue(struct pq *pq, void *item, int priority) {
    assert(pq);
    if ((pq->count + 1) > pq->allocated) {
        if (pq->allocated == 0) {
            pq->allocated = INITIALITEMS;
        } else {
            pq->allocated *= 2;
        }
        pq->queue = (void **)realloc(pq->queue, pq->allocated * sizeof(void *));
        assert(pq->queue);
        pq->priorities = (int *)realloc(pq->priorities, pq->allocated *
                                                            sizeof(int));
        assert(pq->priorities);
    }
    (pq->queue)[pq->count] = item;
    (pq->priorities)[pq->count] = priority;
    (pq->count)++;
}

/* Scan through all the priorities linearly and find lowest. */
void *deletemin(struct pq *pq) {
    int i;
    int lowestElement = 0;
    void *returnVal;
    if (pq->count <= 0) {
        return NULL;
    }
    for (i = 0; i < pq->count; i++) {
        if ((pq->priorities)[i] < (pq->priorities)[lowestElement]) {
            lowestElement = i;
        }
    }
    returnVal = (pq->queue)[lowestElement];
    /* Delete item from queue by swapping final item into place of deleted
      element. */
    if (pq->count > 0) {
        (pq->priorities)[lowestElement] = (pq->priorities)[pq->count - 1];
        (pq->queue)[lowestElement] = (pq->queue)[pq->count - 1];
        (pq->count)--;
    }
    return returnVal;
}

int empty(struct pq *pq) {
    return pq->count == 0;
}

void freePQ(struct pq *pq) {
    if (!pq) {
        return;
    }
    if (pq->allocated > 0) {
        free(pq->queue);
        free(pq->priorities);
    }
    free(pq);
}

struct list {
    void *item;
    struct list *next;
};

struct list *newlist(void *item) {
    struct list *head = (struct list *)malloc(sizeof(struct list));
    assert(head);
    head->item = item;
    head->next = NULL;
    return head;
}

struct list *prependList(struct list *list, void *item) {
    struct list *head = (struct list *)malloc(sizeof(struct list));
    assert(head);
    head->item = item;
    head->next = list;
    return head;
}

void *peekHead(struct list *list) {
    if (!list) {
        return NULL;
    }
    return list->item;
}

void *deleteHead(struct list **list) {
    void *item;
    struct list *next;
    if (!list || !*list) {
        return NULL;
    }
    /* Store values we're interested in before freeing list node. */
    item = (*list)->item;
    next = (*list)->next;
    free(*list);
    *list = next;
    return item;
}

void freeList(struct list *list) {
    struct list *next;
    /* Iterate through list until the end of the list (NULL) is reached. */
    for (next = list; list != NULL; list = next) {
        /* Store next pointer before we free list's space. */
        next = list->next;
        free(list);
    }
}

struct edge;

/* Definition of a graph. */
struct graph {
    int numVertices;
    int numEdges;
    int allocedEdges;
    struct edge **edgeList;
};

/* Definition of an edge. */
struct edge {
    int start;
    int end;
    int cost;
};

struct graph *newGraph(int numVertices) {
    struct graph *g = (struct graph *)malloc(sizeof(struct graph));
    assert(g);
    /* Initialise edges. */
    g->numVertices = numVertices;
    g->numEdges = 0;
    g->allocedEdges = 0;
    g->edgeList = NULL;
    return g;
}

/* Adds an edge to the given graph. */
void addEdge(struct graph *g, int start, int end, int cost) {
    assert(g);
    struct edge *newEdge = NULL;
    /* Check we have enough space for the new edge. */
    if ((g->numEdges + 1) > g->allocedEdges) {
        if (g->allocedEdges == 0) {
            g->allocedEdges = INITIALEDGES;
        } else {
            (g->allocedEdges) *= 2;
        }
        g->edgeList = (struct edge **)realloc(g->edgeList,
                                              sizeof(struct edge *) * g->allocedEdges);
        assert(g->edgeList);
    }

    /* Create the edge */
    newEdge = (struct edge *)malloc(sizeof(struct edge));
    assert(newEdge);
    newEdge->start = start;
    newEdge->end = end;
    newEdge->cost = cost;

    /* Add the edge to the list of edges. */
    g->edgeList[g->numEdges] = newEdge;
    (g->numEdges)++;
}

/* Returns a new graph which is a deep copy of the given graph (which must be
  freed with freeGraph when no longer used). */
struct graph *duplicateGraph(struct graph *g) {
    struct graph *copyGraph = (struct graph *)malloc(sizeof(struct graph));
    assert(copyGraph);
    copyGraph->numVertices = g->numVertices;
    copyGraph->numEdges = g->numEdges;
    copyGraph->allocedEdges = g->allocedEdges;
    copyGraph->edgeList = (struct edge **)malloc(sizeof(struct edge *) * g->allocedEdges);
    assert(copyGraph->edgeList || copyGraph->numEdges == 0);
    int i;
    /* Copy edge list. */
    for (i = 0; i < g->numEdges; i++) {
        struct edge *newEdge = (struct edge *)malloc(sizeof(struct edge));
        assert(newEdge);
        newEdge->start = (g->edgeList)[i]->start;
        newEdge->end = (g->edgeList)[i]->end;
        newEdge->cost = (g->edgeList)[i]->cost;
        (copyGraph->edgeList)[i] = newEdge;
    }
    return copyGraph;
}

/* Frees all memory used by graph. */
void freeGraph(struct graph *g) {
    int i;
    for (i = 0; i < g->numEdges; i++) {
        free((g->edgeList)[i]);
    }
    if (g->edgeList) {
        free(g->edgeList);
    }
    free(g);
}

// A linked list (LL) node to store a queue entry
struct Queue {
    int front, rear, size;
    unsigned capacity;
    int *array;
};

// function to create a queue
// of given capacity.
// It initializes size of queue as 0
struct Queue *createQueue(unsigned capacity) {
    struct Queue *queue = (struct Queue *)malloc(
        sizeof(struct Queue));
    queue->capacity = capacity;
    queue->front = queue->size = 0;

    // This is important, see the enqueue
    queue->rear = capacity - 1;
    queue->array = (int *)malloc(
        queue->capacity * sizeof(int));
    return queue;
}

// Queue is full when size becomes
// equal to the capacity
int isFull(struct Queue *queue) {
    return (queue->size == queue->capacity);
}

// Queue is empty when size is 0
int isEmpty(struct Queue *queue) {
    return (queue->size == 0);
}

// Function to add an item to the queue.
// It changes rear and size
void enqueue2(struct Queue *queue, int item) {
    if (isFull(queue))
        return;
    queue->rear = (queue->rear + 1) % queue->capacity;
    queue->array[queue->rear] = item;
    queue->size = queue->size + 1;
    // printf("%d enqueued to queue\n", item);
}

// Function to remove an item from queue.
// It changes front and size
int dequeue(struct Queue *queue) {
    if (isEmpty(queue))
        return INT_MIN;
    int item = queue->array[queue->front];
    queue->front = (queue->front + 1) % queue->capacity;
    queue->size = queue->size - 1;
    return item;
}

// Function to get front of queue
int front(struct Queue *queue) {
    if (isEmpty(queue))
        return INT_MIN;
    return queue->array[queue->front];
}

// Function to get rear of queue
int rear(struct Queue *queue) {
    if (isEmpty(queue))
        return INT_MIN;
    return queue->array[queue->rear];
}

int min(int a, int b) {
    if (a < b)
        return a;
    else
        return b;
}
struct solution *graphSolve(struct graph *g, enum problemPart part,
                            int numRooms, int startingRoom, int bossRoom, int numShortcuts,
                            int *shortcutStarts, int *shortcutEnds, int numHeartRooms, int *heartRooms) {
    struct solution *solution = (struct solution *)
        malloc(sizeof(struct solution));
    assert(solution);
    if (part == PART_A) {
        /* IMPLEMENT 2A SOLUTION HERE */
        int visited[numRooms];
        int ans[numRooms];
        for (int i = 0; i < numRooms; i++) {
            visited[i] = 0;
            ans[i] = __INT32_MAX__;
        }
        struct Queue *queue = createQueue(numRooms);
        enqueue2(queue, startingRoom);
        ans[startingRoom] = 0;
        while (!isEmpty(queue)) {
            int currentRoom = dequeue(queue);
            if (visited[currentRoom]) {
                continue;
            }
            visited[currentRoom] = 1;
            for (int i = 0; i < g->numEdges; i++) {
                if (g->edgeList[i]->start == currentRoom && !visited[g->edgeList[i]->end]) {
                    enqueue2(queue, g->edgeList[i]->end);
                    ans[g->edgeList[i]->end] = min(ans[g->edgeList[i]->end], ans[currentRoom] + 1);
                }
                if (g->edgeList[i]->end == currentRoom && !visited[g->edgeList[i]->start]) {
                    enqueue2(queue, g->edgeList[i]->start);
                    ans[g->edgeList[i]->start] = min(ans[g->edgeList[i]->start], ans[currentRoom] + 1);
                }
            }
        }
        solution->heartsLost = ans[bossRoom];
        return solution;
    } else if (part == PART_B) {
        /* IMPLEMENT 2B SOLUTION HERE */
        solution->heartsLost = INT_MAX;

        for (int i = 0; i < numShortcuts; i++) {
            struct graph *shortcutGraph = duplicateGraph(g);
            struct edge *shortcutEdge = (struct edge *)malloc(sizeof(struct edge));
            assert(shortcutEdge);
            shortcutEdge->start = shortcutStarts[i];
            shortcutEdge->end = shortcutEnds[i];
            shortcutEdge->cost = 1;
            addEdge(shortcutGraph, shortcutEdge->start, shortcutEdge->end, shortcutEdge->cost);
            struct solution *shortcutSolution = graphSolve(shortcutGraph, PART_A, numRooms, startingRoom, bossRoom, 0, NULL, NULL, numHeartRooms, heartRooms);
            if (shortcutSolution->heartsLost < solution->heartsLost) {
                solution->heartsLost = shortcutSolution->heartsLost;
            }
            freeGraph(shortcutGraph);
        }
        return solution;
    } else {
        /* IMPLEMENT 2C SOLUTION HERE */
        // Bellman-Ford
        long long int ans[numRooms];
        int isHeart[numRooms];
        for (int i = 0; i < numRooms; i++) {
            ans[i] = 1073741824;
            isHeart[i] = 0;
        }
        for (int i = 0; i < numHeartRooms; i++) {
            isHeart[heartRooms[i]] = 1;
        }

        ans[startingRoom] = 0;
        for (int i = 0; i < numRooms - 1; i++) {
            for (int j = 0; j < g->numEdges; j++) {
                int a = g->edgeList[j]->start;
                int b = g->edgeList[j]->end;
                int c = g->edgeList[j]->cost;
                if (isHeart[b])
                    c--;
                ans[b] = min(ans[b], ans[a] + c);
            }
            for (int j = 0; j < g->numEdges; j++) {
                int a = g->edgeList[j]->end;
                int b = g->edgeList[j]->start;
                int c = g->edgeList[j]->cost;
                if (isHeart[b])
                    c--;
                ans[b] = min(ans[b], ans[a] + c);
            }
            // for (int i = 0; i < numRooms; i++) {
            //     printf("%d ", ans[i]);
            // }
            // printf("\n");
        }

        solution->heartsLost = ans[bossRoom];
    }
    return solution;
}
#include <stdio.h>

int main(int argc, char **argv) {
    /* Read the problem in from stdin. */
    struct graphProblem *problem = readProblem(stdin, PART_C);
    /* Find the solution to the problem. */
    struct solution *solution = findSolution(problem, PART_C);

    /* Report solution */
    printf("%d\n", solution->heartsLost);

    freeProblem(problem);
    freeSolution(solution);

    return 0;
}
