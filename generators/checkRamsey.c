#include <stdlib.h>
#include <stdio.h>

#include "sets.c"

#ifndef RAMSEY
#define RAMSEY

typedef struct Graph{
    int n; //number of vertices
    Set* adj;  //adjacency sets
} Graph;


Graph* initGraph(int n){
    Graph* g = malloc(sizeof(Graph));
    g->n = n;
    g->adj = calloc(n, sizeof(Set ));
    return g;
}

void freeGraph(Graph* g){
    free(g->adj);
    free(g);
}

int adjacent(Graph* g, int i, int j){
    return CONTAINS(g->adj[i],j);
}

#define ADJACENT(g,i,j) (CONTAINS(g->adj[i],j))


Graph* dualGraph(Graph* graph){
    Graph* g = malloc(sizeof(Graph));
    int n = graph->n;
    g->n = n;
    g->adj = malloc(n* sizeof(Set  ));

    for (int i = 0; i < n; i++)
    {
        g->adj[i] = WITHOUT(COMPLEMENT(graph->adj[i],n),i);
    }
    
    return g;
}

Graph complGraph(Graph graph){
    int n = graph.n;
    Set* adjacency = calloc(n,sizeof(Set));

    Graph complement = (Graph) {.n = n, .adj = adjacency};

    for (int i = 0; i < n; i++)
    {
        complement.adj[i] = WITHOUT(COMPLEMENT(graph.adj[i],n),i);
    }
    return complement;
}


void complInplace(Graph graph){
    int n = graph.n;

    for (int i = 0; i < n; i++)
    {
        graph.adj[i] = WITHOUT(COMPLEMENT(graph.adj[i],n),i);
    }
}

int numberOfEdges(Graph graph){
    int count = 0;
    for (int v1 = 0; v1 < graph.n; v1++)
    {
        for (int v2 = v1+1; v2 < graph.n; v2++)
        {
            if(CONTAINS(graph.adj[v1],v2)){
                count++;
            }
        }
    }
    return count;
}


int max;
Set besteKliek;

int findClique(Graph g, Set available, Set  kliek, int grootte, int* found, int* best){
    if(EMPTY(available)){
        if(grootte>max){
            max = grootte;
            besteKliek = kliek;
            *found = 1;
            return grootte;
        }
        return 0;
    }

    while(available != S0){

        int v = FIRSTSET(available);

        if(grootte + best[v] <= max){
            return 0;
        }

        ADD(kliek,v);
        REMOVE(available, v);
        Set Nv = g.adj[v];
        findClique(g,available & Nv, kliek, grootte+1, found, best);
        REMOVE(kliek, v);

        if(*found){
            return grootte;
        }

    }
}

//Find size of largest clique in g
int ostergard(Graph g){
    int n = g.n;
    int* found = malloc(sizeof(int));
    int* best = calloc(n,sizeof(int));
    max = 0;
    besteKliek = S0;

    for(int i = n-1; i>=0; i--){
        *found = 0;
        int k = findClique(g, FROMRANGE(i) & g.adj[i], BIT(i), S1 , found, best);
        best[i] = max;
    }

    free(found);
    int res = best[0]; //=max
    free(best);
    return res;
}

int isRamsey(Graph* g, int k1, int k2){
    if (ostergard(*g)>=k1){
        return 0;
    }
    if(ostergard(complGraph(*g))>=k2){
        return 0;
    }
    return 1;    
}


#endif /*RAMSEY*/
