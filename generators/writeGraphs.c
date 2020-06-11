#include <stdio.h>
#include "checkRamsey.c"

void addEdge(Graph* g, int v0, int v1){
    ADD(g->adj[v0],v1);
    ADD(g->adj[v1],v0);
}

void removeEdge(Graph* g, int v0, int v1){
    REMOVE(g->adj[v0],v1);
    REMOVE(g->adj[v1],v0);
}

void writeG6(Graph g){
    if (g.n<63){
        printf("%c",g.n+63);
    }else{
        printf("%c",126);
        unsigned int n = g.n;
        printf("%c",(char) (n>>12)+63);
        printf("%c",(char) ((n>>6)&RANGE_D(6))+63);
        printf("%c",(char) (n&RANGE_D(6))+63);
    }
    
    unsigned char next = 0;
    int x = 0; int y = 1;
    while(x!=-1){
        for (int j = 5; j >=0; j--)
        {
            if(x==(-1)){
                continue;
            }
            if(CONTAINS(g.adj[x],y)){
                next |= (((unsigned char) 1)<<j);
            }
            
            if(x<y-1){
                x++;
            }else if(y<g.n-1){
                y++;
                x = 0;
            }else{
                x = -1;
            }
            
        }
        printf("%c",(next+63));
        next = (char) 0;
    }
    printf("\n");
}

int checkGraph(Graph g){
    for (int i = 0; i < g.n; i++)
    {
        for (int j = 0; j < g.n; j++)
        {
            if(CONTAINS(g.adj[i],j)!=S0 && CONTAINS(g.adj[j],i)==S0){
                printf("Not symmetric %d %d ", i, j);
                return 0;
            }

            if(CONTAINS(g.adj[j],i)!=S0 && CONTAINS(g.adj[i],j)==S0){
                printf("Not symmetric %d %d ", i, j);
                return 0;
            }
        }
        
    }
    return 1;
}