/**
 * genBlock: a program to generate block-circulant Ramsey colorings for different kinds of avoided subgraphs.
 * 
 * Author: Steven Van Overberghe, Dec 2019
 * 
 * ---------------------------------
 * 
 * Copyright (c) 2020
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include "sets.c"

#ifndef GEN_BLOCK
#define GEN_BLOCK
#define GEN_CYC
#include "writeGraphs.c"


#ifndef SEARCH_SINGLE
    #define SEARCH_SINGLE 0
#endif

#ifndef NULCOUNTS
    #define NULCOUNTS 3
#endif

#define MAX_CLIQUE 64

//search splitting parameters
#ifndef MODULO_SPLIT
    #define MODULO_SPLIT 1
    #define SPLIT_DEPTH 0
    #define MODULO 0
#else
#ifndef SPLIT_DEPTH
        #define SPLIT_DEPTH 6
    #endif
    #ifndef MODULO
        #define MODULO 1
    #endif
#endif

#ifndef CHECK_GRAPHS
    #define CHECK_GRAPHS 1
#endif

long long splitCount;

int n;
int m;
int k; //number of blocks

short* blockSize;
short* blockStart;
short* halfBlock;

short* rowX;

int colorCount = 1;
int* avoidSize;
Set*** distances;
Graph* graphs;

Set* reverseDistances;


Set half;

int count = 0;
int found = 0;

#define SHIFT_START(color, block) (distances[color][block]>>blockStart[block])
#define ADJ(index, color) (ROTATE(distances[color], index, n))
#define ADJ_BLOCK(index, color, block) ((ROTATE(SHIFT_START(color,block),index,blockSize[block]))<<blockStart[block])
#define ADJ_UP(index, color) ((distances[color]<<index)&RANGE(n-1))



int calcRow(int x){
    for (int i = 0; i < k-1; i++)
    {
        if(x<blockStart[i+1]){
            return i;
        }
    }
    return k-1;
    
}

int findRow(int x){
    return rowX[x];
    // return x/(n/k);
}

//Assuming left to right and top to bottom fill
Set adjacency(int x, int color, int column){
    Set result = S0;
    int row = findRow(x);
    for (int i = 0; i < column+1; i++){
        result |= (ROTATE_D(distances[color][row][i],x-blockStart[row],blockSize[i])<<blockStart[i]);
    }
    return result;
}



int checkKClique(Set set, int color, int depth, int column){
    if(depth==0){
        return 1;
    }
    if(depth==1){
        return set!=S0;
    }
    if(POPCOUNT(set)<depth){
        return 0;
    }
    int next;
    ITERATE_SET(set,next)
        Set new = set & adjacency(next, color, column);
        if(checkKClique(new,color,depth-1, column)){
            return 1;
        }

    }
    return 0;
}

int checkCycle(Set searchSet, Set used, int color, int length, int column, int v0){
    if(length==2){
        return !EMPTY(adjacency(v0,color,column)&searchSet);
    }
    int next = 0;
    ITERATE_SET(searchSet,next)
        if(checkCycle(adjacency(next,color, column)&(~used), WITH(used,next),color,length-1,column, v0)){
            return 1;
        }
    }
    return 0;
}

int checkCycle2(Set searchSet, Set used, int color, int length, int target, int center, int column){
    if(length==2){
        return EMPTY(adjacency(target,color,column)&searchSet)?0:1;
    }
    int next = 0;
    ITERATE_SET(searchSet,next)
        if(checkCycle2(adjacency(center,color,column)&adjacency(next,color,column)&(~used), WITH(used,next),color,length-1, target,center,column)){
            return 1;
        }
    }
    return 0;
}

//Not so efficient
int checkCycleInSet(int color, int length, int column){
    for (int i = 0; i < k; i++){
        int center = blockStart[i];
        int v0 = 0;
        Set searchset = adjacency(center, color,column);
        ITERATE_SET(searchset,v0)
            if(checkCycle2(adjacency(center,color,column)&adjacency(v0,color,column)&FROMRANGE(v0),WITH(BIT(center),v0),color,length,v0,center,column)==1){
                return 1;
            }
        }
    }
    return 0;
}

int checkKmE(Set searchSet, int color, int size, int column){
    if(size==0){
        return 1;
    }
    if(size==1){
        return searchSet!=S0;
    }
    if(POPCOUNT(searchSet)<size){
        return 0;
    }
    if(size==2){
        return 1;
    }

    int next;
    ITERATE_SET(searchSet,next)
        //adjacent to all coming ones
        if(checkKmE(adjacency(next,color,column)&searchSet,color,size-1, column)==1){
            return 1;
        }

        //not adjacent to the some vertex
        Set copy = searchSet & (~adjacency(next,color,column));
        int next2;
        ITERATE_SET(copy,next2)
            if(checkKClique(adjacency(next,color,column)&adjacency(next2,color,column)&searchSet,color,size-2, column)){
                return 1;
            }
        }
    }
    return 0;
}

int checkKmE_arc(int x, int color, int size, int row, int column){
    int v0 = blockStart[row];
    Set Nvo = adjacency(v0,color,column);
    Set Nx = adjacency(blockStart[column]+x,color,column);
    

    if(checkKmE(Nvo & Nx,color,size-2,column)){
        return 1;
    }

    //neighbours of either x or 0
    Set available = (Nvo ^ Nx)^(BIT(blockStart[column]+x)|BIT(v0));
    
    int y;
    ITERATE_SET(available,y)
        if(checkKClique(Nvo & Nx & adjacency(y,color,column),color,size-3,column)==1){
            return 1;
        }
    }

    return 0;
}

int checkKm2E(int color, int size, int column){
    Set all = RANGE(n-1);

    int v0 = 0;
    ITERATE_SET(all, v0)
        Set nietBuren = (~adjacency(v0,color,column))&(~RANGE(v0))&RANGE(n-1);
        int v1 = 0;
        ITERATE_SET(nietBuren,v1);
            Set buren = adjacency(v1,color,column)&adjacency(v0,color,column);
            int v2 = 0;
            ITERATE_SET(buren,v2)
                Set over = (~adjacency(v2,color,column))&(~RANGE(v2)) & adjacency(v1,color,column)&adjacency(v0,color,column)&RANGE(n-1);
                int v3 = 0;
                ITERATE_SET(over,v3)
                    Set doorsnede = adjacency(v0,color,column)&adjacency(v1,color,column)
                                    &adjacency(v2,color,column)&adjacency(v3,color,column);
                    if(checkKClique(doorsnede,color,size-4,column)==1){
                        return 1;
                    }
                }
            }
        }

    }
    return 0;
}

int trueSize(int color){
    return (avoidSize[color]) % MAX_CLIQUE;
}



void checkRamsey(){
    for (int i = 0; i < colorCount; i++){
        //check code for errors
        if(ostergard(graphs[i])>=trueSize(i)){
            printf("Error in color %d\n",i);
            exit(2);
        }
    }
}


int representative(Set set, int size){
    for (int a = 1; a < size; a++)
    {
        if(ROTATE(set,a,size)>set){
            return 0;
        }
    }
    return 1;
}

int canRepresentative(Set set, int size){
    Set max = RANGE_D(size);
    for (int i = 1; i < size; i++)
    {
        if(set < ((set << i) & max)){
            return 0;
        }
    }
    return 1;
}


int firstFound = 1;

//try all legal color options for distance x
void try_distance(int x, int column, int row){
     if(x==0 && column>1 && row == 0){
        if(distances[0][column-1][column-1]<distances[0][column-2][column-2]){
            return;
        }
    }

    //found Ramsey graph
    if(column == k){
        count++;
        found = 1;

        for (int i = 0; i < colorCount; i++){
            graphs[i].n=n;
            for (int a = 0; a < n; a++)
            {
                graphs[i].adj[a] = adjacency(a,i,k-1);
            }
            
        }

        for (int i = 0; i < colorCount-1; i++){
            writeG6(graphs[i]);
        }
        fflush(stdout);

        if(CHECK_GRAPHS){
            checkRamsey();
        }
        
        return;
    }

    if(column>0 && row==0){
        if(canRepresentative(reverseDistances[column],x)==0){
            return;
        }
    }

    if(x==0 && row==1){
        if(representative(reverseDistances[column],blockSize[column]) == 0){
            return;
        }
    }

    if(column!=row && x==blockSize[column]){
        try_distance(0,column, row+1);
        return;
    }

    if(row==column && x==0){
        try_distance(x+1,column, row);
        return;
    }

    if(column==row && x==halfBlock[column]+1){
            try_distance(0,column+1,0);
            return;
    }

   if(MODULO_SPLIT!=1 && x==SPLIT_DEPTH){
       splitCount++;
       if(splitCount%(MODULO_SPLIT)!=MODULO){
           return;
       }
   }

    int nextX = x+1;

    for (int i = 0; i < colorCount; i++){

        if(SEARCH_SINGLE && found){
            return;
        }


            //Add x and n-x to the current block-string
            ADD(distances[i][row][column], x);
            ADD(distances[i][column][row], (blockSize[column] - x)%blockSize[column]);

            if(i==0 && row==0 && column>0){
                ADD(reverseDistances[column],blockSize[column]-x-1);
            }

            if(avoidSize[i] > 4 * MAX_CLIQUE){
                if(checkCycleInSet(i, avoidSize[i] - 4 * MAX_CLIQUE - 1,column) == 0){
                    try_distance(nextX, column, row);
                }

            }else if(avoidSize[i] > 3 * MAX_CLIQUE){
                if(
                        (checkKmE_arc(x, i, avoidSize[i] - 3 * MAX_CLIQUE,row,column) == 0) &&
                        (checkKm2E(i, avoidSize[i] - 3 * MAX_CLIQUE,column) == 0)){
                    try_distance(nextX, column,row);
                }
            }else if(avoidSize[i] > 2 * MAX_CLIQUE){
                int v0 = blockStart[row];
                if(checkCycle(WITHOUT(adjacency(x+blockStart[column],i,column),v0) , WITH(BIT(x+blockStart[column]),v0), i, avoidSize[i] - 2 * MAX_CLIQUE - 1,column, v0) == 0){
                    try_distance(nextX, column, row);
                }

            }else if(avoidSize[i] > MAX_CLIQUE){
                if(checkKmE_arc(x, i, avoidSize[i] - MAX_CLIQUE, row, column) == 0){
                        try_distance(nextX, column, row);
                    }
                
            }else{
                int v0 = blockStart[row];
  
                if(checkKClique(adjacency(v0,i,column)&adjacency(blockStart[column]+x,i,column) , i, avoidSize[i] - 2,column) == 0){

                    try_distance(nextX,column, row);
                }
            }

            REMOVE(distances[i][row][column], x);
            REMOVE(distances[i][column][row], (blockSize[column] - x)%blockSize[column]);

            if(i == 0 && row==0 && column>0){
                REMOVE(reverseDistances[column],blockSize[column]-x-1);
            }
        }
 

}

void searchRamsey(){
    if(n%k!=0){
        count = 0;
        return;
    }
    
    fprintf(stderr, "k=%d\n",k);
    fprintf(stderr, "n=%d\n",n);
    // for (int i = 0; i < colorCount; i++)
    // {
    //     fprintf(stderr, "a%d=%d\n",i,avoidSize[i]);
    // }

    distances = calloc(colorCount,sizeof(Set*));
    for (int i = 0; i < colorCount; i++)
    {
        distances[i] = calloc(k,sizeof(Set*));
        for (int j = 0; j < k; j++)
        {
            distances[i][j] = calloc(k,sizeof(Set));
        }
        
    }
    
    blockSize = calloc(k,sizeof(short));
    blockStart = calloc(k,sizeof(short));
    halfBlock = calloc(k,sizeof(short));
    rowX = calloc(n,sizeof(short));

    reverseDistances = calloc(k,sizeof(Set));

    for (int i = 0; i < k; i++){
        blockSize[i] = n/k;
        halfBlock[i] = n/k/2;
    }
    blockStart[0] = 0;
    for (int i = 1; i < k; i++)
    {
        blockStart[i] = i*n/k;
    }

    for (int i = 0; i < n; i++)
    {
        rowX[i] = calcRow(i);
    }

    graphs = calloc(colorCount, sizeof(Graph));
    for (int i = 0; i < colorCount; i++)
    {
        graphs[i].n = 0;
        graphs[i].adj = calloc(n+1,sizeof(Set));
    }
    

    count = 0;
    found = 0;
    splitCount = 0;
    
    try_distance(1, 0, 0);

    free(distances);
}


#include <sys/times.h>
#include <unistd.h>
#define time_factor sysconf(_SC_CLK_TCK)

int main(int argc, char const *argv[])
{
    if (argc < 5) {
        fprintf(stderr,"Provide number of vertices, k, and the Ramsey parameters in increasing order.\n");
        return 2;
    }
    colorCount = argc - 3;
    avoidSize = calloc(colorCount, sizeof(int));
    for (int i = 0; i < colorCount; i++){
        switch (argv[3+i][0]){
            case 'C':
            {
                avoidSize[i] = ((int) strtol(&argv[3 + i][1], NULL, 10)) + 2 * MAX_CLIQUE;
                break;
            }
            case 'E':
            {
                avoidSize[i] = ((int) strtol(&argv[3 + i][1], NULL, 10)) + MAX_CLIQUE;
                break;
            }
            case 'J':
            {
                avoidSize[i] = ((int) strtol(&argv[3 + i][1], NULL, 10)) + MAX_CLIQUE;
                break;
            }
            case 'F':
            {
                avoidSize[i] = ((int) strtol(&argv[3 + i][1], NULL, 10)) + 3 * MAX_CLIQUE;
                break;
            }
            case 'W':
            {
                avoidSize[i] = ((int) strtol(&argv[3 + i][1], NULL, 10)) + 4 * MAX_CLIQUE;
                break;
            }
            default:
            {
                avoidSize[i] = (int) strtol(argv[3 + i], NULL, 10);
                break;
            }
        }

    }

    n = (int) strtol(argv[1],NULL,10) -1;
    k = (int) strtol(argv[2],NULL,10);

    if(MODULO_SPLIT>1){
        fprintf(stderr,"Modulo: %d/%d\n",MODULO,MODULO_SPLIT);
    }

    struct tms TMS;
    unsigned int starttime = 0;

    int nulCount = 0;
    while(nulCount < NULCOUNTS && n<SETSIZE){
        n++;

        if(n%k!=0){
            continue;
        }
        
        searchRamsey();
        
        if(count==0){
            nulCount ++;
        }else{
            nulCount = 0;
        }


        if(count==1){
            fprintf(stderr,"count %d: %d graph\n",n, count);
        }else{
            fprintf(stderr,"count %d: %d graphs\n",n, count);
        }
        times(&TMS);
        unsigned int savetime = starttime + (unsigned int) TMS.tms_utime;
        fprintf(stderr, "CPU time: %.1f seconds.\n", (double) savetime / time_factor);
    }
    return 0;
}



#endif /*GEN_BLOCK*/
