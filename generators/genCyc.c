/**
 * genCyc: a program to generate circulant Ramsey colorings for different kinds of avoided subgraphs.
 * 
 * Author: Steven Van Overberghe
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
#include "checkRamsey.c"

#ifndef GEN_CYC
#define GEN_CYC

//Search for one Ramsey graph vs generate all
#ifndef SEARCH_SINGLE
    #define SEARCH_SINGLE 1
#endif

//Stop after a number of consecutive zero-results
#ifndef NULCOUNTS
    #define NULCOUNTS 4
#endif

#define MAX_CLIQUE 64

//search-splitting parameters
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


int DoPrintCirc = 0;
int DoPrintG6 = 0;
long long splitCount = 0;

int n;
int m;

int colorCount;
int* avoidSize;
Set* distances;
Set* restrictions;
int* restrictionCount;

int count = 0;
int found = 0;

//Neighbourhood of v_index in a certain color
#define ADJ(index, color) (ROTATE(distances[color], index, n))
#define ADJ_UP(index, color) ((distances[color]<<index)&RANGE(n-1))


/**
 * Check if the set contains a k-clique in the given color.
 **/
int checkKClique(Set set, int color, int depth){
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

        if(checkKClique(set & ADJ_UP(next,color),color,depth-1)){
            return 1;
        }

    }
    return 0;
}

int checkCycle(Set searchSet, Set used, int color, int length, int target){
    if(length==2){
        return EMPTY(ADJ(target,color)&searchSet)?0:1;
    }
    int next = 0;
    ITERATE_SET(searchSet,next)
        if(checkCycle(ADJ(next,color)&(~used), WITH(used,next),color,length-1, target)){
            return 1;
        }
    }
    return 0;
}

//check if there is a cycle in the neighbourhood of 0
int checkCycle2(Set searchSet, Set used, int color, int length, int target){
    if(length==2){
        return EMPTY(ADJ(target,color)&searchSet)?0:1;
    }
    int next = 0;
    ITERATE_SET(searchSet,next)
        if(checkCycle2(distances[color]&ADJ(next,color)&(~used), WITH(used,next),color,length-1, target)){
            return 1;
        }
    }
    return 0;
}

//Check if a cycle of specified length is present in the searchset. (including 0)
int checkCycleInSet(Set searchset, int color, int length){
    int v0 = 0;
    ITERATE_SET(searchset,v0)
        if(checkCycle2(searchset&ADJ_UP(v0,color),WITH(S1,v0),color,length,v0)==1){
            return 1;
        }
    }
    return 0;
}


//check if G contains a clique with one edge short
int checkKmE(Set searchSet, int color, int size){
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
        
        //not adjacent to the some vertex
        Set copy = searchSet & (~ADJ_UP(next,color));
        int next2;
        ITERATE_SET(copy,next2)
            if(checkKClique(ADJ_UP(next,color)&ADJ(next2,color)&searchSet,color,size-2)){
                return 1;
            }
        }

        //adjacent to all coming ones
        if(checkKmE(ADJ_UP(next,color)&searchSet,color,size-1)==1){
            return 1;
        }

    }
    return 0;
}


//Checks if 0 and x are in a Kn-e of certain size.
int checkKmE_arc(int x, int color, int size){

    if(checkKmE(distances[color]&ADJ(x,color),color,size-2)){
        return 1;
    }

    
    //neighbours of either x or 0
    Set available = (distances[color] ^ ADJ(x,color))^(BIT(x)|BIT(0));
    int y;
    ITERATE_SET(available,y)
        if(checkKClique(distances[color]&ADJ(x,color)&ADJ(y,color),color,size-3)){
            return 1;
        }
    }

    return 0;
}

//Checks if exactly 2 non-adjacent edges in a Kn are missing. Therefore only finds induced subgraphs.
int checkKm2E(int color, int size){
    Set all = RANGE(n-1);

    int v0 = 0;
    ITERATE_SET(all, v0)
        Set nietBuren = (WITHOUT(~distances[color],0)<<v0)&RANGE(n-1);
        int v1 = 0;
        ITERATE_SET(nietBuren,v1);
            Set buren = ADJ(v1,color)&ADJ(v0,color);
            int v2 = 0;
            ITERATE_SET(buren,v2)
                Set over = (WITHOUT(~distances[color],0)<<v2) & ADJ(v1,color)&ADJ(v0,color);
                int v3 = 0;
                ITERATE_SET(over,v3)
                    Set doorsnede = ADJ(v0,color)&ADJ(v1,color)&ADJ(v2,color)&ADJ(v3,color);
                    if(checkKClique(doorsnede,color,size-4)==1){
                        return 1;
                    }
                }
            }
        }

    }
    return 0;
}



void printCirculants(Set* dists, int colors){
    for (int i = 0; i < colors; i++){
        fprintf(stderr,"Color %d: (",i);
        int first = 1;
        for (int j = 0; j < 8*sizeof(Set); j++){
            if(CONTAINS(dists[i], j)){
                fprintf(stderr,first?"%d":",%d",j);
                first = 0;
            }
        }
        fprintf(stderr,")\n");
    }  
}

Graph* fillCircularGraph(Set N_v0, int n){
    Graph* g = initGraph(n);
    for (int i = 0; i < n; i++)
    {
        g->adj[i] = ROTATE(N_v0,i,n);
    }
    return g;
}

int trueSize(int color){
    return (avoidSize[color]) % MAX_CLIQUE;
}

void checkRamsey(int k, int n){
    for (int i = 0; i < k; i++){
        //check code for errors
        Graph* g1 = fillCircularGraph(distances[i], n);
        if(ostergard(*g1)>=trueSize(i)){
            printf("Error detected in color %d: K%d \n",i,ostergard(*g1));
            printCirculants(distances,colorCount);
            exit(2);
        }
        freeGraph(g1);
    }
}


//try all legal color options for distance x
void try_distance(short x, Set used){
    if(x==m+1){
        if(DoPrintCirc){
            printCirculants(distances,colorCount);
        }
        if(CHECK_GRAPHS){
            checkRamsey(colorCount,n);
        }
        count++;
        found = 1;
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

        if(!CONTAINS(restrictions[i], x)){

            //Avoid color permutations
            if(i>0 && avoidSize[i] == avoidSize[i - 1] && distances[i - 1] == S0){
                continue;
            }

            //Add x and n-x to the used distances
            ADD(distances[i], x);
            ADD(distances[i], n - x);

            Set oldRestrictions = restrictions[i];
            


            if(avoidSize[i] > 4 * MAX_CLIQUE){
                if(checkCycleInSet(ADJ(0,i),i,avoidSize[i] - 4 * MAX_CLIQUE -1)==0){
                    try_distance(nextX, WITH(used, x));
                }
            }else if(avoidSize[i] > 3 * MAX_CLIQUE){
                if(
                        (checkKmE_arc(x, i, avoidSize[i] - 3 * MAX_CLIQUE) == 0) &&
                        (checkKm2E(i, avoidSize[i] - 3 * MAX_CLIQUE) == 0)){
                    try_distance(nextX, WITH(used, x));
                }
            }else if(avoidSize[i] > 2 * MAX_CLIQUE){
                if(checkCycle(WITHOUT(ADJ(x,i),0) , WITH(BIT(x),0), i, avoidSize[i] - 2 * MAX_CLIQUE - 1,0) == 0){
                    try_distance(nextX, WITH(used, x));
                }

            }else if(avoidSize[i] > MAX_CLIQUE){
                if(checkKmE_arc(x, i, avoidSize[i] - MAX_CLIQUE) == 0){
                    try_distance(nextX, WITH(used, x));
                }
            }else{
                if(checkKClique(distances[i]&ADJ(x,i) , i, avoidSize[i] - 2) == 0){
                    try_distance(nextX, WITH(used, x));
                }
            }
            
            REMOVE(distances[i], x);
            REMOVE(distances[i], n - x);
            restrictions[i] = oldRestrictions;
            
        }
    }
    
}

void searchRamsey(){
    distances = calloc(colorCount, sizeof(Set));
    restrictions = calloc(colorCount, sizeof(Set));

    if( n%3==0){
        //add n/3 as impossible in every 3-avoiding color
        for (int i = 0; i < colorCount; i++)
        {
            if(avoidSize[i] == 3){
                ADD(restrictions[i], n / 3);
            }
        }
        
    }
    // 1 should only be considered in non-equivalent colors
    for (int i = 1; i < colorCount; i++)
    {
        if(avoidSize[i - 1] == avoidSize[i]){
            ADD(restrictions[i],1);
        }
    }   
    
    count = 0;
    found = 0;
    splitCount = 0;
    
    try_distance(1, S0);

    free(distances);
    free(restrictions);
}

#include <sys/times.h>
#include <unistd.h>
#define time_factor sysconf(_SC_CLK_TCK)

int main(int argc, char const *argv[])
{
    if (argc < 4) {
        fprintf(stderr, "Provide number of vertices and the Ramsey parameters in increasing order.");
        return 2;
    }

    int colorsNr = 0;
    for (int i = 1; i < argc; i++)
    {
        if (argv[i][0]!='-'){
            colorsNr ++;
        }else if(argv[i][1]=='p'){
            DoPrintCirc = 1;
        }else if(argv[i][1]=='g'){
            DoPrintG6 = 1;
        }
    }
    
    
    colorCount = colorsNr - 1;
    avoidSize = calloc(colorCount, sizeof(int));
    int i = 0;
    for (int j = 2; j < argc; j++){
        switch (argv[j][0]){
        case '-':
        {
            break;
        }
        case 'C':
        {
            avoidSize[i] = ((unsigned int) strtol(&argv[j][1], NULL, 10)) + 2 * MAX_CLIQUE;
            i++;
            break;
        }
        case 'E':
        {
            avoidSize[i] = ((int) strtol(&argv[j][1], NULL, 10)) + MAX_CLIQUE;
            i++;
            break;
        }
        case 'J':
        {
            avoidSize[i] = ((int) strtol(&argv[j][1], NULL, 10)) + MAX_CLIQUE;
            i++;
            break;
        }
        case 'F':
        {
            avoidSize[i] = ((int) strtol(&argv[j][1], NULL, 10)) + 3 * MAX_CLIQUE;
            i++;
            break;
        }
        case 'W':
        {
            avoidSize[i] = ((int) strtol(&argv[j][1], NULL, 10)) + 4 * MAX_CLIQUE;
            i++;
            break;
        }
        default:
        {
            avoidSize[i] = (unsigned int) strtol(argv[j], NULL, 10);
            i++;
            break;
        }
        }
        
    }
    
    
    n = (int) strtol(argv[1],NULL,10) -1;
    m = n/2;

    if(MODULO_SPLIT<1){
        fprintf(stderr,"MODULO SPLIT should be at least 1. \n");
        exit(1);
    }else if (MODULO >= MODULO_SPLIT){
         fprintf(stderr,"Use MODULO < MODULO_SPLIT.\n");
    }
    

    if(MODULO_SPLIT>1){
        fprintf(stderr,"Modulo: %d/%d\n",MODULO,MODULO_SPLIT);
    }

    struct tms TMS;
    unsigned int starttime = 0;

    int nulCount = 0;
    while(nulCount < NULCOUNTS && n<SETSIZE){
        n++;
        m = n/2;
        
        searchRamsey();
        if(count==0){
            nulCount++;
        }else{
            nulCount = 0;
        }

        fprintf(stderr,"%d: %d graphs\n",n, count);
        times(&TMS);
        unsigned int savetime = starttime + (unsigned int) TMS.tms_utime;
        fprintf(stderr, "CPU time: %.1f seconds.\n", (double) savetime / time_factor);
    }

    free(avoidSize);
    return 0;
}



#endif /*GEN_CYC*/
