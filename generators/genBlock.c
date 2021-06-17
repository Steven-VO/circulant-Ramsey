/**
 * genBlock: a program to generate block-circulant Ramsey colorings for different kinds of avoided subgraphs.
 * 
 * Author: Steven Van Overberghe, Dec 2019
 * 
 * ---------------------------------
 * 
 * Copyright (c) 2021
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
    #define MODULO 0
#endif

#ifndef CHECK_GRAPHS
    #define CHECK_GRAPHS 1
#endif

#ifndef SYM_2C
    #define SYM_2C 0
#endif

long long splitCount;
long long callCount = 0LL;

int n; //number of vertices
int m; //size of each block
int k; //number of blocks

short blockSize;
short* blockStart;
short halfBlock;

short* rowX;

int colorCount = 1;
int* avoidSize;
Set** adjs;
Set*** distances;
Graph* graphs;


Set half;

int count = 0;
int found = 0;

#define ADJ(index, color) (adjs[color][index])

int calcRow(int x){
    for (int i = 0; i < k-1; i++){
        if(x<blockStart[i+1]){
            return i;
        }
    }
    return k-1; 
}

int findRow(int x){
    return rowX[x];
}



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
        Set new = set & ADJ(next, color);
        if(checkKClique(new,color,depth-1)){
            return 1;
        }

    }
    return 0;
}

int checkCycle(Set searchSet, Set used, int color, int length, int v0){
    if(length==2){
        return !EMPTY(ADJ(v0,color)&searchSet);
    }
    int next = 0;
    ITERATE_SET(searchSet,next)
        if(checkCycle(ADJ(next,color)&(~used), WITH(used,next),color,length-1, v0)){
            return 1;
        }
    }
    return 0;
}

int checkCycle2(Set searchSet, Set used, int color, int length, int target, int center){
    if(length==2){
        return EMPTY(ADJ(target,color)&searchSet)?0:1;
    }
    int next = 0;
    ITERATE_SET(searchSet,next)
        if(checkCycle2(ADJ(center,color)&ADJ(next,color)&(~used), WITH(used,next),color,length-1, target,center)){
            return 1;
        }
    }
    return 0;
}

//Not so efficient
int checkCycleInSet(int color, int length){
    for (int i = 0; i < k; i++){
        int center = blockStart[i];
        int v0 = 0;
        Set searchset = ADJ(center, color);
        ITERATE_SET(searchset,v0)
            if(checkCycle2(ADJ(center,color)&ADJ(v0,color)&FROMRANGE(v0),WITH(BIT(center),v0),color,length,v0,center)==1){
                return 1;
            }
        }
    }
    return 0;
}

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
        //adjacent to all coming ones
        if(checkKmE(ADJ(next,color)&searchSet,color,size-1)==1){
            return 1;
        }

        //not adjacent to the some vertex
        Set copy = searchSet & (~ADJ(next,color));
        int next2;
        ITERATE_SET(copy,next2)
            if(checkKClique(ADJ(next,color)&ADJ(next2,color)&searchSet,color,size-2)){
                return 1;
            }
        }
    }
    return 0;
}

int checkKmE_arc(int x, int color, int size, int row, int column){
    int v0 = blockStart[row];
    Set Nvo = ADJ(v0,color);
    Set Nx = ADJ(blockStart[column]+x,color);
    

    if(checkKmE(Nvo & Nx,color,size-2)){
        return 1;
    }

    //neighbours of either x or 0
    Set available = (Nvo ^ Nx)^(BIT(blockStart[column]+x)|BIT(v0));
    
    int y;
    ITERATE_SET(available,y)
        if(checkKClique(Nvo & Nx & ADJ(y,color),color,size-3)==1){
            return 1;
        }
    }

    return 0;
}

int checkKm2E(int color, int size){
    Set all = RANGE(n-1);

    int v0 = 0;
    ITERATE_SET(all, v0)
        Set nietBuren = (~ADJ(v0,color))&(~RANGE(v0))&RANGE(n-1);
        int v1 = 0;
        ITERATE_SET(nietBuren,v1);
            Set buren = ADJ(v1,color)&ADJ(v0,color);
            int v2 = 0;
            ITERATE_SET(buren,v2)
                Set over = (~ADJ(v2,color))&(~RANGE(v2)) & ADJ(v1,color)&ADJ(v0,color)&RANGE(n-1);
                int v3 = 0;
                ITERATE_SET(over,v3)
                    Set doorsnede = ADJ(v0,color)&ADJ(v1,color)
                                    &ADJ(v2,color)&ADJ(v3,color);
                    if(checkKClique(doorsnede,color,size-4)==1){
                        return 1;
                    }
                }
            }
        }

    }
    return 0;
}

int biPar[5][2] = {{0,0},{0,0}};

int searchBipartiteClass(Set possibleAdd, Set remaining, int size2, int size2Left, int desiredSize1, int desiredSize2, int color){
    int size1;
    if(size2>=desiredSize2){
        size1 = POPCOUNT(remaining);
        if(size1<desiredSize1){
            return 0;
        }else{
            return 1;
        }
    }

    if(desiredSize1!=1 && remaining==S1){
        return 0;
    }
    size1 = POPCOUNT(remaining);
    if(size1<desiredSize1){
        return 0;
    }

    Set nextRemaining = remaining;
    int y;
    ITERATE_SET(possibleAdd, y)
        size2Left--;
        if(size2+1+size2Left<desiredSize2){
            return 0;
        }
        nextRemaining = remaining & ADJ(y,color);
        if (searchBipartiteClass(possibleAdd,nextRemaining,size2+1,size2Left,desiredSize1,desiredSize2,color)==1){
            return 1;
        }
        
    }

    return 0;
}


int checkBipartite(int x, int color, int size1, int size2, int column, int row){
    Set neighX = ADJ(blockStart[column] + x,color);
    int setsize2 = POPCOUNT(ADJ(blockStart[row],color));

    int result1 = searchBipartiteClass(WITHOUT(ADJ(blockStart[row],color),blockStart[column]+x),neighX,1,setsize2-1,size1,size2,color);

    if(size1!=size2){
        int result2 = searchBipartiteClass(WITHOUT(ADJ(blockStart[row],color),blockStart[column] +x),neighX,1,setsize2-1,size2,size1,color);
        return (result1 || result2);
    }
    return result1;
}

//-------------------------------------------

int unitsNr;
short* units = NULL;

Set* fixUnits = NULL;

void calcUnits(int n){
    unitsNr = 0;
    for (int i = 2; i <= n/2; i++){
        int coPrime = 1;
        for (int j = 2; j <= i; j++){
            if((i%j==0) && (n%j==0)){
                coPrime =0;
                break;
            }
        }
        if(coPrime){
            unitsNr++;
        }
    }
    
    units = calloc(unitsNr,sizeof(short));
    int t = 0;
    for (int i = 2; i <= n/2; i++){
        int coPrime = 1;
        for (int j = 2; j <= i; j++)
        {
            if(i%j==0 && n%j==0){
                coPrime =0;
                break;
            }
        }
        if(coPrime){
            units[t] = i;
            t++;
        }
    }
}

Set multiplySet(Set s, int times, int max){
    Set result = S0;
    int el;
    ITERATE_SET(s, el)
        ADD(result, (el*times)%max);
    }
    return result;
}

int isMultRepresentative0(Set s, Set base){
    for (int index = 0; index < unitsNr; index++){
        Set mult = multiplySet(s,units[index],blockSize);
        if(mult < base){
            return 0;
        }else if(mult == base){
            //automorphisms to keep track of
            ADD(fixUnits[0],index);
        }
    }
    return 1;
}

//Check if s can be made smaller than base by multiplication
//Or smaller than the previous with inherited auts
int isMultRepresentative(Set s, int keepTrackIndex){
    for (int index = 0; index < unitsNr; index++)
    {
        int toMult = units[index];
        Set mult = multiplySet(s,toMult,blockSize);
        if(mult < distances[0][0][0]){
            return 0;
        }
        int blockIndex = 0;
        while(blockIndex < keepTrackIndex && CONTAINS(fixUnits[blockIndex],index)){
            if(mult < distances[0][blockIndex+1][blockIndex+1]){
                //multiplication would give smaller sequence
                return 0;
            }else if(blockIndex == keepTrackIndex - 1 && mult == distances[0][blockIndex+1][blockIndex+1]){
                //automorphisms to keep track of
                ADD(fixUnits[keepTrackIndex],index);
            }
            blockIndex++;
        }
    }
    return 1;
}

Set multRepresentativeSet(Set s){
    Set min = s;
    for (int index = 0; index < unitsNr; index++){
        Set mult = multiplySet(s,units[index],blockSize);
        if(mult < min){
            min = mult;
        }
    }
    return min;
}


//-----------------------------------------------------
//Cyclical canonicity

//Indices represent translation size of cyclical automorphisms of word:
// blocks[0][i1-2][i2-2]
Set** cyclAuts;

int representative(Set set, int size){
    for (int a = 1; a < size; a++){
        if(ROTATE(set,a,size)>set){
            return 0;
        }
    }
    return 1;
}

int representativeS(Set set, int size, Set indices){
    int a;
    ITERATE_SET(indices, a)
        if(ROTATE(set,a,size)>set){
            return 0;
        }
    }
    return 1;
}

int representativeSave(Set set, int size, int i1, int i2){
    for (int a = 1; a < size; a++){
        Set rotation = ROTATE(set,a,size);
        if(rotation > set){
            return 0;
        }else if(rotation == set){
            ADD(cyclAuts[i1][i2],a);
        }
    }
    return 1;
}

int representativeSaveS(Set set, int size, int i1, int i2, Set indices){
    int a;
    ITERATE_SET(indices, a)
        Set rotation = ROTATE(set,a,size);
        if(rotation > set){
            return 0;
        }else if(rotation == set){
            ADD(cyclAuts[i1][i2],a);
        }
    }
    return 1;
}

int canRepresentative(Set set, int size){
    Set max = RANGE_D(blockSize);
    for (int i = 1; i < size; i++)
    // for (int i = size-1; i > 0; i--)
    {
        if(set < ((set << i) & max )){
            return 0;
        }
    }
    return 1;
}

Set representativeSet(Set set, int size){
    Set max = set;
    for (int a = 1; a < size; a++){
        Set rot = ROTATE(set,a,size);
        if(rot > max){
            max = rot;
        }
    }
    return max;
}

//--------------------------------------------------

//assumes column and column+1 are completely filled
int aboveDiagCanonical(int column){
    int row = 0;
    while(row < column){
        if(distances[0][row][column] > distances[0][row][column+1]){
            return 0;
        }else if(distances[0][row][column] == distances[0][row][column+1]){
            row++;
        }else{
            return 1;
        }
    }
    return 1;
}

//------------------------------------


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



void addNewEdge(int b0, int b1, int x, int color){
    int x0 = b0*blockSize;
    int x1 = b1*blockSize;

    //Split the computation to save the modulo
    int next = x1+x;
    for (int y = 0; y < blockSize-x; y++)
    {
        ADD(adjs[color][x0+y],next);
        ADD(adjs[color][next],x0+y);
        next++;
    }
    next = x1;
    for (int y = blockSize-x; y < blockSize; y++)
    {
        int next = x1+x+y-blockSize;
        ADD(adjs[color][x0+y],next);
        ADD(adjs[color][next],x0+y);
        next++;
    }
}

void removeNewEdge(int b0, int b1, int x, int color){
    int x0 = b0*blockSize;
    int x1 = b1*blockSize;

    int next = x1+x;
    for (int y = 0; y < blockSize-x; y++)
    {
        REMOVE(adjs[color][x0+y],next);
        REMOVE(adjs[color][next],x0+y);
        next++;
    }
    next = x1;
    for (int y = blockSize-x; y < blockSize; y++)
    {
        int next = x1+x+y-blockSize;
        REMOVE(adjs[color][x0+y],next);
        REMOVE(adjs[color][next],x0+y);
        next++;
    }
}


long long pruneCount = 0LL;

//try all legal color options for distance x
void try_distance(int x, int column, int row, char flipRemains){
    callCount++;

    //found Ramsey graph
    if(column == k){
        count++;
        found = 1;

        for (int i = 0; i < colorCount; i++){
            graphs[i].n=n;
            for (int a = 0; a < n; a++){
                graphs[i].adj[a] = ADJ(a,i);
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

    if(column>0 && row==0 && x < blockSize){
        if(canRepresentative(distances[0][row][column],x) == 0){
            return;
        }
    }

    //Checks if the current column is cyclicaly canonical
    if(x == blockSize){
        if(column == 1){
            //cyclic representative
            if(representative(distances[0][row][column],blockSize) == 0){
                return;
            }
            //flipping all blocks
            Set representative = representativeSet(multiplySet(distances[0][row][column],blockSize-1,blockSize),blockSize);
            if(distances[0][row][column] > representative){
                return;
            }else if(distances[0][row][column] == representative){
                flipRemains = 1;
            }else{
                flipRemains = 0;
            }
        }else if(row == 0){
            cyclAuts[column - 2][row] = S0;
            if(representativeSave(distances[0][row][column],blockSize,column - 2,row) == 0){
                return;
            }
        }
        else if(row < column - 1){
            cyclAuts[column - 2][row] = S0;
            if( !EMPTY(cyclAuts[column - 2][row-1]) && 
                representativeSaveS(distances[0][row][column],blockSize,
                    column - 2,row,cyclAuts[column - 2][row-1]) == 0){
                return;
            }
        }else{
            if( !EMPTY(cyclAuts[column - 2][row-1]) && 
                representativeS(distances[0][row][column],blockSize,cyclAuts[column - 2][row-1]) == 0){
                return;
            }
        }
    }

    if(flipRemains && x==blockSize && column>1 && row==0){
        Set representative = representativeSet(multiplySet(distances[0][row][column],blockSize-1,blockSize),blockSize);
        if(distances[0][row][column] > representative){
            return;
        }else if(distances[0][row][column] == representative){
            flipRemains = 1;
        }else{
            flipRemains = 0;
        }
    }

    //For 2-coloured numbers with symmetric parameters: first colour should be smaller than second
    //We accept the graph with the smallest minimal block.  If compl(G) already has a smaller block only taking the first block in consideration,
    //then the minimal representation of compl(G) will also be smaller.  Since compl(compl(G))=G, at least one of both will be accepted.
    if(SYM_2C && row == column && x==halfBlock+1){
        if(distances[0][0][0] > multRepresentativeSet(distances[1][row][row])){
            return;
        }
    }


    //Blocks on the diagonal must be non-decreasing
    if(row==column && x==halfBlock+1 && column>0){
        if(distances[0][row][row] < distances[0][row-1][row-1]){
            return;
        }else if(column>1 && distances[0][row][row] == distances[0][row-1][row-1]){
            //If two consecutive diagonal blocks are equal, we should look at the blocks above to break ties.
            if(aboveDiagCanonical(column-1) == 0){
                return;
            }
        }
    }
    //Pre-check on non-decreasingness
    if(row==column && x<halfBlock && x>1 && column>0){
        Set currentMask = RANGE_D(m) ^ RANGE_D(m-x+1);
        if((distances[0][row][row] &  currentMask) < (distances[0][row-1][row-1] & currentMask)){
            pruneCount++;
            return;
        }
    }

    //check multiplicative canonicity
    if(row == column && x==halfBlock+1 && column==0){
        fixUnits[0] = S0;
        if(isMultRepresentative0(distances[0][column][column],distances[0][column][column]) == 0){
            return;
        }
    }
    
    if(row == column && x==halfBlock+1 && column>0){

        if(distances[0][column][column] == distances[0][column-1][column-1]){
            //Don't do the same work twice
            fixUnits[column] = fixUnits[column-1];
        }else{
            fixUnits[column] = S0;
            if(isMultRepresentative(distances[0][column][column],column)==0){
                return;
            }
            //TODO: remaining multiplicative automorphisms at the end, should be checked further.
            //This will not change the runtime, however.
        }
        
    }


    if(MODULO_SPLIT != 1 && column==1 && row == 0 && x==blockSize){
       splitCount++;
       if(splitCount%(MODULO_SPLIT)!=MODULO){
           return;
       }
    }

    if(column!=row && x==blockSize){
        try_distance(0,column, row+1, flipRemains);
        return;
    }

    if(row==column && x==0){
        try_distance(x+1,column, row,flipRemains);
        return;
    }

    if(column==row && x==halfBlock+1){
        try_distance(0,column+1,0,flipRemains);
        return;
    }

    int nextX = x+1;

    for (int i = 0; i < colorCount; i++){

        if(SEARCH_SINGLE && found){
            return;
        }

            int newDist;
            if(row==column){
                //Add x and n-x to the current block-string
                ADD(distances[i][row][column], x);
                ADD(distances[i][row][column], blockSize - x);
                addNewEdge(row, column, x, i);
                newDist = x;
            }else{
                ADD(distances[i][row][column], blockSize-x-1);
                addNewEdge(row, column, blockSize-x-1, i);
                newDist = blockSize - x - 1;
            }
            

            if(avoidSize[i] > 5 * MAX_CLIQUE){
                if(checkBipartite(newDist,i,biPar[i][0],biPar[i][1],column,row)==0){
                    try_distance(nextX, column, row, flipRemains);
                }
            } else 
            if(avoidSize[i] > 4 * MAX_CLIQUE){
                if(checkCycleInSet(i, avoidSize[i] - 4 * MAX_CLIQUE - 1) == 0){
                    try_distance(nextX, column, row, flipRemains);
                }

            }else if(avoidSize[i] > 3 * MAX_CLIQUE){
                if(
                        (checkKmE_arc(newDist, i, avoidSize[i] - 3 * MAX_CLIQUE,row,column) == 0) &&
                        (checkKm2E(i, avoidSize[i] - 3 * MAX_CLIQUE) == 0)){
                    try_distance(nextX, column,row, flipRemains);
                }
            }else if(avoidSize[i] > 2 * MAX_CLIQUE){
                int v0 = blockStart[row];
                if(checkCycle(WITHOUT(ADJ(newDist+blockStart[column],i),v0) , WITH(BIT(newDist+blockStart[column]),v0), i, avoidSize[i] - 2 * MAX_CLIQUE - 1, v0) == 0){
                    try_distance(nextX, column, row, flipRemains);
                }

            }else if(avoidSize[i] > MAX_CLIQUE){
                if(checkKmE_arc(newDist, i, avoidSize[i] - MAX_CLIQUE, row, column) == 0){
                        try_distance(nextX, column, row, flipRemains);
                    }
                
            }else{
                int v0 = blockStart[row];
  
                if(checkKClique(ADJ(v0,i)&ADJ(blockStart[column]+newDist,i) , i, avoidSize[i] - 2) == 0){

                    try_distance(nextX,column, row, flipRemains);
                }
            }

            if(row==column){
                REMOVE(distances[i][row][column], x);
                REMOVE(distances[i][row][column], blockSize - x);
                removeNewEdge(row, column, x, i);
            }else{
                REMOVE(distances[i][row][column], blockSize-x-1);
                removeNewEdge(row, column, blockSize-x-1, i);
            }
        }
 

}

//Set-up all necessary parameters
void searchRamsey(){
    if(n%k!=0){
        count = 0;
        return;
    }
    
    fprintf(stderr, "n=%d\n",n);
    fprintf(stderr, "k=%d\n",k);

    adjs = calloc(colorCount,sizeof(Set*));
    for (int i = 0; i < colorCount; i++)
    {
        adjs[i] = calloc(n,sizeof(Set));
    }

    distances = calloc(colorCount,sizeof(Set*));
    for (int i = 0; i < colorCount; i++)
    {
        distances[i] = calloc(k,sizeof(Set*));
        for (int j = 0; j < k; j++)
        {
            distances[i][j] = calloc(k,sizeof(Set));
        }
        
    }
    
    blockSize = n/k;
    m = blockSize;
    blockStart = calloc(k,sizeof(short));
    halfBlock = blockSize/2;
    rowX = calloc(n,sizeof(short));


    blockStart[0] = 0;
    for (int i = 1; i < k; i++){
        blockStart[i] = i*n/k;
    }

    for (int i = 0; i < n; i++)
    {
        rowX[i] = calcRow(i);
    }

    graphs = calloc(colorCount, sizeof(Graph));
    for (int i = 0; i < colorCount; i++){
        graphs[i].n = 0;
        graphs[i].adj = calloc(n+1,sizeof(Set));
    }


    cyclAuts = calloc(k-2, sizeof(Set*));
    for (int i = 2; i < k; i++){
        cyclAuts[i-2] = calloc(i-1, sizeof(Set));
    }
    

    calcUnits(blockSize);
    fixUnits = calloc(k,sizeof(Set));
    

    count = 0;
    found = 0;
    splitCount = 0;
    callCount = 0LL;
    pruneCount = 0LL;
    
    try_distance(1, 0, 0, 1);
    
    for (int i = 0; i < colorCount; i++){
        free(adjs[i]);
    }
    free(adjs);

    
    for (int i = 0; i < colorCount; i++){
        for (int j = 0; j < k; j++){
            free(distances[i][j]);
        }
        free(distances[i]);
        
    }
    free(distances);


    free(fixUnits);
    free(units);
}


#include <sys/times.h>
#include <unistd.h>
#define time_factor sysconf(_SC_CLK_TCK)

int main(int argc, char const *argv[])
{
    if (argc < 5) {
        fprintf(stderr,"Provide number of vertices, k, and the Ramsey parameters.\n");
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
            case 'B':
            {
                //ugly
                int avoid1 = ((int) strtol(&argv[3+i][1], NULL, 10));
                int avoid2 = ((int) strtol(&argv[3+i][2], NULL, 10));
                avoid1 = avoid2<10?(avoid1-avoid2)/10:(avoid1-avoid2)/100;
                biPar[i][0] = avoid1;
                biPar[i][1] = avoid2;
                avoidSize[i] = avoid1+avoid2 + 5 * MAX_CLIQUE;
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


        fprintf(stderr, "calls: %lld\n",callCount);

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
