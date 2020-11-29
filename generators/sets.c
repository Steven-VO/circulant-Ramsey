#include<stdio.h>

#ifndef SETS_H
#define SETS_H

#ifndef SETSIZE
    #define SETSIZE 64
#endif

#if SETSIZE==128
    typedef unsigned __int128 Set;
#else
    typedef unsigned long long Set;
#endif

#define S1 ((Set) (1ULL))
#define S0 ((Set) (0ULL))
#define CONTAINS(set,i) ((set) & (S1<<(i)))
#define ADD(set,i) ((set) |= (S1<<(i)))
#define REMOVE(set,i) ((set) &= ~(S1<<(i)))
#define WITH(set,i) ((set) | (S1<<(i)))
#define WITHOUT(set,i) ((set) & ~(S1<<(i)))
#define BIT(i) (S1<<(i))
#define EMPTY(set) ((set)==S0)

#define RANGE_D(i) (BIT(i)-1)
#define XRANGE(min,max) (RANGE_D(max) & (~RANGE_D(min)))
#define RANGE(i) (i>SETSIZE-2 ? (~S0) : (BIT(i+1)-1))
#define FROMRANGE(i) (~RANGE_D(i))
#define COMPLEMENT(set,max) (~(set) & RANGE_D(max))

    //dangerous and unchecked
#define ROTATE_D(set, amount, n) (((set) << (amount) | (set) >> (n - (amount))) & RANGE(n-1))

    //safer
#define ROTATE(set, amount, n) ( (((set) << (amount)) | (((set) >> (n - (amount))) & RANGE((amount)-1))) & RANGE(n-1) )

#define POPCOUNT_LL(set) (__builtin_popcountll(set))

#if SETSIZE==128
    #define POPCOUNT(set) (POPCOUNT_LL(set) + POPCOUNT_LL(set>>64))
#else
    #define POPCOUNT(set) POPCOUNT_LL(set)
#endif

#if SETSIZE==128
   #define FIRSTSET(set) ((unsigned long long)(set) == S0 ? 64+__builtin_ctzll(set>>64) : __builtin_ctzll(set))
#else
    #define FIRSTSET(set) (__builtin_ctzll(set))
#endif

//changes the set!
#define ITERATE_SET(set, index) while((set)!=S0){index = FIRSTSET(set); set^=BIT(index);


void printint(int a){
    fprintf(stderr,"%d\n",a);
}

void printSet(Set a){
    for (int i = 8*sizeof(Set)-1;i>=0; i--)
    {
        fprintf(stderr,"%d",(unsigned int)((a & BIT(i) )>> i));
    }
    fprintf(stderr,"\n");
    
}

void printReverseSet(Set a){
    for (int i = 0; i < 8*sizeof(Set); i++)
    {
        fprintf(stderr,"%d",(unsigned int)((a & BIT(i) )>> i));
    }
    fprintf(stderr,"\n");
}

#endif /*SETS_H*/