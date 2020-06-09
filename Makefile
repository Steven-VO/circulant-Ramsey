

circulant:
	rm -rf ./genCyc
	gcc -O3 ./generators/genCyc.c -o genCyc -march=native -mpopcnt

allCirculant:
	rm -rf ./genCyc
	gcc -O3 ./generators/genCyc.c -o genCyc -march=native -mpopcnt -DSEARCH_SINGLE=0 -DNULCOUNTS=4 -DSETSIZE=64

largeCirculant:
	rm -rf ./genCyc
	gcc -O3 ./generators/genCyc.c -o genCyc -march=native -mpopcnt -DSEARCH_SINGLE=1 -DNULCOUNTS=4 -DSETSIZE=128
