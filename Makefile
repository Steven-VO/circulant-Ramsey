default: circulant block

circulant:
	rm -rf ./genCyc
	gcc -O3 ./generators/genCyc.c -o genCyc -march=native -mpopcnt

allCirculant:
	rm -rf ./genCyc
	gcc -O3 ./generators/genCyc.c -o genCyc -march=native -mpopcnt -DSEARCH_SINGLE=0 -DNULCOUNTS=4 -DSETSIZE=64

largeCirculant:
	rm -rf ./genCyc
	gcc -O3 ./generators/genCyc.c -o genCyc -march=native -mpopcnt -DSEARCH_SINGLE=1 -DNULCOUNTS=4 -DSETSIZE=128

block:
	rm -rf ./genBlock
	gcc -O3 ./generators/genBlock.c -o genBlock -march=native -mpopcnt -DSEARCH_SINGLE=1 -DNULCOUNTS=3 -DSETSIZE=64

allBlock:
	rm -rf ./genBlock
	gcc -O3 ./generators/genBlock.c -o genBlock -march=native -mpopcnt -DSEARCH_SINGLE=0 -DNULCOUNTS=3 -DSETSIZE=64

largeBlock:
	rm -rf ./genBlock
	gcc -O3 ./generators/genBlock.c -o genBlock -march=native -mpopcnt -DSEARCH_SINGLE=1 -DNULCOUNTS=3 -DSETSIZE=128

clean:
	rm -rf ./genCyc
	rm -rf ./genBlock
