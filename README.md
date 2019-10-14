# Thanos

As graphs become larger and more complex, it is becoming nearly impossible to process them without graph partitioning. Graph partitioning creates many subgraphs which can be processed in parallel thus delivering high-speed computation results. However, graph partitioning is a difficult task. In this open-source package, we introduce Thanos, a fast graph partitioning tool which uses the cross-decomposition algorithm that iteratively partitions a graph. It also produces balanced loads of partitions. The algorithm is well suited for parallel GPU programming which leads to fast and high-quality graph partitioning solutions. Experimental results show that we have achieved a 30x speedup and 35% better edge cut reduction compared to the CPU version of METIS on average.

## Getting Started

```
git clone https://github.com/dannyk0104/Thanos
cd Thanos
make
./thanos /path-to-tsv-file
```

If you want to change the partition size,  
open Thanos.cu and change the

```
#define PARTITION 16
```

to whatever partition size you want to change.

When you run the tool, it will print out the two results on the terminal.  
The first one is the result of random partition before running cross-decomposition.  
The second one is the result of running cross-decomposition.

## Input File Format

Thanos takes tsv files as an input.  
It assumes the vertex numbering starts from 1 not 0.  
The tsv file should be formatted as follows.

```
destination source weight
```

The file must be sorted with source first and the destinations for the source must be sorted too.
