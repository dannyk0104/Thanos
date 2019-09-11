# Thanos

As size of graphs is getting complex and larger, processing such graphs become practically impossible to work on CPUs.  
Rencently, GPUs are proven to be fit for such purpose.  
However, utilizing a GPU to process a large graph is still considered to be inefficient and slow.  
Processing them with subgraphs by partitioning has become important for many applications in areas of computing.  
When a graph is partitioned into multiple sub-graphs, those sub-graphs should be in equal sizes that fit into each GPU with maximum size to maximize GPU utilization.  
Also, it is crucial to reduce number of memory access outside of each sub-graph to avoid global memory access.  
Not only load balancing and quality of partition is important but also partitioning time is very important.  
To achieve all those goals, we introduce Cross-Decomposition algorithm that iteratively partitions a graph.  
The algorithm suits very well for parallal GPU programming which leads to fast and high quality graph partitioning.

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
