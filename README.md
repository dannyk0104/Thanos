# Thanos

Thanos is a fast graph partitioning tool that uses cross-decomposition algorithm.  
Thanos runs on NVIDIA GPUs.

## Getting Started

```
git clone https://github.com/dannyk0104/Thanos
cd Thanos
make
./thanos /path-to-tsv-file
```

If you want to change the partition size,  
open Thanos.cu and change the  
'#define PARTITION 16'  
to whatever partition size you want to change.

When you run the tool, it will print out the two results on the terminal.  
The first one is the result of random partition before running cross-decomposition.  
The second one is the result of running cross-decomposition.

## Input File Format

Thanos takes tsv files as an input.  
It assumes the vertex numbering starts from 1 not 0.  
The tsv file should be formatted as follows.  
destination source weight  
The file must be sorted with source first and the destinations for the source must be sorted too.
