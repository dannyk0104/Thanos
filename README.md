# Thanos

Thanos is a fast graph partitioning tool that uses cross-decomposition algorithm.
Thanos runs on NVIDIA GPUs.

## Getting Started

```
git clone https://github.com/dannyk0104/Thanos
cd Thanos
make
```

If you want to change the partition size, open Thanos.cu and change the
'#define PARTITION 16' to whatever partition size you want to change.

When you run the tool, it will print out the two results on the terminal.
The first one is the result of random partition before running cross-decomposition.
The second one is the result of running cross-decomposition.
