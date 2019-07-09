#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include <stdlib.h>
#include <cstdint>
#include <random>
#include <algorithm>
#include <string>
#include <sstream>
#include <iterator>
#include <map>
#include <string.h>
#include <ctime>
#include <chrono>
#include <math.h>
#include <assert.h>
#include <cstring>

using namespace std::chrono;
using namespace std;

/* change the number if you want to partition with different size*/
#define PARTITION 4


/* !Kernel code that runs the cross-decomposition algorithm
*/
__global__ void CrossDecomposition_kernel(  uint8_t* orig_P, uint8_t* new_P, 
                                            uint64_t* cardi, uint64_t* new_cardi,
                                            const uint64_t num_nodes, const float h, uint64_t capacity,
                                            uint32_t* coo_row, uint32_t* coo_col, 
                                            uint64_t* row_ptr, bool is_divisible)
{
    uint64_t idx = blockIdx.x*blockDim.x+threadIdx.x;
	if(idx < num_nodes){		
        uint64_t cur_node = idx;
        uint64_t connected_and_in_curpart[PARTITION] = {0};
        uint64_t degree = (row_ptr[cur_node+1] - row_ptr[cur_node]);

        for(int j=0; j<degree; j++){
            uint64_t cur_col_ptr = row_ptr[cur_node] + j;
            uint64_t cur_col = coo_col[cur_col_ptr];
            for (int cur_part = 0; cur_part < PARTITION; cur_part++) {
                if (orig_P[cur_col] == cur_part)
                    connected_and_in_curpart[cur_part] += 1;
            }
        }
        
        float cost[PARTITION] = {0};
        for (int i = 0; i < PARTITION; i++){
            cost[i] = h* connected_and_in_curpart[i] + 
                        (1-h)*(num_nodes - (cardi[i] + degree - connected_and_in_curpart[i]));
        } 

        //initialize arg_sort array
        uint8_t arg_sort[PARTITION];
        for (uint8_t i = 0; i < PARTITION; i++)
            arg_sort[i] = i;

        for (int i = 0; i < PARTITION-1; i++){       
            for (int j = 0; j < PARTITION-i-1; j++){
                if (cost[j] < cost[j+1]){
                    float temp = cost[j];
                    cost[j] = cost[j+1];
                    cost[j+1] = temp;
                    int temp2 = arg_sort[j];
                    arg_sort[j] = arg_sort[j+1];
                    arg_sort[j+1] = temp2;
                }
            }
        }

        unsigned long long int old_size;
        for (int i = 0; i < PARTITION; i++){
            if(!is_divisible){
                if(arg_sort[i]==PARTITION-1){
                    old_size = atomicAdd((unsigned long long int*)&new_cardi[arg_sort[i]],(unsigned long long int) 1);
                    if(old_size >= capacity+(num_nodes%PARTITION)){
                        old_size = atomicSub((unsigned int*)&new_cardi[arg_sort[i]],(unsigned int) 1);
                    }
                    else{
                        new_P[cur_node] = arg_sort[i];
                        break;
                    }
                }
                else{
                    old_size = atomicAdd((unsigned long long int*)&new_cardi[arg_sort[i]],(unsigned long long int) 1);
                    if(old_size >= capacity){
                        old_size = atomicSub((unsigned int*)&new_cardi[arg_sort[i]],(unsigned int) 1);
                    }
                    else{
                        new_P[cur_node] = arg_sort[i];
                        break;
                    }
                }
            }
            else{
                old_size = atomicAdd((unsigned long long int*)&new_cardi[arg_sort[i]],(unsigned long long int) 1);
                if(old_size >= capacity){
                    old_size = atomicSub((unsigned int*)&new_cardi[arg_sort[i]],(unsigned int) 1);
                }
                else{
                    new_P[cur_node] = arg_sort[i];
                    break;
                }
            }
        }
    }
}

/* !Kernel to count the number of edges to evaluate the quality of partition
*/
__global__ void evalEdges(uint8_t* P, 
                          uint32_t* coo_row, uint32_t* coo_col, 
                          uint64_t* edges_per_part, uint64_t coo_size)
{
    uint64_t idx = blockIdx.x*blockDim.x+threadIdx.x;
    if(idx<coo_size){
        uint32_t src = coo_row[idx];
        uint32_t dest = coo_col[idx];
        if(src != dest){
            uint8_t src_part = P[src];
            uint8_t dest_part = P[dest];
            atomicAdd( (unsigned long long int*)&edges_per_part[src_part*PARTITION+dest_part], (unsigned long long int) 1);
        }
    }
    return;
} 

/* !Simple function to check CUDA runtime error
*/
void checkCuda(cudaError_t result) {
    if (result != cudaSuccess) {
        fprintf(stderr, "CUDA Runtime Error: %s\n",
            cudaGetErrorString(result));
        assert(result == cudaSuccess);
    }
}

/* !Thanos class definition
*/
class Thanos{
    private:
        float h = 0.9;
        bool is_divisible = true; 
        uint64_t num_nodes, edgecount, coo_size, capacity, rem;
        /*  following vectors are used to read files and create COO + CSR row pointer.
            They will be copied to GPU memories   */
        std::vector<uint32_t> edge_vec_src, edge_vec_dest; 
        std::vector<uint64_t> row_ptrs; 

        uint8_t *row_P, *col_P; //row and column partition arrays
        uint32_t* coo_row, *coo_col; //COO format
        uint64_t* row_ptr; //row pointer of CSR format
        uint64_t *row_cardi, *row_new_cardi, *col_cardi, *col_new_cardi; //cardinality arrays
        uint64_t* edges_per_part; //this var is for evaluating the partition quality

        cudaError_t err;
        
        void readGraph_tsv(const char *filename);
        void readGraph_DARPA_CSR(   const char *filename,
                                    std::vector<uint32_t> &edge_vec_src,
                                    std::vector<uint32_t> &edge_vec_dest,
                                    std::vector<uint64_t> &row_ptrs,
                                    uint64_t &edgecount, uint64_t &nodecount);
        void AllocateGPUMem();
        void initMem();
        void initParts(uint8_t *P, uint64_t *cardi, const uint64_t num_nodes);
        void CrossDecomposition();
        void evaluatePartition();
        void printEdgesPerPar(uint64_t *edges_per_part);

    public:
        Thanos(){};
        Thanos(const char *filename);
        ~Thanos();

};

/* !Constructor for Thanos
    Construct Thanos with tsv file will run everyting for you
*/
Thanos::Thanos(const char *filename){
    readGraph_tsv(filename);
    AllocateGPUMem();
    initMem();
    evaluatePartition();
    CrossDecomposition();
    evaluatePartition();
}

/* !Destructor for Thanos. 
    Deallocates all the GPU memories
*/
Thanos::~Thanos(){
    cudaFree(coo_row);
    cudaFree(coo_col);
    cudaFree(row_ptr);
    cudaFree(row_P);
    cudaFree(col_P);
    cudaFree(row_cardi);
    cudaFree(row_new_cardi);
    cudaFree(col_cardi);
    cudaFree(col_new_cardi);
    cudaFree(edges_per_part);
}

/* !Host function to call the cross-decomposition kernel.
    You can change the boundary of for loop to control the number of iterations
*/
void Thanos::CrossDecomposition(){
    dim3 dimGrid(ceil(((float)num_nodes)/1024),1,1);
    dim3 dimBlock(1024,1,1);

    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    for(int i=0; i<3; i++){
        CrossDecomposition_kernel<<<dimGrid, dimBlock>>>(   row_P, col_P, col_cardi, col_new_cardi,
                                                            num_nodes, h, capacity,
                                                            coo_row, coo_col, 
                                                            row_ptr, is_divisible);

        CrossDecomposition_kernel<<<dimGrid, dimBlock>>>(   col_P, row_P, row_cardi, row_new_cardi,
                                                            num_nodes, h, capacity,
                                                            coo_row, coo_col, 
                                                            row_ptr, is_divisible);

        checkCuda(cudaMemcpy(row_cardi, row_new_cardi, PARTITION*sizeof(uint64_t), cudaMemcpyHostToHost));
        std::fill(row_new_cardi, row_new_cardi + PARTITION, 0);
        checkCuda(cudaMemcpy(col_cardi, col_new_cardi, PARTITION*sizeof(uint64_t), cudaMemcpyHostToHost));
        std::fill(col_new_cardi, col_new_cardi + PARTITION, 0);
    }
    
    err = cudaGetLastError();
    if (err != cudaSuccess)
        printf("Error1: %s\n", cudaGetErrorString(err));
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>( t2 - t1 ).count();
    cout<<"RUN Time in Sec(only kernel): "<<duration*pow(10,-6)<<endl;  
}

/* !Host function to call evalEdges kernel. After kernel call,
    it calls printEdgesPerPar function to show the result on the terminal
*/
void Thanos::evaluatePartition(){
    std::fill(edges_per_part, edges_per_part + PARTITION*PARTITION, 0);
    dim3 dimGrid0(ceil(((float)coo_size)/1024),1,1);
    dim3 dimBlock0(1024,1,1); //1024
    evalEdges<<<dimGrid0, dimBlock0>>>(row_P,coo_row, coo_col,
        edges_per_part, coo_size);
    err = cudaGetLastError();
    if (err != cudaSuccess)
        printf("Error1: %s\n", cudaGetErrorString(err));
    checkCuda(cudaDeviceSynchronize());
    printEdgesPerPar(edges_per_part);
}


/* !Function to allocate all the memories required for Thanos on GPU.
*/
void Thanos::AllocateGPUMem(){
    
    checkCuda(cudaMallocManaged((void**)&coo_row, coo_size*sizeof(uint32_t)));
    checkCuda(cudaMallocManaged((void**)&coo_col, coo_size*sizeof(uint32_t)));
    checkCuda(cudaMallocManaged((void**)&row_ptr, row_ptrs.size()*sizeof(uint64_t)));

    checkCuda(cudaMallocManaged((void**)&row_P, num_nodes*sizeof(uint8_t)));
    checkCuda(cudaMallocManaged((void**)&col_P, num_nodes*sizeof(uint8_t)));

    checkCuda(cudaMallocManaged((void**)&row_cardi, PARTITION*sizeof(uint64_t)));
    checkCuda(cudaMallocManaged((void**)&row_new_cardi, PARTITION*sizeof(uint64_t)));
    checkCuda(cudaMallocManaged((void**)&col_cardi, PARTITION*sizeof(uint64_t)));
    checkCuda(cudaMallocManaged((void**)&col_new_cardi, PARTITION*sizeof(uint64_t)));

    checkCuda(cudaMallocManaged((void**)&edges_per_part, PARTITION*PARTITION*sizeof(uint64_t)));
    
}

/* !Initialize the memories that are allocated in function Allocate GPU Mem
*/
void Thanos::initMem(){
    checkCuda(cudaMemcpy(coo_row, &edge_vec_src[0], coo_size*sizeof(uint32_t), cudaMemcpyHostToHost));
    checkCuda(cudaMemcpy(coo_col, &edge_vec_dest[0], coo_size*sizeof(uint32_t), cudaMemcpyHostToHost));
    checkCuda(cudaMemcpy(row_ptr, &row_ptrs[0], row_ptrs.size()*sizeof(uint64_t), cudaMemcpyHostToHost));
    std::fill(row_cardi, row_cardi + PARTITION, 0);
    std::fill(row_new_cardi, row_new_cardi + PARTITION, 0);
    std::fill(col_cardi, col_cardi + PARTITION, 0);
    std::fill(col_new_cardi, col_new_cardi + PARTITION, 0);
    initParts(row_P, row_cardi, num_nodes);
    initParts(col_P, col_cardi, num_nodes);
    
}

/* !Function to read the graph and updates
    variables capacity and remainder(rem)
*/
void Thanos::readGraph_tsv(const char *filename){
    readGraph_DARPA_CSR(filename, edge_vec_src, edge_vec_dest, row_ptrs, edgecount, num_nodes);
    coo_size = edge_vec_src.size();
    cout<<"Reading Graph Done, num_nodes = "<<num_nodes<<" COO size: "<<coo_size<<endl;

    capacity = floor(num_nodes/PARTITION);
    rem = num_nodes%PARTITION;
    if(rem!=0){
        cout<<"Size of Each Partition is: "<<capacity<<endl;
        cout<<"Size of Last Partition is: "<<capacity+rem<<endl;
        is_divisible = false;
    }
    else{
        cout<<"N is divided perfectly"<<endl;
        cout<<"Size of Each Partition is: "<<capacity<<endl;
    }
}


/* !Helper function to read the TSV file of graph. 
    The TSV file must be in format of
    destination src weight
    The file has to be sorted with src first and
    destination also has to be sorted. 
*/
void Thanos::readGraph_DARPA_CSR(const char *filename,
                         std::vector<uint32_t> &edge_vec_src,
                         std::vector<uint32_t> &edge_vec_dest,
                         std::vector<uint64_t> &row_ptrs,
                         uint64_t &edgecount, uint64_t &nodecount)
{
    int key, val, weight;
    std::ifstream ss(filename);
    std::vector<std::pair<int, long long int>> temp_row_ptrs_vec;
    edgecount = 0;
    nodecount = 0;
    int prevkey = -1;
    if (ss.is_open() && ss.good()){
        while (ss >> val){
            ss >> key;
            ss >> weight;
            nodecount = std::max<int>(nodecount, key);
            key--;
            val--;
            //if(key < val) {
            if (prevkey != key){
                prevkey = key;
                edge_vec_src.push_back(key);
                edge_vec_dest.push_back(key);
                temp_row_ptrs_vec.push_back(std::pair<uint32_t, uint64_t>(key, edgecount));
                edgecount++;
            }
            edge_vec_src.push_back(key);
            edge_vec_dest.push_back(val);
            edgecount++;
            //}
        }
        ss.close();
    }

    uint64_t *temp_row_ptrs = new uint64_t[nodecount + 1];
    std::fill(temp_row_ptrs, temp_row_ptrs + nodecount, -1);
    temp_row_ptrs[nodecount] = edgecount;

    std::vector<std::pair<int, long long int>>::iterator begin = temp_row_ptrs_vec.begin();
    std::vector<std::pair<int, long long int>>::iterator end = temp_row_ptrs_vec.end();

    for (std::vector<std::pair<int, long long int>>::iterator it = begin; it != end; ++it)
        temp_row_ptrs[it->first] = it->second;

    long long int cur_val = edgecount;
    for (int i = nodecount; i >= 0; i--){
        long long int val = temp_row_ptrs[i];
        if (val < 0)
            temp_row_ptrs[i] = cur_val;
        else
            cur_val = val;
    }
    row_ptrs.insert(row_ptrs.begin(), temp_row_ptrs, temp_row_ptrs + nodecount + 1);
    delete[] temp_row_ptrs;
}

/* !Initialize the partition with uniform distribution and
    update the cardinality array
*/
void Thanos::initParts( uint8_t *P, uint64_t *cardi,
                        const uint64_t num_nodes)
{
    random_device rd;
    mt19937 mt(rd());
    uniform_int_distribution<int> dist(0, PARTITION - 1);
    for (size_t i = 0; i < num_nodes; i++){
        uint64_t part = dist(mt);
        P[i] = part;
        cardi[part]++;
    }
}

/* !This function nicely outputs the result of evaluation to the terminal
    if you want to see the number of edges in each partition and between partitions,
    uncomment the couts
*/
void Thanos::printEdgesPerPar(uint64_t *edges_per_part)
{
    uint32_t total_internal_edges = 0, total_external_edges = 0;
    cout << "********************************************************" << endl;
    map<pair<uint8_t, uint8_t>, bool> track;
    for (uint8_t i = 0; i < PARTITION; i++){
        for (uint8_t j = 0; j < PARTITION; j++){
            if (i == j){
                // cout << "Internal Edges for Partition " << (int)i << " :" << (edges_per_part[i * PARTITION + j]) / 2 << endl;
                total_internal_edges += edges_per_part[i * PARTITION + j] / 2;
            }
            else if (track.find(make_pair(i, j)) == track.end()){
                // cout << "Between "
                //  << "PARTITION " << (int)i << " and " << (int)j << " Edges: " << edges_per_part[i * PARTITION + j] << endl;
                total_external_edges += (edges_per_part[i * PARTITION + j]);
                track[make_pair(i, j)] = true;
                track[make_pair(j, i)] = true;
            }
        }
    }
    // cout << "--------------------------------------------------------" << endl;
    cout << "Total Internal Edges: " << total_internal_edges << endl;
    cout << "Total External Edges: " << total_external_edges << endl;
    cout << "********************************************************" << endl;
    return;
}






int main(int argc, char **argv){
    /*simply construct Thanos with providing the path for the tsv file. 
    It will run everything for you*/
    Thanos th(argv[1]);
    return 0;
}

