#ifndef CUDA_UTILS_CUH
#define CUDA_UTILS_CUH

#include <cuda_runtime.h>

struct LaunchConfig 
{
    int gridSize;
    int blockSize;
};

class CudaUtils 
{
public:
    template <typename KernelFunc>
    static LaunchConfig getOptimalConfig(KernelFunc kernel, int numElements) {
        int minGridSize;
        int blockSize;
        
        cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, kernel, 0, numElements);
        
        int gridSize = (numElements + blockSize - 1) / blockSize;

        return {gridSize, blockSize};
    }
};

#endif // CUDA_UTILS_CUH