#ifndef CUDA_UTILS_CUH
#define CUDA_UTILS_CUH

#include <cuda_runtime.h>
#include <algorithm>

struct LaunchConfig 
{
    int gridSize;
    int blockSize;
};

class CudaUtils 
{
public:

    static void setSMMultiplier(int multiplier)
    {
        if (multiplier > 0)
        {
            int deviceId;
            cudaGetDevice(&deviceId);
            cudaDeviceProp prop;
            cudaGetDeviceProperties(&prop, deviceId);
            
            _getStoredSMLimit() = prop.multiProcessorCount * multiplier;
        }
        else
        {
            _getStoredSMLimit() = 0;
        }
    }

    template <typename KernelFunc>
    static LaunchConfig getOptimalConfig(KernelFunc kernel, int numElements) 
    {
        int minGridSize;
        int blockSize;

        cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, kernel, 0, numElements);
        int optimalGrids = (numElements + blockSize - 1) / blockSize;
        int finalGridSize = optimalGrids;

        int smLimit = _getStoredSMLimit();

        if (smLimit > 0)
        {
            finalGridSize = std::min(optimalGrids, smLimit);
        }

        return {finalGridSize, blockSize};
    }

private:
    static int& _getStoredSMLimit() {
        static int limit = 0;
        return limit;
    }
};

#endif // CUDA_UTILS_CUH