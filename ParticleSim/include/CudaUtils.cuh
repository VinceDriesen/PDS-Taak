#ifndef CUDA_UTILS_CUH
#define CUDA_UTILS_CUH

#include <cuda_runtime.h>
#include <algorithm>
#include <cstdio>
#include <string>

struct LaunchConfig 
{
    int gridSize;
    int blockSize;
};

class CudaUtils 
{
public:
    struct Timer 
    {
        cudaEvent_t _start, _stop;

        Timer() {
            cudaEventCreate(&_start);
            cudaEventCreate(&_stop);
        }

        ~Timer() {
            cudaEventDestroy(_start);
            cudaEventDestroy(_stop);
        }

        void start(cudaStream_t stream = 0) {
            cudaEventRecord(_start, stream);
        }

        void stop(cudaStream_t stream = 0) {
            cudaEventRecord(_stop, stream);
        }

        float elapsed() {
            cudaEventSynchronize(_stop);
            float milliseconds = 0;
            cudaEventElapsedTime(&milliseconds, _start, _stop);
            return milliseconds;
        }
    };

    static int getSMCount()
    {
        int deviceId;
        cudaGetDevice(&deviceId);
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, deviceId);
        return prop.multiProcessorCount;
    }

    static std::string getDeviceName()
    {
        int deviceId;
        cudaGetDevice(&deviceId);
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, deviceId);
        return std::string(prop.name);
    }

    static void setGridLimit(int limit)
    {
        _getStoredGridLimit() = limit;
    }

    template <typename KernelFunc>
    static LaunchConfig getOptimalConfig(KernelFunc kernel, int numElements, const char* kernelName = "Unknown Kernel") 
    {
        int minGridSize;
        int blockSize;

        cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, kernel, 0, numElements);
        
        int optimalGrids = (numElements + blockSize - 1) / blockSize;
        int finalGridSize = optimalGrids;
        int limit = _getStoredGridLimit();

        if (limit > 0)
        {
            finalGridSize = std::min(optimalGrids, limit);
        }

        printf("Kernel: %s | Grid: %d | Block: %d\n", kernelName, finalGridSize, blockSize);

        return {finalGridSize, blockSize};
    }

private:
    static int& _getStoredGridLimit() {
        static int limit = 0;
        return limit;
    }
};

#endif // CUDA_UTILS_CUH