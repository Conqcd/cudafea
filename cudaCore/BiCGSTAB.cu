#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include<thrust/host_vector.h>
#include<thrust/device_vector.h>
#include<thrust/copy.h>
#include<thrust/sort.h>
#include<thrust/execution_policy.h>


void Try()
{
	
	dim3 blockSize(16 ,16);
	dim3 threadSize(16, 16);
	
	
	// GPU_Tracer << <blockSize, threadSize >> > (samples_per_pixel, model.dev_image,
	//                                            background, model, depth, screen.camera,infinity);
    // HANDLE_ERROR(cudaMalloc((void**)&dev_ids, sizeof(tinyobj::index_t) * 3 * numtri));

	// HANDLE_ERROR(
	// 	cudaMemcpy(model.image, model.dev_image, sizeof(DoubleColor) * model.height * model.width,cudaMemcpyDeviceToHost
	// 	));

}