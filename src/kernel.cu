__global__
void filterScore(float *R, int m, int n, float *threshold, int* occIdx, float* occScore, int *nOcc)
{
  	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
	int idx = j*m + i; // column-major storage

  	if (i < m && j < n && R[idx] >= threshold[i]) { // FIXME: consistency with other filters (threshold)
		int resPos = atomicAdd(nOcc, 1);
		occScore[resPos] = R[idx];
		occIdx[resPos] = idx;
	}
}

void kernel_wrapper(float *d_R, int m, int n, float *d_threshold, int* d_occIdx, float* d_occScore, int *d_nOcc)
{
	dim3 threadsPerBlock(32, 32);
	dim3 numBlocks((m + threadsPerBlock.x - 1) / threadsPerBlock.x,
                       (n + threadsPerBlock.y - 1) / threadsPerBlock.y);

	filterScore<<<numBlocks, threadsPerBlock>>>(d_R, m, n, d_threshold, d_occIdx, d_occScore, d_nOcc);
}
