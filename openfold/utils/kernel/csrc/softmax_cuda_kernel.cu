// Copyright 2021 AlQuraishi Laboratory
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// modified from fastfold/model/fastnn/kernel/cuda_native/csrc/softmax_cuda_kernel.cu

#include <math_constants.h>
#include <torch/extension.h>
#include <c10/cuda/CUDAGuard.h>

#include <iostream>

#include "ATen/ATen.h"
#include "ATen/cuda/CUDAContext.h"
#include "compat.h"

#define CHECK_CUDA(x) TORCH_CHECK(x.is_cuda(), #x " must be a CUDA tensor")
#define CHECK_CONTIGUOUS(x) TORCH_CHECK(x.is_contiguous(), #x " must be contiguous")
#define CHECK_INPUT(x) \
    CHECK_CUDA(x);     \
    CHECK_CONTIGUOUS(x)

__inline__ __device__ float WarpAllReduceMax(float val) {
    for (int mask = 1; mask < 32; mask *= 2) {
        val = max(val, __shfl_xor_sync(0xffffffff, val, mask));
    }
    return val;
}

__inline__ __device__ float WarpAllReduceSum(float val) {
    for (int mask = 1; mask < 32; mask *= 2) {
        val += __shfl_xor_sync(0xffffffff, val, mask);
    }
    return val;
}


template<typename T>
__global__ void attn_softmax_inplace_(
    T *input, 
    long long rows, int cols
) {
    int threadidx_x = threadIdx.x / 32;
    int threadidx_y = threadIdx.x % 32;
    long long row_offset = (long long)(blockIdx.x * 4 + threadidx_x);
    int cols_per_thread = (cols + 31) / 32;
    int cols_this_thread = cols_per_thread;

    int last_y = (cols / cols_per_thread);

    if (threadidx_y == last_y) {
        cols_this_thread = cols - cols_per_thread * last_y;
    }
    else if (threadidx_y > last_y) {
        cols_this_thread = 0;
    }

    float buf[32];

    int lane_id = threadidx_y;

    if (row_offset < rows) {
        T *row_input = input + row_offset * cols;
        T *row_output = row_input;

        #pragma unroll
        for (int i = 0; i < cols_this_thread; i++) {
            int idx = lane_id * cols_per_thread + i;
            buf[i] = static_cast<float>(row_input[idx]);
        }

        float thread_max = -1 * CUDART_INF_F;
        #pragma unroll
        for (int i = 0; i < cols_this_thread; i++) {
            thread_max = max(thread_max, buf[i]);
        }

        float warp_max = WarpAllReduceMax(thread_max);

        float thread_sum = 0.f;
        #pragma unroll
        for (int i = 0; i < cols_this_thread; i++) {
            buf[i] = __expf(buf[i] - warp_max);
            thread_sum += buf[i];
        }

        float warp_sum = WarpAllReduceSum(thread_sum);
        #pragma unroll
        for (int i = 0; i < cols_this_thread; i++) {
            row_output[lane_id * cols_per_thread + i] =
                static_cast<T>(__fdividef(buf[i], warp_sum));
        }
    }
}


void attn_softmax_inplace_forward_(
    at::Tensor input, 
    long long rows, int cols
) {
    CHECK_INPUT(input);
    const at::cuda::OptionalCUDAGuard device_guard(device_of(input));

    int grid = (rows + 3) / 4;
    dim3 block(128);

    if (input.dtype() == torch::kFloat32) {
        attn_softmax_inplace_<float><<<grid, block>>>(
            (float *)input.data_ptr(),
            rows, cols
        );
    } 
    else {
        attn_softmax_inplace_<at::BFloat16><<<grid, block>>>(
            (at::BFloat16 *)input.data_ptr(), 
            rows, cols
        );
    }
}


template<typename T>
__global__ void attn_softmax_inplace_grad_(
    T *output,
    T *d_ov,
    T *values,
    long long rows, 
    int cols_output,
    int cols_values
) {
    int threadidx_x = threadIdx.x / 32;
    int threadidx_y = threadIdx.x % 32;
    long long row_offset = (long long)(blockIdx.x * 4 + threadidx_x);
    int cols_per_thread = (cols_output + 31) / 32;
    int cols_this_thread = cols_per_thread;
    int rows_values = cols_output;
    // values are set to the beginning of the current 
    // rows_values x cols_values leaf matrix
    long long value_row_offset = row_offset - row_offset % rows_values;
    int last_y = (cols_output / cols_per_thread);

    if (threadidx_y == last_y) {
        cols_this_thread = cols_output - cols_per_thread * last_y;
    }
    else if (threadidx_y > last_y) {
        cols_this_thread = 0;
    }

    float y_buf[32];
    float dy_buf[32];

    int lane_id = threadidx_y;

    if (row_offset < rows) {
        T *row_output = output + row_offset * cols_output;
        T *row_d_ov = d_ov + row_offset * cols_values;
        T *row_values = values + value_row_offset * cols_values;

        float thread_max = -1 * CUDART_INF_F;

        // Compute a chunk of the output gradient on the fly
        int value_row_idx = 0;
        int value_idx = 0;
        #pragma unroll
        for (int i = 0; i < cols_this_thread; i++) {
            T sum = 0.;
            #pragma unroll
            for (int j = 0; j < cols_values; j++) {
                value_row_idx = ((lane_id * cols_per_thread) + i);
                value_idx = value_row_idx * cols_values + j;
                sum += row_d_ov[j] * row_values[value_idx];
            }
            dy_buf[i] = static_cast<float>(sum);
        }

        #pragma unroll
        for (int i = 0; i < cols_this_thread; i++) {
            y_buf[i] = static_cast<float>(row_output[lane_id * cols_per_thread + i]);
        }

        float thread_sum = 0.;

        #pragma unroll
        for (int i = 0; i < cols_this_thread; i++) {
            thread_sum += y_buf[i] * dy_buf[i];
        }

        float warp_sum = WarpAllReduceSum(thread_sum);

        #pragma unroll
        for (int i = 0; i < cols_this_thread; i++) {
            row_output[lane_id * cols_per_thread + i] = static_cast<T>(
                    (dy_buf[i] - warp_sum) * y_buf[i]
            );
        }
    }
}


void attn_softmax_inplace_backward_(
    at::Tensor output,
    at::Tensor d_ov, 
    at::Tensor values,
    long long rows, 
    int cols_output,
    int cols_values
) {
    CHECK_INPUT(output);
    CHECK_INPUT(d_ov);
    CHECK_INPUT(values);
    const at::cuda::OptionalCUDAGuard device_guard(device_of(output));

    int grid = (rows + 3) / 4;
    dim3 block(128);

    if (output.dtype() == torch::kFloat32) {
        attn_softmax_inplace_grad_<float><<<grid, block>>>(
            (float *)output.data_ptr(),
            (float *)d_ov.data_ptr(), 
            (float *)values.data_ptr(),
            rows, cols_output, cols_values
        );
    } else {
        attn_softmax_inplace_grad_<at::BFloat16><<<grid, block>>>(
            (at::BFloat16 *)output.data_ptr(),
            (at::BFloat16 *)d_ov.data_ptr(), 
            (at::BFloat16 *)values.data_ptr(),
            rows, cols_output, cols_values
        );
    }
}
