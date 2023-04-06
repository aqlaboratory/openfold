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

// modified from fastfold/model/fastnn/kernel/cuda_native/csrc/softmax_cuda.cpp

#include <torch/extension.h>

void attn_softmax_inplace_forward_(
    at::Tensor input, 
    long long rows, int cols
);
void attn_softmax_inplace_backward_(
    at::Tensor output, 
    at::Tensor d_ov,
    at::Tensor values,
    long long rows, 
    int cols_output,
    int cols_values
);


PYBIND11_MODULE(TORCH_EXTENSION_NAME, m) {
    m.def(
        "forward_", 
        &attn_softmax_inplace_forward_, 
        "Softmax forward (CUDA)"
    );
    m.def(
        "backward_", 
        &attn_softmax_inplace_backward_, 
        "Softmax backward (CUDA)"
    );
}
