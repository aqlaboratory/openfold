#include <torch/extension.h>

// modified from fastfold/model/fastnn/kernel/cuda_native/csrc/softmax_cuda.cpp

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
