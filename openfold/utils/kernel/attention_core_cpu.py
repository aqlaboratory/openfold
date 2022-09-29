import torch
from openfold.utils.tensor_utils import permute_final_dims

def forward_(attention_logits, rows, cols):
    chunk_dim = -4
    chunk_size = 256

    for al_c in torch.split(attention_logits, chunk_size, dim=chunk_dim):
        tmp = torch.nn.functional.softmax(al_c, dim=-1)
        al_c.copy_(tmp)

def backward_(attention_logits, #output
              grad_output, #d_ov
              v, #values
              rows,
              cols_output,
              cols_values):

    seq_dim = -4 # [1, 5120, 8, 384, 384] #右から3番目以降の次元で連続して分割可能な次元
    chunk_size = 5120

    for al_c, go_c, v_c in zip(torch.split(attention_logits, chunk_size, dim=seq_dim),
                               torch.split(grad_output, chunk_size, dim=seq_dim),
                               torch.split(v, chunk_size, dim=seq_dim)):
        v_t = permute_final_dims(v_c, (1, 0))
        dy_buf = torch.matmul(go_c, v_t)
        tmp_sum = torch.sum(al_c * dy_buf, dim=-1)
        al_c.mul_(dy_buf.sub_(tmp_sum.unsqueeze(-1)))
