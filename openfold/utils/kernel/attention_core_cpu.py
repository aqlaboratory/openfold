import torch
from openfold.utils.tensor_utils import permute_final_dims

def flat(t, keep_dims):
    return t.reshape([-1] + list(t.shape[-keep_dims:]))

def forward_(attention_logits, rows, cols):
    chunk_size = 512 * 1024

    for al_c in torch.split(flat(attention_logits, keep_dims=1), chunk_size, dim=0):
        tmp = torch.nn.functional.softmax(al_c, dim=-1)
        al_c.copy_(tmp)

def backward_(attention_logits, #output
              grad_output, #d_ov
              v, #values
              rows,
              cols_output,
              cols_values):
    chunk_size = 2048

    for al_c, go_c, v_c in zip(torch.split(flat(attention_logits, keep_dims=2), chunk_size, dim=0),
                               torch.split(flat(grad_output, keep_dims=2), chunk_size, dim=0),
                               torch.split(flat(v, keep_dims=2), chunk_size, dim=0)):
        v_t = permute_final_dims(v_c, (1, 0))
        dy_buf = torch.matmul(go_c, v_t)
        tmp_sum = torch.sum(al_c * dy_buf, dim=-1)
        al_c.mul_(dy_buf.sub_(tmp_sum.unsqueeze(-1)))
