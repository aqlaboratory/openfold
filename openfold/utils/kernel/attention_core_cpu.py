import torch
from openfold.utils.tensor_utils import permute_final_dims

def forward_(attention_logits, rows, cols):
    col_max, _ = torch.max(attention_logits, dim=-1)
    attention_logits.sub_(col_max.unsqueeze(-1))
    torch.exp(attention_logits, out=attention_logits)
    col_sum = torch.sum(attention_logits, dim=-1)
    attention_logits.div_(col_sum.unsqueeze(-1))

def backward_(attention_logits, #output
              grad_output, #d_ov
              v, #values
              rows,
              cols_output,
              cols_values):

    seq_dim = -4
    chunk_size = 128
    for al_c, go_c, v_c in zip(torch.split(attention_logits, chunk_size, dim=seq_dim),
                               torch.split(grad_output, chunk_size, dim=seq_dim),
                               torch.split(v, chunk_size, dim=seq_dim)):
        v_t = permute_final_dims(v_c, (1, 0))
        dy_buf = torch.matmul(go_c, v_t)
        tmp_sum = torch.sum(al_c * dy_buf, dim=-1)
        al_c.mul_(dy_buf.sub_(tmp_sum.unsqueeze(-1)))
