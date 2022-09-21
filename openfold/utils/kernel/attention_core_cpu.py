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
    v_t = permute_final_dims(v, (1, 0))
    dy_buf = torch.matmul(grad_output, v_t)
    tmp_sum = torch.sum(attention_logits * dy_buf, dim=-1)
    attention_logits.mul_(dy_buf.sub_(tmp_sum.unsqueeze(-1)))
