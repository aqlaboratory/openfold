import copy
import os

import torch
import deepspeed


local_rank = int(os.getenv('LOCAL_RANK', '0'))
world_size = int(os.getenv('WORLD_SIZE', '1'))

class Model(torch.nn.Module):
    def __init__(self):
        super().__init__()

        self.ml = torch.nn.ModuleList()
        for _ in range(4000):
            self.ml.append(torch.nn.Linear(500, 500))

    def forward(self, batch):
        for i, l in enumerate(self.ml):
            # print(f"{i}: {l.weight.device}")
            batch = l(batch)

        return batch


class DummyDataset(torch.utils.data.Dataset):
    def __init__(self):
        self.batch = torch.rand(500, 500) 

    def __getitem__(self, idx):
        return copy.deepcopy(self.batch)

    def __len__(self):
        return 1000

dd = DummyDataset()
dl = torch.utils.data.DataLoader(dd)
example = next(iter(dl)).to(f"cuda:{local_rank}")

model = Model()
model = model.to(f"cuda:{local_rank}")

model = deepspeed.init_inference(
    model,
    mp_size=world_size,
    checkpoint=None,
    replace_method=None,
    #replace_method="auto"
)

out = model(example)
#if not torch.distributed.is_initialized() or torch.distributed.get_rank() == 0:
#    print(out)
