# Copyright 2021 AlQuraishi Laboratory
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import argparse
import json


parser = argparse.ArgumentParser(description='''Outputs a DeepSpeed 
                                                configuration file to 
                                                stdout''')

parser.add_argument("--gradient_clipping", type=float, default=None,
                    help="Value of gradient clipping")
p = parser.add_argument_group("Optimizer")
p.add_argument("--optimizer", default=None,
               help='''Choice of optimizer. Choose between "Adam" or 
                       "OneBitAdam"''')
p.add_argument("--lr", dest="lr", type=float, default=1e-3,
               help="The learning rate")
p.add_argument("--freeze_step", type=int, default=100,
               help='''Number of warm-up steps before 1-bit compression 
                       activates. Applies only when --optimizer is 
                       OneBitAdam''')
p.add_argument("--cuda_aware", action="store_true", default=False, 
               help='''Indicates that the underlying MPI library supports 
                       CUDA-Aware communication. Applies only when 
                       --optimizer is OneBitAdam''')
p.add_argument("--comm_backend_name", type=str, default="nccl",
               help='''Communication implementation for OneBitAdam. Choose
                       from nccl and mpi''')
p.add_argument("--eps", type=float, default=1e-8,
               help="Adam epsilon parameter")


sched = parser.add_argument_group("Scheduler")
sched.add_argument(
    "--scheduler", type=str, default=None,
    help='''The LR scheduler. Choose from "LRRangeTest", "OneCycle", WarmupLR, 
            and WarmupDecayLR". Documentation for each can be found here: 
            deepspeed.readthedocs.io/en/latest/schedulers.html'''
)
                    
range_test = sched.add_argument_group("LRRangeTest")
range_test.add_argument(
    "--lr_range_test_min_lr", type=float, default=1e-04
)
range_test.add_argument(
    "--lr_range_test_step_size", type=int, default=2000
)
range_test.add_argument(
    "--lr_range_test_step_rate", type=float, default=1.0
)
range_test.add_argument(
    "--lr_range_test_staircase", type=bool, default=False
)

cycle = sched.add_argument_group("OneCycle")
cycle.add_argument(
    "--cycle_min_lr", type=float, default=1e-06
)
cycle.add_argument(
    "--cycle_max_lr", type=float, default=1e-03
)
cycle.add_argument(
    "--cycle_decay_lr_rate", type=float, default=0
)
cycle.add_argument(
    "--cycle_first_step_size", type=int, default=2000
)
cycle.add_argument(
    "--cycle_second_step_size", type=int, default=None
)
cycle.add_argument(
    "--cycle_first_stair_count", type=int, default=0       
)
cycle.add_argument(
    "--cycle_second_stair_count", type=int, default=0     
)
cycle.add_argument(
    "--cycle_decay_step_size", type=int, default=0
)
cycle.add_argument(
    "--cycle_momentum", type=bool, default=True
)
cycle.add_argument(
    "--cycle_min_mom", type=float, default=0.8
)
cycle.add_argument(
    "--cycle_max_mom", type=float, default=0.9
)
cycle.add_argument(
    "--cycle_decay_mom_rate", type=float, default=0
)

warmup = sched.add_argument_group("WarmupLR")
warmup.add_argument(
    "--warmup_min_lr", type=float, default=0.
)
warmup.add_argument(
    "--warmup_max_lr", type=float, default=0.001
)
warmup.add_argument(
    "--warmup_num_steps", type=int, default=1000
)

warmup_decay = sched.add_argument_group("WarmupDecayLR")
warmup_decay.add_argument(
    "--warmup_decay_total_num_steps", type=int, default=1e05
)
warmup_decay.add_argument(
    "--warmup_decay_min_lr", type=float, default=0.
)
warmup_decay.add_argument(
    "--warmup_decay_max_lr", type=float, default=0.001
)
warmup_decay.add_argument(
    "--warmup_decay_num_steps", type=int, default=1000
)


p = parser.add_argument_group("Half-precision training (fp16)")
p.add_argument("--fp16", dest="fp16", action="store_true", default=False,
               help="""Whether to train in 16-bit/mixed-precision mode. 
                       Mutually exclusive with --amp""")

p = parser.add_argument_group("Half-precision training (bfloat16)")
p.add_argument("--bfloat16", dest="bfloat16", action="store_true",
               default=False,
               help="""Whether to train in 16-bit bfloat16 mode. Mutually
                       exclusive with --amp and --fp16. Requires hardware
                       support""")

p = parser.add_argument_group("AMP")
p.add_argument("--amp", action="store_true", default=False,
                help="""Whether to enable AMP training. Mutually exclusive with 
                      --fp16""")
p.add_argument("--opt_level", action="store_true", default=False,
                help="""AMP optimization level. One of "O0", "O1", "O2", or 
                        "O3".""")

p = parser.add_argument_group("Activation checkpointing")
p.add_argument("--partition_activations", action="store_true", 
                default=False,
               help="Activation checkpointing")
p.add_argument("--cpu_checkpointing", action="store_true", default=False,
               help="Offload activation checkpoints to CPU")
p.add_argument("--profile", action="store_true", 
                default=False,
               help="Whether to profile activation checkpointing")


p = parser.add_argument_group("ZeRO optimization")
p.add_argument("--zero_stage", type=int, default=2,
               help="ZeRO optimizer stage")
p.add_argument("--allgather_partitions", action="store_true", 
               default=False,
               help='''Allgather collective vs. broadcast collectives 
                          for parameter gathering''')
p.add_argument("--allgather_bucket_size", type=int, default=1e9,
               help="Number of elements allgathered at one time")
p.add_argument("--overlap_comm", action="store_true", default=False,
               help='''Whether to overlap gradient reduction and backward 
                          pass''')
p.add_argument("--reduce_scatter", action="store_true", default=False,
               help="Use reduce to average gradients")
p.add_argument("--reduce_bucket_size", type=int, default=1e9,
               help="Number of elements reduced at one time")
p.add_argument("--offload_optimizer", action="store_true", default=False,
               help='''Offload optimizer state to CPU. Valid only when 
                          --stage is 2 or 3''')
p.add_argument("--pin_memory", action="store_true", default=False,
                help="Speeds up offloaded throughput at the cost of memory")

p = parser.add_argument_group("Flops profiler")
p.add_argument("--flops_profiler", action="store_true", default=False,
               help="Whether to enable the DeepSpeed Flops Profiler")
p.add_argument("--profile_step", type=int, default=1,
               help='''The global training step at which to run the flops
                       profiler. Has no effect unless --flops_profiler is 
                       given''')
p.add_argument("--module_depth", type=int, default=-1,
               help='''Depth to which aggregated module info is printed. Has 
                       no effect unless --flops_profiler is given''')
p.add_argument("--top_modules", type=int, default=3,
               help='''Number of top modules to print in the aggregated 
                       profile. Has no effect unless --flops_profiler is
                       given''')
p.add_argument("--detailed_flops_profile", action="store_true", 
                default=False,
               help='''Whether the flops_profiler should be detailed. Has
                       no effect unless --flops_profiler is given''')
        
args = parser.parse_args()

d = {}

# Optimizer settings
if(args.optimizer is not None):
    optimizer = {}
    optimizer["type"] = args.optimizer
    params = {}
    params["lr"] = args.lr
    params["eps"] = args.eps
    
    if(args.optimizer == "OneBitAdam"):
        params["freeze_step"] = args.freeze_step
        params["cuda_aware"] = args.cuda_aware
        params["comm_backend_name"] = args.comm_backend_name
    
    optimizer["params"] = params
    d["optimizer"] = optimizer

# LR scheduler
if(args.scheduler is not None):
    scheduler = {}
    scheduler["type"] = args.scheduler
    params = {}
    if(args.scheduler == "LRRangeTest"):
        params["lr_range_test_min_lr"] = args.lr_range_test_min_lr
        params["lr_range_test_step_size"] = args.lr_range_test_step_size
        params["lr_range_test_step_rate"] = args.lr_range_test_step_rate
        params["lr_range_test_staircase"] = args.lr_range_test_staircase
    elif(args.scheduler == "OneCycle"):
        params["cycle_min_lr"] = args.cycle_min_lr
        params["cycle_max_lr"] = args.cycle_max_lr
        params["decay_lr_rate"] = args.cycle_decay_lr_rate
        params["cycle_first_step_size"] = args.cycle_first_step_size
        params["cycle_second_step_size"] = args.cycle_second_step_size
        params["cycle_first_stair_count"] = args.cycle_first_stair_count
        params["cycle_second_stair_count"] = args.cycle_second_stair_count
        params["cycle_momentum"] = args.cycle_momentum
        params["cycle_min_mom"] = args.cycle_min_mom
        params["cycle_max_mom"] = args.cycle_max_mom
        params["decay_mom_rate"] = args.cycle_decay_mom_rate
    elif(args.scheduler == "WarmupLR"):
        params["warmup_min_lr"] = args.warmup_min_lr
        params["warmup_max_lr"] = args.warmup_max_lr
        params["warmup_num_steps"] = args.warmup_num_steps
    elif(args.scheduler == "WarmupDecayLR"):
        params["total_num_steps"] = args.warmup_decay_total_num_steps
        params["warmup_min_lr"] = args.warmup_decay_min_lr
        params["warmup_max_lr"] = args.warmup_decay_max_lr
    else:
        raise ValueError("Invalid scheduler")

    scheduler["params"] = params
    d["scheduler"] = scheduler

# 16-bit training
if(sum([args.amp, args.fp16, args.bfloat16]) > 1):
    raise ValueError("Only one of --fp16, --amp, or --bfloat16 can be enabled")

if(args.amp):
    amp = {}
    amp["enabled"] = True
    amp["pin_memory"] = args.opt_level
    d["amp"] = amp
elif(args.fp16):
    fp16 = {}
    fp16["enabled"] = args.fp16
    d["fp16"] = fp16
elif(args.bfloat16):
    bfloat16 = {}
    bfloat16["enabled"] = args.bfloat16
    d["bfloat16"] = bfloat16

# Activation checkpointing
ac = {}
ac["partition_activations"] = args.partition_activations
ac["cpu_checkpointing"] = args.cpu_checkpointing
ac["profile"] = args.profile
d["activation_checkpointing"] = ac

# ZeRO optimization
zo = {}
zo["stage"] = args.zero_stage
zo["allgather_partitions"] = args.allgather_partitions
zo["allgather_bucket_size"] = args.allgather_bucket_size
zo["reduce_bucket_size"] = args.reduce_bucket_size
zo["overlap_comm"] = args.overlap_comm
zo["reduce_scatter"] = args.reduce_scatter

if(args.offload_optimizer):
    oo = {}
    oo["device"] = "cpu"
    oo["pin_memory"] = args.pin_memory
    zo["offload_optimizer"] = oo

d["zero_optimization"] = zo

# Flops Profiler
flops_profiler = {}
flops_profiler["enabled"] = args.flops_profiler
flops_profiler["profile_step"] = args.profile_step
flops_profiler["module_depth"] = args.module_depth
flops_profiler["top_modules"] = args.top_modules
flops_profiler["detailed"] = args.detailed_flops_profile
d ["flops_profiler"] = flops_profiler

if(args.gradient_clipping):
    d["gradient_clipping"] = args.gradient_clipping

print(json.dumps(d, indent=2))
