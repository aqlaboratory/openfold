#!/usr/bin/env python

# This script extracts fp32 consolidated weights from a zero 2 and 3 DeepSpeed checkpoints. It gets
# copied into the top level checkpoint dir, so the user can easily do the conversion at any point in
# the future. Once extracted, the weights don't require DeepSpeed and can be used in any
# application.
#
# example: python zero_to_fp32.py . pytorch_model.bin

import argparse
import torch
import glob
import math
import os
from collections import OrderedDict
import re

# while this script doesn't use deepspeed to recover data, since the checkpoints are pickled with
# DeepSpeed data structures it has to be available in the current python environment.
import deepspeed
from deepspeed.utils import logger

debug = 0

# load to cpu
device = torch.device('cpu')


def get_model_state_file(checkpoint_dir, zero_stage):
    if not os.path.isdir(checkpoint_dir):
        raise FileNotFoundError(f"Directory '{checkpoint_dir}' doesn't exist")

    # there should be only one file
    if zero_stage == 2:
        file = os.path.join(checkpoint_dir, "mp_rank_00_model_states.pt")
    elif zero_stage == 3:
        file = os.path.join(checkpoint_dir, "zero_pp_rank_0_mp_rank_00_model_states.pt")

    if not os.path.exists(file):
        raise FileNotFoundError(f"can't find model states file at '{file}'")

    return file


def get_optim_files(checkpoint_dir):
    # XXX: need to test that this simple glob rule works for multi-node setup too
    optim_files = sorted(glob.glob(os.path.join(checkpoint_dir, "*_optim_states.pt")))

    if len(optim_files) == 0:
        raise FileNotFoundError(
            f"can't find '*_optim_states.pt' files in directory '{checkpoint_dir}'")

    return optim_files


def parse_model_state(file):
    state_dict = torch.load(file, map_location=device)

    if "buffer_names" not in state_dict:
        raise ValueError(f"{file} is not a model state checkpoint")
    buffer_names = state_dict["buffer_names"]
    if debug:
        print("Found buffers:", buffer_names)

    # recover just the buffers while restoring them to fp32 if they were saved in fp16
    buffers = {
        k: v.float()
        for k,
        v in state_dict["module"].items() if k in buffer_names
    }
    return buffers


def parse_optim_states(files, ds_checkpoint_dir):

    total_files = len(files)
    state_dicts = []
    for f in files:
        state_dicts.append(torch.load(f, map_location=device))

    if not "zero_stage" in state_dicts[0]['optimizer_state_dict']:
        raise ValueError(f"{files[0]} is not a zero checkpoint")
    zero_stage = state_dicts[0]['optimizer_state_dict']["zero_stage"]
    world_size = state_dicts[0]['optimizer_state_dict']["partition_count"]
    param_shapes = state_dicts[0]["param_shapes"]
    # For ZeRO-2 each param group can have different partition_count as data parallelism for expert
    # parameters can be different from data parallelism for non-expert parameters. So we can just
    # use the max of the partition_count to get the dp world_size.

    if type(world_size) is list:
        world_size = max(world_size)

    if world_size != total_files:
        raise ValueError(
            f"Expected {world_size} of '*_optim_states.pt' under '{ds_checkpoint_dir}' but found {total_files} files. "
            "Possibly due to an overwrite of an old checkpoint, or a checkpoint didn't get saved by one or more processes."
        )

    # the groups are named differently in each stage
    if zero_stage == 2:
        fp32_groups_key = "single_partition_of_fp32_groups"
    elif zero_stage == 3:
        fp32_groups_key = "fp32_flat_groups"
    else:
        raise ValueError(f"unknown zero stage {zero_stage}")

    if zero_stage == 2:
        fp32_flat_groups = [
            state_dicts[i]['optimizer_state_dict'][fp32_groups_key]
            for i in range(len(state_dicts))
        ]
    elif zero_stage == 3:
        # if there is more than one param group, there will be multiple flattened tensors - one
        # flattened tensor per group - for simplicity merge them into a single tensor
        #
        # XXX: could make the script more memory efficient for when there are multiple groups - it
        # will require matching the sub-lists of param_shapes for each param group flattened tensor

        fp32_flat_groups = [
            torch.cat(state_dicts[i]['optimizer_state_dict'][fp32_groups_key],
                      0) for i in range(len(state_dicts))
        ]

    return zero_stage, world_size, param_shapes, fp32_flat_groups


def _get_fp32_state_dict_from_zero_checkpoint(ds_checkpoint_dir):
    """
    Returns fp32 state_dict reconstructed from ds checkpoint

    Args:
        - ``ds_checkpoint_dir``: path to the deepspeed checkpoint folder (where the optimizer files are)

    """
    print(f"Processing zero checkpoint '{ds_checkpoint_dir}'")

    optim_files = get_optim_files(ds_checkpoint_dir)
    zero_stage, world_size, param_shapes, fp32_flat_groups = parse_optim_states(optim_files, ds_checkpoint_dir)
    print(
        f"Detected checkpoint of type zero stage {zero_stage}, world_size: {world_size}")

    model_file = get_model_state_file(ds_checkpoint_dir, zero_stage)
    buffers = parse_model_state(model_file)

    if zero_stage == 2:
        return _get_fp32_state_dict_from_zero2_checkpoint(world_size,
                                                          param_shapes,
                                                          fp32_flat_groups,
                                                          buffers)
    elif zero_stage == 3:
        return _get_fp32_state_dict_from_zero3_checkpoint(world_size,
                                                          param_shapes,
                                                          fp32_flat_groups,
                                                          buffers)


def _get_fp32_state_dict_from_zero2_checkpoint(world_size,
                                               param_shapes,
                                               fp32_flat_groups,
                                               buffers):

    # Reconstruction protocol:
    #
    # XXX: document this

    if debug:
        for i in range(world_size):
            for j in range(len(fp32_flat_groups[0])):
                print(f"fp32_flat_groups[{i}][{j}].shape={fp32_flat_groups[i][j].shape}")

    # XXX: memory usage doubles here (zero2)
    num_param_groups = len(fp32_flat_groups[0])
    merged_single_partition_of_fp32_groups = []
    for i in range(num_param_groups):
        merged_partitions = [sd[i] for sd in fp32_flat_groups]
        full_single_fp32_vector = torch.cat(merged_partitions, 0)
        merged_single_partition_of_fp32_groups.append(full_single_fp32_vector)
    avail_numel = sum([
        full_single_fp32_vector.numel()
        for full_single_fp32_vector in merged_single_partition_of_fp32_groups
    ])

    if debug:
        wanted_params = sum([len(shapes) for shapes in param_shapes])
        wanted_numel = sum(
            [sum(shape.numel() for shape in shapes.values()) for shapes in param_shapes])
        # not asserting if there is a mismatch due to possible padding
        print(f"Have {avail_numel} numels to process.")
        print(f"Need {wanted_numel} numels in {wanted_params} params.")

    state_dict = OrderedDict()

    # buffers
    state_dict.update(buffers)
    if debug:
        print(f"added {len(buffers)} buffers")

    # params
    # XXX: for huge models that can't fit into the host's RAM we will have to recode this to support
    # out-of-core computing solution
    total_numel = 0
    total_params = 0
    for shapes, full_single_fp32_vector in zip(param_shapes, merged_single_partition_of_fp32_groups):
        offset = 0
        avail_numel = full_single_fp32_vector.numel()
        for name, shape in shapes.items():

            unpartitioned_numel = shape.numel()
            total_numel += unpartitioned_numel
            total_params += 1

            if debug:
                print(
                    f"{name} full shape: {shape} unpartitioned numel {unpartitioned_numel} "
                )
            state_dict[name] = full_single_fp32_vector.narrow(
                0,
                offset,
                unpartitioned_numel).view(shape)
            offset += unpartitioned_numel

        # Z2 started to align to 2*world_size to improve nccl performance. Therefore both offset and
        # avail_numel can differ by anywhere between 0..2*world_size. Due to two unrelated complex
        # paddings performed in the code it's almost impossible to predict the exact numbers w/o the
        # live optimizer object, so we are checking that the numbers are within the right range
        align_to = 2 * world_size

        def zero2_align(x):
            return align_to * math.ceil(x / align_to)

        if debug:
            print(f"original offset={offset}, avail_numel={avail_numel}")

        offset = zero2_align(offset)
        avail_numel = zero2_align(avail_numel)

        if debug:
            print(f"aligned  offset={offset}, avail_numel={avail_numel}")

        # Sanity check
        if offset != avail_numel:
            raise ValueError(
                f"consumed {offset} numels out of {avail_numel} - something is wrong")

    print(
        f"Reconstructed fp32 state dict with {total_params} params {total_numel} elements"
    )

    return state_dict


def zero3_partitioned_param_info(unpartitioned_numel, world_size):
    remainder = unpartitioned_numel % world_size
    padding_numel = (world_size - remainder) if remainder else 0
    partitioned_numel = math.ceil(unpartitioned_numel / world_size)
    return partitioned_numel, padding_numel


def _get_fp32_state_dict_from_zero3_checkpoint(world_size,
                                               param_shapes,
                                               fp32_flat_groups,
                                               buffers):

    # Reconstruction protocol: For zero3 we need to zip the partitions together at boundary of each
    # param, re-consolidating each param, while dealing with padding if any

    avail_numel = fp32_flat_groups[0].numel() * world_size
    # merge list of dicts, preserving order
    param_shapes = {k: v for d in param_shapes for k, v in d.items()}

    if debug:
        for i in range(world_size):
            print(f"fp32_flat_groups[{i}].shape={fp32_flat_groups[i].shape}")

        wanted_params = len(param_shapes)
        wanted_numel = sum(shape.numel() for shape in param_shapes.values())
        # not asserting if there is a mismatch due to possible padding
        print(f"Have {avail_numel} numels to process.")
        print(f"Need {wanted_numel} numels in {wanted_params} params.")

    state_dict = OrderedDict()

    # buffers
    state_dict.update(buffers)
    if debug:
        print(f"added {len(buffers)} buffers")

    # params
    # XXX: for huge models that can't fit into the host's RAM we will have to recode this to support
    # out-of-core computing solution
    offset = 0
    total_numel = 0
    total_params = 0
    for name, shape in param_shapes.items():

        unpartitioned_numel = shape.numel()
        total_numel += unpartitioned_numel
        total_params += 1

        partitioned_numel, partitioned_padding_numel = zero3_partitioned_param_info(unpartitioned_numel, world_size)

        if debug:
            print(
                f"{total_params} {name} full shape: {shape} partition0 numel={partitioned_numel} partitioned_padding_numel={partitioned_padding_numel}"
            )

        # XXX: memory usage doubles here
        state_dict[name] = torch.cat(
            tuple(fp32_flat_groups[i].narrow(0,
                                             offset,
                                             partitioned_numel)
                  for i in range(world_size)),
            0).narrow(0,
                      0,
                      unpartitioned_numel).view(shape)
        offset += partitioned_numel

    offset *= world_size

    # Sanity check
    if offset != avail_numel:
        raise ValueError(
            f"consumed {offset} numels out of {avail_numel} - something is wrong")

    print(
        f"Reconstructed fp32 state dict with {total_params} params {total_numel} elements"
    )

    return state_dict


def get_fp32_state_dict_from_zero_checkpoint(checkpoint_dir, tag=None):
    """
    Convert ZeRO 2 or 3 checkpoint into a single fp32 consolidated state_dict that can be loaded with
    ``load_state_dict()`` and used for training without DeepSpeed or shared with others, for example
    via a model hub.

    Args:
        - ``checkpoint_dir``: path to the desired checkpoint folder
        - ``tag``: checkpoint tag used as a unique identifier for checkpoint. If not provided will attempt to load tag in 'latest' file. e.g., ``global_step14``

    Returns:
        - pytorch ``state_dict``

    Note: this approach may not work if your application doesn't have sufficient free CPU memory and
    you may need to use the offline approach using the ``zero_to_fp32.py`` script that is saved with
    the checkpoint.

    A typical usage might be ::

        from deepspeed.utils.zero_to_fp32 import get_fp32_state_dict_from_zero_checkpoint
        # do the training and checkpoint saving
        state_dict = get_fp32_state_dict_from_zero_checkpoint(checkpoint_dir) # already on cpu
        model = model.cpu() # move to cpu
        model.load_state_dict(state_dict)
        # submit to model hub or save the model to share with others

    In this example the ``model`` will no longer be usable in the deepspeed context of the same
    application. i.e. you will need to re-initialize the deepspeed engine, since
    ``model.load_state_dict(state_dict)`` will remove all the deepspeed magic from it.

    If you want it all done for you, use ``load_state_dict_from_zero_checkpoint`` instead.

    """
    if tag is None:
        latest_path = os.path.join(checkpoint_dir, 'latest')
        if os.path.isfile(latest_path):
            with open(latest_path, 'r') as fd:
                tag = fd.read().strip()
        else:
            raise ValueError(f"Unable to find 'latest' file at {latest_path}")

    ds_checkpoint_dir = os.path.join(checkpoint_dir, tag)

    if not os.path.isdir(ds_checkpoint_dir):
        raise FileNotFoundError(f"Directory '{ds_checkpoint_dir}' doesn't exist")

    return _get_fp32_state_dict_from_zero_checkpoint(ds_checkpoint_dir)


def convert_zero_checkpoint_to_fp32_state_dict(checkpoint_dir, output_file, tag=None):
    """
    Convert ZeRO 2 or 3 checkpoint into a single fp32 consolidated ``state_dict`` file that can be
    loaded with ``torch.load(file)`` + ``load_state_dict()`` and used for training without DeepSpeed.

    Args:
        - ``checkpoint_dir``: path to the desired checkpoint folder. (one that contains the tag-folder, like ``global_step14``)
        - ``output_file``: path to the pytorch fp32 state_dict output file (e.g. path/pytorch_model.bin)
        - ``tag``: checkpoint tag used as a unique identifier for checkpoint. If not provided will attempt to load tag in the file named ``latest`` in the checkpoint folder, e.g., ``global_step14``
    """

    state_dict = get_fp32_state_dict_from_zero_checkpoint(checkpoint_dir, tag)
    print(f"Saving fp32 state dict to {output_file}")
    torch.save(state_dict, output_file)


def load_state_dict_from_zero_checkpoint(model, checkpoint_dir, tag=None):
    """
    1. Put the provided model to cpu
    2. Convert ZeRO 2 or 3 checkpoint into a single fp32 consolidated ``state_dict``
    3. Load it into the provided model

    Args:
        - ``model``: the model object to update
        - ``checkpoint_dir``: path to the desired checkpoint folder. (one that contains the tag-folder, like ``global_step14``)
        - ``tag``: checkpoint tag used as a unique identifier for checkpoint. If not provided will attempt to load tag in the file named ``latest`` in the checkpoint folder, e.g., ``global_step14``

    Returns:
        - ``model`: modified model

    Make sure you have plenty of CPU memory available before you call this function. If you don't
    have enough use the ``zero_to_fp32.py`` utility to do the conversion. You will find it
    conveniently placed for you in the checkpoint folder.

    A typical usage might be ::

        from deepspeed.utils.zero_to_fp32 import load_state_dict_from_zero_checkpoint
        model = load_state_dict_from_zero_checkpoint(trainer.model, checkpoint_dir)
        # submit to model hub or save the model to share with others

    Note, that once this was run, the ``model`` will no longer be usable in the deepspeed context
    of the same application. i.e. you will need to re-initialize the deepspeed engine, since
    ``model.load_state_dict(state_dict)`` will remove all the deepspeed magic from it.

    """
    logger.info(f"Extracting fp32 weights")
    state_dict = get_fp32_state_dict_from_zero_checkpoint(checkpoint_dir, tag)

    logger.info(f"Overwriting model with fp32 weights")
    model = model.cpu()
    model.load_state_dict(state_dict, strict=False)

    return model

def get_global_step_from_zero_checkpoint(checkpoint_dir):
    global_step = -1
    latest_path = os.path.join(checkpoint_dir, 'latest')
    if os.path.isfile(latest_path):
        with open(latest_path, 'r') as fd:
            tag = fd.read().strip()
            match = re.match(r"global_step([0-9]+)", tag)
            global_step = int(match.group(1))
    else:
            raise ValueError(f"Unable to find 'latest' file at {latest_path}")
    return global_step

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "checkpoint_dir",
        type=str,
        help="path to the desired checkpoint folder, e.g., path/checkpoint-12")
    parser.add_argument(
        "output_file",
        type=str,
        help=
        "path to the pytorch fp32 state_dict output file (e.g. path/checkpoint-12/pytorch_model.bin)"
    )
    parser.add_argument("-d", "--debug", action='store_true', help="enable debug")
    args = parser.parse_args()

    debug = args.debug

    convert_zero_checkpoint_to_fp32_state_dict(args.checkpoint_dir, args.output_file)
