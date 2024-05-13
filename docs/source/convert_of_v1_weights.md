## Weights Renaming 

As part of the [OpenFold v2 update](https://github.com/aqlaboratory/openfold/releases/tag/v2.0.0) with the integration of multimer prediction, certain model layers of the AlphaFold model were renamed. For example.

`module.model.template_angle_embedder.*` is now referred to as 
`module.model.template_embedder.template_single_embedder.*`

If you have some checkpoints that were trained using OpenFold v1 or older, and now want to resume training on OpenFold v2, you may need to convert your checkpoints.

## FAQ

### Do I need to convert my checkpoints / model weights?

If you want to run inference or resume training from a checkpoint that was trained with OpenFold V1, you will need to convert your checkpoint.

If you want load model weights only, without starting from a specific time step, then you should not need to convert your checkpoints. The training of the model will begin from `step=0` in this case. To do so, you'll need both the `--resume_from_ckpt` and `--resume_model_weights_only`  flags. This example allows you train starting from the pre-trained openfold weights:

```bash
$ python3 $OPENFOLD_DIR/train_openfold.py test_data_epoch/mmcifs test_data_epoch/alignments test_data_epoch/template_mmcifs $OUTPUT_DIR 2021-09-30 \
	...
	--resume_from_ckpt openfold/resources/openfold_params/finetuning_2.pt \
	--resume_model_weights_only

```

### How do I convert my checkpoints? 

Use [`scripts/convert_v1_to_v2_weights.py`](https://github.com/aqlaboratory/openfold/blob/main/scripts/convert_v1_to_v2_weights.py) e.g.

	`python scripts/convert_v1_to_v2_weights.py checkpoints/6-209.ckpt checkpoints/6-209.ckpt.converted`

Then, to resume training, set the following flags:

```bash
$ python3 $OPENFOLD_DIR/train_openfold.py test_data_epoch/mmcifs test_data_epoch/alignments test_data_epoch/template_mmcifs $OUTPUT_DIR 2021-09-30 \
	...
	--resume_from_ckpt checkpoints/6-209.ckpt.converted
```