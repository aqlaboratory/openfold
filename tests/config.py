import ml_collections as mlc


monomer_consts = mlc.ConfigDict(
    {
        "model": "model_1_ptm",  # monomer:model_1_ptm, multimer: model_1_multimer_v3
        "is_multimer": False,  # monomer: False, multimer: True
        "chunk_size": 4,
        "batch_size": 2,
        "n_res": 22,
        "n_seq": 13,
        "n_templ": 3,
        "n_extra": 17,
        "n_heads_extra_msa": 8,
        "eps": 5e-4,
        # For compatibility with DeepMind's pretrained weights, it's easiest for
        # everyone if these take their real values.
        "c_m": 256,
        "c_z": 128,
        "c_s": 384,
        "c_t": 64,
        "c_e": 64,
        "msa_logits": 23,  # monomer: 23, multimer: 22
        "template_mmcif_dir": None  # Set for test_multimer_datamodule
    }
)

multimer_consts = mlc.ConfigDict(
    {
        "model": "model_1_multimer_v3",  # monomer:model_1_ptm, multimer: model_1_multimer_v3
        "is_multimer": True,  # monomer: False, multimer: True
        "chunk_size": 4,
        "batch_size": 2,
        "n_res": 22,
        "n_seq": 13,
        "n_templ": 3,
        "n_extra": 17,
        "n_heads_extra_msa": 8,
        "eps": 5e-4,
        # For compatibility with DeepMind's pretrained weights, it's easiest for
        # everyone if these take their real values.
        "c_m": 256,
        "c_z": 128,
        "c_s": 384,
        "c_t": 64,
        "c_e": 64,
        "msa_logits": 22,  # monomer: 23, multimer: 22
        "template_mmcif_dir": None  # Set for test_multimer_datamodule
    }
)

consts = monomer_consts 

config = mlc.ConfigDict(
    {
        "data": {
            "common": {
                "masked_msa": {
                    "profile_prob": 0.1,
                    "same_prob": 0.1,
                    "uniform_prob": 0.1,
                },
            }
        }
    }
)
