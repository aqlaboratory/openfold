import ml_collections as mlc

consts = mlc.ConfigDict(
    {
        "batch_size": 2,
        "n_res": 11,
        "n_seq": 13,
        "n_templ": 3,
        "n_extra": 17,
        "eps": 5e-4,
        # For compatibility with DeepMind's pretrained weights, it's easiest for
        # everyone if these take their real values.
        "c_m": 256,
        "c_z": 128,
        "c_s": 384,
        "c_t": 64,
        "c_e": 64,
    }
)
