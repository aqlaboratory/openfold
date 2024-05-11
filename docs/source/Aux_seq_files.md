# Auxiliary Sequence Files for OpenFold Training

The training dataset of OpenFold is very large. The `pdb` directory alone contains 185,000 mmcifs;  each chain for has multiple sequence alignment files and mmcif files. 

OpenFold introduces a few new file structures for faster access to alignments and mmcif data.

This documentation will explain the benefits of having the condensed file structure, and explain the contents of each of the files.

## Default alignment file structure 

One way to store mmcifs and alignments files would be to have a directory for each mmcif chain. 

For example, consider two protein as a case study
```
- OpenProteinSet
  └── mmcifs 
	 ├── 3lrm.cif
	 └── 6kwc.cif
	 ...
```

In the `alignments` directory, [PDB:6KWC](https://www.rcsb.org/structure/6KWC) is a monomer with one chain, and thus would have one alignment direcotry. [PDB:3LRM](https://www.rcsb.org/structure/3lrm), a homotetramer, would have one alignment directory for each of its four chains.
```
- OpenProteinSet 
  └── alignments 
	  └── 3lrm_A 
		  ├── bfd_uniclust_hits.a3m
		  ├── mgnify_hits.a3m
		  ├── pdb70_hits.hhr
		  └── uniref90_hits.a3m 
	  └── 3lrm_B 
		  ├── bfd_uniclust_hits.a3m
		  ├── mgnify_hits.a3m
		  ├── pdb70_hits.hhr
		  └── uniref90_hits.a3m 
	  └── 3lrm_C 
		  ├── bfd_uniclust_hits.a3m
		  ├── mgnify_hits.a3m
		  ├── pdb70_hits.hhr
		  └── uniref90_hits.a3m 
	  └── 3lrm_D 
		  ├── bfd_uniclust_hits.a3m
		  ├── mgnify_hits.a3m
		  ├── pdb70_hits.hhr
		  └── uniref90_hits.a3m 
	  └── 6kwc_A 
		  ├── bfd_uniclust_hits.a3m
		  ├── mgnify_hits.a3m
		  ├── pdb70_hits.hhr
		  └── uniref90_hits.a3m 
		...
```

In practice, the IO overhead of having one directory per protein chain makes accessing the alignments slow. 

## OpenFold DB file structure 

 Here we describe a new filesystem that can be used by OpenFold for more efficient access of alignment file and index file contents

All together, the file directory would look like:
```
- OpenProteinSet 
  ├── duplicate_pdb_chains.txt
  └── pdb
	  ├── mmcif_cache.json 
	  └── mmcifs 
		  ├── 3lrm.cif
		  └── 6kwc.cif
	  └── alignment_db
		  ├── alignment_db_0.db 
          ├── alignment_db_1.db
          ...
          ├── alignment_db_9.db
		  └── alignment_db.index 
```

We will describe each of the file types here.

### Alignments db files and index files

To speed up access of MSAs, OpenFold has an alternate alignments storage procedure. Instead of storing dedicated files for each single alignment, we consolidate large sets of alignments to single files referred to as _alignments_db's_. This can reduce I/O overhead and in practice we recommend using around 10 `alignments_db_x.db` files to store the total training set of OpenFold. During training, OpenFold can access each alignment using byte index pointers that are stored in a separate index file (`alignments_db.index`). The alignments for the `3LRM` and `6KWC` examples would be recorded in the index file as follows:

```alignments_db.index
{
	...
	"3lrm_A": {
        "db": "alignment_db_0.db",
        "files": [
            ["bfd_uniclust_hits.a3m", 212896478938, 1680200],
            ["mgnify_hits.a3m", 212893696883, 2782055],
            ["pdb70_hits.hhr", 212898159138, 614978],
            ["uniref90_hits.a3m", 212898774116, 6165789]
        ]
    },
    "6kwc_A": {
        "db": "alignment_db_1.db",
        "files": [
            ["bfd_uniclust_hits.a3m", 415618723280, 380289],
            ["mgnify_hits.a3m", 415618556077, 167203],
            ["pdb70_hits.hhr", 415619103569, 148672],
            ["uniref90_hits.a3m", 415617547852, 1008225]
        ]
    }
	...
}
```

For each entry, the corresponding `alignment_db` file and the byte start location and number of bytes to read the respective alignments are given. For example, the alignment information in `bfd_uniclust_hits.a3m` for chain `3lrm_A` can be found in the database file `alignment_db_0.db`, starting at byte location `212896478938` and reading in the next `1680200` bytes.

### Chain cache files and mmCIF cache files

Information from the mmcif files can be parsed in advance to create a `chain_cache.json` or a `mmcif_cache.json`. For OpenFold, the `chain_cache.json` is used to sample chains for training, and the `mmcif_cache.json` is used to prefilter templates. 

Here's what the chain_cache.json entry looks like for our examples:

```chain_cache.json
{
	...
	    "3lrm_A": {
        "release_date": "2010-06-30",
        "seq": "MFAFYFLTACISLKGVFGVSPSYNGLGLTPQMGWDNWNTFACDVSEQLLLDTADRISDLGLKDMGYKYIILDDCWSSGRDSDGFLVADEQKFPNGMGHVADHLHNNSFLFGMYSSAGEYTCAGYPGSLGREEEDAQFFANNRVDYLKYANCYNKGQFGTPEISYHRYKAMSDALNKTGRPVFYSLCNWGQDLTFYWGSGIANSWRMSGDVTAEFTRPDSRCPCDGDEYDCKYAGFHCSIMNILNKAAPMGQNAGVGGWNDLDNLEVGVGNLTDDEEKAHFSMWAMVKSPLIIGANVNNLKASSYSIYSQASVIAINQDSNGIPATRVWRYYVSDTDEYGQGEIQMWSGPLDNGDQVVALLNGGSVSRPMNTTLEEIFFDSNLGSKKLTSTWDIYDLWANRVDNSTASAILGRNKTATGILYNATEQSYKDGLSKNDTRLFGQKIGSLSPNAILNTTVPAHGIAFYRLRPSSDYKDDDDK",
        "resolution": 2.7,
        "cluster_size": 6
    },
    "3lrm_B": {
        "release_date": "2010-06-30",
        "seq": "MFAFYFLTACISLKGVFGVSPSYNGLGLTPQMGWDNWNTFACDVSEQLLLDTADRISDLGLKDMGYKYIILDDCWSSGRDSDGFLVADEQKFPNGMGHVADHLHNNSFLFGMYSSAGEYTCAGYPGSLGREEEDAQFFANNRVDYLKYANCYNKGQFGTPEISYHRYKAMSDALNKTGRPVFYSLCNWGQDLTFYWGSGIANSWRMSGDVTAEFTRPDSRCPCDGDEYDCKYAGFHCSIMNILNKAAPMGQNAGVGGWNDLDNLEVGVGNLTDDEEKAHFSMWAMVKSPLIIGANVNNLKASSYSIYSQASVIAINQDSNGIPATRVWRYYVSDTDEYGQGEIQMWSGPLDNGDQVVALLNGGSVSRPMNTTLEEIFFDSNLGSKKLTSTWDIYDLWANRVDNSTASAILGRNKTATGILYNATEQSYKDGLSKNDTRLFGQKIGSLSPNAILNTTVPAHGIAFYRLRPSSDYKDDDDK",
        "resolution": 2.7,
        "cluster_size": 6
    },
    "3lrm_C": {
        "release_date": "2010-06-30",
        "seq": "MFAFYFLTACISLKGVFGVSPSYNGLGLTPQMGWDNWNTFACDVSEQLLLDTADRISDLGLKDMGYKYIILDDCWSSGRDSDGFLVADEQKFPNGMGHVADHLHNNSFLFGMYSSAGEYTCAGYPGSLGREEEDAQFFANNRVDYLKYANCYNKGQFGTPEISYHRYKAMSDALNKTGRPVFYSLCNWGQDLTFYWGSGIANSWRMSGDVTAEFTRPDSRCPCDGDEYDCKYAGFHCSIMNILNKAAPMGQNAGVGGWNDLDNLEVGVGNLTDDEEKAHFSMWAMVKSPLIIGANVNNLKASSYSIYSQASVIAINQDSNGIPATRVWRYYVSDTDEYGQGEIQMWSGPLDNGDQVVALLNGGSVSRPMNTTLEEIFFDSNLGSKKLTSTWDIYDLWANRVDNSTASAILGRNKTATGILYNATEQSYKDGLSKNDTRLFGQKIGSLSPNAILNTTVPAHGIAFYRLRPSSDYKDDDDK",
        "resolution": 2.7,
        "cluster_size": 6
    },
    "3lrm_D": {
        "release_date": "2010-06-30",
        "seq": "MFAFYFLTACISLKGVFGVSPSYNGLGLTPQMGWDNWNTFACDVSEQLLLDTADRISDLGLKDMGYKYIILDDCWSSGRDSDGFLVADEQKFPNGMGHVADHLHNNSFLFGMYSSAGEYTCAGYPGSLGREEEDAQFFANNRVDYLKYANCYNKGQFGTPEISYHRYKAMSDALNKTGRPVFYSLCNWGQDLTFYWGSGIANSWRMSGDVTAEFTRPDSRCPCDGDEYDCKYAGFHCSIMNILNKAAPMGQNAGVGGWNDLDNLEVGVGNLTDDEEKAHFSMWAMVKSPLIIGANVNNLKASSYSIYSQASVIAINQDSNGIPATRVWRYYVSDTDEYGQGEIQMWSGPLDNGDQVVALLNGGSVSRPMNTTLEEIFFDSNLGSKKLTSTWDIYDLWANRVDNSTASAILGRNKTATGILYNATEQSYKDGLSKNDTRLFGQKIGSLSPNAILNTTVPAHGIAFYRLRPSSDYKDDDDK",
        "resolution": 2.7,
        "cluster_size": 6
    },
	"6kwc_A": {
        "release_date": "2021-01-27",
        "seq": "GSTIQPGTGYNNGYFYSYWNDGHGGVTYTNGPGGQFSVNWSNSGEFVGGKGWQPGTKNKVINFSGSYNPNGNSYLSVYGWSRNPLIEYYIVENFGTYNPSTGATKLGEVTSDGSVYDIYRTQRVNQPSIIGTATFYQYWSVRRNHRSSGSVNTANHFNAWAQQGLTLGTMDYQIVAVQGYFSSGSASITVS",
        "resolution": 1.297,
        "cluster_size": 195
    },
	...
}
```

The mmcif_cache.json file would contain similar information, but condensed by mmcif id, e.g.

```mmcif_cache.json
{
    "3lrm": {
        "release_date": "2010-06-30",
        "chain_ids": [
            "A",
            "B",
            "C",
            "D"
        ],
        "seqs": [
            "MFAFYFLTACISLKGVFGVSPSYNGLGLTPQMGWDNWNTFACDVSEQLLLDTADRISDLGLKDMGYKYIILDDCWSSGRDSDGFLVADEQKFPNGMGHVADHLHNNSFLFGMYSSAGEYTCAGYPGSLGREEEDAQFFANNRVDYLKYANCYNKGQFGTPEISYHRYKAMSDALNKTGRPVFYSLCNWGQDLTFYWGSGIANSWRMSGDVTAEFTRPDSRCPCDGDEYDCKYAGFHCSIMNILNKAAPMGQNAGVGGWNDLDNLEVGVGNLTDDEEKAHFSMWAMVKSPLIIGANVNNLKASSYSIYSQASVIAINQDSNGIPATRVWRYYVSDTDEYGQGEIQMWSGPLDNGDQVVALLNGGSVSRPMNTTLEEIFFDSNLGSKKLTSTWDIYDLWANRVDNSTASAILGRNKTATGILYNATEQSYKDGLSKNDTRLFGQKIGSLSPNAILNTTVPAHGIAFYRLRPSSDYKDDDDK",
            "MFAFYFLTACISLKGVFGVSPSYNGLGLTPQMGWDNWNTFACDVSEQLLLDTADRISDLGLKDMGYKYIILDDCWSSGRDSDGFLVADEQKFPNGMGHVADHLHNNSFLFGMYSSAGEYTCAGYPGSLGREEEDAQFFANNRVDYLKYANCYNKGQFGTPEISYHRYKAMSDALNKTGRPVFYSLCNWGQDLTFYWGSGIANSWRMSGDVTAEFTRPDSRCPCDGDEYDCKYAGFHCSIMNILNKAAPMGQNAGVGGWNDLDNLEVGVGNLTDDEEKAHFSMWAMVKSPLIIGANVNNLKASSYSIYSQASVIAINQDSNGIPATRVWRYYVSDTDEYGQGEIQMWSGPLDNGDQVVALLNGGSVSRPMNTTLEEIFFDSNLGSKKLTSTWDIYDLWANRVDNSTASAILGRNKTATGILYNATEQSYKDGLSKNDTRLFGQKIGSLSPNAILNTTVPAHGIAFYRLRPSSDYKDDDDK",
            "MFAFYFLTACISLKGVFGVSPSYNGLGLTPQMGWDNWNTFACDVSEQLLLDTADRISDLGLKDMGYKYIILDDCWSSGRDSDGFLVADEQKFPNGMGHVADHLHNNSFLFGMYSSAGEYTCAGYPGSLGREEEDAQFFANNRVDYLKYANCYNKGQFGTPEISYHRYKAMSDALNKTGRPVFYSLCNWGQDLTFYWGSGIANSWRMSGDVTAEFTRPDSRCPCDGDEYDCKYAGFHCSIMNILNKAAPMGQNAGVGGWNDLDNLEVGVGNLTDDEEKAHFSMWAMVKSPLIIGANVNNLKASSYSIYSQASVIAINQDSNGIPATRVWRYYVSDTDEYGQGEIQMWSGPLDNGDQVVALLNGGSVSRPMNTTLEEIFFDSNLGSKKLTSTWDIYDLWANRVDNSTASAILGRNKTATGILYNATEQSYKDGLSKNDTRLFGQKIGSLSPNAILNTTVPAHGIAFYRLRPSSDYKDDDDK",
            "MFAFYFLTACISLKGVFGVSPSYNGLGLTPQMGWDNWNTFACDVSEQLLLDTADRISDLGLKDMGYKYIILDDCWSSGRDSDGFLVADEQKFPNGMGHVADHLHNNSFLFGMYSSAGEYTCAGYPGSLGREEEDAQFFANNRVDYLKYANCYNKGQFGTPEISYHRYKAMSDALNKTGRPVFYSLCNWGQDLTFYWGSGIANSWRMSGDVTAEFTRPDSRCPCDGDEYDCKYAGFHCSIMNILNKAAPMGQNAGVGGWNDLDNLEVGVGNLTDDEEKAHFSMWAMVKSPLIIGANVNNLKASSYSIYSQASVIAINQDSNGIPATRVWRYYVSDTDEYGQGEIQMWSGPLDNGDQVVALLNGGSVSRPMNTTLEEIFFDSNLGSKKLTSTWDIYDLWANRVDNSTASAILGRNKTATGILYNATEQSYKDGLSKNDTRLFGQKIGSLSPNAILNTTVPAHGIAFYRLRPSSDYKDDDDK"
        ],
        "no_chains": 4,
        "resolution": 2.7
    },
	    "6kwc": {
        "release_date": "2021-01-27",
        "chain_ids": [
            "A"
        ],
        "seqs": [
            "GSTIQPGTGYNNGYFYSYWNDGHGGVTYTNGPGGQFSVNWSNSGEFVGGKGWQPGTKNKVINFSGSYNPNGNSYLSVYGWSRNPLIEYYIVENFGTYNPSTGATKLGEVTSDGSVYDIYRTQRVNQPSIIGTATFYQYWSVRRNHRSSGSVNTANHFNAWAQQGLTLGTMDYQIVAVQGYFSSGSASITVS"
        ],
        "no_chains": 1,
        "resolution": 1.297
    },
    ...
}
```


### Duplicate pdb chain files 

Duplicate chains occur across pdb entries. Some of these chains are the homomeric units of a multimer, others are subunits that are shared across different protein.

To reduce storage overhead of creating / storing identical data for duplicate entries, we have a duplicate chain file. Each line stores the all chains that are identical. Our `6kwc` and `3lrm` examples would be stored as follows.

```duplicate_pdb_chains.txt
...
6kwc_A
3lrm_A 3lrm_B 3lrm_C 3lrm_D
...
```


