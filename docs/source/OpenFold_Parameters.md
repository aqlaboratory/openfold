# Notes on OpenFold Training and Parameters

For OpenFold model parameters, v. 06_22.

## Training details

OpenFold was trained using OpenFold on 44 A100s using the training schedule from Table 4 in
the AlphaFold supplement. AlphaFold was used as the pre-distillation model.
Training data is hosted publicly in the "OpenFold Training Data" RODA repository.

To improve model diversity, we forked training after the initial training phase
and finetuned an additonal branch without templates.

## Parameter files

Parameter files fall into the following categories:

    initial_training.pt:
        OpenFold at the end of the initial training phase.
    finetuning_x.pt:
        Checkpoints in chronological order corresponding to peaks in the
        validation LDDT-Ca during the finetuning phase. Roughly evenly spaced
        across the 45 finetuning epochs.

        NOTE: finetuning_1.pt, which was included in a previous release, has
        been deprecated.
    finetuning_no_templ_x.pt
        Checkpoints in chronological order corresponding to peaks during an
        additional finetuning phase also starting from the 'initial_training.pt'
        checkpoint but with templates disabled.
    finetuning_no_templ_ptm_x.pt
        Checkpoints in chronological order corresponding to peaks during the
        pTM training phase of the `no_templ` branch. Models in this category
        include the pTM module and comprise the most recent of the checkpoints
        in said branch.
    finetuning_ptm_x.pt:
        Checkpoints in chronological order corresponding to peaks in the pTM
        training phase of the mainline branch. Models in this category include
        the pTM module and comprise the most recent of the checkpoints in said
        branch.

Average validation LDDT-Ca scores for each of the checkpoints are listed below.
The validation set contains approximately 180 chains drawn from CAMEO over a
three-month period at the end of 2021.

    initial_training: 0.9088
    finetuning_2: 0.9061
    finetuning_3: 0.9075
    finetuning_4: 0.9059
    finetuning_5: 0.9054
    finetuning_no_templ_1: 0.9014
    finetuning_no_templ_2: 0.9032
    finetuning_no_templ_ptm_1: 0.9025
    finetuning_ptm_1: 0.9075
    finetuning_ptm_2: 0.9097