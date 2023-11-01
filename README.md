# 「富岳」実装版OpenFold (OpenFold for Fugaku)

## ソースコードと必要ファイルの取得

「富岳」実装版OpenFoldのソースコードをクローンする
```shell
git clone https://github.com/RIKEN-RCCS/OpenFold-for-Fugaku.git
cd OpenFold-for-Fugaku
```

`openfold/resources/`ディレクトリに`stereo_chemical_props.txt`をダウンロードする
```shell
wget -P openfold/resources https://git.scicore.unibas.ch/schwede/openstructure/-/raw/7102c63615b64735c4941278d92b554ec94415f8/modules/mol/alg/src/stereo_chemical_props.txt
```

## データセットの取得
データセットを保存する任意のディレクトリ(以降`DATADIR`と記載)を用意しておく。
aria2をインストールしたのち、学習済みモデルおよび推論や学習に必要なデータセットをダウンロードする。

```shell
git clone https://github.com/spack/spack.git
SPACK_ROOT=`pwd`/spack
source $SPACK_ROOT/share/spack/setup-env.sh
spack install aria2
spack load aria2

DATADIR="データセットを保存するディレクトリ"
bash scripts/download_alphafold_params.sh $DATADIR
bash scripts/download_alphafold_dbs.sh $DATADIR reduced_dbs
```

## 初期設定
`scripts/setenv`の以下の環境変数を設定する

- `PREFIX`: インストールディレクトリ (必須)
- `OPENFOLDDIR`: このREADME.mdがあるディレクトリ (必須) 
- `DATADIR`: データセットをダウンロードしたディレクトリ (必須)
- `PJSUB_OPT`: ジョブ投入時の追加オプション (任意)
    - シングルアカウントの場合は`-g <グループID>`を指定する
    - 必要におうじて`-x PJM_LLIO_GFSCACHE=<ボリューム>`を指定する

## インストール
インストールジョブを投入するスクリプトを実行する
```
./scripts/install_fugaku.sh
```
ジョブの出力は`install_openfold.*_0.out`と`install_openfold.*_1.out`に書き込まれます。ジョブの終了後、それらのファイルの最後に`Finished`が出力されていることを確認してください。


## 推論
### 対象配列が少数である場合
inferenceディレクトリのシングルノード版スクリプトで前処理と推論を連続して行います。詳細は[こちら](inference/README.md)

### 対象配列が多数である場合
1. preproc_fugakuディレクトリのマルチノード版スクリプトで前処理を行います。詳細は[こちら](preproc_fugaku/README.md)
2. inferenceディレクトリのスクリプトで推論を行います。詳細は[こちら](inference/README.md)

## Copyright notice
Copyright 2023 RIKEN & Fujitsu Limited

This software includes AlphaFold and OpenFold that are distributed in the Apache License 2.0.
However, DeepMind's pretrained parameters fall under the CC BY 4.0 license, a copy of which is downloaded to `openfold/resources/params` by the installation script. Note that the latter replaces the original, more restrictive CC BY-NC 4.0 license as of January 2022.

## Changes from the original OpenFold
- Add code to preprocess multiple sequences using multiple nodes in Fugaku
- Add code to infer multiple sequences using multiple nodes in Fugaku
- Port the memory-efficient attention module to CPUs
- Improve processing efficiency of batch matrix multiplications (BMMs) by making the input tensors to BMMs contiguous
- Speed up reading of mmcif data by pickling and lz4 compression beforehand
