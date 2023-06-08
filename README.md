# 「富岳」実装版OpenFold (Fugaku implementation for OpenFold)

## ソースコードと必要ファイルの取得

「富岳」実装版OpenFoldのソースコードをクローン
```shell
git clone <TBD>
cd openfold
```

`openfold/resources/`ディレクトリに`stereo_chemical_props.txt`をダウンロード
```shell
wget -P openfold/resources https://git.scicore.unibas.ch/schwede/openstructure/-/raw/7102c63615b64735c4941278d92b554ec94415f8/modules/mol/alg/src/stereo_chemical_props.txt
```

## データセットの取得
学習済みモデル、推論や学習に必要なデータセットをダウンロード
```
bash scripts/download_data.sh <データディレクトリ>
bash scripts/download_alphafold_dbs.sh <データディレクトリ>
```


## 初期設定
`scripts/setenv`の以下の環境変数を設定する

- `PREFIX`: インストールディレクトリ (必須)
- `OPENFOLDDIR`: このREADME.mdがあるディレクトリ (必須) 
- `DATADIR`: データディレクトリ (必須)
- `PJSUB_OPT`: pjsub時の追加オプション (任意)
    - シングルアカウントの場合は`-g <グループID>`を指定する
    - 必要におうじて`-x PJM_LLIO_GFSCACHE=<ボリューム>`を指定する

## インストール (PREFIXにビルド済み実行環境を指定した場合は不要)
インストール用バッチジョブを投入するスクリプトを実行
```
./scripts/install_fugaku.sh
```
ジョブの出力は`install_openfold.*_0.out`と`install_openfold.*_1.out`に書き込まれる。それらのファイルの最後に`Finished`が出力されていることを確認する。


## 前処理と推論の実行
### 少数の配列の場合
inferenceディレクトリのシングルノード版スクリプトで前処理と推論を連続して行う。詳細は[こちら](inference/README.md)

### 多数の配列の場合
1. preproc_fugakuディレクトリのマルチノード版スクリプトで前処理を行う。詳細は[こちら](preproc_fugaku/README.md)
2. inferenceディレクトリのスクリプトで推論を行う。詳細は[こちら](inference/README.md)

## 学習の実行
1. trainディレクトリのスクリプトで学習を行う。詳細は[こちら](train/README.md)

## [任意] mmcifの変換
mmcifを事前にパースしpickle化およびlz4圧縮しておくことで、特に複数ノード学習での速度低下を防ぐことができる。
1. convert_mmcifディレクトリに移動
2. run.shのInputDirとOutputDirを設定
3. ジョブ投入
```
pjsub -g {グループID} run.sh
```
4. 学習および推論ディレクトリのparametersに記載されているmmcifのディレクトリを変換したものに置き換える


## Copyright notice

This software includes AlphaFold and OpenFold that are distributed in the Apache License 2.0.
However, DeepMind's pretrained parameters fall under the CC BY 4.0 license, a copy of which is downloaded to `openfold/resources/params` by the installation script. Note that the latter
replaces the original, more restrictive CC BY-NC 4.0 license as of January 2022.

## Changes from OpenFold
- TBD
