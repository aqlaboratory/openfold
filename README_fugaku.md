# OpenFold for Fugaku

## コードの準備 (/vol0006/ra000018/data/OpenFold/src/openfold_v1.0.0をコピーしてきた際は不要)
1. OpenFoldをgit cloneもしくはコピーし、そのディレクトリに移動
2. DeepSpeedCPUを`third_party`以下にgit cloneもしくはコピー
3. stereo_chemical_props.txtを`openfold/openfold/resources/stereo_chemical_props.txt`に配置
    - 共有ディレクトリの`misc/stereo_chemical_props.txt`に以下のURLからダウンロードしたものが置いてある
        - https://git.scicore.unibas.ch/schwede/openstructure/-/raw/7102c63615b64735c4941278d92b554ec94415f8/modules/mol/alg/src/stereo_chemical_props.txt

## 初期設定
`scripts/setenv`を編集し、以下の変数に値を設定

- PREFIX: インストール先ディレクトリ (必須)
    - ビルド済み実行環境を使う場合は `/vol0006/ra000018/data/OpenFold/local`
- OPENFOLDDIR: このREADMEがあるディレクトリ (必須) 
- DATADIR: データディレクトリ (必須)
    - 現在はこのディレクトリ `/vol0006/ra000018/data/OpenFold/data`
- PJSUB_OPT: pjsub時の追加オプション (任意)
    - シングルアカウントの場合は`-g $(グループID)`を指定する
    - 場合によっては`-x PJM_LLIO_GFSCACHE=/vol0006`が必要

## インストール (PREFIXにビルド済み実行環境を指定した場合は不要)
- インストール用バッチジョブを投入するスクリプトを実行
```
./scripts/install_fugaku.sh
```

## 前処理と推論
### 少数の配列の場合
inferenceディレクトリのシングルノード版スクリプトで前処理と推論を同時に行う。詳細は[こちら](inference/README.md)

### 多数の配列の場合
1. preproc_fugakuディレクトリのマルチノード版スクリプトで前処理を行う。詳細は[こちら](preproc_fugaku/README.md)
2. inferenceディレクトリのスクリプトで推論を行う。詳細は[こちら](inference/README.md)


## 学習の実行
1. trainディレクトリに移動
2. parametersファイルを設定 (デフォルトではsampleディレクトリの101M_Aで事前学習)。データセットにPDBの1/100を使用する場合はparameters_rand1891をparametersに上書き
3. ジョブ投入
```
./Submit_training <NodeShape> <TimeLimit>

# 例: 2x2x2ノードで2時間制限で実行する場合
# ./Submit_training 2x2x2 2:0:0
```
4. logディレクトリに出力される結果を確認

## 前処理の実行
1. preprocディレクトリに移動
2. ジョブ投入
```
./Submit_preproc <TimeLimit> <InputDir> <OutputDir>

# 例: 1時間制限で、inputディレクトリ中の全てのfastaを前処理し、outputディレクトリに出力する場合
./Submit_preproc 1:0:0 input output
```

## mmcifの変換 (任意)
mmcifをpythonでパースする処理が遅いので、事前にパースしpickle化およびlz4圧縮しておくことができる。
1. convert_mmcifディレクトリに移動
2. run.shのInputDirとOutputDirを設定
3. ジョブ投入
```
pjsub -g {グループID} run.sh
```
4. 学習および推論ディレクトリのparametersに記載されているmmcifのディレクトリを変換したものに置き換える
