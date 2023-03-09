# 富士通最適化版 推論 (inference/relaxation) スクリプト

## 事前準備
1. 前処理を実行する
2. [README](../README.md)を参考にmmcifキャッシュを作成する

## 使用方法

### シングルノード版

1. 

### マルチノード版

1. `filter_fasta.sh`を用いて入力のfastaファイルから処理するシーケンスをフィルタする
  - `./filter_fasta.sh $InputFasta $OutputFasta $AlignmentDir [$PredictionDir [$NGList]]`
    - `InputFasta`: 入力のfastaファイルのパス
    - `OutputFasta`: 出力されるfastaファイルのパス
    - `AlignmentDir`: 前処理で出力されたalignmentディレクトリのパス
    - [任意] `$PredictionDir`: 2回目以降のジョブで、推論結果が含まれるディレクトリのパス
    - [任意] `$NGList`: 2回目以降のジョブで、失敗したシーケンス名のリストのパス

1. `inference/parameters_multi`の以下の必須項目を設定する
  - `MMCIFCache`: 事前準備で作成したmmcifキャッシュのパス
  - `InputFasta`: (フィルタ済みの) 入力シーケンスのfastaファイルのパス
  - `AlignmentDir`: 前処理で出力されたalignmentディレクトリのパス
  - `OutputDir`: 出力ディレクトリのパス
  - `Timeout`: 入力シーケンスごとのタイムアウト時間 [秒]

1. 必要があれば`inference/parameters_multi`の以下の項目を変更する
  - `--jax_param_path`: 使用するAlphaFold2パラメータ
  - `max_template_date`: 指定した日付以前のタンパク質構造をテンプレートとして使用する

1. ノード数と制限時間を決める
  - ノード数: 入力シーケンス数以下の数
  - 制限時間: 任意の時間。`estimate_time.awk`を用いて処理時間を推定し、(ノード数)×(制限時間)がおよそ推定処理時間となるように設定してもよい。

1. `Submit_inference_multi`により推論のジョブを投入する
  - `Submit_inference_multi {$NumNodes|$NodeShape} $TimeLimit`
  - `NumNodes`: ノード数
  - `NodeShape`: 3次元ノード形状 (XxYxZ)
  - `TimeLimit`: ジョブの時間制限
  - 例
    - 16ノードで1時間の時間制限で実行する場合: `./Submit_inference_multi 16 1:00:00`
    - 2x3x8のノード形状で6時間の時間制限で実行する場合: `./Submit_inference_multi 2x3x8 6:00:00`

1. 実行結果は`log/ノード数/Submit_inference_multi.*`に出力される
  - 各プロセスの出力は`output/0/out.1.*`に書き出される

1. 未処理のシーケンスがある場合は再度ジョブを実行する
  - 実行に失敗したシーケンスのリスト (NGリスト) を作成
    - `find_ng.sh $LogDir > $NGList`により正常に推論できなかったシーケンス名のリストを生成する
  - 再度fastaファイルをフィルタし、`inference/parameters_multi`の`InputFasta`を変更
  - ジョブを投入


