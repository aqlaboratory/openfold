# 富士通最適化版 推論 (inference/relaxation) スクリプト

## 事前準備
1. [README_org.md](../README_org.md#training)を参考にmmcif_cacheとchain_data_cacheを作成する
2. [前処理ディレクトリ](../preproc_fugaku)のスクリプトを使用し前処理を行う (マルチノード版の場合は必須)

## 使用方法

### シングルノード版

1. `inference/parameters`の以下の必須項目を設定する
    - `MMCIFCache`: 事前準備で作成したmmcifキャッシュのパス
    - `InputFastaDir`: 入力シーケンスのfastaファイルが含まれるディレクトリのパス
    - `AlignmentDir`: 前処理で出力されたalignmentディレクトリのパス。このスクリプトで前処理も行う場合には指定しない。
    - `OutputDir`: 出力ディレクトリのパス。`$LOGDIR`とした場合はログディレクトリとなる。

1. 必要があれば`inference/parameters`の以下の項目を変更する
    - `--jax_param_path`: 使用するAlphaFold2パラメータ
    - `max_template_date`: 指定した日付以前のタンパク質構造をテンプレートとして使用する

1. `Submit_inference`により推論のジョブ(1ノード)を投入する
    - `./Submit_inference $TimeLimit`
      - `TimeLimit`: ジョブの時間制限
    - 例
      - 1時間の時間制限で実行する場合: `./Submit_inference 1:00:00`

1. 実行結果は`log/ノード数/Submit_inference.*`に出力される
    - プロセスの出力は`output/0/out.1.0`に書き出される

### マルチノード版

1. `inference/parameters_multi`の以下の必須項目を設定する
    - `MMCIFCache`: 事前準備で作成したmmcifキャッシュのパス
    - `InputFasta`: (フィルタ済みの) 入力シーケンスのfastaファイルのパス
    - `AlignmentDir`: 前処理で出力されたalignmentディレクトリのパス
    - `AlignmentLogDir`: 前処理で出力されたログディレクトリのパス
    - `OutputDir`: 出力ディレクトリのパス。`$LOGDIR`とした場合はジョブ毎に生成されるログディレクトリとなる。
    - `Timeout`: 入力シーケンスごとのタイムアウト時間 [秒]

2. 必要があれば`inference/parameters_multi`の以下の項目を変更する
    - `--jax_param_path`: 使用するAlphaFold2パラメータ
    - `max_template_date`: 指定した日付以前のタンパク質構造をテンプレートとして使用する
    - `--ignore_timeout_chain_history`: 指定したjobid以前のジョブでタイムアウトで失敗した履歴を無視し、再実行する
    - `--ignore_failed_chain_history`: 指定したjobid以前のジョブでメモリ不足などで失敗した履歴を無視し、再実行する

3. ノード数と制限時間を決める
    - ノード数: 入力シーケンス数以下の数
    - 制限時間: 任意の時間。`estimate_time.awk`を用いて処理時間を推定し、(ノード数)×(制限時間)がおよそ推定処理時間となるように設定してもよい。
      - `./estimate_time.awk $InputFasta`

4. `Submit_inference_multi`により推論のジョブを投入する
    - `./Submit_inference_multi {$NumNodes|$NodeShape} $TimeLimit`
    - `NumNodes`: ノード数
    - `NodeShape`: 3次元ノード形状 (XxYxZ)
    - `TimeLimit`: ジョブの時間制限
    - 例
      - 16ノードで1時間の時間制限で実行する場合: `./Submit_inference_multi 16 1:00:00`
      - 2x3x8のノード形状で6時間の時間制限で実行する場合: `./Submit_inference_multi 2x3x8 6:00:00`

5. 実行ログは`log/ノード数/Submit_inference_multi.*`に出力されるので、必要に応じて確認する
    - 各プロセスの出力は`output/0/out.1.*`に書き出される
    - シーケンス毎の処理結果は`$OutputDir/result/*.csv`に出力される

6. 未処理のシーケンスがある場合は`4.ジョブの投入`を再度行う
