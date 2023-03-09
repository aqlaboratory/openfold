# 富士通最適化版 推論 (inference/relaxation) スクリプト

## 事前準備
1. 前処理を実行する
2. [README](../README.md)を参考にmmcifキャッシュを作成する

## 使用方法

### シングルノード版

1. 

### マルチノード版

1. `inference/parameters_multi`の以下の必須項目を設定する
  - `MMCIFCache`: 事前準備で作成したmmcifキャッシュのパス
  - `InputFasta`: 入力シーケンスのfastaファイルのパス
  - `AlignmentDir`: 前処理で出力されたalignmentディレクトリのパス
  - `OutputDir`: 出力ディレクトリのパス
  - `Timeout`: 入力シーケンスごとのタイムアウト時間 [秒]

2. 必要があれば`inference/parameters_multi`の以下の項目を変更する
  - `--jax_param_path`: 使用するAlphaFold2パラメータ
  - `max_template_date`: 指定した日付以前のタンパク質構造をテンプレートとして使用する

3. ノード数と制限時間を決める
  - ノード数: 入力シーケンス数以下の数
  - 制限時間: 任意の時間。`estimate_time.awk`を用いて処理時間を推定し、(ノード数)×(制限時間)がおよそ推定処理時間となるように設定してもよい。

4. `Submit_inference_multi`により推論のジョブを投入する
  - 16ノードで1時間の時間制限で実行する場合: `./Submit_inference_multi 16 1:00:00`
  - 2x3x8のノード形状で6時間の時間制限で実行する場合: `./Submit_inference_multi 2x3x8 6:00:00`

5. 実行の状況は`log/ノード数/Submit_inference_multi.*`に出力されるログを確認する
  - 各プロセスの出力は`output/0/out.1.*`に書き出される

6. ジョブ終了後に推論の結果を調査する
