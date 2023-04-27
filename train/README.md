# 富士通最適化版 学習 (training) スクリプト

## 事前準備
1. 前処理を実行する
2. [README](../README.md)を参考にmmcif_cacheとchain_data_cacheを作成する

## 使用方法

1. `parameters`の以下の必須項目を設定する
    - ...

1. 必要があれば`parameters`の以下の項目を変更する
    - ...

1. `Submit_training`により学習のジョブを投入する
    - `Submit_training <NodeShape> <TimeLimit>`
    - `TimeLimit`: ジョブの時間制限
    - 例
      - 2x2x2ノードで1時間の時間制限で実行する場合: `./Submit_training 2x2x2 1:00:00`

1. 実行結果は`log/<ノード数>/Submit_training.*`に出力される
    - 各プロセスの出力は`output/0/out.1.*`に書き出される
