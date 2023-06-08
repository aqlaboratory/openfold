# 富士通最適化版 学習 (training) スクリプト

## 準備
1. 前処理を実行する
2. [README](../README.md)を参考にmmcif_cacheとchain_data_cacheを作成する

### mmcifの変換 (任意)
事前にmmcifをパースしpickle化およびlz4圧縮します。特に複数ノードによる学習での速度低下を防ぐことができます。
1. convert_mmcifディレクトリに移動
2. run.shのInputDirとOutputDirを設定
3. ジョブ投入
```
pjsub -g <グループID> run.sh
```
4. 学習や推論ディレクトリのparametersに記載されているmmcifのディレクトリを変換したものに置き換える

## 使用方法

1. `parameters`の以下の必須項目を設定する
    - ...

1. 必要におうじて`parameters`の以下の項目を変更する
    - ...

1. `Submit_training`により学習のジョブを投入する
    - `./Submit_training <NodeShape> <TimeLimit>`
      - `NodeShape` : ジョブを実行するノードの3次元形状
      - `TimeLimit` : ジョブの最長実行時間
    - 例
      - 2x2x2ノードで1時間の時間制限で実行する場合: `./Submit_training 2x2x2 1:00:00`

1. 実行結果は`log/<ノード数>/Submit_training.*`に出力される
    - 各プロセスの出力は`output/0/out.1.*`に書き出される
