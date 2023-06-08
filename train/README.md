# 富士通最適化版 学習 (training) スクリプト

## 準備
1. [README_org.md](../README_org.md#training)を参考にmmcif_cacheとchain_data_cacheを作成する
2. [前処理ディレクトリ](../preproc_fugaku)のスクリプトを使用し前処理を行う


### mmcifの変換 (任意)
事前にmmcifをパースしpickle化およびlz4圧縮します。特に複数ノードによる学習での速度低下を防ぐことができます。
1. [convert_mmcifディレクトリ](../convert_mmcif)に移動する
2. run.shの`OutputDir`に変換後の出力先ディレクトリを設定する
3. ジョブを投入する
```
pjsub -g <グループID> run.sh
```
4. trainやinferenceディレクトリの`parameters`ファイルを編集し、`MMCIFDir`を変換後のmmcifディレクトリに置き換える

## 使用方法

1. `parameters`の以下の項目を設定する
    - `AlignmentDir` : 前処理で出力されたalignmentのディレクトリ
    - `MMCIFCache` : 準備1.で生成したmmcif_cache.jsonのパス
    - `ChainDataCache` : 準備1.で生成したchain_data_cache.jsonのパス

1. `Submit_training`により学習のジョブを投入する
    - `./Submit_training <NodeShape> <TimeLimit>`
      - `NodeShape` : ジョブを実行するノードの3次元形状
      - `TimeLimit` : ジョブの最長実行時間
    - 例
      - 2x2x2ノードで1時間の時間制限で実行する場合: `./Submit_training 2x2x2 1:00:00`

1. 実行結果は`log/<ノード数>/Submit_training.*`に出力される
    - 各プロセスの出力は`output/0/out.1.*`に書き出される
