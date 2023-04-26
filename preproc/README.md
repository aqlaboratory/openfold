# 通常版 前処理（MSA/テンプレート検索）スクリプト

fastaファイルに対して前処理を行うスクリプトです。複数ノードによるfastaファイル単位での並列化に対応していますが、高速に処理する必要がある場合には処理単位での並列化を行う[preproc_fugaku](../preproc_fugaku)を使用してください。

## 使用方法

1. 出力先ディレクトリが存在しない場合は作成する
2. `./Submit_preproc`でジョブを投入する

```
./Submit_preproc <Number of nodes> <Time limit> <Input dir> <Output dir> [<Input file list>]

# 例: 2ノード、1時間制限で、example/input内のfastaファイルを前処理してexample/alignmentに出力する
./Submit_preproc 2 1:00:00 example/input example/alignment
```
    * Input file list: Input dir内で処理対象とするファイル名を1行1ファイル名で記載したファイル。指定がなければInput dirの全てのファイルが処理対象となる。

