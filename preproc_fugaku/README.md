# 富士通研究所最適化版OpenFold前処理（MSA/テンプレート検索）スクリプト

## 動作原理

TBD

## 使用方法

1. [`scripts/setenv`](../scripts/setenv)の以下の環境変数をセットする
    * `$PREFIX`: OpenFoldパッケージ・ライブラリ・実行ファイル等のパス
    * `$OPENFOLDDIR`: OpenFoldリポジトリのパス
    * `$DATADIR`: 各種データベースを含むディレクトリのパス
    * `$PJSUB_OPT`: `pjsub`に渡されるオプション

2. `$DATADIR`に[AlphaFold2のディレクトリ構造](https://github.com/deepmind/alphafold/tree/v2.3.1#genetic-databases)と同様にMGnify, PDB70, small BFD, Uniref90を配置する

3. 入力する配列を1つのFASTAファイルに結合する
    * pdb_mmcifや複数の.cifファイルを前処理する場合、[`data_dir_to_fasta.py`](..//scripts/data_dir_to_fasta.py)を使用する

4. [`submit.sh`](submit.sh)の以下の環境変数をセットする
    * `$InputFile`: 3. のFASTAファイルのパス
    * `$OutputDir`: 前処理結果を出力するディレクトリのパス（ジョブ投入時に存在しない場合、ジョブ実行時に自動作成）
    * `DoStaging`: Pythonモジュール、実行ファイル、データベースについてLLIO transfer（富岳）またはメモリステージング（それ以外）を行うかどうか。富岳では各SIO（87 GiB/ノード）が使用するデータベース全体を保持する必要があるため、十分なノード数（例えば>=12ノード以上）を使用する場合に有効にする
    * `$NumNodes`: 各ジョブのノード数
    * `$JobTime_*`: 各ジョブの制限時間
    * `$Timeout_*`: 各ジョブ中の検索ツールの入力配列ごとの制限時間
    * `$NumProcs*`, `$NumThreads*`: 各検索ツールのノード内プロセス数、プロセス内スレッド数

5. `submit.sh`を実行してジョブを投入する

## トラブルシューティング

### LLIO transfer時に容量不足のエラーが表示される

各データベースとPythonモジュール等の容量の合計がジョブの使用できるSIOの最大容量を超える場合はLLIO transferが失敗します。失敗した場合でも前処理は実行できますが、速度が低下する場合があります。ある程度の個数（>10配列）の前処理を行う場合はノード数 `$NumNodes`を増やすことを検討してください。

### mpi4pyが存在せず実行に失敗する

本最適化実装ではオリジナルのOpenFoldでは使用しないmpi4pyを使用します。オリジナルのOpenFold用に構築した環境を使用する場合は、[インストールスクリプト](../scripts/install_fugaku_others.sh)を参考にしてmpi4pyを追加で導入してください。

### `$OutputDir`内に拡張子が`.temp`のファイルが生成される

`.temp`ファイルはMSA/テンプレートファイルを書きだす際に一時的に生成されるファイルで、ジョブやプロセスが強制終了されると削除されずに残る場合があります。削除してください。
