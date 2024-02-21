# 富士通最適化版 前処理（MSA/テンプレート検索）スクリプト

## 出力フォーマットについて
各データベースについて検索を実行した場合、以下のファイルが出力されます。

### `SubDirectorySize=0`の場合
```
（$OutputDirのパス）
└──（タンパク質名）
　   ├── mgnify_hits.a3m    : MGnifyから検索したMSA
　   ├── pdb70_hits.hhr     : PDB70から検索したテンプレート
　   ├── small_bfd_hits.a3m : small BFDから検索したMSA（$ConvertSmallBFDToA3M=1の場合）
　   ├── small_bfd_hits.sto : small BFDから検索したMSA（$ConvertSmallBFDToA3M=0の場合）
　   └── uniref90_hits.a3m  : UniRef90から検索したMSA
```

### `SubDirectorySize>0`の場合
```
（$OutputDirのパス）
└──（サブディレクトリ）
　 └──（タンパク質名）
　 　   ├── mgnify_hits.a3m    : MGnifyから検索したMSA
　 　   ├── pdb70_hits.hhr     : PDB70から検索したテンプレート
　 　   ├── small_bfd_hits.a3m : small BFDから検索したMSA（$ConvertSmallBFDToA3M=1の場合）
　 　   ├── small_bfd_hits.sto : small BFDから検索したMSA（$ConvertSmallBFDToA3M=0の場合）
　 　   └── uniref90_hits.a3m  : UniRef90から検索したMSA
```

`small_bfd_hits.a3m`と`small_bfd_hits.sto`には以下の違いがあります。

* 一般に`a3m`のほうがファイルサイズが小さくなります。
* `$MaxHits_small_bfd`を設定している場合、ヒット件数がその値に制限されます。制限しない場合、一部のタンパク質ではファイルサイズが非常に大きくなる場合があります。

## 使用方法

1. [`scripts/setenv`](../scripts/setenv)の以下の環境変数をセットする
    * `$PREFIX`: OpenFoldパッケージ・ライブラリ・実行ファイル等のパス
    * `$OPENFOLDDIR`: OpenFoldリポジトリのパス
    * `$DATADIR`: 各種データベースを含むディレクトリのパス
    * `$PJSUB_OPT`: `pjsub`に渡されるオプション
    * 以下の[`preproc_fugaku`](.)用設定（デフォルト値は富岳のもの）
	* `OPENFOLD_MAX_MEM`: ノードあたりの最大使用メモリサイズ
	* `OPENFOLD_PBS_RESOURCE`: （PBSジョブスケジューラを用いる場合のみ）リソース名
	* `TMPDIR`: テンポラリファイル用ディレクトリ

2. `$DATADIR`に[AlphaFold2のディレクトリ構造](https://github.com/deepmind/alphafold/tree/v2.3.1#genetic-databases)と同様にMGnify, PDB70, small BFD, Uniref90を配置する

3. 入力する配列を1つのFASTAファイルに結合する
    * pdb_mmcifや複数の`.cif`ファイルを前処理する場合、[`data_dir_to_fasta.py`](..//scripts/data_dir_to_fasta.py)を使用する

4. [`submit.sh`](submit.sh)の以下の環境変数をセットする
    * `$InputFile`: 3. のFASTAファイルのパス
    * `$OutputDir`: 前処理結果を出力するディレクトリのパス（ジョブ投入時に存在しない場合、ジョブ実行時に自動作成）
    * `$TempDir`: 一部の前処理結果ファイルがジョブ実行中に一時的に保存されるディレクトリのパス。GBレベルのファイルが作成される可能性があるため、`/tmp`や`/worktmp`は非推奨
    * `$DoStaging`: Pythonモジュール、実行ファイル、データベースについてLLIO transfer（富岳）またはメモリステージング（それ以外）を行うかどうか。富岳では各SIO（87 GiB/ノード）が使用するデータベース全体を保持する必要があるため、十分なノード数（例えば>=12ノード以上）を使用する場合に有効にする
    * `$ConvertSmallBFDToA3M`: small BFDについて、出力ファイルをstoからa3mに変換する
    * `$CreateDirOnDemand`: 各配列の出力ディレクトリについて、ジョブ開始時ではなく前処理ファイルを書き出す直前に作成する。配列数が多くファイル操作に時間がかかる場合に有効
        * 後述の`$SubDirectorySize`と同時に使用する場合、サブディレクトリも書き出す直前に作成する
    * `$SubDirectorySize`: 入力FASTAファイルの先頭から指定された配列の個数ごとに"0"から始まるサブディレクトリに分割して出力を行う。0の場合はサブディレクトリを作成しない
    * `$LimitMaxMem*`: プロセスあたりの最大メモリ量を（ノードメモリ量）/（ノード内プロセス数）に制限するかどうか
    * `$NumNodes`: 各ジョブのノード数
    * `$JobTime_*`: 各ジョブの制限時間
    * `$Timeout_*`: 各ジョブ中の検索ツールの入力配列ごとの制限時間
    * `$NumProcs*`, `$NumThreads*`: 各検索ツールのノード内プロセス数、プロセス内スレッド数
    * `$StreamSTOSize`: MSAとしてこのサイズを超えるstoファイルが出力された場合、a3mへの変換をストリーミングで行う実装を使用する。小さな値を設定するほど各プロセスのメモリ消費量が小さくなるが、変換速度は低下する
    * `$MaxHits_*`: 各ジョブの最大MSAヒット件数。`$MaxHits_*=""`の場合、ヒット件数を制限しない

5. `submit.sh`を実行してジョブを投入する

6. 一部のファイルの前処理が完了せずにジョブが終了する場合、`$NumNodes`等を調整して再度ジョブを投入する
    * `$InputFile`の内容が同じ場合、前処理未完了の配列のみが自動で検出されて前処理が行われる

## トラブルシューティング

### LLIO transfer時に容量不足のエラーが表示される

各データベースとPythonモジュール等の容量の合計がジョブの使用できるSIOの最大容量を超える場合はLLIO transferが失敗します。失敗した場合でも前処理は実行できますが、速度が低下する場合があります。ある程度の個数（>10配列）の前処理を行う場合はノード数 `$NumNodes`を増やすことを検討してください。

### mpi4pyが存在せず実行に失敗する

本最適化実装ではオリジナルのOpenFoldでは使用しないmpi4pyを使用します。オリジナルのOpenFold用に構築した環境を使用する場合は、[インストールスクリプト](../scripts/install_fugaku_others.sh)を参考にしてmpi4pyを追加で導入してください。

### `$OutputDir`内に拡張子が`.temp`のファイルが生成される

`.temp`ファイルはMSA/テンプレートファイルを書きだす際に一時的に生成されるファイルで、ジョブやプロセスが強制終了されると削除されずに残る場合があります。削除してください。

### `$OutputDir`に生成されるディレクトリ名が`$InputFile`の配列のラベル名と異なる
現在の実装では同じラベル名の複数の配列が入力されることを許容していますが、出力パスの衝突を避けるために末尾に適当な文字列が追加されます。

### PBS等の一般的なスパコン・クラスタ上でPDB70についての検索ジョブが失敗する
PDB70を用いるテンプレート検索はUniref90のMSA検索結果を使用するため、Uniref90のジョブが先に完了している必要があります。しかし、現在の実装では富岳の`pjsub`によるステップジョブにのみ対応しているため、他のジョブスケジューラの場合は`Do_uniref90`と`Do_pdb70`を同時に有効にするとこれら2つのジョブが同時に実行され、検索処理が失敗する場合があります。先にUniref90のジョブが完了したことを確認したうえでPDB70のジョブを投入するか、ジョブスケジューラに実装されている方法でジョブの依存関係を設定してください。
