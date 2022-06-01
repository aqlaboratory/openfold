# OpenFold for Fugaku

1. OpenFoldをgit cloneもしくはコピーし、そのディレクトリに移動
2. DeepSpeedCPUを`third_party`以下にgit cloneもしくはコピー
3. stereo_chemical_props.txtを`openfold/openfold/resources/stereo_chemical_props.txt`に配置
    - 共有ディレクトリの`misc/stereo_chemical_props.txt`に以下のURLからダウンロードしたものが置いてある
        - https://git.scicore.unibas.ch/schwede/openstructure/-/raw/7102c63615b64735c4941278d92b554ec94415f8/modules/mol/alg/src/stereo_chemical_props.txt
4. `scripts/setenv`を編集しパスを設定
    - PREFIX: インストール先ディレクトリ (必須)
    - OPENFOLDDIR: このREADMEがあるディレクトリ (必須)
    - DATADIR: データディレクトリ (必須)
        - 現在はこのディレクトリ `/vol0003/share/hp210263/OpenFold/data`
    - PJSUB_OPT: pjsub時の追加オプション (任意)
        - シングルアカウントの場合は`-g $(グループID)`を指定する
6. インストール用バッチジョブを投入するスクリプトを実行
```
./scripts/install_fugaku.sh
```
