#

1. OpenFoldをgit cloneもしくはコピーし、そのディレクトリに移動
2. DeepSpeedCPUをthird_party以下にgit clone
3. stereo_chemical_props.txtのダウンロード
    - 以下のURLからファイルをダウンロードし、`openfold/openfold/resources/stereo_chemical_props.txt`に配置
        https://git.scicore.unibas.ch/schwede/openstructure/-/raw/7102c63615b64735c4941278d92b554ec94415f8/modules/mol/alg/src/stereo_chemical_props.txt
4. `scripts/setenv`を編集しパスを設定
5. インストールを実行
    ```
./scripts/install_fugaku.sh`を実行
    ```