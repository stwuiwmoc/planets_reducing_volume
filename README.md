# combined_minimize.py

- 1024x1024の曲面形状測定結果に対して差分が最小になるような軸外し切り取り放物面を、3パラメーターの非線形最適化と行列計算で導出します
- zernike多項式の第10項までで表される特徴的な鏡面形状に対しても同様の導出ができます。
- 研究自体の背景・学術的成果に関してはP-CG18-P05.pdfをご覧ください。

# DEMO
1024x1024の鏡面形状測定値(.csv)と10x36の駆動行列(.csv)を用いると、最適化結果を以下のように図示します
![wh4_r00_id1_f1100](https://user-images.githubusercontent.com/82031414/120265031-12dc0400-c2da-11eb-9cb4-99ad706557d2.png)

# Features
- 放物面のパラメータ調整による最適化と駆動行列による最適化を、選択または組み合わせることができます。
- IPythonコンソールでinputの指示に従って入力するだけで、測定結果/zernike多項式の切り替え、最適化を行う内容などを選べます。
- scipyの効率的な最適化関数によって、一つの目標形状に対して30分程度の計算時間で終了します。

# Requirement

- numpy      : 1.19.2
- scipy      : 1.6.1
- pandas     : 1.2.3
- matplotlib : 3.3.4
- PIL        : 8.1.2
- proper     : 3.2.4

# Usage

1. 同階層ディレクトリに、csv形式で曲面形状測定結果と駆動行列を配置します
    - 駆動行列はmake_opration_matrix.pyと_Fxxフォルダ内の有限要素法解析データを用いて作成できます。追加ライブラリは必要ありません。
2. ファイルを実行すると、IPythonコンソールで入力が促されます。
3. Input code type
    - 一つの曲面形状について最適化するか、zernike多項式第10項までの各項全てに対して最適化するかを選択します。
    - Input data type : 上記でInputを選択した場合、何に対して最適化するかを決定します。
4. filter select
    - 曲面形状に対して適用する平滑化フィルタの種類を選択します。
6. Input filter size or sigma
    - 平滑化フィルタの平滑化サイズ（ガウシアンフィルタの場合はσ）を入力します。
7. Ideal Surface Minimize
    - 3パラメーターの調整による最適化を行うかどうかを選択します。
8. Warping Harness Minimize
    - 駆動行列による最適化を行うかどうかを選択します。
9. 最適化が開始し、図が /mkfokder/combined_minimize/ に出力されます。
