# HSC_photometry
Photometry of HSC for the specific target (Ra, Dec)

＊COIASの画像中から変光星が検出された時に、その変光星の測光を行う。

＊COIASの画像は、移動天体探しのため同じ領域を4枚以上撮影された画像を用いているが、COIASに実装されていない同領域のHSC-PDR3の画像も解析対象とする。

主な環境(conda env exportの結果）

python=3.13.1
astropy=7.0.0
photutils=2.2.0
matplotlib-base=3.9.3

＃使い方

./create_target.py 　#処理するfitsファイルの自動記入と測光したい天体の座標を記入

実行後 targets.txt ができあがる。

./run_all_photometry.sh でバッチ処理。hsc_photometry_centroid.pyが実行される。

hsc_photometry_centroid.pyの説明

ra,decをX,Y座標に変換。X,Y座標を中心に±25pixel(50X50pixel)の領域で光度重心を測定

光度重心のX,Y座標に対してaperture測光。測光半径は12pixel

誤差はhdu[3].dataのvarianceを利用

出力X,Y座標はds9表記
