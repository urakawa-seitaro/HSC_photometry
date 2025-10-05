# HSC_photometry
Photometry of HSC for the specific target (Ra, Dec)

＊COIASの画像中から変光星が検出された時に、その変光星の測光を行う。

＊COIASの画像は、移動天体探しのため同じ領域を4枚以上撮影された画像を用いているが、COIASに実装されていない同領域のHSC-PDR3の画像も解析対象とする。

主な環境(conda env exportの結果）

python=3.13.1
astropy=7.0.0
photutils=2.2.0
matplotlib-base=3.9.3
