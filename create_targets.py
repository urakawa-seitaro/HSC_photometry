#!/Users/coias/opt/anaconda3/envs/astropy/bin/python3
# create_targets.py

import os
import glob
import sys

def create_targets_file():
    """
    カレントディレクトリ内の全てのFITSファイルと指定された単一座標を使用して、
    targets.txt ファイルを作成します。
    """
    
    # ユーザーから測光対象の座標を入力として受け取る
    print("--- targets.txt 作成 ---")
    
    # RAの入力
    ra = input("測光を行う天体の赤経 (RA) を入力してください (例: 123.456 または 10h01m02s): ").strip()
    if not ra:
        print("エラー: RAが入力されませんでした。処理を中断します。")
        sys.exit(1)
        
    # Decの入力
    dec = input("測光を行う天体の赤緯 (Dec) を入力してください (例: +10.123 または -05d30m00s): ").strip()
    if not dec:
        print("エラー: Decが入力されませんでした。処理を中断します。")
        sys.exit(1)
        
    # FITSファイルリストの取得
    # '*' はワイルドカードで、カレントディレクトリの全ての .fits ファイルを検索します。
    fits_files = glob.glob("*.fits")
    
    if not fits_files:
        print("警告: カレントディレクトリにFITSファイル (*.fits) が見つかりませんでした。")
        sys.exit(0)
        
    output_filename = "targets.txt"
    
    # ファイル書き込み
    with open(output_filename, 'w') as f:
        f.write("# FITS_FILE_NAME\tRA (deg or h:m:s)\tDec (deg or d:m:s)\n")
        
        for filename in sorted(fits_files):
            # ファイル名、RA、Dec をタブ区切りで書き込む
            line = f"{filename}\t{ra}\t{dec}\n"
            f.write(line)

    print("\n--- 作成完了 ---")
    print(f"FITSファイル {len(fits_files)} 件に対応する '{output_filename}' を作成しました。")
    print(f"共通の座標: RA={ra}, Dec={dec}")
    
    # 作成されたファイルの最初の数行を表示
    print("\n[targets.txt の内容 (先頭5行)]")
    with open(output_filename, 'r') as f:
        for i, line in enumerate(f):
            print(line.strip())
            if i >= 4:
                break

if __name__ == "__main__":
    create_targets_file()