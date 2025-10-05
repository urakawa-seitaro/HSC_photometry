#!/bin/bash

# 実行するPythonスクリプト名
PYTHON_SCRIPT="hsc_photometry_centroid.py"

# 入力リストファイル名
TARGETS_FILE="targets.txt"

# ----------------------------------------------------
# 実行前のチェックと権限付与 (★★★ 修正箇所 ★★★)
# ----------------------------------------------------

# Pythonスクリプトの存在確認
if [ ! -f "$PYTHON_SCRIPT" ]; then
    echo "エラー: Pythonスクリプト '$PYTHON_SCRIPT' が見つかりません。"
    exit 1
fi

# Pythonスクリプトに実行権限を付与
chmod +x "$PYTHON_SCRIPT"
echo "Pythonスクリプトに実行権限を付与しました: $PYTHON_SCRIPT"

# ターゲットファイルの存在確認
if [ ! -f "$TARGETS_FILE" ]; then
    echo "エラー: ターゲットリストファイル '$TARGETS_FILE' が見つかりません。"
    echo "ファイル名, RA, Decを記載した '$TARGETS_FILE' を作成してください。"
    exit 1
fi

echo "--- 測光処理を開始します ---"

# targets.txt を行ごとに読み込み、コメント行 (#) と空行をスキップ
grep -vE '^\s*#|^\s*$' "$TARGETS_FILE" | while IFS=$'\t ' read -r fits_file ra dec; do
    
    # FITSファイルの存在確認
    if [ ! -f "$fits_file" ]; then
        echo "警告: FITSファイル '$fits_file' が見つかりません。スキップします。"
        continue
    fi
    
    echo "--------------------------------------------------------"
    echo "実行中: $fits_file (RA: $ra, Dec: $dec)"
    
    # Pythonプログラムの実行: 'python' コマンドを使わず、直接実行する形式に変更 (シバンのおかげ)
    # ./[スクリプト名] [FITSファイル] [RA] [DEC]
    ./"$PYTHON_SCRIPT" "$fits_file" "$ra" "$dec" # ★呼び出し形式を修正★
    
    # 実行結果の確認
    if [ $? -ne 0 ]; then
        echo "致命的なエラー: $fits_file の処理中にエラーが発生しました。スクリプトを中断します。"
        exit 1
    fi
    
done

echo "--------------------------------------------------------"
echo "--- 全ての測光処理が完了しました。 ---"
echo "結果は 'photometry_results_centroid.txt' に出力されています。"