#!/Users/coias/opt/anaconda3/envs/astropy/bin/python3
# hsc_photometry_centroid.py

import argparse
import os
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord, Angle
from astropy import units as u
from astropy.wcs import WCS
from astropy.time import Time
from photutils.aperture import CircularAperture, aperture_photometry
from photutils.centroids import centroid_2dg # ★光度重心計算のために追加★

# ------------------------------------------------------------------------------
# ヘルパー関数 (変更なし)
# ------------------------------------------------------------------------------
def get_hsc_calibration_info(fits_file_path):
    """
    FITSファイルから観測時刻(UT)、フィルター名、FLUXMAG0を取得する。
    """
    try:
        with fits.open(fits_file_path) as hdul:
            primary_header = hdul[0].header
            
            if 'TIME-OBS' in primary_header:
                t = Time(primary_header['TIME-OBS'], format='isot', scale='utc')
            elif 'DATE-AVG' in primary_header:
                t = Time(primary_header['DATE-AVG'], format='isot', scale='tai')
            elif 'TIME-MID' in primary_header:
                t = Time(primary_header['TIME-MID'], format='isot', scale='utc')
            else:
                raise KeyError("時刻情報ヘッダ ('TIME-OBS', 'DATE-AVG', 'TIME-MID') が見つかりません。")
            
            time_obj_utc = t.utc 
            
            try:
                filter_name = primary_header['FILTER'].strip()
            except KeyError:
                raise KeyError("フィルター情報ヘッダ ('FILTER') が見つかりません。")

            try:
                fluxmag0 = float(primary_header['FLUXMAG0'])
            except KeyError:
                try:
                    # FLUXMAG0導出ロジック (省略せず記載)
                    entry_hdu_index = primary_header['AR_HDU'] - 1
                    entry_hdu = hdul[entry_hdu_index]
                    photocalib_id = primary_header['PHOTOCALIB_ID']
                    
                    (pointer_entry,) = entry_hdu.data[entry_hdu.data['id'] == photocalib_id]
                    
                    photocalib_hdu_index = entry_hdu_index + pointer_entry['cat.archive']
                    photocalib_hdu = hdul[photocalib_hdu_index]
                    start_row, end_row = pointer_entry['row0'], pointer_entry['row0'] + pointer_entry['nrows']
                    
                    (photocalib_data,) = photocalib_hdu.data[start_row:end_row]
                    calibration_mean = photocalib_data['calibrationMean']
                    
                    fluxmag0 = (1.0e23 * 10**(48.6 / -2.5) * 1.0e9) / calibration_mean
                    
                except Exception as inner_e:
                    raise KeyError(f"FLUXMAG0の導出に必要な拡張情報が不足しています。詳細: {inner_e}")
            
            exptime = float(primary_header.get('EXPTIME', 1.0))
            
            return {'time_obj_utc': time_obj_utc, 'fluxmag0': fluxmag0, 'exptime': exptime, 'filter_name': filter_name}
            
    except Exception as e:
        print(f"エラー: FITSファイルからのキャリブレーション情報取得に失敗しました。詳細: {e}")
        return None

def parse_ra_dec(ra_str, dec_str):
    """RA, Decの文字列をastropy.coordinates.SkyCoordオブジェクトに変換します。"""
    try:
        ra = float(ra_str) * u.degree
        dec = float(dec_str) * u.degree
        coord = SkyCoord(ra, dec, frame='icrs')
    except ValueError:
        try:
            ra = Angle(ra_str, unit=u.hour)
            dec = Angle(dec_str, unit=u.degree)
            coord = SkyCoord(ra, dec, frame='icrs')
        except Exception:
            raise ValueError("RA/Decの入力形式が不正です。度数または時分秒形式で指定してください。")
            
    return coord

# ------------------------------------------------------------------------------
# Main関数
# ------------------------------------------------------------------------------
def main():
    # 1. 引数の設定 (変更なし)
    parser = argparse.ArgumentParser(
        description="HSC FITS画像から指定座標の星の光度重心を求め、再測光します。"
    )
    parser.add_argument("fits_file", help="解析するFITS画像ファイル名")
    parser.add_argument("ra", help="天体の赤経 (RA)。度数または時分秒形式")
    parser.add_argument("dec", help="天体の赤緯 (Dec)。度数または度分秒形式")
    args = parser.parse_args()

    # 2. 座標のパース (変更なし)
    try:
        target_coord = parse_ra_dec(args.ra, args.dec)
    except ValueError as e:
        print(f"エラー: {e}")
        return

    # 3. FITSファイルの読み込みとデータ準備
    try:
        hdul = fits.open(args.fits_file)
        
        # 科学画像データ (SCI: HDU 1) と WCS情報
        image_data = hdul[1].data.astype(np.float64)
        w = WCS(hdul[1].header)
        
        # 分散データ (VAR: HDU 3)
        variance_data = hdul[3].data.astype(np.float64)

    except Exception as e:
        print(f"エラー: FITSファイルの読み込みまたは必要なHDUの取得に失敗しました。詳細: {e}")
        if 'hdul' in locals() and hdul: hdul.close()
        return

    calibration_info = get_hsc_calibration_info(args.fits_file)
    if calibration_info is None:
        hdul.close()
        return
    
    # 4. 天体座標を初期ピクセル座標に変換
    try:
        # SkyCoordオブジェクト全体を渡す (前回エラー対応)
        pixel_coord = w.world_to_pixel(target_coord) 
    except Exception as e:
        print(f"エラー: 世界座標からピクセル座標への変換に失敗しました。詳細: {e}")
        hdul.close()
        return

    # 最初の中心座標 (floatまたは0D array)
    initial_x, initial_y = pixel_coord

    # 5. 光度重心の計算 (Centroiding)
    
    # 切り出し範囲の設定: 中心から各方向25ピクセル -> 計50x50
    HALF_SIZE = 25
    
    # 切り出し範囲のインデックス (Numpyは [Ymin:Ymax, Xmin:Xmax] の順)
    x_min = max(0, int(initial_x - HALF_SIZE))
    x_max = min(image_data.shape[1], int(initial_x + HALF_SIZE))
    y_min = max(0, int(initial_y - HALF_SIZE))
    y_max = min(image_data.shape[0], int(initial_y + HALF_SIZE))
    
    # 画像と分散データの切り出し
    cutout_data = image_data[y_min:y_max, x_min:x_max]
    cutout_variance = variance_data[y_min:y_max, x_min:x_max]
    
    # 光度重心の計算
    try:
        # cutout_data内の重心 (相対座標) を計算
        centroid_x_rel, centroid_y_rel = centroid_2dg(cutout_data)
        
        # 絶対座標に変換: (切り出し範囲の開始点) + (相対座標)
        centroid_x = centroid_x_rel + x_min
        centroid_y = centroid_y_rel + y_min
        
        # 求めた重心を出力用変数に代入
        final_x_center = centroid_x
        final_y_center = centroid_y
        
    except Exception as e:
        print(f"警告: 光度重心の計算に失敗しました。初期座標 ({initial_x:.3f}, {initial_y:.3f}) で測光を続行します。詳細: {e}")
        final_x_center = initial_x
        final_y_center = initial_y

    # 6. 求めた光度重心に対する測光の実行
    R_APERTURE = 12.0
    
    # 測光アパーチャを光度重心の座標に設定
    apertures = CircularAperture((final_x_center, final_y_center), r=R_APERTURE)
    
    # 分散データ内の負の値をゼロにクリッピングしてから平方根を取る (RuntimeWarning対策)
    safe_variance = np.maximum(variance_data, 0.0)
    
    # 測光の実行
    phot_table = aperture_photometry(image_data, apertures, error=np.sqrt(safe_variance))
    
    flux = phot_table['aperture_sum'][0]
    flux_err = phot_table['aperture_sum_err'][0]

    # 7. 等級と等級誤差の計算 (NameError対策済み)
    
    magnitude = np.nan
    magnitude_error = np.nan
    
    try:
        flux_mag_0 = calibration_info['fluxmag0']
        
        if flux > 0:
            magnitude = 2.5 * np.log10(flux_mag_0 / flux)
            magnitude_error = 1.0857 * (flux_err / flux)
        else:
            print("警告: 測光フラックスが非正のため、等級計算をスキップしました (出力は NaN)。")

    except Exception as e:
        print(f"警告: 等級計算中に予期せぬエラーが発生しました。詳細: {e}")
        
    # 8. 出力情報の準備
    observation_time_utc = calibration_info['time_obj_utc'].isot 
    filter_name = calibration_info['filter_name']
    
    output_ra = target_coord.ra.to_string(unit=u.hour, sep=':', precision=2, pad=True)
    output_dec = target_coord.dec.to_string(unit=u.degree, sep=':', precision=1, alwayssign=True)
    
    # 9. 結果のファイル出力 (★★★ 光度重心座標を追加 ★★★)
    output_filename = "photometry_results_centroid.txt" # ファイル名を変更して区別
    
    # 出力文字列の作成: 日時(UT), RA, Dec, X(初期), Y(初期), X(重心), Y(重心), フィルター, 等級, 等級誤差
    output_line = (
        f"{observation_time_utc}\t{output_ra}\t{output_dec}\t"
        f"{initial_x + 1 :.3f}\t{initial_y + 1:.3f}\t"
        f"{final_x_center + 1 :.3f}\t{final_y_center + 1 :.3f}\t" # ★光度重心を追加★
        f"{filter_name}\t{magnitude:.3f}\t{magnitude_error:.3f}\n"
    )

    file_exists = os.path.exists(output_filename)
    with open(output_filename, 'a') as f:
        if not file_exists:
            # ヘッダーを修正
            f.write("#日時(UT)\tRA(h:m:s)\tDec(d:m:s)\tX_initial\tY_initial\tX_centroid\tY_centroid\tフィルター\t等級(mag)\t等級誤差(mag_err)\n")
        f.write(output_line)

    print(f"\n測光が完了しました。")
    print(f"  初期座標 (X, Y): {initial_x:.3f}, {initial_y:.3f}")
    print(f"  光度重心 (X, Y): {final_x_center:.3f}, {final_y_center:.3f}")
    print(f"  測光等級: {magnitude:.3f}")
    print(f"  等級誤差: {magnitude_error:.3f}")
    print(f"結果は '{output_filename}' に出力されました。")
    
    hdul.close()

if __name__ == "__main__":
    main()