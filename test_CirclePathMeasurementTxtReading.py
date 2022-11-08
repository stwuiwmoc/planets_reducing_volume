# %%

import importlib

import matplotlib.pyplot as plt
import numpy as np

import planets_optimize_myclass as pom

if __name__ == "__main__":
    importlib.reload(pom)

    CONSTS = pom.Constants(
        physical_radius=925e-3,
        ignore_radius=25e-3,
        pixel_number=1024,
        zernike_max_degree=10,
        offset_height_percent=2)

    original_filepath = "raw_data/220117xrmEAi.v5.60.hei.txt"
    deformed_filepath = "raw_data/220117xrmFEi.v5.60.hei.txt"

    diff = pom.CirclePathMeasurementTxtReading(
        Constants=CONSTS,
        original_txt_fpath=original_filepath,
        deformed_txt_fpath=deformed_filepath
    )

    # 変形前後の差分のプロット
    fig1 = plt.figure()
    gs1 = fig1.add_gridspec(1, 1)

    ax11 = fig1.add_subplot(gs1[0, 0])
    ax11.plot(diff.df_diff["degree"], diff.df_diff["height"])

    ax11.set_xlabel("degree")
    ax11.set_ylabel("height [m]")
    ax11.grid()

    fig1.tight_layout()

    # 変形前、後それぞれの高さと角度のプロット
    fig2 = plt.figure(figsize=(10, 10))
    gs2 = fig2.add_gridspec(2, 1)

    ax21 = fig2.add_subplot(gs2[0, 0])
    ax21.plot(
        diff.df_raw_original["degree"].iloc[:30],
        diff.df_raw_original["height"].iloc[:30],
        color="blue")
    ax21.plot(
        diff.df_raw_original["degree"].iloc[-30:],
        diff.df_raw_original["height"].iloc[-30:],
        color="red")
    ax21.grid()
    ax21.set_xlim(-100, -50)

    ax21 = fig2.add_subplot(gs2[1, 0])
    ax21.plot(
        diff.df_raw_deformed["degree"].iloc[:30],
        diff.df_raw_deformed["height"].iloc[:30],
        color="blue")
    ax21.plot(
        diff.df_raw_deformed["degree"].iloc[-30:],
        diff.df_raw_deformed["height"].iloc[-30:],
        color="red")
    ax21.grid()
    ax21.set_xlim(-100, -50)

    # 変形前の円環パスx, y座標のプロット
    # 色はデータのインデックス小さい順（高さではない）
    fig3 = plt.figure(figsize=(10, 10))
    fig3.suptitle(original_filepath)
    gs3 = fig3.add_gridspec(2, 2)

    ax31 = fig3.add_subplot(gs3[0, 0])
    ax31.scatter(
        diff.df_raw_original["x"],
        diff.df_raw_original["y"],
        c=np.arange(len(diff.df_raw_original["x"])))
    ax31.set_xlabel("x [m]")
    ax31.set_ylabel("y [m]")

    # 円環パスの開始位置と終了位置が被っている所を拡大
    ax32 = fig3.add_subplot(gs3[0, 1])
    ax32.scatter(
        diff.df_raw_original["x"],
        diff.df_raw_original["y"],
        c=np.arange(len(diff.df_raw_original["x"])))
    ax32.set_xlim(0.15, 0.25)
    ax32.set_ylim(-0.9, -0.8)
    ax32.grid()
    ax32.set_xlabel("x [m]")
    ax32.set_ylabel("y [m]")

    # ax32でプロットした測定点のx座標とz座標で棒グラフ作成
    bar_width = 0.001
    ax33 = fig3.add_subplot(gs3[1, 1])

    # 円環パスの開始位置付近
    # 座標がほぼ同じ点だと棒グラフが完全に被ってしまうので棒の太さ分だけずらす
    ax33.bar(
        x=diff.df_raw_original["x"].iloc[:5] + bar_width,
        height=diff.df_raw_original["height"].iloc[:5],
        width=bar_width,
        color="purple",
        label="start_side")

    # 円環パスの終了位置付近
    ax33.bar(
        x=diff.df_raw_original["x"].iloc[-5:],
        height=diff.df_raw_original["height"].iloc[-5:],
        width=bar_width,
        color="yellow",
        label="end_side")

    ax33.set_xlim(0.15, 0.25)
    ax33.legend()
    ax33.set_xlabel("x [m]")
    ax33.set_ylabel("z [m]")

    fig3.tight_layout()
