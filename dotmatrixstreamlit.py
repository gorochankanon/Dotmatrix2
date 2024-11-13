import streamlit as st
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from io import StringIO


# Dot Matrixプロットを行う関数
def dotmatrix(f1, f2, win):
    # シーケンスの読み込み
    record1=next(SeqIO.parse(f1, "fasta"))
    record2=next(SeqIO.parse(f2, "fasta"))

    seq1=record1.seq
    seq2=record2.seq

    # 配列1と配列2の部分配列の長さ
    len1=len(seq1) - win + 1
    len2=len(seq2) - win + 1

    # 画像の幅と高さを設定
    width=500
    height=500

    image=np.zeros((height, width))  # 実際に描く行列

    # ハッシュテーブルを初期化
    hash={}

    # seq1の部分配列をハッシュに登録
    for x in range(len1):
        sub1=seq1[x:x + win]  # 部分配列を取得
        if sub1 not in hash:
            hash[sub1]=[]  # 初めての部分配列なら空リストを作成
        hash[sub1].append(x)  # その部分配列の位置xをリストに追加

    # seq2の部分配列を調べて行列に反映
    for y in range(len2):
        sub2=seq2[y:y + win]  # seq2の部分配列を取得
        py=int(y / len2 * height)  # yを画像の位置pyにスケーリング
        if sub2 in hash:
            for x in hash[sub2]:
                px=int(x / len1 * width)  # xを画像の位置pxにスケーリング
                image[py,px]=1  # [py, px]にドットを描画
    # 行列を画像として表示
    plt.imshow(image,extent=(1,len1,len2,1),cmap="Grays")
    #plt.show() 
    st.pyplot(plt) # 書 き 加える：Streamlit上にMatplotlibを表示

# タイトル
st.title("Dot Matrix")

# 配列ファイルのアップローダー
file1=st.sidebar.file_uploader("Sequence file 1:")
file2=st.sidebar.file_uploader("Sequence file 2:")
win=st.sidebar.slider("Window size:",4,100,10)  # ウィンドウサイズスライダー
if file1 and file2: #2つのファイルがアップロードされていれば
    with StringIO(file1.getvalue().decode("utf-8")) as f1,\
        StringIO(file2.getvalue().decode("utf-8")) as f2:
        dotmatrix(f1,f2,win)
