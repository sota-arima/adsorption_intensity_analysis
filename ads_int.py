import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import MDAnalysis as mda
import os

#chane working directory
current_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(current_dir)
print(os.getcwd())

#解析データの読み込み
data = mda.Universe('Au_c-PEG2k_200ns_ads_int_no_centering.gro', 'Au_c-PEG2k_200ns_ads_int_no_centering.xtc')
surface_bottom = float(input('金表面のz座標(下側)を入力 (Å)'))
surface_upper = float(input('金表面のz座標(上側)を入力 (Å)'))

#PEGのO原子を選択
O_atoms = data.select_atoms('resname POL and name O*')

#吸着判定距離
z_threshold = float(input('吸着とみなす距離を入力 (Å)'))

#フレーム数を取得
num_frames = len(data.trajectory)
print(num_frames)

#Ads_Intのリスト
Ads_Int_list = []
Time_list = []

print(len(O_atoms))
#position = data.atoms.positions
#print(position)
#最初のフレームに移動
for n in range(num_frames):
    data.trajectory[n]
    
    #z座標の抽出
    O_z = O_atoms.positions[:,2]

    #PEGのO原子と金表面とのz座標距離を計算
    z_distances_bottom = np.abs(O_z - surface_bottom)
    z_distances_upper = np.abs(O_z - surface_upper)

    #吸着しているO原子の数
    ads_O = np.sum((z_distances_bottom < z_threshold) | (z_distances_upper < z_threshold))
    total_O = len(O_z)

    #吸着強度の計算
    ads_int = ads_O/total_O
    Ads_Int_list.append(ads_int)
    Time_list.append(n)
    #print(f'吸着強度{ads_int}')
    
#吸着強度のデータフレームの読み込みと確認
Ads_Int_df = pd.DataFrame({
    'adsorption_intensity':Ads_Int_list
})

print(Ads_Int_df.head())
print(Ads_Int_df.tail())
print(Ads_Int_df.columns)

#CSV出力
output_flie = os.path.join(current_dir,'ads_int_Au_c-PEG2k.csv')
Ads_Int_df.to_csv(output_flie,index=False)
print(f'解析結果を{output_flie}に保存しました')

#吸着強度の時間変化のプロット
plt.plot(Ads_Int_df.index/1000, Ads_Int_df['adsorption_intensity'],color='red')
plt.xlabel('Time (ns)')
plt.ylabel('Adsorption Intensity')
plt.ylim(-0.05,1.0)
plt.show()

#吸着強度の分布のプロット
bins =30
hist,bin_edges = np.histogram(Ads_Int_df['adsorption_intensity'],bins=bins)
bin_width = bin_edges[1]-bin_edges[0]
bin_centers = (bin_edges[:-1]+bin_edges[1:])/2
plt.bar(bin_centers,hist,color='red',width=bin_width,edgecolor='black')
plt.xlabel('Adsorption Intensity')
plt.ylabel('Count')
plt.xlim(-0.05,1.0)
plt.show()

#吸着確率の計算
ads_frame_num = (Ads_Int_df['adsorption_intensity'] > 0.25).sum()
print(f'吸着フレーム数{ads_frame_num}')
print(f'全フレーム数{num_frames}')
print(f'吸着確率{ads_frame_num/num_frames}')
