import os
import numpy as np
import matplotlib.pyplot as plt

# -------- 提取函数（保持不变） -------- #
def extract_brate_bler_from_file(filename):
    brate_list = []
    bler_list = []
    with open(filename, encoding="utf-8") as f:
        lines = f.readlines()
    in_dl = False
    brate_idx = None
    bler_idx = None
    for line in lines:
        line = line.strip()
        if "---DL---" in line:
            in_dl = True
            brate_idx = None
            bler_idx = None
            continue
        if "---UL---" in line:
            in_dl = False
            continue
        if in_dl:
            if "brate" in line and "(%)" in line:
                headers = line.split()
                for i, h in enumerate(headers):
                    if h == "brate" and brate_idx is None:
                        brate_idx = i
                    if h == "(%)" and bler_idx is None:
                        bler_idx = i
                continue
            if brate_idx is not None and bler_idx is not None:
                cols = line.split()
                if len(cols) > max(brate_idx, bler_idx):
                    brate_raw = cols[brate_idx]
                    bler_raw = cols[bler_idx]
                    # brate 统一 Mbps
                    if brate_raw.lower().endswith("m"):
                        try:
                            brate = float(brate_raw[:-1])
                        except:
                            brate = None
                    elif brate_raw.lower().endswith("k"):
                        try:
                            brate = float(brate_raw[:-1]) / 1000.0
                        except:
                            brate = None
                    else:
                        try:
                            brate = float(brate_raw)
                        except:
                            brate = None
                    if brate is not None:
                        brate_list.append(brate)
                    # bler 百分比转小数
                    if "%" in bler_raw:
                        try:
                            bler = float(bler_raw.replace("%", "")) / 100.0
                            bler_list.append(bler)
                        except:
                            pass
    return np.array(brate_list), np.array(bler_list)

# -------- 读取文件 -------- #
file_groups = {
    "MCS8": ["with.txt", "without.txt"],
    "MCS10": ["10_with.txt", "10_without.txt"]
}

# 提取数据
group_labels = []
bler_with = []
bler_without = []
brate_with = []
brate_without = []

for mcs, files in file_groups.items():
    group_labels.append(mcs)
    for f in files:
        if os.path.exists(f):
            brate, bler = extract_brate_bler_from_file(f)
            print(f"{f}: brate={len(brate)}, bler={len(bler)}")
            if f.endswith("with.txt") or f.endswith("10_with.txt"):
    # 是 with predistortion
             bler_with.append(bler)
             brate_with.append(brate)
            else:
             bler_without.append(bler)
             brate_without.append(brate)

# -------- 绘图配置（论文风格）-------- #
plt.rcParams.update({
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 11,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 10,
    'figure.dpi': 300,
    'figure.figsize': (6, 3),
    'axes.grid': True
})

# 横坐标位置
positions_with = [x - 0.15 for x in range(1, len(group_labels)+1)]
positions_without = [x + 0.15 for x in range(1, len(group_labels)+1)]

# -------- 画 BLER -------- #
plt.figure()
b1 = plt.boxplot(bler_with, positions=positions_with, widths=0.25,
                 patch_artist=True, boxprops=dict(facecolor='skyblue'),
                 medianprops=dict(color='black'))

b2 = plt.boxplot(bler_without, positions=positions_without, widths=0.25,
                 patch_artist=True, boxprops=dict(facecolor='salmon'),
                 medianprops=dict(color='black'))

plt.xticks(ticks=range(1, len(group_labels)+1), labels=group_labels, fontsize=15)
plt.ylabel("BLER", fontsize=15)
plt.title("BLER grouped by MCS", fontsize=15)
plt.legend([b1["boxes"][0], b2["boxes"][0]],
           ["With Predistortion", "Without Predistortion"],
           loc="upper right", fontsize=10,prop={'weight': 'bold'})
plt.tight_layout()
plt.tick_params(axis='y', labelsize=11)

plt.savefig("bler_grouped_boxplot.pdf")
plt.show()

# -------- 画 Bitrate -------- #
plt.figure()
b1 = plt.boxplot(brate_with, positions=positions_with, widths=0.25,
                 patch_artist=True, boxprops=dict(facecolor='skyblue'),
                 medianprops=dict(color='black'))

b2 = plt.boxplot(brate_without, positions=positions_without, widths=0.25,
                 patch_artist=True, boxprops=dict(facecolor='salmon'),
                 medianprops=dict(color='black'))

plt.xticks(ticks=range(1, len(group_labels)+1), labels=group_labels, fontsize=15)
plt.ylabel("Throughput (Mbps)", fontsize=15)
plt.title("Throughput grouped by MCS", fontsize=15)
plt.legend([b1["boxes"][0], b2["boxes"][0]],
           ["With Predistortion", "Without Predistortion"], fontsize=10,prop={'weight': 'bold'})
plt.ylim(4, 9)
plt.tick_params(axis='y', labelsize=15)

plt.tight_layout()
plt.savefig("brate_grouped_boxplot.pdf")
plt.show()
