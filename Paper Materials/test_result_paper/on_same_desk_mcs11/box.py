import os
import matplotlib.pyplot as plt
import numpy as np
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42


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
                    if "%" in bler_raw:
                        try:
                            bler = float(bler_raw.replace("%", "")) / 100.0
                            bler_list.append(bler)
                        except:
                            pass
    return np.array(brate_list), np.array(bler_list)

# 配置文件组
file_groups = {
    "MCS11": ["12_with.txt", "12_without.txt"],
    "MCS13": ["13_with.txt", "13_without.txt"],
    "MCS15": ["14_with.txt", "14_without.txt"]
}

# 设置风格
plt.rcParams.update({
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 11,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 10,
    'figure.dpi': 300,
    'figure.figsize': (6, 4.5),
    'axes.grid': True
})

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
            if "without" in os.path.basename(f):  # 精确判断
                bler_without.append(bler)
                brate_without.append(brate)
            else:
                bler_with.append(bler)
                brate_with.append(brate)

# 设置横坐标位置
positions_with = [x - 0.15 for x in range(1, len(group_labels)+1)]
positions_without = [x + 0.15 for x in range(1, len(group_labels)+1)]
plt.rcParams.update({'font.size': 22})

# ✅ BLER grouped boxplot
plt.figure(figsize=(7, 4))

b1 = plt.boxplot(bler_with, positions=positions_with, widths=0.25,
                 patch_artist=True, boxprops=dict(facecolor='skyblue'),
                 medianprops=dict(color='black'))
b2 = plt.boxplot(bler_without, positions=positions_without, widths=0.25,
                 patch_artist=True, boxprops=dict(facecolor='salmon'),
                 medianprops=dict(color='black'))

plt.xticks(ticks=range(1, len(group_labels)+1), labels=group_labels, fontsize=18)
plt.ylabel("BLER", fontsize=18)
plt.ylim(0, 0.5)  # 放大 y 轴
#plt.title("BLER grouped by MCS", fontsize=15)
plt.legend([b1["boxes"][0], b2["boxes"][0]],
           ["With SAPE", "Without SAPE"],
           loc="upper right", prop={'size': 15, 'weight': 'bold'})
plt.tick_params(axis='y', labelsize=18)
plt.tight_layout()
plt.savefig("bler_boxplot.pdf")
#plt.show()

plt.figure(figsize=(7, 4))
b1 = plt.boxplot(brate_with, positions=positions_with, widths=0.25,
                 patch_artist=True, boxprops=dict(facecolor='skyblue'),
                 medianprops=dict(color='black'))
b2 = plt.boxplot(brate_without, positions=positions_without, widths=0.25,
                 patch_artist=True, boxprops=dict(facecolor='salmon'),
                 medianprops=dict(color='black'))

plt.xticks(ticks=range(1, len(group_labels)+1), labels=group_labels, fontsize=18)
plt.ylabel("Throughput (Mbps)", fontsize=18)
#plt.ylim(4, 9)  # 放大 y 轴
#plt.title("Throughput grouped by MCS", fontsize=15)
plt.legend([b1["boxes"][0], b2["boxes"][0]],
           ["With SAPE", "Without SAPE"],
           loc="upper left", prop={'size': 14, 'weight': 'bold'})
plt.tick_params(axis='y', labelsize=15)

plt.tight_layout()
plt.savefig("brate_boxplot.pdf")
#plt.show()
