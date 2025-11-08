import os
import matplotlib.pyplot as plt
import numpy as np

# Provided function to extract brate and bler from file
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
                    # brate 统一Mbps
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
                    # bler 百分数转小数
                    if "%" in bler_raw:
                        try:
                            bler = float(bler_raw.replace("%", "")) / 100.0
                            bler_list.append(bler)
                        except:
                            pass
    return np.array(brate_list), np.array(bler_list)

# File mapping: MCS values and "with"/"without"
file_groups = {
    "11": ["12_with.txt", "12_without.txt"],
    "13": ["13_with.txt", "13_without.txt"],
    "15": ["14_with.txt", "14_without.txt"]
}
# 设置论文风格全局字体
plt.rcParams.update({
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 11,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 10,
    'figure.dpi': 300,
    'figure.figsize': (6, 3),  # ✅ 调整图尺寸
    'axes.grid': True
})

# ✅ 修改你的横轴标签
custom_labels = [
    "MCS11_with_predistortion",
    "MCS11_without_predistortion",
    "MCS13_with_predistortion",
    "MCS13_without_predistortion",
    "MCS15_with_predistortion",
    "MCS15_without_predistortion"
]
# ✅ 新的 lists，用于存放成功提取的数据和对应标签
valid_labels = []
valid_bler_data = []
valid_brate_data = []

for mcs, files in file_groups.items():
    for f in files:
        if os.path.exists(f):
            brate, bler = extract_brate_bler_from_file(f)
            label = f"{mcs}_{'with' if f.endswith('_with.txt') else 'without'}"
            print(f"File: {f}, Label: {label}, brate: {len(brate)}, bler: {len(bler)}")
            if len(brate) > 0 and len(bler) > 0:
                valid_labels.append(label)
                valid_brate_data.append(brate)
                valid_bler_data.append(bler)
            else:
                print(f"⚠️ 跳过空数据文件: {f}")
# ⚠️ 确保数据和标签数量一致
assert len(valid_bler_data) == len(custom_labels), "标签数量和数据不一致！"

# ✅ BLER Boxplot
plt.figure()
plt.boxplot(valid_bler_data, patch_artist=True)
plt.xticks(ticks=range(1, len(custom_labels)+1), labels=custom_labels, rotation=30, ha='right')
plt.ylabel("BLER")
plt.title("Block Error Rate across MCS and Predistortion Configurations")
plt.tight_layout()
plt.savefig("bler_boxplot.pdf")  # 可换成 .png
plt.show()

# ✅ Bitrate Boxplot
plt.figure()
plt.boxplot(valid_brate_data, patch_artist=True)
plt.xticks(ticks=range(1, len(custom_labels)+1), labels=custom_labels, rotation=30, ha='right')
plt.ylabel("Throughput (Mbps)")
plt.title("Throughput across MCS and Predistortion Configurations")
plt.tight_layout()
plt.savefig("brate_boxplot.pdf")
plt.show()

