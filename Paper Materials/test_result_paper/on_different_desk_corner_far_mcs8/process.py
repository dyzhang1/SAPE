import matplotlib.pyplot as plt
import numpy as np

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

def plot_cdf(data, label):
    sorted_data = np.sort(data)
    cdf = np.arange(1, len(sorted_data)+1) / len(sorted_data)
    plt.plot(sorted_data, cdf, label=label)

# 读取两组文件
brate_with, bler_with = extract_brate_bler_from_file("with_mcs10")
brate_without, bler_without = extract_brate_bler_from_file("without_mcs10")

# 画brate CDF
plt.figure(figsize=(7, 5))
plot_cdf(brate_with, "With Predistortion")
plot_cdf(brate_without, "Without Predistortion")
plt.xlabel("DL Throughput (Mbps)")
plt.ylabel("CDF")
plt.title("CDF of DL Throughput (brate)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

# 画block error rate CDF
plt.figure(figsize=(7, 5))
plot_cdf(bler_with, "With Predistortion")
plot_cdf(bler_without, "Without Predistortion")
plt.xlabel("Block Error Rate (BLER)")
plt.ylabel("CDF")
plt.title("CDF of Block Error Rate (BLER)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

