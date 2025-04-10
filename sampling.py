import numpy as np
import json
import argparse
import matplotlib.pyplot as plt
from scipy.stats.qmc import LatinHypercube
from scipy.spatial.distance import pdist

def main():
    parser = argparse.ArgumentParser(description="LHS 采样")
    parser.add_argument("--num_samples", type=int, required=True, help="采样数量")
    parser.add_argument("--atoms_type", type=int, required=True, help="原子种类数")
    parser.add_argument("--DIS_min", type=float, required=True, help="DIS 最小值")
    parser.add_argument("--DIS_max", type=float, required=True, help="DIS 最大值")
    parser.add_argument("--step", type=float, required=True, help="DIS 采样步长")
    parser.add_argument("--volume_min", type=float, required=False, help="体积最小值", default=100)
    parser.add_argument("--volume_max", type=float, required=False, help="体积最大值", default=300)
    parser.add_argument("--area_min", type=float, required=False, help="面积最小值", default=50)
    parser.add_argument("--area_max", type=float, required=False, help="面积最大值", default=200)
    parser.add_argument("--height_min", type=float, required=False, help="层高最小值", default=14)
    parser.add_argument("--height_max", type=float, required=False, help="层高最大值", default=20)

    args = parser.parse_args()

    num_DIS_variables = args.atoms_type ** 2  # 计算 DIS 维度
    num_total_variables = num_DIS_variables + 3  # DIS + 体积 + 面积 + 层高

    DIS_possible = np.arange(args.DIS_min, args.DIS_max + args.step, args.step)

    # 生成 LHS 采样
    sampler = LatinHypercube(d=num_total_variables)
    samples = sampler.random(n=args.num_samples)

    # 变换到实际数值范围
    indices = (samples[:, :num_DIS_variables] * len(DIS_possible)).astype(int)
    indices = np.clip(indices, 0, len(DIS_possible) - 1)
    DIS_samples = DIS_possible[indices]

    volume_samples = np.round(args.volume_min + (args.volume_max - args.volume_min) * samples[:, num_DIS_variables]).astype(int)
    area_samples = np.round(args.area_min + (args.area_max - args.area_min) * samples[:, num_DIS_variables + 1]).astype(int)
    height_samples = np.round(args.height_min + (args.height_max - args.height_min) * samples[:, num_DIS_variables + 2]).astype(int)

    all_samples = np.hstack((DIS_samples, volume_samples[:, None], area_samples[:, None], height_samples[:, None]))

    # 保存 JSON 文件
    samples_list = np.around(all_samples, 2).tolist()
    with open("LHS_samples.json", "w") as f:
        json.dump(samples_list, f, indent=4)

    # 评估采样质量
    #evaluate_sampling(all_samples, num_DIS_variables)

def evaluate_sampling(samples, num_DIS_variables):
    """ 评估 LHS 采样质量 """
    print("\n=== 采样质量评估 ===")
    print("均值:", np.mean(samples, axis=0))
    print("标准差:", np.std(samples, axis=0))

    # 计算最近邻距离
    distances = pdist(samples)
    min_dist = np.min(distances)
    print(f"最小数据点距离: {min_dist:.4f}")

    # 绘制直方图
    fig, axes = plt.subplots(1, 4, figsize=(16, 4))

    # DIS 变量（仅绘制第一个）
    axes[0].hist(samples[:, 0], bins=20, alpha=0.7, color='b')
    axes[0].set_title("DIS 分布")

    # 体积
    axes[1].hist(samples[:, num_DIS_variables], bins=20, alpha=0.7, color='r')
    axes[1].set_title("体积分布")

    # 面积
    axes[2].hist(samples[:, num_DIS_variables + 1], bins=20, alpha=0.7, color='g')
    axes[2].set_title("面积分布")

    # 层高
    axes[3].hist(samples[:, num_DIS_variables + 2], bins=10, alpha=0.7, color='orange')
    axes[3].set_title("层高分布")

    plt.tight_layout()
    plt.show()
    plt.save('./sampling.png')
if __name__ == "__main__":
    main()

