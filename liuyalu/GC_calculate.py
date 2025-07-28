"""
计算 FASTA 文件中所有序列的 GC 含量GC% = (G + C) / (A + T + G + C) * 100
用法如下：
在terminal中输入：python gc_content.py dna_example.fasta
"""
# 导入pathlib中的Path类能够获取fasta文件的路径
import sys
from pathlib import Path

def calculate_gc(fasta_path):
    # 定义一个空字符串，用于存储读取的序列
    seq = ""
    # 打开FASTA文件
    with open(fasta_path) as fh:
        # 遍历文件中的每一行
        for line in fh:
            # 如果行以">"开头，则跳过该行
            if line.startswith(">"):
                continue
            # 将行中的空格去除，并将序列转换为大写，然后添加到seq字符串中
            seq += line.strip().upper()

    # 统计序列中G、C、A、T的数量
    g = seq.count("G")
    c = seq.count("C")
    a = seq.count("A")
    t = seq.count("T")
    # 计算序列中总的碱基数
    total = a + t + g + c

    # 如果序列中没有检测到有效碱基序列，则抛出异常
    if total == 0:
        raise ValueError("FASTA 文件中没有检测到有效碱基序列！")

    # 计算序列中G和C的比例
    gc_ratio = (g + c) / total
    # 将比例转换为百分比
    gc_percent = gc_ratio * 100
    # 返回GC百分比和G、C、A、T的数量以及总的碱基数
    return gc_percent, g, c, a, t, total

def main():
    # 检查命令行参数个数是否为2，如果不是，则打印用法并退出程序
    if len(sys.argv) != 2:
        print("用法: python gc_content.py <fasta_file>")
        sys.exit(1)

    # 获取命令行参数中的fasta文件路径
    # sys.argv[1]是命令行输入的第一个参数，然后将fasta文件路径转换为Path对象
    fasta_file = Path(sys.argv[1])
    # 检查文件是否存在，如果不存在，则打印错误信息并退出程序
    if not fasta_file.exists():
        print(f"错误：文件 {fasta_file} 不存在")
        sys.exit(1)

    # 调用calculate_gc函数计算fasta文件中的GC含量
    gc, g, c, a, t, total = calculate_gc(fasta_file)
    # 打印文件路径
    print(f"文件: {fasta_file}")
    # 打印A、T、G、C的碱基数和总碱基数
    print(f"A: {a}, T: {t}, G: {g}, C: {c}, 总碱基数: {total}")
    # 打印GC含量，保留两位小数
    print(f"GC含量: {gc:.2f}%")

# 运行主函数
if __name__ == "__main__":
    main()