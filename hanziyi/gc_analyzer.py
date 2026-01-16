"""
计算DNA序列中G和C碱基的比例
1.将序列转化为大写
2.检查序列中是否包含非AGCT字符，如果包含，则删除
3.计算GC含量，即(G+C)的长度除以总长度
4.使用try...except来处理异常
"""
def gc_analyzer(seq):
    try:
        #转化为大写
        seq = seq.upper().strip()
        if not seq:
            return "错误：序列为空"
        #检查是否只包含AGCT
        valid_chars = set("ATCG")
        for char in seq:
            if char not in valid_chars:
                raise ValueError(f"序列包含非法字符: '{char}'")
        #计算GC含量
        g_count = seq.count('G')
        c_count = seq.count('C')
        total = len(seq)
        if total == 0:
            return 0.0
        gc_percentage = (g_count+c_count)/total * 100
        return gc_percentage
    except ValueError as e:
        return f"输入序列错误: {e}"
    except Exception as e:
        return f"发生未知错误: {e}"

#测试用例
# if __name__ == "__main__":
#     # 测试用例
#     test_sequences = [
#         "ATCGATCG",  # 正常序列，GC含量50%
#         "AAAAA",  # 正常序列，GC含量0%
#         "GGGGCCCC",  # 正常序列，GC含量100%
#         "ATCGN",  # 包含非法字符
#         "atcgatcg",  # 小写字母
#         "",  # 空序列
#         "A T C G",  # 包含空格
#     ]
#
#     for seq in test_sequences:
#         result = gc_analyzer(seq)
#         print(f"序列: '{seq}'")
#         print(f"GC含量: {result}")


"""
统计k-mer频率
1.检查序列长度和k,确保k≤序列长度
2.遍历序列，提取每个长度为k的子串
3.使用字典来统计每个字符串出现的次数
4.处理异常
"""

def kmers_count(seq,k):
    try:
        # 转化为大写
        seq = seq.upper().strip()
        # 序列是否为空
        if not seq:
            return {"错误": "序列为空"}
        if k <= 0:
            return {"错误": f"k值必须为正整数"}
        if k > len(seq):
            return {"错误": f"k值({k})大于序列长度({len(seq)})"}
        # 检查是否只包含ATCG字符
        valid_chars = set('ATCG')
        for char in seq:
            if char not in valid_chars:
                raise ValueError(f"序列包含非法字符: '{char}'")

        kmer_counts = {}
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i + k]
            kmer_counts[kmer] = kmer_counts.get(kmer, 0) + 1
        # 转化为频率
        total_kmers = len(seq) - k + 1
        kmer_freq = {kmer: format(count / total_kmers,".4f") for kmer, count in kmer_counts.items()}
        return kmer_freq
    except ValueError as e:
        return {"错误": f"输入序列错误: {e}"}
    except Exception as e:
        return {"错误": f"发生未知错误: {e}"}


# 测试k-mer统计
if __name__ == "__main__":
    # 测试用例
    print("测试k-mer频率统计:")
    test_cases = [
        ("ATCGATCGATCG", 3),  # 正常序列，k=3
        ("AAAAA", 2),  # 重复序列
        ("ATCG", 5),  # k大于序列长度
        ("ATCGN", 2),  # 包含非法字符
        ("", 2),  # 空序列
        ("ATCGATCG", -1),  # k为负值
    ]

    for seq, k in test_cases:
        result = kmers_count(seq, k)
        print(f"序列: '{seq}'")
        print(f"k-mer频率k={k}: {result}")

    # 更复杂的例子
    complex_seq = "ATCGATCGATCGATCG"
    for k in [2, 3, 4]:
        result = kmers_count(complex_seq, k)
        print(f"序列: '{complex_seq}'")
        print(f"k-mer频率k={k}: {result}")
