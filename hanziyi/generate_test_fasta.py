
import sys
import random

def generate_test_fasta(num_sequences=10,min_len=50,max_len=200):
    """生成测试fasta数据"""
    bases = ['A','G','C','T']
    for i in range(num_sequences):
        #输出header
        sys.stdout.write(f'>seq{i}_len_{random.randint(min_len, max_len)}\n')
        #生成随机序列
        seq_length = random.randint(min_len,max_len)  #随机选一个长度
        seq = ''.join(random.choices(bases,k=seq_length))
        # 按FASTA格式每行80字符输出
        for j in range(0, len(seq), 80):
            sys.stdout.write(seq[j:j + 80] + '\n')
if __name__ == '__main__':
    generate_test_fasta()