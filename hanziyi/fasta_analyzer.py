"""
1.读取FASTA文件（多行序列转换为单行序列）。
2.如果指定了输入文件，则从文件读取；否则从sys.stdin读取。
3.对每条序列，计算长度和GC含量。
4.根据给定的过滤条件（--min, --max, --gcmin, --gcmax）进行过滤。
5.将过滤后的序列信息（如ID，长度，GC含量）输出到stdout或指定的文件，同时将日志信息（如总序列数，过滤后的序列数）输出到stderr。
"""
import random
from typing import Dict, List, Tuple, Generator
import os
import sys
import argparse
from gc_analyzer import gc_analyzer
def read_fasta(input_stream)->Generator[Tuple[str,str],None,None]:
    """
    从输入流读取FASTA格式数据
    :param input_stream:输入流（文件对象或sys.stdin）
    :return:生成器，每次生成(header, sequence)元组
    """
    header = None
    sequence_lines = []
    for line in input_stream:
        line = line.strip()
        if not line:  #跳过空行
            continue
        if line.startswith('>'):
            if header is not None:
                #返回一个包含两个元素的元组 (header, sequence)
                yield header,''.join(sequence_lines)
                #重置用于下一个序列
                sequence_lines = []
            header = line[1:].split()[0]  #去掉‘>',取第一部分作为ID
        else:
            sequence_lines.append(line.upper())
    #输出最后一个序列
    if header is not None and sequence_lines:
        yield header,''.join(sequence_lines)

def analyze_fasta(input_stream, min_len=0, max_len=float('inf'),
                  gc_min=0, gc_max=100, log_stream=sys.stderr):
    """
    分析FASTA序列，应用过滤条件
    :param input_stream:输入流
    :param min_len:最小长度过滤
    :param max_len:最大长度过滤
    :param gc_min:最小GC含量过滤
    :param gc_max:最大GC含量过滤
    :param log_stream:日志输出流
    :return:
    filtered_sequences: 过滤后的序列列表[(header, sequence, length, gc_content)]
    stats: 统计信息字典
    """
    total_sequences = 0
    filtered_sequences = []
    log_stream.write(f"开始处理序列，过滤条件: "
                     f"长度 {min_len}-{max_len if max_len != float('inf') else 'inf'}, "
                     f"GC含量 {gc_min}%-{gc_max}%\n")
    for header,sequence in read_fasta(input_stream):
        total_sequences += 1
        #计算长度和GC含量
        length = len(sequence)
        gc_content = gc_analyzer(sequence)
        #应用过滤条件
        if min_len<=length<=max_len and gc_min<=gc_content<=gc_max:
            filtered_sequences.append((header, sequence, length, gc_content))
    # 计算统计信息
    stats = {
        'total_sequences': total_sequences,
        'filtered_sequences': len(filtered_sequences),
        'lengths': [seq[2] for seq in filtered_sequences],
        'gc_contents': [seq[3] for seq in filtered_sequences]
    }
    # 输出日志
    log_stream.write(f"处理完成: "
                     f"总序列数={total_sequences}, "
                     f"通过过滤={len(filtered_sequences)}\n")

    return filtered_sequences, stats

def print_sequence_distribution(sequences,output_stream):
    """
    输出序列长度分布信息
    参数:
        sequences: 序列列表[(header, sequence, length, gc_content)]
        output_stream: 输出流
    """
    if not sequences:
        output_stream.write("没有符合条件的序列\n")
        return
    # 提取长度和GC含量
    lengths = [seq[2] for seq in sequences]
    gc_contents = [seq[3] for seq in sequences]
    #计算统计信息
    max_len = max(lengths)
    min_len = min(lengths)
    avg_len = sum(lengths)/len(lengths)

    max_gc = max(gc_contents)
    min_gc = min(gc_contents)
    avg_gc = sum(gc_contents)/len(gc_contents)
    # 长度分布直方图数据
    output_stream.write("序列长度和GC含量分布统计\n")
    output_stream.write(f"序列总数: {len(sequences)}\n")
    output_stream.write(f"长度范围: {min_len} - {max_len} bp (平均: {avg_len:.1f} bp)\n")
    output_stream.write(f"GC含量范围: {min_gc:.1f}% - {max_gc:.1f}% (平均: {avg_gc:.2f}%)\n\n")

    # 长度分布直方图
    output_stream.write("长度分布:\n")
    #计算直方图
    if max_len-min_len>0:
        bin_size = max(1,(max_len-min_len)//10)
        bins = {}
        for length in lengths:
            bin_num = length//bin_size
            bins[bin_num] = bins.get(bin_num,0)+1
        # 输出直方图
        for bin_num in sorted(bins.keys()):
            start = bin_num * bin_size
            end = (bin_num + 1) * bin_size - 1
            count = bins[bin_num]
            bar = '#' * (count * 50 // len(sequences)) if sequences else ''
            output_stream.write(f"{start:6d}-{end:6d} bp: {count:4d} {bar}\n")
    output_stream.write("\n")

def output_results(sequences,output_stream,output_format='table'):
    """
    输出过滤序列
    :param sequences:序列列表
    :param output_stream:输出流
    :param output_format:输出格式 ('table' 或 'fasta')
    """
    if output_format == 'table':
        output_stream.write("\n过滤后的序列列表:\n")
        output_stream.write(f"{'ID':<30} {'长度':>8} {'GC含量':>8} {'N%':>6}\n")
        for header, sequence, length, gc_content in sequences:
            n_percent = sequence.count('N') / length * 100 if length > 0 else 0
            # 截断过长的header
            display_header = header[:27] + "..." if len(header) > 30 else header
            output_stream.write(f"{display_header:<30} {length:>8} {gc_content:>7.1f}% {n_percent:>5.1f}%\n")

    elif output_format == 'fasta':
        for header, sequence, length, gc_content in sequences:
            output_stream.write(f">{header} length={length} gc={gc_content:.2f}%\n")
            # FASTA格式每行通常不超过80个字符
            for i in range(0, len(sequence), 80):
                output_stream.write(sequence[i:i + 80] + "\n")

def main():
    parser = argparse.ArgumentParser(
        description='FASTA序列分析工具，支持长度和GC含量过滤',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=""
    )
    # 输入参数
    parser.add_argument('input', nargs='?', default='-',
                        help='输入FASTA文件（使用"-"或省略表示从标准输入读取）')

    # 过滤参数
    parser.add_argument('--min', type=int, default=0,
                        help='最小序列长度（默认: 0）')
    parser.add_argument('--max', type=int, default=0,
                        help='最大序列长度（0表示无限制，默认: 0）')
    parser.add_argument('--gcmin', type=float, default=0.0,
                        help='最小GC含量百分比（默认: 0.0）')
    parser.add_argument('--gcmax', type=float, default=100.0,
                        help='最大GC含量百分比（默认: 100.0）')

    # 输出参数
    parser.add_argument('--stdout', action='store_true',
                        help='将结果输出到标准输出（默认行为，与--out互斥）')
    parser.add_argument('--out', type=str,
                        help='将结果输出到指定文件（FASTA格式）')
    parser.add_argument('--format', choices=['table', 'fasta'], default='table',
                        help='输出格式：table（表格）或fasta（FASTA格式，默认: table）')
    parser.add_argument('--verbose', '-v', action='count', default=0,
                        help='详细输出级别（-v: 标准，-vv: 详细）')

    args = parser.parse_args()
    #设置日志级别
    log_level = args.verbose
    # 处理max参数：0表示无限制
    max_len = float('inf') if args.max == 0 else args.max
    # 检查输出参数冲突
    if args.out and args.stdout:
        sys.stderr.write("不能同时指定 --out 和 --stdout 参数\n")
        sys.exit(1)
    try:
        #确定输入源
        if args.input=='-' or args.input is None:
            sys.stderr.write("从标准输入读取数据\n")
            input_stream = sys.stdin
        else:
            if not os.path.exists(args.input):
                sys.stderr.write(f"输入文件不存在: {args.input}\n")
                sys.exit(1)
            sys.stderr.write(f"从文件读取数据: {args.input}\n")
            input_stream = open(args.input, 'r')
        # 分析FASTA序列
        filtered_sequences, stats = analyze_fasta(
            input_stream,
            min_len=args.min,
            max_len=max_len,
            gc_min=args.gcmin,
            gc_max=args.gcmax,
            log_stream=sys.stderr
        )
        # 关闭输入流（如果是文件）
        if input_stream is not sys.stdin:
            input_stream.close()
        # 确定输出流
        if args.out:
            output_stream = open(args.out, 'w')
            sys.stderr.write(f"结果将输出到文件: {args.out}\n")
        else:
            output_stream = sys.stdout
            sys.stderr.write("结果将输出到标准输出\n")
        #输出序列分布统计
        if log_level>=1:
            print_sequence_distribution(filtered_sequences,sys.stderr)
        #输出过滤后的序列
        output_results(filtered_sequences,output_stream,args.format)
        #关闭输出流
        if output_stream is not sys.stdout:
            output_stream.close()
            sys.stderr.write(f"结果已保存到: {args.out}\n")
        # 最终统计信息
        sys.stderr.write(f"\n分析完成。"
                         f"总序列数: {stats['total_sequences']}, "
                         f"通过过滤: {stats['filtered_sequences']}\n")
    except KeyboardInterrupt:
        sys.stderr.write("\n 用户中断操作\n")
        sys.exit(130)
    except Exception as e:
        sys.stderr.write(f"处理过程中发生错误: {str(e)}\n")
        sys.exit(1)


if __name__=="__main__":
    main()