# FASTA Tools

一个用于FASTA序列文件处理和分析的Python工具包。

## 功能特性

- 支持从文件或标准输入读取FASTA格式序列
- 序列长度过滤（最小/最大长度）
- GC含量过滤（最小/最大GC含量）
- 序列质量验证（非法字符检测）
- 多种输出格式（表格、FASTA、统计信息）
- 详细的统计信息输出
- 支持管道操作
- 完整的日志系统

## 安装

### 从源码安装

```bash
git clone https://github.com/Muriel-bit/BestPractice4Bioinformatics.git
cd feature/gc-analyzer
pip install -e 

python generate_test_fasta.py | python fasta_analyzer.py --min 100 
--max 150 --gcmin 30 --gcmax 70 -vv --stdout
```