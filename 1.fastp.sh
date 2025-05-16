#!/bin/bash
#SBATCH --job-name=fastp         # 任务名称
#SBATCH --output=fastp.log       # 标准输出文件
#SBATCH --error=fastp.error.txt         # 错误输出文件
#SBATCH --ntasks=1                # 任务数
#SBATCH --cpus-per-task=10         # 每个任务的 CPU 核心数
#SBATCH --mem=30G            # 设置任务总内存为 30 GB
#SBATCH --time=100:00:00           # 运行时间上限（HH:MM:SS）此处100小时
#SBATCH --partition=student        # 分区名称，指定teacher分区或者student分区，不同分区有不同的节点

wd=/osmgfs10000/home/wub/RNA-seq/  #该项目工作区设置
cd /osmgfs10000/home/wub/RNA-seq/   #移动至工作区
db=/osmgfs10000/home/wub/db/         #参考数据库位置
soft=/osmgfs10000/home/wub/conda/    #软件位置
source activate   #激活基本环境
cd ${wd}           #再次移动至工作区，这一步是确认的步骤，可不加

conda activate fastp  #激活软件
mkdir -p temp/QC result/QC


    time tail -n+2 result/mydata.txt|cut -f1|rush -j 4 \
      "fastp -i Seq/{}_1.fq.gz -I Seq/{}_2.fq.gz \
        -j temp/QC/{}_fastp.json -h temp/QC/{}_fastp.html \
        -o temp/QC/{}_1.fastq  -O temp/QC/{}_2.fastq \
        > temp/QC/{}.log 2>&1"

#{}_1.fq.gz是输入需要质控的文件，{}是mydata中第一列的样本名称，_1.fq.gz是样本序列文件的其他部分，不同公司这一部分有区别，根据实际情况选择
    
	# 质控后结果汇总
    echo -e "SampleID\tRaw\tClean" > temp/fastp
    for i in `tail -n+2 result/mydata.txt|cut -f1`;do
        echo -e -n "$i\t" >> temp/fastp
        grep 'total reads' temp/QC/${i}.log|uniq|cut -f2 -d ':'|tr '\n' '\t' >> temp/fastp
        echo "" >> temp/fastp
        done
    sed -i 's/ //g;s/\t$//' temp/fastp
	
	 # 按mydata排序
    awk 'BEGIN{FS=OFS="\t"}NR==FNR{a[$1]=$0}NR>FNR{print a[$1]}' temp/fastp result/mydata.txt \
      > result/QC/fastp.txt
    cat result/QC/fastp.txt
