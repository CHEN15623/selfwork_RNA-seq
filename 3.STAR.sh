#!/bin/bash
#SBATCH --job-name=STAR         # 任务名称
#SBATCH --output=STAR.log       # 标准输出文件
#SBATCH --error=STAR.error.txt         # 错误输出文件
#SBATCH --ntasks=1                # 任务数
#SBATCH --cpus-per-task=30         # 每个任务的 CPU 核心数
#SBATCH --mem=50G            # 设置任务总内存为 30 GB
#SBATCH --time=100:00:00           # 运行时间上限（HH:MM:SS）此处100小时
#SBATCH --partition=student        # 分区名称，指定teacher分区或者student分区，不同分区有不同的节点

wd=/osmgfs10000/home/wub/RNA-seq/  #该项目工作区设置
cd /osmgfs10000/home/wub/RNA-seq/   #移动至工作区
db=/osmgfs10000/home/wub/db/         #参考数据库位置
soft=/osmgfs10000/home/wub/conda/    #软件位置
source activate   #激活基本环境
cd ${wd}           #再次移动至工作区，这一步是确认的步骤，可不加
conda activate STAR



mkdir -p temp/STAR

tail -n +2 result/mydata.txt | cut -f1 | rush -j 2  \
"STAR  --runThreadN 10 \
--genomeDir /osmgfs10000/home/wub/db/STAR/GRCm39_index \
--readFilesIn temp/QC/{}_1.fastq temp/QC/{}_2.fastq \
--outFileNamePrefix temp/STAR/{1}_data \
--outSAMtype BAM SortedByCoordinate \
> temp/STAR/{}.log 2>&1"