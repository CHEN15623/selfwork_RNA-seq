#!/bin/bash
#SBATCH --job-name=featureCounts        # 任务名称
#SBATCH --output=featureCounts.log       # 标准输出文件
#SBATCH --error=featureCounts.error.txt         # 错误输出文件
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
conda activate subread



mkdir -p result/featureCounts

# 运行计数
featureCounts -T 20 -p -t exon -g gene_id \
-a ${db}/STAR/Mus_musculus.GRCm39.113.gtf \
-o result/featureCounts/gene_counts.txt \
temp/STAR/*Aligned.sortedByCoord.out.bam

# 提示完成
echo "featureCounts 运行完成于 $(date)"