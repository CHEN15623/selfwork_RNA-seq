
# 转录组数据分析流程
# -----1、下机数据质控，使用FASTQC或fastp，fastp融合了低质量过滤，去接头，裁切，功能强大便捷
# -----2、宿主基因组比较，使用HISAT2或STAR
# -----3、基因定量，使用RSEM或featureCounts或HTSeq-count



# 环境搭建,此分析流程在conda 24.11.2版本进行
# 创建工作目录
mkdir -p /osmgfs10000/home/wub/RNA-seq/ 
mkdir -p /osmgfs10000/home/wub/RNA-seq/temp Seq result
mkdir -p /osmgfs10000/home/wub/db/
# seq目录存放下机双端数据  temp目录临时文件存放 result结果文件存放
# 编辑mydata.txt文件，存放于result目录，mydata格式在实例数据中
wd= /osmgfs10000/home/wub/RNA-seq/
db= /osmgfs10000/home/wub/db/
# fastp安装
conda create -n fastp
conda activate fastp
conda install -c bioconda fastp
fastp -v
# 本流程版本fastp 0.23.4

# 安装rush进行数据并行处理
# rush安装应写入conda系统路径，保证在不同环境中都能调用rush，也可以直接讲rush写入系统路径，哪怕不用congda也能进行并行处理
# rush v0.5.4

# 1、数据质控
cd ${wd}
mkdir -p temp/QC result/QC
conda activate fastp
    time tail -n+2 result/mydata.txt|cut -f1|rush -j 4 \
      "fastp -i Seq/{}_1.fq.gz -I Seq/{}_2.fq.gz \
        -j temp/QC/{}_fastp.json -h temp/QC/{}_fastp.html \
        -o temp/QC/{}_1.fastq  -O temp/QC/{}_2.fastq \
        > temp/QC/{}.log 2>&1"

# mydata.txt是所有样品信息的汇总 -j 是并行处理的样本数，根据目前服务器资源设置
#  {}_1.fq.gz是输入需要质控的文件，{}是mydata中第一列的样本名称，_1.fq.gz是样本序列文件的其他部分，不同公司这一部分有区别，根据实际情况选择
    
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


# 2、转录组需要使用STAR比对参考基因组
conda create -n STAR #创建新环境
conda activate STAR
conda install -c bioconda star  #安装STAR
# 下载参考数据库，以小鼠为例
mkdir -p db/STAR
cd db/STAR
wget -c https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
wget -c https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz
wget -c https://ftp.ensembl.org/pub/release-113/gtf/mus_musculus/Mus_musculus.GRCm39.113.gtf.gz
gunzip Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
gunzip Mus_musculus.GRCm39.cdna.all.fa.gz
gunzip Mus_musculus.GRCm39.113.gtf.gz

# 构建STAR索引，索引时间较长，占用内存巨大，不建议使用命令行运行
STAR --runMode genomeGenerate \
--runThreadN 10 \
--genomeDir /osmgfs10000/home/wub/db/STAR/GRCm39_index \
--genomeFastaFiles /osmgfs10000/home/wub/db/STAR/Mus_musculus.GRCm39.dna.primary_assembly.fa \
--sjdbGTFfile /osmgfs10000/home/wub/db/STAR/Mus_musculus.GRCm39.113.gtf \
--sjdbOverhang 99

# STAR比对参考基因组
cd ${wd}           
conda activate STAR
mkdir -p temp/STAR

tail -n +2 result/mydata.txt | cut -f1 | rush -j 2  \
"STAR  --runThreadN 10 \
--genomeDir /osmgfs10000/home/wub/db/STAR/GRCm39_index \
--readFilesIn temp/QC/{}_1.fastq temp/QC/{}_2.fastq \
--outFileNamePrefix temp/STAR/{1}_data \
--outSAMtype BAM SortedByCoordinate \
> temp/STAR/{}.log 2>&1"

# 3、基因定量，此流程使用featureCounts
#   featureCounts依赖subread
#   下载subread
conda create subread
conda activate subread
conda install -c bioconda subread
cd ${wd}           
mkdir -p result/featureCounts

# 运行计数
featureCounts -T 20 -p -t exon -g gene_id \
-a ${db}/STAR/Mus_musculus.GRCm39.113.gtf \
-o result/featureCounts/gene_counts.txt \
temp/STAR/*Aligned.sortedByCoord.out.bam

# 提示完成
echo "featureCounts 运行完成于 $(date)"
# 4、差异基因分析
# 数据处理生成用于 DESeq2 的分组表
cut -f1,2 result/mydata.txt | tail -n +2 > samples_tmp.txt
awk 'BEGIN{print "sample\tcondition"} {print $0}' samples_tmp.txt > samples.txt
rm samples_tmp.txt
# 简化结果，准备 gene_counts_simple.txt
cut -f1,7- result/featureCounts/gene_counts.txt > result/featureCounts/gene_counts_simple.txt

# 加载R包，建议新建环境或用windows端运行，Linux系统R包下载编码易出错，需要多试几次。
conda activate r_env

# 进入R交互式界面
R
# 安装需要的包，安装不成功的话换其他方法，还是不行的话使用PC版R，由于包依赖比较多，版本不兼容较常见，本流程使用的r版本为老版本R version 4.3.3
install.packages("tidyverse")
install.packages("BiocManager")
library(BiocManager)
BiocManager::install("DESeq2")
# 退出
q()

# 上述安装失败可选择备用方法
conda activate r
conda install -y r-tidyverse r-BiocManager
R
library(BiocManager)
BiocManager::install("DESeq2")
q()

cd${wd}
Rscript deseq2_analysis.R #deseq2_analysis.R是Deseq2的R分析脚本，具体参数调整自行修改
# 接下来使用windows版本的R进行分析，富集分析依赖于enrich.R脚本，该版本使用的为R4.5.0


#  --以上内容仅供参考，所有脚本所依赖的软件版本更替可能会导致依赖项不兼容，具体参考各软件对应官网--  #
