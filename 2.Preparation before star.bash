#转录组需要使用STAR比对参考基因组
conda create -n STAR #创建新环境
conda activate STAR
conda install -c bioconda star  #安装STAR
#下载参考数据库，以小鼠为例
mkdir -p db/STAR
cd db/STAR
wget -c https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
wget -c https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz
wget -c https://ftp.ensembl.org/pub/release-113/gtf/mus_musculus/Mus_musculus.GRCm39.113.gtf.gz
gunzip Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
gunzip Mus_musculus.GRCm39.cdna.all.fa.gz
gunzip Mus_musculus.GRCm39.113.gtf.gz

#构建STAR索引
STAR --runMode genomeGenerate \
--runThreadN 10 \
--genomeDir /osmgfs10000/home/wub/db/STAR/GRCm39_index \
--genomeFastaFiles /osmgfs10000/home/wub/db/STAR/Mus_musculus.GRCm39.dna.primary_assembly.fa \
--sjdbGTFfile /osmgfs10000/home/wub/db/STAR/Mus_musculus.GRCm39.113.gtf \
--sjdbOverhang 99

