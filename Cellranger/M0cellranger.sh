##usage：M0cellranger.sh Seqdata_folder count
##Seqdata_folder,包含测序文件夹的目录，每个文件夹的测序文件来自同一个样本
##count,文件名中的sample名用_符号分割的数量，测序文件名如190065B_Ca2_4_2_R1.fq.gz,则count为2
##cellranger count 命令中--transcriptome需要根据情况选择参考目录
folder=$1;
c=$2;
old=`pwd`;
cd $folder;
#进入Seqdata_folder后,读取每个文件夹
dir=$(ls -l |awk '/^d/ {print $NF}');
for fold in $dir;do
	#先取sample名,然后跑cellranger count
	fl=`ls $fold|head -1`;
	sample="";
	for((i=1;i<=$c;i++));do
		h=`echo $fl|cut -d "_" -f $i`;
		if [ ${#sample} == 0 ];then
			sample=$h;
		else
			sample=$sample"_"$h;
		fi;
	done;
	echo $fold;
	echo $sample;
	output=$old"/M0cellranger/";
	echo $output;
	if [ ! -d $output ]; then mkdir $output; fi;
	cellranger count \
	--id="sample_"$sample \
	--sample=$sample \
	--transcriptome=/home/rlhua/reference/refdata-cellranger-GRCh38-3.0.0 \
	--fastqs=$fold \
	--localcores=16 \
	--expect-cells=6000 \
	--localmem=64;
done;
mv count_* $output;
cd $old;
