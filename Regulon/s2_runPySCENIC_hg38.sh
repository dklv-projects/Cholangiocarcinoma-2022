echo $0
cmddir=$(readlink -f $(dirname $0))
echo $cmddir
echo $1
echo $2
#if [ ! -d $outputdir ];then mkdir -p $outputdir;fi

## input loom
#f_loom_path_scenic=output/expr_sample.loom,the expression data used for regulon identification
#f_loom_path_expr=output/expr_withoutRLP.loom, the expression data to dig
f_loom_path_scenic=$1
f_loom_path_expr=$2

## reference database
cisDB=$cmddir"/cisTarget_databases_hg38"
f_tfs=$cisDB"/hs_hgnc_tfs.txt"
f_motif_path=$cisDB"/motifs-v9-nr.hgnc-m0.001-o0.0.tbl"
f_db_names=$cisDB"/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather"

## outputs
grn_output=s2_scenic.adj.tsv
ctx_output=s2_scenic.reg.tsv
f_pyscenic_output=s2.pyscenic.loom


arboreto_with_multiprocessing.py \
    $f_loom_path_scenic \
    $f_tfs \
    --method grnboost2 \
    --output $grn_output \
    --num_workers 20 \
    --seed 777

pyscenic ctx \
    $grn_output \
    $f_db_names \
    --annotations_fname $f_motif_path \
    --expression_mtx_fname $f_loom_path_scenic \
    --output $ctx_output \
    --num_workers 10
    ## Note:
    ## The reason why I didn't use `--mask_dropouts` parameters:
    ## In R version of SCENIC, mask dropouts is False (default when pySCENIC>0.9.6, current version: 0.10.1). 
    ## Here we match the behavior of R version.

## set TMPDIR to current path, in case of no enough disk space on /tmp/
export TMPDIR=`pwd` 

pyscenic aucell \
    $f_loom_path_expr \
    $ctx_output \
    --output $f_pyscenic_output \
    --num_workers 20 \
    --seed 777

sample_name=s3
threads=8
min_regulon_size=10
s3_postSCENIC=$cmddir"/s3_postSCENIC.py"
python $s3_postSCENIC $f_pyscenic_output $ctx_output $sample_name $threads $min_regulon_size

