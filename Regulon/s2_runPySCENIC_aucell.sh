## inputs
sn=$1
f_loom_path_expr=output/celltype_loom/${sn}.loom
## outputs
ctx_output=output/s2_control.reg.tsv
f_pyscenic_output=output/celltype_loom/s2_${sn}.pyscenic.loom

# ## set TMPDIR to current path, in case of no enough disk space on /tmp/
export TMPDIR=`pwd` 

pyscenic aucell \
    $f_loom_path_expr \
    $ctx_output \
    --output $f_pyscenic_output \
    --num_workers 20 \
    --seed 777

sample_name=output/celltype_loom/s3_${sn}
threads=10
min_regulon_size=10
python s3_postSCENIC.py $f_pyscenic_output $ctx_output $sample_name $threads $min_regulon_size

