mkdir -p logs/slurm

snakemake \
-j250 \
--latency-wait=180 \
--cluster="sbatch -c {threads} \
           --mem={resources.mem_mb}M \
           --time={resources.time} \
           --output=logs/slurm/%j.out \
           --error=logs/slurm/%j.out \
           --mail-user=w.vanrheenen-2@umcutrecht.nl \
           --mail-type=FAIL" \
--cluster-status ./cluster_status.py


# use include to read lines from other snakefile.