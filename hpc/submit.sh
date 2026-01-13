#!/bin/bash
# Submit job to HPC cluster


# User and system configuration
user=emima@dtu.dk
model=XeonE5_2660v3
queue=compute
walltime=72:00
n_cores=12
ram_per_core=2

# Job name with date and time
job_name="job_$(date +'%Y%m%d_%H%M%S')"

echo "Submitting job: $job_name"
bsub -J $job_name \
    -q ${queue} \
    -o err_out_hpc/${job_name}.out \
    -e err_out_hpc/${job_name}.err \
    -n ${n_cores} \
    -R "span[hosts=1]" \
    -R "rusage[mem=${ram_per_core}GB]" \
    -W ${walltime} \
    -R "select[model == ${model}]" \
    job1.sub



