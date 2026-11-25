## Quick guide

To compile code first load dependencies needed:

```bash
module load genoa

module load gompi
```

Then:

```bash
# Compile src code (once)
make cavity_par

# Change to target resoluton directory
cd [path to target resolution directory]

# Submit a job
sbatch run.sh
```