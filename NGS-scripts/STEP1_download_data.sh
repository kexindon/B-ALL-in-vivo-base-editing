ssh kexindon@luria.mit.edu
srun --pty bash
# Lab server: 
cd /net/bmc-lab2/data/lab/sanchezrivera/kexindong/
# How to save sequencing data to the server:
srun rsync -av /net/bmc-pub17/data/bmc/public/Hemann/251030Hem /net/bmc-lab2/data/lab/sanchezrivera/kexindong/