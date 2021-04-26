from gwf import Workflow, AnonymousTarget
from os.path import join
from glob import glob
#gwf = Workflow()
gwf = Workflow(
        defaults={
            'cores': 1,
            'memory': '16g',
            'walltime': '01:00:00'})

#bam_files=glob('~/genomedk/matov/umiseq_analysis/CRUK5Mb/*/*consensus.sort.bam')
bam_files=glob('*consensus.sort.bam')
#realpath C70A06997D_cfdna_N289-54_consensus.sort.bam
#/home/matov/matovanalysis/umiseq_analysis/CRUK5Mb/C70A06997D_cfdna_N289-54_consensus.sort.bam

print(bam_files)

output_dir = "."

bed_file = "NEW_METHOD_hg38_08feb2016_capture_targets.bed"
for bam in bam_files:
    bam_id = bam.split("/")[-1].split(".")[0].replace("-","_")
    output = join(output_dir, f"{bam_id}.txt")
    gwf.target(
        f"UMI_seq_fragment_length_data_{bam_id}",
        inputs=[bam, bed_file],
        outputs=[output],
        walltime="04:00:00",
    ) << """
    python UMI_seq_fragment_length.py {} {} {}
    """.format(
        bam, bed_file, output
    )
