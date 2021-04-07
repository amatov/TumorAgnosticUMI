from gwf import Workflow, AnonymousTarget
from os.path import join
from glob import glob
gwf = Workflow()

output_dir = "/genomedk/matov/umiseq_analysis/CRUK5Mb/"
bams = glob("*.sorted.bam")
bed_file = "NEW_METHOD_hg38_08feb2016_capture_targets.bed"
for bam in bams:
    bam_id = bam.split("/")[-1].split(".")[0]
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
