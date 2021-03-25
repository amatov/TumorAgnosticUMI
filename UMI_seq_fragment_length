import argparse
import pysam
import py2bit
import os
import csv

from utils import regions


class Variation:
    def __init__(self, pos, ref, alt):
        self.pos = pos
        self.ref = ref
        self.alt = alt


def str_set_join(s1, s2):
    return "".join(set(s1 + s2))


class MutationFinder:
    def __init__(
        self,
        bam_file,
        output,
    ):
        # Check if index exists, if not create an index file
        bai_filename = f"{bam_file}.bai"
        if not os.path.exists(bai_filename):
            print(f"No index file found ({bai_filename}), generating...")
            pysam.index(bam_file)

        # Open files for reading and writing
        self._output_fp = open(output, "w")
        self.writer = csv.writer(self._output_fp, delimiter="\t")
        self.bam_file = pysam.AlignmentFile(bam_file, "rb")

    def __del__(self):
        self._output_fp.close()

    def write_finding(self, chrom, pos, ref, alt, length):
        context = self.tb.sequence(
            chrom, pos - self.context_flank, pos + self.context_flank + 1
        )
        assert len(context) == 2 * self.context_flank + 1
        self.writer.writerow(
            [
                chrom,
                str(pos),
                ref,
                alt,
                str(length),
                context,
            ]
        )

    def find_mutations(self, chrom, start, stop, mapq=20):

        for pileupcolumn in self.bam_file.pileup(
            contig=chrom,
            start=start,
            stop=stop,
            truncate=True,
            min_mapping_quality=mapq,
            ignore_overlaps=False,
            flag_require=2,  # proper paired
        ):
            ref_pos = pileupcolumn.reference_pos

            # Go through reads
            mem = {}
            for read in pileupcolumn.pileups:
                if read.is_del or read.is_refskip:
                    continue

                read = read.alignment

                if read.is_duplicate or read.is_secondary or read.is_supplementary:
                    continue

                query_name = read.query_name

                if query_name not in mem:
                    mem[query_name] = (
                        read.reference_start,
                        read.reference_end,
                        read.is_reverse,
                        read.is_read1,
                        read.query_sequence[read.query_position]

                    )
                else:
                    mem_start, mem_end, mem_reverse, mem_is_read1, mem_read_allel = mem[query_name]

                    del mem[query_name]
                    if (
                        mem_start is None
                        or mem_end is None
                        or read.reference_start is None
                        or read.reference_end is None
                        or read.is_reverse == mem_reverse
                    ):
                        continue

                    if read.is_reverse:
                        start = mem_start
                        end = read.reference_end
                        start_is_first = mem_is_read1
                    else:
                        start = read.reference_start
                        end = mem_end
                        start_is_first = not mem_is_read1

                    if start < end:
                        continue

                    read_allel = read.query_sequence[read.query_position]
                    if read_allel not in "ATGC" or read_allel != mem_read_allel:
                        continue

                    length = abs(read.template_length)

                    yield ref_pos, read_allel, length

def main(bam_file, bit_file, bed_file, vcf_file, output_file, is_hg38):
    tmp_file = "{}.tmp".format(output_file)
    finder = MutationFinder(
        bam_file, bit_file, vcf_file, tmp_file, long_vcf_chrom_format=is_hg38
    )
    for chrom, start, stop in regions(bed_file):
        finder.find_mutations(chrom, start, stop)
    os.rename(tmp_file, output_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("bam_file", help="bam file")
    parser.add_argument("bit_file")
    parser.add_argument("bed_file")
    parser.add_argument("vcf_file")
    parser.add_argument("output_file")
    parser.add_argument("hg_version", choices=["hg19", "hg38"])
    args = parser.parse_args()
    is_hg38 = args.hg_version == "hg38"
    main(
        args.bam_file,
        args.bit_file,
        args.bed_file,
        args.vcf_file,
        args.output_file,
        is_hg38,
    )
