import argparse
import pysam
import os
import csv
import numpy as np

class Variation:
    def __init__(self, pos, ref, alt):
        self.pos = pos
        self.ref = ref
        self.alt = alt

def str_set_join(s1, s2):
    return "".join(set(s1 + s2))

class FragmentFinder:
    def __init__(
        self,
        bam_file
    ):
        # Check if index exists, if not create an index file
        bai_filename = f"{bam_file}.bai"
        if not os.path.exists(bai_filename):
            print(f"No index file found ({bai_filename}), generating...")
            pysam.index(bam_file)

        # Open files for reading and writing
        self.bam_file = pysam.AlignmentFile(bam_file, "rb")

    def find_fragments(self, chrom, start, stop, mapq=20):

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

def regions(bed_file):
    with open(bed_file) as fd:
        reader = csv.reader(fd, delimiter = "\t")

        for line in reader:
            yield line[0], int(line[1]), int(line[2]) # yields: chromosome, start position, end positionm

nucleosome_index_map = {"A":0, "T":1, "G":2, "C":3}

def main(bam_file, bed_file, output_file,max_length=700):
    finder = FragmentFinder(bam_file)

    region_lst = list(regions(bed_file))

    tensor = np.zeros((len(region_lst), 4, max_length-1), dtype=np.uint16)

    for index,region in enumerate(region_lst):
        chrom, start, stop = region
        for ref_pos,read_allel,length in finder.find_fragments(chrom, start, stop):
            if 1<=length<=max_length:
                nucle = nucleosome_index_map.get(read_allel)
                tensor[index,nucle,length]+=1

    with open(output_file, "w") as fd:
        writer = csv.writer(fd, delimiter = "\t")
        print(tensor)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("bam_file", help="bam file")
    parser.add_argument("bed_file")
    parser.add_argument("output_file")
    args = parser.parse_args()
    main(
        args.bam_file,
        args.bed_file,
        args.output_file,
    )
