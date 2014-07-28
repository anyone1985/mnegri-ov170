#!/usr/bin/python3 

import argparse
import vcf
import pysam

def filter_file(source, destination):

    reader = vcf.VCFReader(filename=source)

    with open(destination, "w") as handle:
        writer = vcf.VCFWriter(handle, reader)

        for record in reader:
            if all(not sample.called for sample in record):
                continue
            writer.write_record(record)

    final = pysam.tabix_index(destination, preset="vcf", force=True)

    return final


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("source", help="Source file")
    parser.add_argument("destination", help="Destination template")

    options = parser.parse_args()

    destination = (options.destination + ".vcf" 
            if not str(options.destination).endswith(".vcf")
            else options.destination)

    final = filter_file(options.source, destination)

    print("Final saved file: {}".format(final))


if __name__ == "__main__":
    main()

