#!/usr/bin/python3
# -*- coding: utf-8 -*-

#  Copyright 2014 Luca Beltrame <luca.beltrame@marionegri.it>
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.


import argparse

import pandas as pd
from sqlalchemy import create_engine
from sqlalchemy import MetaData, Table, Column, Integer, String, Float
from sqlalchemy import ForeignKey
from sqlalchemy.orm import relationship
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import Session
import vcf

# These work around proper gene association bugs (due to amplicon extending
# past boundaries, or annotation picking things on the wrong strand)

substitutions = {
    "ZEB1-AS1": "ZEB1",      # Different strand
    "RHNO1": "FOXM1",        # Different strand
    # "PLA2G6": np.nan,     # off-target?
    "C3orf72": "FOXL2",      # Different strand
    # "MC1R": np.nan,       # past the promoter
    "ACAA1": "MYD88",        # Different strand
    "VIM-AS1": "VIM",        # Different strand
    "LOC100507424": "FOXM1",  # Wrong annotation?
    "MTOR-AS1": "MTOR",
    "EGFR-AS1": "EGFR",
    "WRAP53": "TP53",
    "EPM2AIP1": "MLH1",
    "C5orf22": "DROSHA",
    "C9orf53": "CDKN2A",
    "LYRM5": "KRAS",
    "N4BP2L1": "BRCA2",
    "RMDN3": "RAD51",
    "NBR2": "BRCA1",
    "CNTD2": "AKT2",
    "HSCB": "CHEK2"
}


Base = declarative_base()


class Variant(Base):

    __tablename__ = "variants"

    variant_id = Column(Integer, primary_key=True)
    sample_id = Column(Integer, ForeignKey("samples.sample_id"))

    chrom = Column(String, nullable=False)
    start = Column(Integer, nullable=False)
    end = Column(Integer, nullable=False)
    gene = Column(String)
    depth = Column(Integer)
    fraction = Column(Float)
    ref = Column(String)
    alt = Column(String)

    def __repr__(self):
        return "<Variant {0}:{1}-{2} - {3}>".format(self.chrom, self.start,
                                                    self.end, self.gene)


class Sample(Base):
    __tablename__ = "samples"
    sample_id = Column(Integer, primary_key=True)

    name = Column(String, nullable=False)
    sample_name = Column(String)  # Tissue bank name
    phenotype = Column(String)

    variants = relationship("Variant", backref="samples")

    def __init__(self, name, phenotype=None, sample_name=None):
        self.name = name
        self.phenotype = phenotype
        self.sample_name = sample_name

    def __repr__(self):
        return "<Sample {0} - {1}>".format(self.name, self.phenotype)


def _correct_names(name: str) -> str:

    """Handle wrongly spelled sample IDs in VCFs."""

    name = name.split("_")[1]

    if len(name) < 5:
        name = "2" + name if len(name) == 4 else "20" + name

    return name


def create_database_samples(reader: vcf.VCFReader) -> dict:

    all_samples = reader.samples

    db_samples = dict()

    for item in all_samples:
        correct_name = _correct_names(item)
        entry = Sample(name=item, phenotype="PS-O" if "primary" in item
                       else "SCR", sample_name=correct_name)
        db_samples[item] = entry

    return db_samples


def create_mutation_database(vcf_file: str, database_file: str):

    engine = create_engine('sqlite:///{}'.format(database_file))
    Base.metadata.create_all(engine)

    session = Session(bind=engine)

    reader = vcf.VCFReader(filename=vcf_file)
    db_samples = create_database_samples(reader)

    recordnum = 0

    for record in reader:

        if recordnum % 1000 == 0:
            print("{} records processed.".format(recordnum))
            session.commit()

        chrom = record.CHROM
        start = record.POS - 1
        end = start + len(record.REF)

        alt = pd.Series(record.ALT)

        genes = [item.split("|")[5] for item in record.INFO["EFF"]]

        if len(genes) > 1:

            if all([item == genes[0] for item in genes]):
                gene = genes[0]
            else:
                genes = [substitutions[gene] if gene in substitutions
                         else gene for gene in genes]
                genes = [gene for gene in genes if gene]
                gene = pd.Series(genes).unique().item()

        else:
            gene = genes[0]

        called_samples = [item for item in record.samples if item.called
                          and item.is_variant]

        for call in called_samples:

            try:
                depth = int(call.data.DP)
            except Exception:
                depth = int(record.INFO["DP"])

            if (hasattr(call.data, "FA") and call.data.FA is not None
                and not isinstance(call.data.FA, list)):

                fraction = call.data.FA

            elif hasattr(call.data, "FREQ") and call.data.FREQ is not None:
                fraction = float(call.data.FREQ)

            elif hasattr(call.data, "AD") and call.data.AD is not None:
                fraction = call.data.AD[1] / depth

            alt_base = ",".join([item.sequence for item in record.ALT])

            variant = Variant(chrom=chrom, start=start, end=end, depth=depth,
                              fraction=fraction, gene=gene, ref=record.REF,
                              alt=alt_base)

            sample_instance = db_samples[call.sample]
            sample_instance.variants.append(variant)
            session.add(sample_instance)
            recordnum += 1

    session.commit()


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("source", help="Source VCF")
    parser.add_argument("destination", help="Destination database")

    options = parser.parse_args()

    create_mutation_database(options.source, options.destination)


if __name__ == '__main__':
    main()