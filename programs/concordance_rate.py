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

from collections import Counter
from pathlib import Path

import vcf


def count_values(element: Path, min_depth: int,
                 min_fraction: float) -> Counter:

    """Calculate concordant / discordant calls for a paired Ft-S-Sd-S sample.

    :return: A Counter instance containing counts for concordant and
    discordant calls.
    """

    with element.open("rb") as handle:
        reader = vcf.Reader(handle)
        contents = list()

        for record in reader:
            if record.FILTER:
                continue
            fractions = list()
            valid_calls = [call for call in record if call.called and
                           call.is_variant]

            if not valid_calls:
                continue

            depth = min([call.data.DP for call in valid_calls])

            if depth < min_depth:
                continue

            if all([hasattr(item.data, "FA") for item in valid_calls]):
                fraction = min((item.data.FA for item in valid_calls))
            elif all([hasattr(item.data, "FREQ") for item in valid_calls]):
                fraction = min((item.data.FREQ for item in valid_calls))

            if fraction >= min_fraction:
                if len(valid_calls) == 2:
                    if len(record.ALT) != 1:
                        # Two genotypes that are different
                        contents.append("discordant")
                    else:
                        contents.append("concordant")
                else:
                        contents.append("discordant")

    return Counter(contents)


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--min-depth", type=int, default=0,
                        help="Minimum depth")
    parser.add_argument("--min-fraction", type=float, default=0.01,
                        help="Minimum fraction")
    parser.add_argument("source",
                        help="Source VCF file (with Ft-S-Sd-S pair)")

    options = parser.parse_args()

    counts = count_values(Path(options.source), options.min_depth,
                          options.min_fraction)


if __name__ == "__main__":
    main()
