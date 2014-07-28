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
import sarge
import sys


QUERIES = {
    "synonymous": "where type='snp' and impact = 'synonymous_coding' and ",
    "nonsynonymous": "where type='snp' and impact = 'non_syn_coding' and ",
    "indel": "where type='indel' and ",
    "ins": "where type='indel' and sub_type='ins' and ",
    "del": "where type='indel' and sub_type='del' and ",
    "other": ("where impact not in ('synonymous_coding', 'non_syn_coding')"
              " and type != 'indel' and "),
    "all": "where "
    }

GEMINI_COMMAND = (
    'gemini query -q "select variant_id from variants {0}"'
    ' --sample-filter "phenotype={1}" {2} --header'
)


def extract_data(database: str, query_type: str,
                 phenotype: int=1) -> pd.DataFrame:

    command = GEMINI_COMMAND.format(QUERIES[query_type] + " in_1kg=0 ",
                                    phenotype, database)
    print("Running", command)
    command = sarge.capture_stdout(command)

    series = pd.read_table(command.stdout, squeeze=True)

    return set(series.values)


def count_data(group1: set, group2: set, group1_name: str="PS-O",
               group2_name: str="SCR") -> dict:

    common_items = len(group1.intersection(group2))
    group1_items = len(group1)
    group2_items = len(group2)

    result = {
        group1_name: group1_items - common_items,
        group2_name: group2_items - common_items,
        "Shared": common_items
        }

    return result


def main():

    valid_choices = list(QUERIES.keys()) + ["complete"]

    parser = argparse.ArgumentParser()

    parser.add_argument("--group1", default="PS-O",
                        help="Name of group 1 [PS-O]")
    parser.add_argument("--group2", default="SCR",
                        help="Name of group 2 [SCR]")
    parser.add_argument("-t", "--type", choices=valid_choices,
                        help="Type of counting to use",
                        default="nonsynonymous")
    parser.add_argument("database", help="Source GEMINI database")

    options = parser.parse_args()

    if options.type != "complete":

        group1_set = extract_data(options.database, options.type, 1)
        group2_set = extract_data(options.database, options.type, 2)

        result = count_data(group1_set, group2_set, options.group1,
                            options.group2)

        group1 = result[options.group1]
        group2 = result[options.group2]
        shared = result["Shared"]

        print("Type:\t{0}".format(options.type))
        print("{0}:\t{1}".format(options.group1, group1),
              "{0}:\t{1}".format(options.group2, group2),
              "Shared:\t{0}".format(shared),
              "Total:\t{0}".format(group1 + group2 + shared),
              sep="\n")
    else:
        df = pd.DataFrame(
            index=[options.group1, options.group2, "Shared", "Total"])

        for count_type in ["nonsynonymous", "synonymous", "indel", "other",
                           "all"]:
            group1_set = extract_data(options.database, count_type, 1)
            group2_set = extract_data(options.database, count_type, 2)
            result = count_data(group1_set, group2_set, options.group1,
                                options.group2)

            group1 = result[options.group1]
            group2 = result[options.group2]
            shared = result["Shared"]

            column = [result[options.group1], result[options.group2],
                      result["Shared"], group1 + group2 + shared]
            df[count_type] = column

        df.to_csv(sys.stdout, sep="\t")


if __name__ == '__main__':
    main()