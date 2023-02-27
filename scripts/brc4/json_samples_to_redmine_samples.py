#!env python3

# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import argparse
import json
from pathlib import Path
from typing import List


class Dataset:

    def __init__(self, ds) -> None:
        self.component = ds["component"]
        self.species = ds["species"]
        self.name = ds["name"]
        self.samples = ds["runs"]

    def __str__(self) -> str:
        lines = []
        for sample in self.samples:
            lines.append(f"{sample['name']} : {','.join(sample['accessions'])}")
        return "\n".join(lines)

    def tab_list(self) -> str:
        """Returns the full description of the samples (component, sepcies, name) in a tab list string"""
        lines = []
        columns_init = [self.component, self.species, self.name]
        for sample in self.samples:
            columns = columns_init + [sample['name']]
            lines.append("\t".join(columns))
        return "\n".join(lines)


class DatasetCollection:

    def __init__(self, json_path: Path) -> None:
        self.datasets: List[Dataset] = []
        with json_path.open('r') as fh:
            for json_dataset in json.load(fh):
                dataset = Dataset(json_dataset)
                self.datasets.append(dataset)

    def __str__(self) -> str:
        lines = []
        for dataset in self.datasets:
            lines.append(str(dataset))
        return "\n".join(lines)

    def print_list(self) -> str:
        lines = []
        for dataset in self.datasets:
            lines.append(dataset.tab_list())
        return "\n".join(lines)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Convert a sample list from a json file into a sample list for Redmine')

    parser.add_argument('--json', type=str, required=True,
                        help='Path to the json file')

    parser.add_argument('--print_list', action='store_true',
                        help='Print a tab list of everything, instead of the default list for Redmine')
    args = parser.parse_args()
    data = DatasetCollection(Path(args.json))
    if args.print_list:
        print(data.print_list())
    else:
        print(data)
