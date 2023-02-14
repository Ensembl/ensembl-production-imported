#!env python3

import argparse
import json
from pathlib import Path


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

    def full(self) -> str:
        lines = []
        cols_init = (self.component, self.species, self.name)
        for sample in self.samples:
            cols = list(cols_init)
            cols.append(sample['name'])
            lines.append("\t".join(cols))
        return "\n".join(lines)


class DatasetCollection:

    def __init__(self, json_path: Path) -> None:
        self.datasets = []
        with json_path.open('r') as fh:
            for json_dataset in json.load(fh):
                dataset = Dataset(json_dataset)
                self.datasets.append(dataset)
    
    def __str__(self) -> str:
        lines = []
        for dataset in self.datasets:
            lines.append(str(dataset))
        return "\n".join(lines)

    def full(self) -> str:
        lines = []
        for dataset in self.datasets:
            lines.append(dataset.full())
        return "\n".join(lines)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert a sample list from a json file into a sample list for Redmine')
    
    parser.add_argument('--json', type=str, required=True,
                        help='Path to the json file')
    
    parser.add_argument('--full', action='store_true',
                        help='Print a tab list of everything')
    args = parser.parse_args()
    data = DatasetCollection(Path(args.json))
    if args.full:
        print(data.full())
    else:
        print(data)