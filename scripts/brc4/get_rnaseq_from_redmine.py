#!env python3

from typing import Dict, List, Tuple
import argparse
import json
import re
from unidecode import unidecode
from pathlib import Path
from redminelib import Redmine
from redminelib.resources import Issue as RedmineIssue
from redminelib.managers import ResourceManager as RedmineResource

url = 'https://redmine.apidb.org'
default_fields = dict(
    status_name='Data Processing (EBI)',
    cf_17="Data Processing (EBI)",
)
insdc_pattern = r'^GC[AF]_\d{9}(\.\d+)?$'
accession_api_url = "https://www.ebi.ac.uk/ena/browser/api/xml/%s"
veupathdb_id = 1976


def load_abbrevs(path: str) -> List[str]:
    """
    Load a list of organism abbrevs from a file. Expected to be one per line.

    Args:
        path: Path to the organism abbrevs file.

    Returns:
        A list of all organism_abbrevs.

    """
    if not path:
        print("Warning: I don't have a list of older abbrevs to compare with.")
        return []
    abbrevs = []
    with open(path, "r") as abbr_file:
        for line in abbr_file:
            line = line.rstrip()
            if line:
                fields = line.split("\t")
                if len(fields) == 1:
                    abbrevs.append(line)
                else:
                    raise Exception(
                        "Can't load current abbrevs from a multicolumn string")
    return abbrevs


def retrieve_rnaseq_datasets(redmine: Redmine, output_dir_path: str, build: int = None,
                             abbrevs_file: str = None) -> None:
    """
    Get RNA-Seq metadata from Redmine, store them in json files.
    Each issue/dataset is stored as one file in the output dir.

    Args:
        redmine: A connected Redmine object.
        output_dir_path: Directory where the dataset files are to be stored.
        build: BRC build number.
        abbrevs_file: Path to a list of organism_abbrevs that are already in use.
    """

    all_abbrevs = load_abbrevs(abbrevs_file)

    issues = get_issues(redmine, "RNA-seq", build)
    if not issues:
        print("No files to create")
        return
    
    # Create the output dir
    output_dir = Path(output_dir_path)
    output_dir.mkdir(exist_ok=True)
    
    # Write all datasets in files
    all_datasets = []
    used_names = []

    problems = []
    ok_datasets = []
    warn_abbrevs = []
    
    for issue in issues:
        dataset, problem = parse_dataset(issue)
        
        if problem:
            problems.append({"issue": issue, "desc": problem})
            continue

        try:
            component = dataset["component"]
            organism = dataset["species"]
            dataset_name = dataset["name"]
            
            if dataset_name in used_names:
                problems.append(
                    {"issue": issue, "desc": f"Dataset name already used: {dataset_name}"}
                )
                continue
            else:
                used_names.append(dataset_name)
            
            if abbrevs_file and organism not in all_abbrevs:
                warn_abbrevs.append({"issue": issue, "desc": organism})
                
            ok_datasets.append({"issue": issue, "desc": organism})
            
            # Create directory
            dataset_dir = output_dir / component
            dataset_dir.mkdir(exist_ok=True)
            
            # Create file
            file_name = organism + "_" + dataset_name + ".json"
            dataset_file = dataset_dir / file_name
            with open(dataset_file, "w") as f:
                json.dump([dataset], f, indent=True)
        except Exception as error:
            problems.append({"issue": issue, "desc": str(error)})
            pass
        all_datasets.append(dataset)

    print("%d issues total" % len(issues))
    print_summaries(problems, "issues with problems")
    print_summaries(
        warn_abbrevs,
        "issues using unknown organism_abbrevs (maybe new ones). Those are still imported"
    )
    print_summaries(ok_datasets, "datasets imported correctly")

    # Create a single merged file as well
    merged_file = Path(output_dir) / "all.json"
    with open(merged_file, "w") as f:
        json.dump(all_datasets, f, indent=True)

       
def print_summaries(summaries: Dict, description: str) -> None:
    """Print a summary of various counts.

    This will print one line for each issue in the dict, with its description, the issue id
    and the issue subject.

    Args:
        summaries: Dict with 2 keys:
            issue: A Redmine Issue object.
            desc: A description for that issue.
    """
    desc_length = 64
    
    if summaries:
        print()
        print(f"{len(summaries)} {description}:")
        for summary in summaries:
            desc = summary["desc"]
            issue = summary["issue"]
            print(f"\t{desc:{desc_length}}\t{issue.id}\t({issue.subject})")
    
 
def parse_dataset(issue: RedmineIssue) -> Tuple[Dict, str]:
    """
    Extract RNA-Seq dataset metadata from a Redmine issue.

    Args:
        issue: A Redmine issue.

    Returns:
        A tuple of 2 objects:
            datasets: A dict representing a dataset, with the following keys:
                component: String for the BRC component DB.
                species: String for the organism abbrev.
                name: String for the internal dataset name.
            problem: A string description if there was a parsing problem
            (empty string otherwise).
    """
    customs = get_custom_fields(issue)
    dataset = {
        "component": "",
        "species": "",
        "name": "",
        "runs": [],
    }
    problem = ""

    dataset["component"] = get_custom_value(customs, "Component DB")
    dataset["species"] = get_custom_value(customs, "Organism Abbreviation").strip()
    dataset["name"] = get_custom_value(customs, "Internal dataset name").strip()

    if not dataset["species"]:
        problem = "Missing Organism Abbreviation"
    elif not check_organism_abbrev(dataset["species"]):
        problem = f"Wrong Organism Abbreviation format: '{dataset['species']}'"
    elif not dataset["name"]:
        problem = "Missing Internal dataset name"
    else:
        dataset["name"] = normalize_name(dataset["name"])
    
    # Get samples/runs
    samples_str = get_custom_value(customs, "Sample Names")
    try:
        samples = parse_samples(samples_str)
        
        if not samples:
            problem = "Missing Samples"
        
        dataset["runs"] = samples
    except Exception as e:
        problem = str(e)
    
    return dataset, problem


def check_organism_abbrev(name: str) -> bool:
    """Check the organism_abbrevs string format to avoid special characters.

    Args:
        name: organism_abbrev to check.
    
    Returns:
        True if the organism_abbrev format is correct.
        False otherwise.
    """
    return not re.search(r'[ \/\(\)#\[\]:]', name)


def normalize_name(old_name: str) -> str:
    """Remove special characters from an organism_abbrev, keep ascii only.

    Args:
        old_name: the organism_abbrev to format.
    
    Returns:
        The formatted organism_abbrev.
    """
    
    # Remove any diacritics
    name = old_name.strip()
    name = unidecode(name)
    name = re.sub(r"[ /]", "_", name)
    name = re.sub(r"[;:.,()\[\]{}]", "", name)
    name = re.sub(r"\+", "_plus_", name)
    name = re.sub(r"\*", "_star_", name)
    name = re.sub(r"%", "pc_", name)
    name = re.sub(r"_+", "_", name)
    if re.search(r"[^A-Za-z0-9_.-]", name):
        print("WARNING: name contains special characters: %s (%s)" % (old_name, name))
        return
    
    return name


def parse_samples(sample_str: str) -> List[Dict]:
    """Parse a list of samples from a Redmine task.

    Args:
        sample_str: The value of the field 'Sample Names' from an RNA-Seq Redmine task.
    
    Returns:
        A list of samples dicts, with the following keys:
            name: the name of the sample.
            accessions: a list of string representing the SRA accessions for that sample.
    """
    samples = []
    
    # Parse each line
    lines = sample_str.split("\n")

    sample_names = dict()
    for line in lines:
        line = line.strip()
        if line == "":
            continue

        # Get sample_name -> accessions
        parts = line.split(":")
        if len(parts) > 2:
            end = parts[-1]
            start = ":".join(parts[:-1])
            parts = [start, end]
        
        if len(parts) == 2:
            sample_name = parts[0].strip()
            
            if sample_name in sample_names:
                raise Exception("Several samples have the same name '%s'" % sample_name)
            else:
                sample_names[sample_name] = True
            
            accessions_str = parts[1].strip()
            accessions = [x.strip() for x in accessions_str.split(",")]
            
            if not validate_accessions(accessions):
                if validate_accessions(sample_name.split(",")):
                    raise Exception("Sample name and accessions are switched?")
                else:
                    raise Exception(f"Invalid accession among '{accessions}'")
            
            sample = {
                "name": normalize_name(sample_name),
                "accessions": accessions
            }
            samples.append(sample)
        else:
            raise Exception("Sample line doesn't have 2 parts: '%s'" % line)
    
    return samples


def validate_accessions(accessions: str) -> bool:
    """Check SRA accessions format, to make sure we get proper ones.

    Args:
        accessions: a list of strings to check

    Return:
        True if all strings are proper SRA accessions.
        False if at least one is not a proper SRA accession.
    """
    if "" in accessions:
        return False
    for acc in accessions:
        if not re.search(r'^[SE]R[RSXP]\d+$', acc):
            return False
    return True


def get_custom_fields(issue: RedmineIssue) -> Dict:
    """Put all Redmine custom fields in a dict instead of an array.

    Args:
        issue: A Redmine issue.
    
    Returns:
        A dict where each key is a custom field.
    """
    
    cfs = {}
    for c in issue.custom_fields:
        cfs[c["name"]] = c
    return cfs


def get_custom_value(customs: Dict, key: str) -> str:
    """Retrieve a custom value from a customs dict.

    Args:
        customs: Dict of customs values gotten from get_custom_fields.
        key: Key to extract the value from the custom dict.
    
    Returns:
        A single value.
        Throws an exception if there are more than 1 value.
        If there is no such key in the dict, return an empty string.
    """
   
    try:
        value = customs[key]["value"]
        if isinstance(value, list):
            if len(value) == 1:
                value = value[0]
            elif len(value) > 1:
                raise Exception("More than 1 values for key %s" % (key))
        return value
    except KeyError:
        print("No field %s" % (key))
        return ""
    

def get_issues(redmine: Redmine, datatype: str, build: int = None) -> List[RedmineIssue]:
    """Retrieve all issue for new genomes, be they with or without gene sets.

    Args:
        redmine: A Redmine connected object.
        datatype: What datatype to use to filter the issues.
        build: The BRC build to use to filter.

    Returns:
        A list of Redmine issues.
    """
    
    other_fields = {"cf_94": datatype}
    if build:
        version_id = get_version_id(redmine, build)
        other_fields["fixed_version_id"] = version_id

    return list(get_ebi_issues(redmine, other_fields))


def get_version_id(redmine: Redmine, build: int) -> int:
    """Given a build number, get the Redmine version id for it.

    Args:
        redmine: A Redmine connected object.
        build: The BRC build to use to filter.

    Returns:
        The version id from Redmine for that build.
    """
    versions = redmine.version.filter(project_id=veupathdb_id)
    version_name = "Build " + str(build)
    version_id = [version.id for version in versions if version.name == version_name]
    return version_id

   
def get_ebi_issues(redmine, other_fields=dict()) -> RedmineResource:
    """Get EBI issues from Redmine, add other fields if provided.

    Args:
        redmine: A Redmine connected object.
        other_fields: A dict of fields to provide to filter the issues.

    Returns:
        A Redmine resource set.
    """
    # Other fields replace the keys that already exist in default_fields
    search_fields = {**default_fields, **other_fields}
    
    return redmine.issue.filter(**search_fields)
    

def main():
    parser = argparse.ArgumentParser(description='Retrieve metadata from Redmine')
    
    parser.add_argument('--key', type=str, required=True,
                        help='Redmine authentification key')
    parser.add_argument('--output_dir', type=str, required=True,
                        help='Output_dir')
    # Choice
    parser.add_argument('--get', choices=['rnaseq', 'dnaseq'], required=True,
                        help='Get rnaseq, or dnaseq issues')
    # Optional
    parser.add_argument('--build', type=int,
                        help='Restrict to a given build')
    parser.add_argument('--current_abbrevs', type=str,
                        help='File that contains the list of current organism_abbrevs')
    args = parser.parse_args()
    
    # Start Redmine API
    redmine = Redmine(url, key=args.key)
    
    # Choose which data to retrieve
    if args.get == 'rnaseq':
        retrieve_rnaseq_datasets(redmine, args.output_dir, args.build, args.current_abbrevs)
    elif args.get == 'dnaseq':
        # TODO
        # retrieve_dnaseq_datasets(redmine, args.output_dir, args.build, args.current_abbrevs)
        print("Not yet implemented")
    else:
        print("Need to say what data you want to --get: rnaseq? dnaseq?")


if __name__ == "__main__":
    main()
