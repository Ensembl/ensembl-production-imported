#!env python3

from redminelib import Redmine
import argparse
import os, json, re, time
import requests
import xml.etree.ElementTree as ET
from unidecode import unidecode
from pathlib import Path
 
url = 'https://redmine.apidb.org'
default_fields = dict(
        status_name = 'Data Processing (EBI)',
        cf_17 = "Data Processing (EBI)",
        )
insdc_pattern = "^GC[AF]_\d{9}(\.\d+)?$"
accession_api_url = "https://www.ebi.ac.uk/ena/browser/api/xml/%s"
veupathdb_id = 1976

def load_abbrevs(path):
    """
    Load a list of organism abbrevs from a file. Expected to be one per line
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
                    raise Exception("Can't load current abbrevs from a multicolumn string")
    return abbrevs

def retrieve_rnaseq_datasets(redmine, output_dir_path, build=None, abbrevs_file=None):
    """
    Get RNA-Seq metadata from Redmine, store them in json files.
    Each issue/dataset is stored as one file in the output dir
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
            problems.append({"issue" : issue, "desc" : problem})
            continue

        try:
            component = dataset["component"]
            organism = dataset["species"]
            dataset_name = dataset["name"]
            
            if dataset_name in used_names:
                problems.append({"issue" : issue, "desc" : "Dataset name already used: {dataset_name}"})
                continue
            else:
                used_names.append(dataset_name)
            
            if abbrevs_file and organism not in all_abbrevs:
                warn_abbrevs.append({"issue" : issue, "desc" : organism})
                

            ok_datasets.append({"issue" : issue, "desc" : organism})
            
            # Create directory
            dataset_dir = output_dir / component
            dataset_dir.mkdir(exist_ok=True)
            
            # Create file
            file_name = organism + "_" + dataset_name + ".json"
            dataset_file = dataset_dir / file_name
            with open(dataset_file, "w") as f:
                json.dump([dataset], f, indent=True)
        except Exception as error:
            problems.append({"issue" : issue, "desc" : str(error)})
            pass
        all_datasets.append(dataset)

    print("%d issues total" % len(issues))
    print_summaries(problems, "issues with problems")
    print_summaries(warn_abbrevs, "issues using unknown organism_abbrevs (maybe new ones). Those are still imported")
    print_summaries(ok_datasets, "datasets imported correctly")

    # Create a single merged file as well
    merged_file = Path(output_dir) / "all.json"
    with open(merged_file, "w") as f:
        json.dump(all_datasets, f, indent=True)
        
def print_summaries(summaries, description):
    desc_length = 64
    
    if summaries:
        print()
        print(f"{len(summaries)} {description}:")
        for summary in summaries:
            desc = summary["desc"]
            issue = summary["issue"]
            print(f"\t{desc:{desc_length}}\t{issue.id}\t({issue.subject})")
    
 
def parse_dataset(issue):
    """
    Extract RNA-Seq dataset metadata from a Redmine issue
    Return a nested dict
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

    #print("\n%s\tParsing issue %s (%s)" % (dataset["component"], issue.id, issue.subject))
    
    failed = False
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

def check_organism_abbrev(name):
    """Limited check: do not accept organism_abbrevs with special characters"""
    return not re.search("[ \/\(\)#\[\]:]", name)
    
def normalize_name(old_name):
    """Remove special characters, keep ascii only"""
    
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

def parse_samples(sample_str):
    samples = []
    
    # Parse each line
    lines = sample_str.split("\n")

    sample_names = dict()
    for line in lines:
        line = line.strip()
        if line == "": continue

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
                    raise Exception(f"Sample name and accessions are switched?")
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

def validate_accessions(accessions):
    """
    Check the accessions, to make sure we get proper ones
    """
    if "" in accessions: return False
    for acc in accessions:
        if not re.search("^[SE]R[RSXP]\d+$", acc):
            return False
    return True

def get_custom_fields(issue):
    """
    Put all Redmine custom fields in a dict instead of an array
    Return a dict
    """
    
    cfs = {}
    for c in issue.custom_fields:
        cfs[c["name"]] = c
    return cfs

def get_custom_value(customs, key):
   
    try:
        value = customs[key]["value"]
        if isinstance(value, list):
            if len(value) == 1:
                value = value[0]
            elif len(components) > 1:
                raise Exception("More than 1 values for key %s" % (key))
        return value
    except:
        print("No field %s" % (key))
        return ""
    

def get_issues(redmine, datatype, build=None):
    """
    Retrieve all issue for new genomes, be they with or without gene sets
    Return a Redmine ResourceSet
    """
    
    other_fields = { "cf_94" : datatype }
    if build:
        version_id = get_version_id(redmine, build)
        other_fields["fixed_version_id"] = version_id

    return list(get_ebi_issues(redmine, other_fields))

def get_version_id(redmine, build):
    """
    Given a build number, get the version id for it
    """
    versions = redmine.version.filter(project_id=veupathdb_id)
    version_name = "Build " + str(build)
    version_id = [version.id for version in versions if version.name == version_name]
    return version_id
    
def get_ebi_issues(redmine, other_fields=dict()):
    """
    Get EBI issues from Redmine, add other fields if provided
    Return a Redmine ResourceSet
    """

    # Other fields replace the keys that already exist in default_fields
    search_fields = { **default_fields, **other_fields }
    
    return redmine.issue.filter(**search_fields)
    

def main():
    # Parse command line arguments
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
        retrieve_dnaseq_datasets(redmine, args.output_dir, args.build, args.current_abbrevs)
    else:
        print("Need to say what data you want to --get: rnaseq? dnaseq?")

if __name__ == "__main__":
    main()
