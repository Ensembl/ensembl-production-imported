### Bulk Handover wrapper:

*Wrapper script [Bulk_CoreDB_Handover_Beta.sh](https://github.com/Ensembl/ensembl-production-imported/blob/main/scripts/utility/Bulk_CoreDB_Handover_Beta.sh) facilitates processing multiple ensembl core database handovers to Ensembl Automation.*

- The wrapper integrates basic info at runtime and makes indpenedent calls to the ensembl-prodinf-tools [handover client](https://github.com/Ensembl/ensembl-prodinf-tools/blob/main/src/scripts/handover_client.py)
- The wrapper requires some basic information, such as division and release channel (i.e. Main site or beta/MVP)
- Allows for 'dry runs' before submission to check details are accurate and error free. 
- Standard output to terminal can be shown with or without colourised text.
- Required input: flat text file of core database names listed one per line.
- Captures to log file(s) (per core). Recording such info as handover-token-ID. 

**Usage**:

```
Bulk_CoreDB_Handover_Beta <Core DB host> <INPUT_CORES_LIST> <main or beta ?> <Division> <'Short HO context description'>
```
E.g: 
```
Bulk_CoreDB_Handover_Beta dummy-host-1 Cores.list.txt main metazoa 'Metazoa handover E112'
````
