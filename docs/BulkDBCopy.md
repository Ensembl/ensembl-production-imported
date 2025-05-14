### Bulk Database copy wrapper:

*Wrapper script [Bulk_CoreDB_Copy.sh](https://github.com/Ensembl/ensembl-production-imported/blob/main/scripts/utility/Bulk_CoreDB_Copy.sh) facilitates processing multiple database copy jobs.*

- Integrates basic info such as MYSQL source and target hosts, core databases and makes indpenedent calls to [dbcopy_client](https://github.com/Ensembl/ensembl-prodinf-tools/blob/main/src/scripts/dbcopy_client.py)
- Wrapper allows for perforning dry runs before submission to check details are accurate and error free. 
- Standard output to terminal can be shown with or without colourised text.
- Required input: Sorted flat textfiles of core database names one per line. e.g:
- Not suitable for creating database copys on the same MYSQL host, only between independent hosts.
- Can perform **submission or new jobs**, or **listing existing db_copy jobs**.

*Source_cores.txt*:
```
prefix_apis_melifera_core_56_112_1
prefix_drosophila_melanogaster_core_56_112_1
prefix_homo_sapiens_core_56_112_1
```
*Target_cores.txt*:
```
apis_melifera_core_112_1
drosophila_melanogaster_core_112_1
homo_sapiens_core_112_1
```
*Usage*:

```
Bulk_DB_Copy <SOURCE host> <Source Core(s) list file> <TARGET host> <Target Core(s) list file> <Operation: submit | list>>
```
E.g. 
```
Bulk_DB_Copy dummy-source-host-1 Source_cores.txt dummy-target-host-1 Target_cores.txt submit
```

