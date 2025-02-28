# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# Code inpired by
# https://github.com/Ensembl/ensembl-refget/blob/main/pipeline/indexer/create_indexdb.py
# Thanks, Arne :)

import argparse
import sys

from datetime import datetime, UTC

import tkrzw


def get_args():
    parser = argparse.ArgumentParser()

    # ( Copied from https://github.com/Ensembl/ensembl-refget/blob/main/pipeline/indexer/create_indexdb.py )
    # If you use /dev/shm on Slurm, the RAM to hold the DB will be accounted to
    # your process. You should allow for approx. 20GB per 100M entries.
    parser.add_argument("--dbfile",
        help= 'Database file to query.',
        required=True
    )
    parser.add_argument("--batch",
        help= 'Default batch size to query for. Default is 1k',
        default= 1_000,
        required=False,
        type=int
    )
    #
    args = parser.parse_args()
    return args


## MAIN ##
def main():
    args = get_args()

    _start = datetime.now(UTC)
    print(f"Opening DB {args.dbfile} with ({_start}: {_start - _start})", file=sys.stderr)
    db = tkrzw.DBM()
    db.Open(args.dbfile, False,
            dbm="HashDBM",
            no_wait=True).OrDie()

    _opened = datetime.now(UTC)
    print(f"DB open OK. Quering ({_opened}: {_opened - _start})", file=sys.stderr)

    # queadding uniparc data
    queried_cnt = 0

    # TODO:
    queries = list()
    cnt = 0
    for cnt, line in enumerate(sys.stdin, start = 1):
        key = line.strip()
        queries.append(key)
    queried_cnt = cnt
    
    res = db.GetMulti(*queries)
    print(res)
    print(len(res))


    _queried = datetime.now(UTC)

    # Closes the database.
    print(f"Queried {queried_cnt} times. Closing DB ({_queried}: {_queried - _start})", file=sys.stderr)
    db.Close().OrDie()
    _closed = datetime.now(UTC)

    _closed = datetime.now(UTC)
    print(f"DB close OK ({_closed}: {_closed - _start})", file=sys.stderr)


if __name__ == "__main__":
    main()

