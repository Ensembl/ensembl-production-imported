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
from datetime import datetime, UTC
import sys

import tkrzw


def get_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--dbfile", help="Database file to query.", required=True)
    parser.add_argument(
        "--batch",
        help="Default batch size to query for. Default is 5k",
        default=5_000,
        required=False,
        type=int,
    )

    args = parser.parse_args()
    return args


def dump(data, stream):
    if not data:
        return 0, 0
    found, non_unique = len(data), 0

    for item in sorted(data.items(), key=lambda i: i[0]):
        decoded = list(map(lambda b: b.decode("utf-8"), item))
        if decoded[1].count("\t") > 0:
            non_unique += 1
        print("\t".join(decoded), file=stream)

    return found, non_unique


## MAIN ##
def main():
    args = get_args()

    _start = datetime.now(UTC)
    print(f"Opening DB {args.dbfile} ({_start}: {_start - _start})", file=sys.stderr)
    db = tkrzw.DBM()
    db.Open(args.dbfile, False, dbm="HashDBM", no_wait=True).OrDie()

    _opened = datetime.now(UTC)
    print(f"DB open OK. Quering ({_opened}: {_opened - _start})", file=sys.stderr)

    # queadding uniparc data
    queried_cnt = 0

    queries = []
    cnt, found, non_unique = 0, 0, 0
    for cnt, line in enumerate(sys.stdin, start=1):
        key = line.strip().upper()
        queries.append(key)
        if cnt % args.batch == 0:
            res = db.GetMulti(*queries)
            queries = []
            _found, _non_unique = dump(res, sys.stdout)
            found += _found
            non_unique += _non_unique
            _batch = datetime.now(UTC)
            print(
                f"Queried {cnt} times, found {found}, non unique {non_unique} ({_batch}: {_batch - _start})",
                file=sys.stderr,
            )
    queried_cnt = cnt

    # process last batch
    res = db.GetMulti(*queries)
    queries = []
    _found, _non_unique = dump(res, sys.stdout)
    found += _found
    non_unique += _non_unique

    _queried = datetime.now(UTC)

    # Closes the database.
    print(
        f"Queried {queried_cnt} times, found {found}, non unique {non_unique}. Closing DB ({_queried}: {_queried - _start})",
        file=sys.stderr,
    )
    db.Close().OrDie()
    _closed = datetime.now(UTC)

    _closed = datetime.now(UTC)
    print(f"DB close OK ({_closed}: {_closed - _start})", file=sys.stderr)


if __name__ == "__main__":
    main()
