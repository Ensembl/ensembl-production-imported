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
        help=('Database file. Will be created if it does not exist.'
              ' It is recommended to create the DB in /dev/shm if available,'
              ' then copy it to permanent storage.' 
              ),
        required=True
    )
    # This is using a file hash database. This specifies the number of buckets
    # that the hash will have. Ideally, this number should be larger than the
    # number of keys inserted. If there are more entries, keys will hash to the
    # same bucket and be stored in a linked list. Eventually, performance will
    # deteriorate, so this should be chosen to be large enough.
    # If a DB grows over time to exceed the initial value, it should be rebuilt.
    # It's possible to rebuild an existing database. This is not currently part
    # of this script.
    # The dbsize option is only applied when creating a database. When opening
    # an existing DB, it is ignored.
    parser.add_argument("--dbsize",
        help=(
            'Tunes the number of hash buckets for the DB. This should ideally be'
            ' about 20%% more than the number of entries you expect.'
            ' Default is 1.2 billion.'
        ),
        default=1_200_000_000,
        required=False,
        type=int
    )
    #
    args = parser.parse_args()
    return args


def add_uniparc_data(db, stream, _start):
    _batch_start = datetime.now(UTC)
    cnt, collisions = 0, 0
    for cnt, line in enumerate(stream, start = 1):
        uniparc_id, md5u = line.strip().split() # NB: split on WS
        md5u = md5u.encode("utf-8")
        status, prev = db.SetAndGet(md5u, uniparc_id.encode("utf-8"), True) # overwrite: True
        if prev:
            collisions += 1
            extended = "\t".join([ prev.decode("utf-8"), uniparc_id ]).encode("utf-8")
            db.Set(md5u, extended, True) # overwrite: True
            _collision = datetime.now(UTC)
            print(f"Collision for '{md5u} : { extended }' ({_collision}: {_collision - _start})", file=sys.stderr)
        if cnt % 100_000_000 == 0:
            _info = datetime.now(UTC)
            print(f"Loaded {cnt} records ({collisions} collisions) ({_info}: {_info - _start})", file=sys.stderr)
    return cnt, collisions

## MAIN ##
def main():
    args = get_args()

    _start = datetime.now(UTC)
    print(f"Opening DB {args.dbfile} with {args.dbsize} buckets ({_start}: {_start - _start})", file=sys.stderr)
    db = tkrzw.DBM()
    # This is a Tkrzw file hash DB. Open as writeable.
    # The DB supports compression. This is not enabled because it saves a few
    # percent disk space but is almost half as fast
    # see https://dbmx.net/tkrzw/api-python/tkrzw.html
    db.Open(args.dbfile, True, # writable
            dbm="HashDBM",
            no_wait=True,
            truncate=True,
            sync_hard=True,
            offset_width=5, # 2^(2^5) = 2^32 ~ 4.29e9 ?
            align_pow=3,
            update_mode="UPDATE_IN_PLACE",
            num_buckets=args.dbsize).OrDie()

    _opened = datetime.now(UTC)
    print(f"DB open OK. Loading ({_opened}: {_opened - _start})", file=sys.stderr)

    # adding uniparc data
    loaded_cnt, collisions = add_uniparc_data(db, sys.stdin, _start)
    _loaded = datetime.now(UTC)

    # Closes the database.
    print(f"Loaded {loaded_cnt} records ({collisions} collisions). Closing DB ({_loaded}: {_loaded - _start})", file=sys.stderr)
    db.Close().OrDie()
    _closed = datetime.now(UTC)

    _closed = datetime.now(UTC)
    print(f"DB close OK ({_closed}: {_closed - _start})", file=sys.stderr)


if __name__ == "__main__":
    main()

