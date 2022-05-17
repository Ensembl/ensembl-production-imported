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

"""This script creates a database of stable ids to check their uniqueness

How to use:

1. You need to create the (empty) database with --create.
   This will replace any file with that same name.
2. Add stable_ids from EnsEMBL cores on a given server (you can filter by db name with --prefix)
3. Print out a summary with --summary, or the whole list of duplicates with --list_duplicates

Example:
    $ python check_stable_ids.py --db ids.sqlite3 --create
    $ python check_stable_ids.py --db ids.sqlite3 --add \
        --host $HOST --port $PORT --user $USER --password $PASSWORD
    $ python check_stable_ids.py --db ids.sqlite3 --summary
    $ python check_stable_ids.py --db ids.sqlite3 --list_duplicates
    
"""

import argparse
from typing import List
import mysql.connector
from mysql.connector.cursor import MySQLCursor
import os
import errno
from sqlalchemy import Column, Integer, Index, String, ForeignKey, insert, select, text
from sqlalchemy.orm import declarative_base
from sqlalchemy import create_engine
from dataclasses import dataclass


@dataclass
class Feature:
    """Simple feature object"""

    name: str
    biotype: str
    feature: str


class CoreServer(object):
    """Interface to a MySQL server with cores in it
    """

    def __init__(self, host: str, port: str, user: str, password: str):
        """Init the server object and connect to it
        
        Args:
            host: MySQL server host
            port: MySQL server port
            user: MySQL server user
            password: MySQL server password
        """
        
        self.host = host
        self.port = port
        self.user = user
        self.password = password
        self.engine = None
        self.db = None
        self.cores: List[str] = []

        self.connect()

    def connect(self) -> None:
        """Create a connection to the server
        """
        
        self.db = mysql.connector.connect(
            user=self.user,
            passwd=self.password,
            host=self.host,
            port=self.port)

    def _cursor(self) -> MySQLCursor:
        return self.db.cursor()
    
    def get_cores(self, prefix: str, build: int) -> None:
        """Retrieve the list of cores on the server
        
        Args:
            prefix: only use cores that start with this prefix.
            build: only use cores from this build.
        """

        if not self.db:
            self.connect()
        
        query = "SHOW DATABASES LIKE '%_core_%'"

        cursor = self._cursor()
        cursor.execute(query)
        
        for db in cursor:
            if prefix and not db[0].startswith(prefix):
                continue
            if build and f"_core_{build}_" not in db[0]:
                continue
            self.cores.append(db[0])
    
    def get_features(self, core_name: str, feature: str) -> List[Feature]:
        """Retrieve all the stable ids for a feature table in a given core
        
        Args:
            core_name: name of the core database to use
            feature: feature name (gene, transcript, or translation)
        
        Returns:
            A list of Features
        """
        
        self.db.database = core_name
        cursor = self._cursor()
        
        biotype = ''
        query = ''
        if feature in ('gene', 'transcript'):
            query = f"SELECT stable_id, biotype FROM {feature}"
        elif feature == 'translation':
            query = f"SELECT stable_id FROM {feature}"
        else:
            raise Exception(f"Unsupported feature type: {feature}")

        cursor.execute(query)
        
        for row in cursor:
            name = row[0]
            if feature in ('gene', 'transcript'):
                biotype = row[1]
            feat = Feature(name=name, feature=feature, biotype=biotype)
            yield feat
    
    def get_core_metadata(self, core_name: str, key: str) -> list:
        """Retrieve a metadata value from a given core
        
        Args:
            core_name: name of the core database to use
            key: metadata_key to get the metadata_value from
        
        Returns:
            List of all values
        """
        
        self.db.database = core_name
        cursor = self._cursor()
        
        query = "SELECT meta_value FROM meta WHERE meta_key = %s"
        cursor.execute(query, [key])
        
        rows = []
        for row in cursor:
            rows.append(row[0])
        
        return rows


##############################################################################
# Prepare the SQLalchemy schema for the stable_id database
Base = declarative_base()


class Db(Base):
    __tablename__ = 'db'

    db_id = Column(Integer, primary_key=True)
    db_name = Column(String)
    production_name = Column(String)


class StableId(Base):
    __tablename__ = 'stable_id'

    name_id = Column(Integer, primary_key=True)
    db_id = Column(ForeignKey('db.db_id'))
    feature = Column(String)
    name = Column(String)
    biotype = Column(String)
    Index('ix_feat_name', feature, name)


class StableIdDB(object):
    """Representation of an SQLite database of stable ids
    """
    all_features = ('gene', 'transcript', 'translation')
    
    def __init__(self, path: str):
        """Init the database object
        
        Args:
            path: path to the SQLite database file
        """
        self.path = path
        self.engine = None
    
    def create(self) -> None:
        """Create the SQLite database
        
        Note:
            Recreate the file from scratch if it already exists
        """
        
        if os.path.exists(self.path):
            os.remove(self.path)
        
        url = f"sqlite+pysqlite:///{self.path}"
        self.engine = create_engine(url, echo=False, future=True)
        Base.metadata.create_all(self.engine)

    def connect(self) -> None:
        """Connect to an existing SQLite database
        """
        
        if not os.path.exists(self.path):
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), self.path)

        if not self.engine:
            url = f"sqlite+pysqlite:///{self.path}"
            self.engine = create_engine(url, echo=False, future=True)

    def add_stable_ids(self, core_server: CoreServer,
                       features: List[Feature] = all_features) -> None:
        """Get the stable ids for a list of features from the cores in a server and store them in the db
        
        If no features are provided, the features list from all_features is used
        
        Args:
            core_server: a CoreServer object
            features: a list of features to extract the stable_ids
                      from (gene, transcript, translation)
        """
        
        with self.engine.connect() as conn:
            for core in core_server.cores:
                prod_name = core_server.get_core_metadata(core, 'species.production_name')
                print(f"Load data from {prod_name[0]}")
                
                for feature in features:
                    features = core_server.get_features(core, feature)
                    db_id = self._get_db_id(conn, core, prod_name[0])
                    to_insert = [
                        {
                            'db_id': db_id,
                            'feature': feature.feature,
                            'biotype': feature.biotype,
                            'name': feature.name
                        }
                        for feature in features]
                    conn.execute(insert(StableId), to_insert)
                    conn.commit()

    def _get_db_id(self, conn, db_name: str, prod_name: str) -> None:
        """Get the db_id for a given core from the database
        
        Note:
            Add the core db to the database if it is not already in it.
        
        Args:
            conn: a database connection (made from engine.connect)
            db_name: the core name
            prod_name: the production_name from the core
        """

        stmt_get = select(Db).where(Db.db_name == db_name)
        results = conn.execute(stmt_get)

        rows = []
        for row in results:
            rows.append(row)
        
        if len(rows) > 1:
            raise Exception(f"Several rows in database loaded for db_name {db_name}")
        elif len(rows) == 1:
            row = rows[0]
            db_id = row["db_id"]
        else:
            db_stmt = insert(Db).values(db_name=db_name, production_name=prod_name)
            result = conn.execute(db_stmt)
            conn.commit()
            db_id = result.inserted_primary_key[0]
        
        return db_id

    def get_duplicated_ids(self, feature: str) -> list:
        """Retrieve all duplicated stable_ids of the same feature between different core dbs
        
        Args:
            feature: the feature of the stable_ids to compare

        Returns:
            A list of dicts with 3 keys: db1, db2, name
            ordered by the stable_id names, then db1 and db2
        """
        
        query = text("""SELECT db1.db_name, db2.db_name, s1.name,
                   FROM db db1 LEFT JOIN stable_id s1 ON db1.db_id=s1.db_id,
                        db db2 LEFT JOIN stable_id s2 ON db2.db_id=s2.db_id
                   WHERE s1.name_id != s2.name_id
                        AND s1.feature=s2.feature
                        AND s1.feature = :feature
                        AND s1.name = s2.name
                    ORDER BY s1.name, db1.production_name, db2.production_name
                """)
        
        dup_ids = []
        no_reverse = {}
        with self.engine.connect() as conn:
            for record in conn.execute(query, {'feature': feature}):
                (db1, db2, name) = record
                
                id_hash1 = f"{db1}_{db2}_{name}"
                id_hash2 = f"{db2}_{db1}_{name}"
                if id_hash1 in no_reverse or id_hash2 in no_reverse:
                    continue
                else:
                    no_reverse[id_hash1] = 1
                    no_reverse[id_hash2] = 1
                
                dup_id = {'db1': db1, 'db2': db2, 'name': name}
                dup_ids.append(dup_id)
        
        return dup_ids
        
    def get_duplicates_summary(self, feature: str) -> list:
        """Show a summary of all duplicates for a given feature
        
        Args:
            feature: the feature of the stable_ids to compare

        Returns:
            A list of dicts with 3 keys: db1, db2, count
            ordered by db1, db2
        """
        
        query = text("""SELECT db1.db_name, db2.db_name, count(*)
                   FROM db db1 LEFT JOIN stable_id s1 ON db1.db_id=s1.db_id,
                        db db2 LEFT JOIN stable_id s2 ON db2.db_id=s2.db_id
                   WHERE s1.name_id != s2.name_id
                        AND s1.feature=s2.feature
                        AND s1.feature = :feature
                        AND s1.name = s2.name
                    GROUP BY db1.db_id, db2.db_id
                    ORDER BY db1.production_name, db2.production_name
                """)
        
        dup_ids = []
        no_reverse = {}
        with self.engine.connect() as conn:
            for record in conn.execute(query, {'feature': feature}):
                (db1, db2, count) = record
                
                id_hash1 = f"{db1}_{db2}"
                id_hash2 = f"{db2}_{db1}"
                if id_hash1 in no_reverse or id_hash2 in no_reverse:
                    continue
                else:
                    no_reverse[id_hash1] = 1
                    no_reverse[id_hash2] = 1
                
                dup_id = {'db1': db1, 'db2': db2, 'count': count}
                dup_ids.append(dup_id)
        
        return dup_ids

    def get_duplicated_ids_all_features(self) -> list:
        """Retrieve all duplicated stable_ids between all features between different core dbs

        Returns:
            A list of dicts with 5 keys: db1, db2, feat1, feat2, biotype1, biotype2, name
            ordered by the stable_id names, then db1 and db2
        """
        
        query = text("""SELECT db1.db_name, db2.db_name,
                             s1.feature, s2.feature,
                             s1.biotype, s2.biotype,
                              s1.name
                   FROM db db1 LEFT JOIN stable_id s1 ON db1.db_id=s1.db_id,
                        db db2 LEFT JOIN stable_id s2 ON db2.db_id=s2.db_id
                   WHERE s1.name_id != s2.name_id
                        AND s1.name = s2.name
                    ORDER BY s1.name, db1.production_name, db2.production_name,
                             s1.feature, s2.feature
                """)
        
        dup_ids = []
        no_reverse = {}
        with self.engine.connect() as conn:
            for record in conn.execute(query):
                (db1, db2, feat1, feat2, biotype1, biotype2, name) = record
                
                db_key1 = f"{db1}_{feat1}"
                db_key2 = f"{db2}_{feat2}"
                db_key_list = [db_key1, db_key2]
                db_key_list.sort()
                db_key = "_".join(db_key_list) + f"_{name}"
                if db_key in no_reverse:
                    continue
                else:
                    no_reverse[db_key] = 1
                
                dup_id = {
                    'db1': db1,
                    'db2': db2,
                    'feat1': feat1,
                    'feat2': feat2,
                    'biotype1': biotype1,
                    'biotype2': biotype2,
                    'name': name
                }
                dup_ids.append(dup_id)
        
        return dup_ids
        
    def get_duplicates_summary_all_features(self) -> list:
        """Show a summary of all duplicates between all features

        Returns:
            A list of dicts with 5 keys: db1, db2, feature1, feature2, count
            ordered by db1, db2
        """
        
        query = text("""SELECT db1.db_name, db2.db_name, s1.feature, s2.feature, count(*)
                   FROM db db1 LEFT JOIN stable_id s1 ON db1.db_id=s1.db_id,
                        db db2 LEFT JOIN stable_id s2 ON db2.db_id=s2.db_id
                   WHERE s1.name_id != s2.name_id
                        AND s1.name = s2.name
                    GROUP BY db1.db_id, db2.db_id, s1.feature, s2.feature
                    ORDER BY db1.production_name, db2.production_name, s1.feature, s2.feature
                """)
        
        dup_ids = []
        no_reverse = {}
        with self.engine.connect() as conn:
            for record in conn.execute(query):
                (db1, db2, feat1, feat2, count) = record
                
                # Make a unique key to avoid having both ways db1 -> db2 and db2 -> db1
                db_key1 = f"{db1}_{feat1}"
                db_key2 = f"{db2}_{feat2}"
                db_key_list = [db_key1, db_key2]
                db_key_list.sort()
                db_key = "_".join(db_key_list)
                if db_key in no_reverse:
                    continue
                else:
                    no_reverse[db_key] = 1
                
                dup_id = {'db1': db1, 'db2': db2, 'feat1': feat1, 'feat2': feat2, 'count': count}
                dup_ids.append(dup_id)
        
        return dup_ids


def print_values(columns: tuple, objects: list) -> None:

    print("#" + "\t".join(columns))
    for obj in objects:
        line_vals = [str(obj[key]) for key in columns]
        print("\t".join(line_vals))


def main():
    # Parse command line arguments
    desc = 'Create a database of stable_ids from a list of cores, and check their uniqueness'
    parser = argparse.ArgumentParser(description=desc)
    
    parser.add_argument('--db', type=str, required=True, help='The SQLite db to create or use')

    parser.add_argument('--host', type=str, help='Host of the server to use')
    parser.add_argument('--port', type=str, help='Port of the server to use')
    parser.add_argument('--user', type=str, help='User of the server to use')
    parser.add_argument('--password', type=str, help='Password of the server to use')

    parser.add_argument('--create', action='store_true',
                        help='Create the db (replace and reinit if it exists)')
    parser.add_argument('--add', action='store_true',
                        help='Add ids from the cores to the db')
    parser.add_argument('--build', type=int, help='Filter addition by build number')
    parser.add_argument('--summary', action='store_true',
                        help='Get a summary of the duplicates')
    parser.add_argument('--all_summary', action='store_true',
                        help='Get a summary of the duplicates for all features')
    parser.add_argument('--list_duplicates', action='store_true',
                        help='Show a complete list of duplicated ids')
    parser.add_argument('--all_list_duplicates', action='store_true',
                        help='Show a complete list of duplicated ids between all features')

    parser.add_argument('--prefix', type=str, help='Optional prefix to filter cores to use')
    args = parser.parse_args()
    
    # Init the db (no action yet)
    iddb = StableIdDB(args.db)

    # Create the database (empty)
    if args.create:
        iddb.create()

    # Add stable ids from the cores, to the database
    elif args.add:
        core_server = CoreServer(host=args.host, port=args.port, user=args.user,
                                 password=args.password)
        core_server.get_cores(args.prefix, args.build)

        iddb.connect()
        iddb.add_stable_ids(core_server)

        cores = core_server.cores
        
        if cores:
            print(f"{len(cores)} cores, starting with {cores[0]}")
        else:
            print("No cores")
    
    # Print out all the duplicated stable_ids, and between which pairs of cores
    elif args.list_duplicates:
        iddb.connect()
        dup_ids = iddb.get_duplicated_ids('gene')

        columns = ('db1', 'db2', 'name')
        print_values(columns, dup_ids)
    
    # Print out all the duplicated stable_ids, and between which pairs of cores and all features
    elif args.all_list_duplicates:
        iddb.connect()
        dup_ids = iddb.get_duplicated_ids_all_features()

        columns = ('db1', 'db2', 'biotype1', 'biotype2', 'feat1', 'feat2', 'name')
        print_values(columns, dup_ids)
    
    # Print out a summary count of duplicates between pairs of cores
    elif args.summary:
        iddb.connect()
        dup_ids = iddb.get_duplicates_summary('gene')

        columns = ('db1', 'db2', 'count')
        print_values(columns, dup_ids)

    # Print out a summary count of duplicates between pairs of cores for all features
    elif args.all_summary:
        iddb.connect()
        dup_ids = iddb.get_duplicates_summary_all_features()

        columns = ('db1', 'db2', 'feat1', 'feat2', 'count')
        print_values(columns, dup_ids)
    else:
        print("No action performed")


if __name__ == "__main__":
    main()
