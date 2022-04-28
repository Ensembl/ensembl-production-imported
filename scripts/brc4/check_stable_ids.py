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
"""

import argparse
import mysql.connector
from mysql.connector.cursor import MySQLCursor
import os, json, re, time
from sqlalchemy import Column, Integer, Index, String, ForeignKey, insert, select
from sqlalchemy.orm import declarative_base
from sqlalchemy import create_engine

####################################################################################################

class CoreServer(object):
    """Interface to a MySQL server with cores in it
    """

    def __init__(self, host: str, port: str, user: str, password: str):
        """Init the database object
        """
        
        self.host = host
        self.port = port
        self.user = user
        self.password = password
        self.engine = None
        self.db = None
        self.cores = []

        self.connect()

    def connect(self) -> None:
        """Create a connection to the server
        """
        
        self.db = mysql.connector.connect(
                user=self.user,
                passwd=self.password,
                host=self.host, 
                port=self.port) 

    def cursor(self) -> MySQLCursor:
        return self.db.cursor()
    
    def get_cores(self, prefix: str) -> None:
        """Retrieve the list of cores on the server
        
        Args:
            prefix: only use cores that start with this prefix.
            
        """

        if not self.db:
            self.connect()
        
        query = "SHOW DATABASES LIKE '%_core_%'"

        cursor = self.cursor()
        cursor.execute(query)
        
        for db in cursor:
            if prefix and not db[0].startswith(prefix): continue
            self.cores.append(db[0])
    
    def get_stable_ids(self, core_name: str, feature: str) -> list:
        """Retrieve all the stable ids for a feature in a given core
        
        Args:
            core_name: name of the core database to use
            feature: feature name (gene, transcript, or translation)
        
        Returns:
            A list of stable_ids strings
        """
        
        self.db.database = core_name
        cursor = self.db.cursor()
        
        query = f"SELECT stable_id FROM {feature}"
        cursor.execute(query)
        
        for row in cursor:
            yield row[0]
    
    def get_core_metadata(self, core_name: str, key: str) -> list:
        """Retrieve a metadata value from a given core
        
        Args:
            core_name: name of the core database to use
            key: metadata_key to get the metadata_value from
        
        Returns:
            List of all values
        """
        
        self.db.database = core_name
        cursor = self.db.cursor()
        
        query = f"SELECT meta_value FROM meta WHERE meta_key = %s"
        cursor.execute(query, [key])
        
        rows = []
        for row in cursor:
            rows.append(row[0])
        
        return rows

####################################################################################################
# Prepare the SQLalchemy schema
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
    Index('ix_feat_name', feature, name)

class StableIdDB(object):
    """Representation of an SQLite database of stable ids
    """
    
    def __init__(self, path: str, replace=False):
        """Init the database object
        """
        self.path = path
        self.engine = None
    
    def create(self) -> None:
        """Create the SQLite database
        
        Note:
            Recreate the file if it already exists
        """
        
        if os.path.exists(self.path):
            os.remove(self.path)
        
        url = f"sqlite+pysqlite:///{self.path}"
        self.engine = create_engine(url, echo=False, future=True)
        Base.metadata.create_all(self.engine)
    
    def connect(self) -> None:
        """Connect to the SQLite database
        """

        if not self.engine:
            url = f"sqlite+pysqlite:///{self.path}"
            self.engine = create_engine(url, echo=False, future=True)
    
    def add_features(self, core_server: CoreServer, feature: str) -> None:
        """Get the stable ids from the cores in a server and store them in the db
        """
        with self.engine.connect() as conn:
            for core in core_server.cores:
                prod_name = core_server.get_core_metadata(core, 'species.production_name')
                print(f"Load {feature} ids from {core} ({prod_name[0]})")
                
                stable_ids = core_server.get_stable_ids(core, feature)
                db_id = self._get_db_id(conn, core, prod_name[0])
                to_insert = [ { 'db_id': db_id, 'feature': feature, 'name': stable_id } for stable_id in stable_ids ]
                result = conn.execute(insert(StableId), to_insert)
                conn.commit()

    def _get_db_id(self, conn, db_name: str, prod_name: str) -> None:

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
        

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
            description='Create a database of stable_ids from a list of cores, and check their uniqueness')
    
    parser.add_argument('--db', type=str, required=True, help='The SQLite db to create or use')

    parser.add_argument('--host', type=str, help='Host of the server to use')
    parser.add_argument('--port', type=str, help='Port of the server to use')
    parser.add_argument('--user', type=str, help='User of the server to use')
    parser.add_argument('--password', type=str, help='Password of the server to use')

    parser.add_argument('--create', action='store_true', help='Create the db')
    parser.add_argument('--add', action='store_true', help='Add to the db')
    parser.add_argument('--summary', action='store_true', help='Get a summary of the db')

    parser.add_argument('--prefix', type=str, help='Optional prefix to filter cores to use')
    args = parser.parse_args()
    
    # Choose which data to retrieve
    iddb = StableIdDB(args.db)
    if args.create:
        iddb.create()
    elif args.add:
        core_server = CoreServer(host=args.host, port=args.port, user=args.user, password=args.password)
        core_server.get_cores(args.prefix)

        iddb.connect()
        iddb.add_features(core_server, 'gene')
        iddb.add_features(core_server, 'transcript')
        iddb.add_features(core_server, 'translation')

        cores = core_server.cores
        
        if cores:
            print(f"{len(cores)} cores, starting with {cores[0]}")
        else:
            print(f"No cores")
    elif args.summary:
        pass
    else:
        print("No action performed")

if __name__ == "__main__":
    main()

