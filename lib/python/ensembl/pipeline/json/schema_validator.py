#!env python3

import eHive
import json
from jsonschema import validate

class schema_validator(eHive.BaseRunnable):

    def run(self):
        json_file = self.param_required("json_file")
        json_schema_file = self.param_required("json_schema")
        
        json_data = self.get_json(json_file)
        schema = self.get_json(json_schema_file)
        
        validate(instance=json_data, schema=schema)

    def get_json(self, json_path):
        with open(json_path) as json_file:
            data = json.load(json_file)
            return data

