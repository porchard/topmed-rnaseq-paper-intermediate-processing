#!/usr/bin/env python
# coding: utf-8


import json
import argparse
import logging
import os

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s')

parser = argparse.ArgumentParser()
parser.add_argument('--json', type=str, help='JSON file to stage')
parser.add_argument('--stage-using', type=str, nargs='+', default=[])
args = parser.parse_args()

JSON = args.json
STAGE_USING = args.stage_using


logging.info('Loading {}'.format(JSON))

with open(JSON) as f:
    data = json.load(f)


logging.info('Determining staging locations')
staging_locations = {}

for key, value in data.items():
    new_location = [os.getcwd()] + [value[i] for i in STAGE_USING if i in list(value.keys())] + [os.path.basename(key)]
    staging_locations[key] = os.path.join(*new_location)


unique_staging_locations = set(staging_locations.values())
if len(unique_staging_locations) != len(staging_locations):
    raise ValueError('Staging locations are not unique!')



# now write the files
logging.info('Writing new json')
new_data = {staging_locations[k]: v for k, v in data.items()}
print(json.dumps(new_data, indent=4, sort_keys=True))


# cp each file to the new staging location
logging.info('Staging')
to_link = len(data)
count = 0
for key, value in data.items():
    count += 1
    if count % 10000 == 0:
        logging.info('Copying file {}/{}'.format(count, to_link))
    new_location = staging_locations[key]
    # make sure the directory exists
    new_location_dir = os.path.dirname(new_location)
    if not os.path.exists(new_location_dir):
        os.makedirs(new_location_dir)
    # make sure there's not already a file here
    if os.path.exists(new_location):
        raise ValueError('File already exists at {}'.format(new_location))
    os.system('cp {} {}'.format(key, new_location))
logging.info('Done')
