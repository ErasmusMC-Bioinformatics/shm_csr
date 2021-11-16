#!/usr/bin/env bash

set -e
DEFAULT_BASE_IMAGE="$1"

python3 create_container_hash.py $DEFAULT_BASE_IMAGE
mulled-build-files --namespace rhpvorderman build-and-test ./container_hash.tsv --verbose

