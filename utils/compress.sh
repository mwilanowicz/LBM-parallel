#!/bin/bash

find ../outputs -type f ! -name "*.xz" ! -name "metadata.csv" ! -name "config.csv" ! -name ".DS_Store" -exec xz -v {} +
