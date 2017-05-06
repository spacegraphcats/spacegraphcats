#!/bin/bash

set -xe

make acido-search
make 15genome-search

make acido-search > scripts/acidio-search.log
make 15genome-search > scripts/15genome-search.log


if [ -z "$(git status --porcelain)" ]; then
  echo "All tracked files are committed. \n"
else
  echo "Results have changed unexpectedly. Please commit changes to search results. \n"
  git status
  exit 1
fi
