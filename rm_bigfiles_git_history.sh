#!/bin/bash

# Find large files in history
# https://stackoverflow.com/questions/10622179/how-to-find-identify-large-commits-in-git-history
# git rev-list --objects --all |
#  git cat-file --batch-check='%(objecttype) %(objectname) %(objectsize) %(rest)' |
#  sed -n 's/^blob //p' |
#  sort --numeric-sort --key=2 |
#  cut -c 1-12,41- |
#  $(command -v gnumfmt || echo numfmt) --field=2 --to=iec-i --suffix=B --padding=7 --round=nearest


# Found files
# 9496bbaa9f4c   53MiB shinyApp/data/analysis4_DATASETS_COMBINED_log15_light_new_ora.rds
# 1d4124ef055f   58MiB shinyApp/data/a4_data_diffcom_all_filtered.rds
# e6d7641768a9   83MiB GEOmetadb.sqlite

# FILE="GEOmetadb.sqlite"
FILE="shinyApp/data/a4_data_diffcom_all_filtered.rds"

# Removal command
# https://stackoverflow.com/questions/2100907/how-to-remove-delete-a-large-file-from-commit-history-in-the-git-repository
git filter-branch --prune-empty -d /dev/shm/scratch \
  --index-filter "git rm --cached -f --ignore-unmatch $FILE" \
  --tag-name-filter cat -- --all
