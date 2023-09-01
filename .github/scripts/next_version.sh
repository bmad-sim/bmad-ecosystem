#!/bin/bash

COUNTER=0
# Current date as YYYYMMDD
DATE=$(date +%Y%m%d)

# Check last git tag
LAST_TAG=$(git describe --abbrev=0 --tags)

# echo "Last Tag was: $LAST_TAG"
# Check if last tag contains the current date
if [[ $LAST_TAG == *$DATE* ]]; then
    # Get the last number of the tag
    COUNTER=$(echo $LAST_TAG | grep -o '[0-9]*$')
    # Increment the counter
    COUNTER=$((COUNTER+1))
fi

# Assemble the new tag
NEW_TAG="$DATE-$COUNTER"
# Print out the new tag to stdout
echo "$NEW_TAG"
