#!/bin/bash

COUNTER=0
# Current date as YYYYMMDD
DATE=$(date +%Y%m%d)

# Check last git tag
LAST_TAG=$(git tag -l | tail -n 1)

# Split the last tag into the date portion and sequence
LAST_TAG_DATE=$(echo $LAST_TAG | grep -o '^[^-]*')
LAST_TAG_SEQUENCE=$(echo $LAST_TAG | grep -o '[0-9]*$')

# Compare LAST_TAG_DATE with current date if it is greater than or equal to the current date, increment the sequence
# The greater date can happen in the case of UTC timezones
if [[ $LAST_TAG_DATE -ge $DATE ]]; then
    COUNTER=$((LAST_TAG_SEQUENCE+1))
    DATE=$LAST_TAG_DATE
else
    COUNTER=0
fi

# Assemble the new tag
NEW_TAG="$DATE-$COUNTER"
# Print out the new tag to stdout
echo "$NEW_TAG"
