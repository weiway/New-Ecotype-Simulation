#!/bin/sh
# @param 1 The Phylip program to run.
# @param 2 The directory containing the infile.
# @param 3 The debug option.

# Change to the directory containing the infile.
cd "$2"

# Remove the outfile and outtree files if they exist.
if [ -f outfile ]; then
  rm outfile
fi
if [ -f outtree ]; then
  rm outtree
fi

# Check to see if this is running on Debian/Ubuntu.
if [ -f /etc/debian_version ]; then
  COMMAND="phylip $1"
else
  COMMAND="$1"
fi

# Look for the debug option.
if [ "$3" = "true" ]; then
  # Run the Phylip program with debugging on.
  cat input | $COMMAND
else
  # Run the Phylip program with debugging off.
  cat input | $COMMAND > /dev/null
fi
