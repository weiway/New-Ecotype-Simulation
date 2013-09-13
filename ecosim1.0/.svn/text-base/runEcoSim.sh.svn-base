#!/bin/sh

# Find the absolute path to this script, and change to that directory.
SCRIPT_NAME=$(echo "$0" | sed "s/.*\\///g")
if [ "$0" = "./$SCRIPT_NAME" ]; then
  ABS_PATH=$(echo `pwd`/$0 | sed "s/\/$SCRIPT_NAME//g" | sed "s/\/\.//g")
else
  ABS_PATH=$(echo "$0" | sed "s/\/$SCRIPT_NAME//g")
fi
if [ "$ABS_PATH" != "$SCRIPT_NAME" ]; then
  cd "$ABS_PATH"
fi
# Check if the debug argument was passed to this script.
DEBUG=0
for ARG in "$@"; do
  if [ "$ARG" = "--debug" ] || [ "$ARG" = "-d" ]; then
    DEBUG=1
  fi
done
# Run Ecotype Simulation, display debug window if requested.
if [ $DEBUG = "1" ]; then
  if [ -f /usr/bin/gnome-terminal ]; then
    gnome-terminal -t "Ecotype Simulation" -x java -jar bin/EcoSim.jar $@ &
  else
    xterm -T "Ecotype Simulation" -e java -jar bin/EcoSim.jar $@ &
  fi
else
  java -jar bin/EcoSim.jar $@ &
fi
