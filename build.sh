#!/bin/bash

# This script installs radexi in $(HOME) directory.

mkdir -v bin
make -j4
if [[ "${?}" -eq 0 ]]
then
  # Make new tree for future directory
  mkdir -v radexi radexi/data radexi/history
  mv -v -f bin/radexi radexi/radexi

  # Don't overwrite user's saved molecules
  if [[ -d "$HOME/radexi" ]]
  then
    mv -vf ./radexi/radexi $HOME/radexi/
  else
    mv -vf ./radexi $HOME/radexi
  fi

  rm -rv bin radexi
else
  echo "Build error. Check the dependencies on https://github.com/doedre/radexi-cmd."
fi
