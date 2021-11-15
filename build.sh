#!/bin/bash

# This script installs radexi in $(HOME) directory.

make -j4

mkdir radexi && mkdir radexi/bin && mkdir radexi/data && mkdir radexi/history
cp bin/radexi radexi/bin/radexi
mv ./radexi $HOME/radexi
