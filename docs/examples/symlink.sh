#!/bin/bash

# Define source and destination directories
SRC_DIR="../../notebooks"
DEST_DIR="."

# Ensure the destination directory exists
mkdir -p "$DEST_DIR"

# Find and symlink all ex*.ipynb files
for notebook in "$SRC_DIR"/ex*.ipynb; do
    if [ -f "$notebook" ]; then
        ln -sf "$notebook" "$DEST_DIR/$(basename "$notebook")"
        echo "Symlinked: $notebook -> $DEST_DIR/$(basename "$notebook")"
        # echo $notebook
    fi
done
