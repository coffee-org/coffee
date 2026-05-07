#!/bin/bash
# 
# Pre-allocate placeholder streams for coronagraph pipeline
# 
echo "Setting up coronagraph environment streams..."

# Create placeholder shared-memory images (1024x1024 floats)
milk-exec "mem.mk2Dim apostart 1024 1024"
milk-exec "mem.mk2Dim pupa0 1024 1024"

echo "Placeholder streams created successfully."
