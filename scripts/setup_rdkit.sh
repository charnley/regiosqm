#!/bin/bash

cd public

mkdir -p rdkit
cd rdkit

# VERSION=

wget https://unpkg.com/@rdkit/rdkit/Code/MinimalLib/dist/RDKit_minimal.js
wget https://unpkg.com/@rdkit/rdkit/Code/MinimalLib/dist/RDKit_minimal.wasm
