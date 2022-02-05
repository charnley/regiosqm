#!/bin/bash

cd public

VERSION=9.1.0
VERSION=9.4.0

wget https://web.chemdoodle.com/downloads/ChemDoodleWeb-$VERSION.zip

unzip ChemDoodleWeb-$VERSION.zip

mkdir chemdoodleweb

cp -r ChemDoodleWeb-$VERSION/install/* chemdoodleweb/
cp -r ChemDoodleWeb-$VERSION/src/* chemdoodleweb/

rm -r ChemDoodleWeb-$VERSION*
