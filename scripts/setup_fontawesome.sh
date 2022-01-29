#!/bin/bash

cd public

VERSION=5.15.4

wget https://use.fontawesome.com/releases/v$VERSION/fontawesome-free-$VERSION-web.zip

unzip fontawesome-free-$VERSION-web.zip

rm -r fontawesome
mv fontawesome-free-$VERSION-web fontawesome
rm fontawesome-free-$VERSION-web.zip
