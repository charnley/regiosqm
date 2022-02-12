#!/bin/bash

SCRIPT=`realpath $0`
SCRIPTPATH=`dirname $SCRIPT`

wget https://raw.githubusercontent.com/neverpanic/google-font-download/master/google-font-download -O $SCRIPTPATH/download_google_font.sh
