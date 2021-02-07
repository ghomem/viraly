#!/bin/bash

MYDIR=`dirname $0`
MYORIGIN="$1"

cd $MYDIR
bokeh serve --disable-index --allow-websocket-origin=$MYORIGIN viral.py viral2.py viral-long.py viral-staging.py viral-marketing.py viral-simple.py hw.py
