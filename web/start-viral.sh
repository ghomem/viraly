#!/bin/bash

MYDIR=`dirname $0`

cd $MYDIR
bokeh serve --disable-index --allow-websocket-origin=lo.gic.li viral.py viral2.py viral-long.py viral-staging.py viral-marketing.py
