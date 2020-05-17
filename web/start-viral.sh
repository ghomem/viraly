#!/bin/bash

MYDIR=`dirname $0`

cd $MYDIR
bokeh serve --allow-websocket-origin=lo.gic.li  viral.py
