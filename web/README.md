This is a bokeh frontend for the viraly script. It can be launched with:

```
#!/bin/bash

MYDIR=`dirname $0`

cd $MYDIR
bokeh serve --allow-websocket-origin=lo.gic.li  viral.py
```

For production use an nginx (or equivalent reverse proxy) should be put in from of the web application.
