This is a bokeh frontend for the viraly script. It can be launched with:

```
#!/bin/bash

MYDIR=`dirname $0`

cd $MYDIR
bokeh serve --allow-websocket-origin=lo.gic.li viral.py
```

For production use an nginx (or equivalent reverse proxy) should be put in front of the web application.

Q: Why is there a viral-staging.py file here? Do I not know that a branch could be used for staging?

A: It is because it is because it is practical to run both staging and production in parallel on the same bokeh server using two nginx locations, the same repo and the same bokeh process.

```
bokeh serve --allow-websocket-origin=lo.gic.li viral.py viral-staging.py
```
