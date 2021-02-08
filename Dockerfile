FROM ubuntu:20.04

MAINTAINER Your Name "homem.gustavo@gmail.com"

ENV MYORIGIN=lo.gic.li

RUN apt-get update -y && \
    DEBIAN_FRONTEND='noninteractive' apt-get install -y python3-numpy python3-scipy python3-pip python3-jinja2 python3-markupsafe python3-dateutil python3-tornado python3-typing-extensions python3-matplotlib python3-distutils

WORKDIR /app

RUN pip3 install bokeh

COPY ./viraly.py                /app
COPY ./web/viral.py             /app
COPY ./web/viral2.py            /app
COPY ./web/viral-long.py        /app
COPY ./web/viral-staging.py     /app
COPY ./web/viral-marketing.py   /app
COPY ./web/viral-simple.py      /app

EXPOSE 5006

# docker run  -e "MYORIGIN=yourhostname.net" -v -d -p 5006:5006 viraly

CMD bokeh serve --disable-index --allow-websocket-origin=localhost:5006 --allow-websocket-origin=$MYORIGIN viral.py viral2.py viral-long.py viral-staging.py viral-marketing.py viral-simple.py
