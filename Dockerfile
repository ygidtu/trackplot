FROM python:3.11-buster

ENV ROOT_DIR=/opt/trackplot
RUN mkdir $ROOT_DIR
COPY ./ $ROOT_DIR

RUN cd $ROOT_DIR && pip install -e .
ENTRYPOINT ["python", "/opt/trackplot/main.py"]