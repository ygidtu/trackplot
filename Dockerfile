FROM python:3.10-buster

ENV ROOT_DIR=/opt/sashimi
RUN mkdir $ROOT_DIR
COPY ./ $ROOT_DIR

RUN cd $ROOT_DIR && pip install -r requirements.txt
ENTRYPOINT ["python", "/opt/sashimi/main.py"]