FROM python:3.9-buster

ENV ROOT_DIR=/opt/sashimi
RUN mkdir $ROOT_DIR
COPY ./ $ROOT_DIR

RUN cd $ROOT_DIR && pip install pipenv && pipenv install --system
ENTRYPOINT ["python", "/opt/sashimi/main.py"]