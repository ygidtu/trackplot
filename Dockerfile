FROM python:3.9-buster

ENV PYPI=https://pypi.douban.com/simple
ENV ROOT_DIR=/opt/sashimi
RUN mkdir $ROOT_DIR
COPY ./ $ROOT_DIR

RUN cd $ROOT_DIR && pip install -r requirements.txt # -i $PYPI
ENTRYPOINT ["python", "/opt/sashimi/main.py"]