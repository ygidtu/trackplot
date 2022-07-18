FROM python:3.9-buster

ENV ROOT_DIR=/opt/sashimi
RUN mkdir $ROOT_DIR
COPY ./ $ROOT_DIR

RUN cd $ROOT_DIR && pip install jupyterlab && pip install -r requirements.txt -i https://pypi.douban.com/simple
ENTRYPOINT ["python", "/opt/sashimi/main.py"]