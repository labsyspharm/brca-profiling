FROM python:3.9

RUN apt-get update -y
RUN apt-get upgrade -y

RUN pip install seaborn
RUN pip install lpocv

COPY . /app/
