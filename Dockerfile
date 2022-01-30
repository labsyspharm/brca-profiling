FROM python:3.9

RUN pip install seaborn
RUN pip install lpocv

COPY . /app/
