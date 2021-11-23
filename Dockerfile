FROM python:3.7.5

RUN pip install numpy==1.21.4
RUN pip install pandas==1.3.4
RUN pip install scipy==1.7.1
RUN pip install seaborn==0.11.2
RUN pip install sklearn
RUN pip install lpocv

COPY . /app/
WORKDIR /app/src
