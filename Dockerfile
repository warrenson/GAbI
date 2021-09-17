#from tensorflow/tensorflow:2.4.2-gpu-jupyter
FROM tensorflow/tensorflow:2.4.2-jupyter

COPY requirements.txt /tmp
RUN pip install -r /tmp/requirements.txt

