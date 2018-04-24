FROM ubuntu:trusty
RUN apt-get -y update 
RUN apt-get -y install gfortran && apt-get -y install mpich libmpich-dev libopenblas-dev liblapack-dev python-nose wget make python-tk
RUN wget https://bootstrap.pypa.io/get-pip.py && python get-pip.py
RUN pip install tqdm pytest numpy pandas attrs scipy 'matplotlib<1.6.0'
