FROM ubuntu
RUN apt-get -y update 
RUN apt-get -y install gfortran && apt-get -y install mpich libmpich-dev libopenblas-dev liblapack-dev python-nose wget
RUN wget https://bootstrap.pypa.io/get-pip.py && python get-pip.py
RUN pip install numpy pandas attrs scipy
RUN pip install 'matplotlib<1.6.0'
