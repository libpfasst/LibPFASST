# In order to update the image that gets used in the CircleCI testing do:
# 1. Create a new hub.docker.com account to push the newly built docker image to
#    e.g. a libpfasst dockerhub account
# 2. `docker build -t libpfasst/libpfasst:latest .`
#    this will create a new image on your compute
# 3. `docker login`
#    using the credentials from the newly created account
# 4. `docker push libpfasst/libpfasst`
#    this will put the image in the docker hub account so it can be pulled
# Finally, go and update the docker image used in the .circleci/config.yml file

FROM ubuntu:xenial
RUN set -ex \
  && apt-get -y update \
  && apt-get -y install \
         gfortran \
         mpich \
         libmpich-dev \
         libopenblas-dev \
         liblapack-dev \
         python-nose \
         wget \
         make \
         python-tk \
         python-pip \
         python-dev 
RUN set -ex \
  && pip install  \
        pytest 
#        tqdm \
#        numpy \
#        pandas \
#        attrs \
#        scipy \
#        'matplotlib<1.6.0'
