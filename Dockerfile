FROM nvidia/cuda:9.1-devel 

RUN apt-get update && apt-get -y install noweb

RUN apt-get -y install libcurl4-openssl-dev

RUN apt-get -y install git

ENV PATH /home/anaconda/bin:$PATH

RUN apt-get install -y wget bzip2 ca-certificates

RUN wget --quiet https://repo.continuum.io/archive/Anaconda3-5.1.0-Linux-x86_64.sh -O ~/anaconda.sh && \
                 /bin/bash ~/anaconda.sh -b -p /home/anaconda && \                                
                    rm ~/anaconda.sh && \
                    ln -s /home/anaconda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
                    echo ". /home/anaconda/etc/profile.d/conda.sh" >> ~/.bashrc && \
                        echo "conda activate base" >> ~/.bashrc

RUN conda install --yes boto3

RUN pip install awscli --upgrade --user

ENV HOME /home
ENV PYTHONPATH /home/CEO/python

RUN cd /home && git clone -b devel_cuda9.1_python3.6 https://github.com/rconan/CEO.git
RUN cd /home/CEO &&  make all cython

ADD docker_test /home/docker_test

ENTRYPOINT ["/home/anaconda/bin/python"]
