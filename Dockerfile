FROM nvidia/cuda:9.1-devel 

RUN apt-get update && apt-get -y install noweb

RUN apt-get -y install libcurl4-openssl-dev

RUN apt-get -y install git

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH

RUN apt-get install -y wget bzip2 ca-certificates

RUN wget --quiet https://repo.continuum.io/archive/Anaconda3-5.1.0-Linux-x86_64.sh -O ~/anaconda.sh && \
        	 /bin/bash ~/anaconda.sh -b -p /opt/conda && \				      
         	    rm ~/anaconda.sh && \
		    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
	    	    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
		        echo "conda activate base" >> ~/.bashrc

RUN conda install --yes boto3

COPY . /tmp/src

RUN cd /tmp/src && git checkout devel_cuda9.1_python3.6 && make all cython
