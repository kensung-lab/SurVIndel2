FROM ubuntu

WORKDIR /home/

COPY htslib-1.18.tar.bz2 .
COPY build_htslib.sh .

RUN apt-get update --fix-missing
RUN DEBIAN_FRONTEND="noninteractive" apt-get -y install tzdata
RUN apt-get install -y unzip build-essential zlib1g-dev autoconf libbz2-dev liblzma-dev libcurl4-openssl-dev cmake libssl-dev

RUN ./build_htslib.sh

COPY CMakeLists.txt ./
COPY *.h ./
COPY *.cpp ./
ADD libs ./libs

RUN cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo -DNATIVE=ON . && make

RUN apt-get install -y python python-dev curl
RUN curl https://bootstrap.pypa.io/pip/2.7/get-pip.py --output get-pip.py
RUN python get-pip.py

RUN pip install pysam pyfaidx numpy scikit-learn

COPY random_pos_generator.py survindel2.py ./

RUN /bin/echo "#!/bin/bash" > run.sh
RUN /bin/echo "" >> run.sh
RUN /bin/echo "python survindel2.py --threads 8 input/\$1 workdir reference/\$2" >> run.sh
RUN chmod a+x run.sh

ENTRYPOINT ["./run.sh"]
