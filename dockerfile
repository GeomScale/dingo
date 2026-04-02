# Use the Gurobi Docker image as a base
FROM gurobi/python

# Install any additional dependencies for dingo
RUN apt-get update && apt-get install -y \
	cmake \
    	lp-solve \
	git \
	wget \
	vim \
	bzip2 \
	g++ \
	&& rm -rf /var/lib/apt/lists/*

# Install dependencies
RUN apt-get update && apt-get install -y libsuitesparse-dev
RUN pip install sparseqr \
		Cython \
		cobra \ 
		kaleido

WORKDIR /app

# Get dingo
RUN git clone https://github.com/GeomScale/dingo.git && \
	cd dingo && \
	git submodule update --init

# Get boost library
WORKDIR /app/dingo
RUN wget -O boost_1_76_0.tar.bz2 https://boostorg.jfrog.io/artifactory/main/release/1.76.0/source/boost_1_76_0.tar.bz2 && \
	tar xjf boost_1_76_0.tar.bz2 && \
	rm boost_1_76_0.tar.bz2

# Get PySQR
RUN apt-get install libsuitesparse-dev

# Install dingo
RUN pip install matplotlib plotly
RUN ["python", "setup.py", "install", "--user"]