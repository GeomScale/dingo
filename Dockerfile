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

# Get PySQR
RUN apt-get install libsuitesparse-dev

# Install Python dependencies
RUN pip install matplotlib \
	plotly \
	networkx \
	pyoptinterface[highs]

# Get dingo
WORKDIR /workspaces/dingo
COPY . .

# Get submodules
RUN git submodule update --init

# Get lp-solve
RUN wget https://sourceforge.net/projects/lpsolve/files/lpsolve/5.5.2.11/lp_solve_5.5.2.11_source.tar.gz &&\
	tar xzvf lp_solve_5.5.2.11_source.tar.gz &&\
	rm lp_solve_5.5.2.11_source.tar.gz

# Get boost library
RUN wget -O boost_1_76_0.tar.bz2 https://archives.boost.io/release/1.76.0/source/boost_1_76_0.tar.bz2 &&\
	tar xjf boost_1_76_0.tar.bz2 &&\
	rm boost_1_76_0.tar.bz2

# Set environmental variable gurobi license path
ENV GRB_LICENSE_FILE=/opt/gurobi/gurobi.lic

# Install dingo
RUN ["python", "setup.py", "install", "--user"]
