# Use the Gurobi Docker image as a base
FROM gurobi/python

# Install any additional dependencies for dingo
RUN apt-get update && apt-get install -y \
    cmake \
    lp-solve \
    && rm -rf /var/lib/apt/lists/*

# Install dingo
RUN pip install dingo

# Set the working directory
WORKDIR /app

# Copy the current directory contents into the container at /app
ADD . /app

# Run dingo when the container launches
CMD ["python", "setup.py"]