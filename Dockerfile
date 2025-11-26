# Dockerfile for LidarForFuel Python package
FROM mambaorg/micromamba:latest

# Set working directory
WORKDIR /app

# Copy environment file
COPY environment.yml /tmp/environment.yml

# Create conda environment from environment.yml
RUN micromamba env create -f /tmp/environment.yml && \
    micromamba clean -afy

# Set environment variables
ENV MAMBA_DEFAULT_ENV=lidarforfuel
ENV PATH="/opt/conda/envs/lidarforfuel/bin:$PATH"

# Copy Python package
COPY python/ /app/python/
COPY test/ /app/test/

# Set PYTHONPATH to include the app directory
ENV PYTHONPATH="/app:$PYTHONPATH"

# Make sure we're using the conda environment
SHELL ["micromamba", "run", "-n", "lidarforfuel", "/bin/bash", "-c"]

# Default command
CMD ["python", "-c", "import sys; print('LidarForFuel Python package ready'); print(f'Python: {sys.version}')"]

