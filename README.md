## Overview
This application allows you to run specific Python scripts to test algorithms or run a default example. The scripts are designed to work with Docker to provide an isolated environment.

## Dokcer hub link:
https://hub.docker.com/repository/docker/omega111111/discrete_logarithm/general

## Setup and Build
1. **Clone the repository**
2. **Build and run the Docker image:**
```
docker build -t app .
```
## Running the Application
```
docker run -it 
```
3. **Run as simple python3 programm**
Use --test flag to run tests
```
python3 discrete_logarithm.pt --tests
```