#!/bin/bash

docker run -e PATIENTDIR --volume $1:/opt/Patient --volume $2:/opt/OSS-DBS -it --rm ningfei/oss-dbs python3 Launcher_OSS_lite.py
