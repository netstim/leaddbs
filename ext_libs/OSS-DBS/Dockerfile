# Dockerfile to build the OSS DBS environment
#
# Author:
# Max Schroeder <max.schroeder@uni-rostock.de>


###############################################################################
# Stage 1
FROM sfbelaine/oss_dbs:python_latest as user_creation

ARG OSS_UID="1000"
ARG OSS_GID="1000"

ENV OSS_UID=$OSS_UID \
    OSS_GID=$OSS_GID

# Create a new user and group in order to run processes as non-root
RUN groupadd -g ${OSS_GID} OSS-DBS
RUN useradd -m -s /bin/bash -g ${OSS_GID} -u ${OSS_UID} OSS-DBS

RUN chown -R ${OSS_UID}:${OSS_GID} /opt/OSS-DBS

# switch user
USER $OSS_UID

###############################################################################
# Stage 2
FROM user_creation as app

# Copy all from the current source code repository
COPY --chown=OSS-DBS:OSS-DBS . /opt/OSS-DBS

WORKDIR /opt/OSS-DBS/OSS_platform

CMD [ "/bin/bash" ]
