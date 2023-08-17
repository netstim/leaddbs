FROM ningfei/oss-dbs:latest

ARG UNAME=OSS-DBS

ARG UID=1000

ARG GID=1000

# Create user and group
RUN groupadd -g $GID $UNAME && \
    useradd -m -s /bin/bash -u $UID -g $GID $UNAME

# Change owner of /opt/OSS-DBS
RUN chown -R $UID:$GID /opt/OSS-DBS

# Switch user
USER $UNAME
