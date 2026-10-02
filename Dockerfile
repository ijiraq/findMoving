# Pinned to an Ubuntu 24.04 build; :latest moved to 26.04, which changes the iraf and ds9 packages.
FROM quay.io/jupyter/base-notebook@sha256:23fbe16371af3c0fdc833fa0962b42fd07811e152b821ff51d19a8d6e9cc86cf AS deploy
USER root
RUN apt upgrade -y
RUN apt update -y

# system settings and permissions
COPY canfar_src/nofiles.conf /etc/security/limits.d/
COPY canfar_src/nsswitch.conf /etc/
RUN apt-get -y install	sssd-ad
RUN apt-get -y install	sssd-tools
# put the standard start up script into place
COPY canfar_src/startup.sh /skaha/startup.sh
# put the initialization of IRAF into the global setup
COPY canfar_src/iraf.sh /etc/profile.d/
# see https://bugzilla.redhat.com/show_bug.cgi?id=1773148
RUN touch /etc/sudo.conf && echo "Set disable_coredump false" > /etc/sudo.conf

# these are some extra linux things (beyond jupyter needs) we find useful
RUN apt-get install -yq gcc git libx11-dev libxt-dev libcfitsio-dev
RUN apt-get install -yq gfortran emacs pip vim adcli parallel

# note that this container can also run a an xterm on X11 Desktop
RUN apt-get install -yq curl xterm xrdp iraf iraf-dev
# complete th configuration of iraf
RUN ln -s /usr/lib/iraf/bin /usr/lib/iraf/bin.linux
RUN ln -s /usr/lib/iraf/noao/bin /usr/lib/iraf/noao/bin.linux
RUN ln -s /usr/lib/iraf/unix/bin /usr/lib/iraf/unix/bin.linux

# get a good version of ds9
RUN curl -fL https://ds9.si.edu/download/ubuntu-24.04_x86_64/xpa.ubuntu-24.04_x86_64.2.2.1.tar.gz | tar -C /usr/bin -xzf -
RUN curl -fL https://ds9.si.edu/download/ubuntu-24.04_x86_64/ds9.ubuntu-24.04_x86_64.8.8.tar.gz | tar -C /usr/bin -zxf -

# Build python packages in a specialty venv for this project
ARG VENV=/opt/findMoving/astropy
RUN python -m venv ${VENV}
RUN . ${VENV}/bin/activate && python -m pip install canfar cadctap cadcdata vos 
RUN . ${VENV}/bin/activate && python -m pip install matplotlib ephem
RUN . ${VENV}/bin/activate && python -m pip install mp_ephem 
RUN . ${VENV}/bin/activate && python -m pip install ccdproc pyraf pyds9

# and scripts that will initlize the project environment
RUN echo ". ${VENV}/bin/activate" > /etc/profile.d/activate_python_venv.sh
ARG BUILDDIR=/opt/findMoving
RUN mkdir -p ${BUILDDIR}
WORKDIR ${BUILDDIR}
COPY src  ./
ARG iraf=/usr/lib/iraf
ARG IRAFARCH=linux
ARG USER=`whoami`
# Install the project software
RUN . /opt/findMoving/astropy/bin/activate && python -m pip install -r requirements.txt
RUN . /opt/findMoving/astropy/bin/activate && python setup.py install 
ENTRYPOINT ["/skaha/startup.sh"]
