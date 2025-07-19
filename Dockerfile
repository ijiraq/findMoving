FROM ubuntu:latest AS deploy
ARG DEBIAN_FRONTEND=noninteractive
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
# see https://bugzilla.redhat.com/show_bug.cgi?id=1773148
RUN touch /etc/sudo.conf && echo "Set disable_coredump false" > /etc/sudo.conf
RUN apt install -y jupyter-notebook


## see https://bugzilla.redhat.com/show_bug.cgi?id=1773148
RUN apt-get install -yq curl xterm xrdp vim adcli parallel iraf
RUN apt-get install -yq gcc git libx11-dev iraf-dev libxt-dev libcfitsio-dev
RUN apt-get install -yq gfortran emacs pip
# put the initialization of IRAF into the global setup

# get a good version of ds9
RUN curl -L https://ds9.si.edu/download/ubuntu22x86/xpa.ubuntu22x86.2.1.20.tar.gz | tar -C /usr/bin -zxf -
# RUN curl -L https://ds9.si.edu/download/ubuntu22x86/xpa.ubuntu22x86.2.1.20.tar.gz | tar -C /usr/bin -xzf -
RUN curl -L https://ds9.si.edu/download/ubuntu22x86/ds9.ubuntu22x86.8.7b1.tar.gz | tar -C /usr/bin -xzf - 
# RUN curl https://ds9.si.edu/download/ubuntu22x86/ds9.ubuntu22x86.8.6b1.tar.gz  | tar -C /usr/bin -xzf -
WORKDIR /opt
RUN apt-get install -yq python3-dev python3-numpy-dev python3-setuptools cython3 python3-pytest-astropy
RUN apt-get install -yq python3-wxgtk4.0
RUN apt-get install -y python3.12-venv
RUN python3 -m venv /opt/findMoving/astropy
RUN /opt/findMoving/astropy/bin/python3 -m pip install 'astropy>=5.1.0,<6.0.0'
# RUN apt-get install -qy python3-pyraf
COPY canfar_src/iraf.sh /etc/profile.d/
RUN ln -s /usr/lib/iraf/bin /usr/lib/iraf/bin.linux
RUN ln -s /usr/lib/iraf/noao/bin /usr/lib/iraf/noao/bin.linux
RUN ln -s /usr/lib/iraf/unix/bin /usr/lib/iraf/unix/bin.linux

RUN /opt/findMoving/astropy/bin/python -m pip install ccdproc pyraf vos matplotlib ephem pyds9 mp_ephem 
# RUN pip install ossos
COPY canfar_src/findMoving.sh /etc/profile.d/
ARG BUILDDIR=/opt/findMoving
RUN mkdir -p ${BUILDDIR}
# COPY ds9_dist/ds9.unknown.8.3.tar.gz ${BUILDDIR}/
# RUN tar xf ${BUILDDIR}/ds9.unknown.8.3.tar.gz ; mv ds9 /usr/bin/
# COPY ds9_dist/xpa.unknown.2.1.20.tar.gz ${BUILDDIR}/
# RUN tar xf ${BUILDDIR}/xpa.unknown.2.1.20.tar.gz ; mv xpa* /usr/bin/
WORKDIR ${BUILDDIR}
COPY src  ./
ARG iraf=/usr/lib/iraf
ARG IRAFARCH=linux
ARG USER=`whoami`
RUN /opt/findMoving/astropy/bin/python3 -m pip install -r requirements.txt
RUN /opt/findMoving/astropy/bin/python3 setup.py install 
# RUN apt-get install -y automake autoconf libx11-dev zlib1g-dev libxml2-dev libxslt1-dev libxft-dev tcl-dev tk-dev zip
# RUN mkdir ds9
# WORKDIR ${BUILDDIR}/ds9
# RUN curl -L https://github.com/SAOImageDS9/SAOImageDS9/archive/refs/tags/v8.3.tar.gz | tar xzf - 
ENTRYPOINT ["/skaha/startup.sh"]
