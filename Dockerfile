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


## see https://bugzilla.redhat.com/show_bug.cgi?id=1773148
RUN apt-get install -yq curl xterm xrdp vim adcli parallel iraf
RUN apt-get install -yq gcc git libx11-dev iraf-dev libxt-dev libcfitsio-dev
RUN apt-get install -yq gfortran emacs xpa-tools pip
# put the initialization of IRAF into the global setup

# get a good version of ds9
WORKDIR /opt
# RUN curl https://ds9.si.edu/download/ubuntu22x86/ds9.ubuntu22x86.8.5.tar.gz  | tar xzf - ; mv ds9 /usr/bin/ds9

RUN apt-get install -yq python3-wxgtk4.0
RUN pip install 'astropy>=5.1.0,<6.0.0'
RUN apt-get install -qy python3-pyraf
COPY canfar_src/iraf.sh /etc/profile.d/
RUN ln -s /usr/lib/iraf/bin /usr/lib/iraf/bin.linux
RUN ln -s /usr/lib/iraf/noao/bin /usr/lib/iraf/noao/bin.linux
RUN ln -s /usr/lib/iraf/unix/bin /usr/lib/iraf/unix/bin.linux

ARG BUILDDIR=/opt/findMoving
RUN mkdir -p ${BUILDDIR}
WORKDIR ${BUILDDIR}
COPY src  ./
RUN python3 -m pip install -r requirements.txt
RUN python3 setup.py install 

CMD /skaha/startup.sh

