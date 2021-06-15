FROM centos:latest
USER root

# install the SSSD-CLIENT needed by Arcade ?
RUN yum install -y sssd-client acl xterm python3

RUN pip3 install astropy

# system settings and permissions
ADD nofiles.conf /etc/security/limits.d/

## see https://bugzilla.redhat.com/show_bug.cgi?id=1773148
RUN touch /etc/sudo.conf && echo "Set disable_coredump false" > /etc/sudo.conf

# generate missing dbus uuid (issue #47)
RUN dbus-uuidgen --ensure

ADD nsswitch.conf /etc/

ARG BUILDDIR=/opt/findMoving
RUN mkdir -p ${BUILDDIR}
WORKDIR ${BUILDDIR}
COPY src  ./
WORKDIR ${BIULDDIR}/src
RUN python3 setup.py install

CMD [ "bash" ]

