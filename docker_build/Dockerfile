FROM tylerreddy/steric-conflicts:latest

#install some more dependencies needed for GROMACS build:
RUN apt-get update && apt-get install -y \
	cmake

#Download and build / install GROMACS:
RUN wget ftp://ftp.gromacs.org/pub/gromacs/gromacs-5.1.2.tar.gz
RUN tar xfz gromacs-5.1.2.tar.gz && cd gromacs-5.1.2 && mkdir build && cd build && cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON && make && make check && sudo make install && /bin/bash -c "source /usr/local/gromacs/bin/GMXRC" && rm -r /gromacs-5.1.2
