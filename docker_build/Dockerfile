FROM tylerreddy/steric-conflicts:latest

#install some more dependencies needed for GROMACS build:
RUN apt-get update && apt-get install -y \
	cmake

#Download and build / install GROMACS:
RUN wget ftp://ftp.gromacs.org/pub/gromacs/gromacs-5.1.2.tar.gz
RUN tar xfz gromacs-5.1.2.tar.gz && cd gromacs-5.1.2 && mkdir build && cd build && cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON && make && make check && sudo make install && rm -r /gromacs-5.1.2 && rm /gromacs-5.1.2.tar.gz
RUN echo "source /usr/local/gromacs/bin/GMXRC" >> ~/.bashrc
ADD *py /steric_conflict_resolution/
ADD test_data /steric_conflict_resolution_work
RUN pip install nose coveralls

#coveralls currently requires a git repo & unique commit to submit the report, so
#do this artificially:
RUN cd /steric_conflict_resolution && \
    git config --global user.email "dummy@dummy.com" && \
    git config --global user.name "Dummy" && \
    git init && \
    git add * && \
    git commit -m "Dummy."
