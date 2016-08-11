FROM tylerreddy/steric-conflicts:latest

#install some more dependencies needed for GROMACS build:
RUN apt-get update && apt-get install -y \
	cmake

#Download and build / install GROMACS:
# need the specific fixes from Berk Hess to handle lambda-growth mdrun of massive systems (not available in release-2016 branch yet)
RUN git clone https://github.com/gromacs/gromacs.git && cd gromacs && git fetch https://gerrit.gromacs.org/gromacs refs/changes/05/6105/1 && git checkout FETCH_HEAD && mkdir build && cd build && cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON && make && make check && sudo make install && rm -r /gromacs
RUN echo "source /usr/local/gromacs/bin/GMXRC" >> ~/.bashrc
ADD *py /steric_conflict_resolution/
ADD test_data /steric_conflict_resolution_work/test_data
RUN pip install nose coveralls mock

#coveralls currently requires a git repo & unique commit to submit the report, so
#do this artificially:
RUN cd /steric_conflict_resolution && \
    git config --global user.email "dummy@dummy.com" && \
    git config --global user.name "Dummy" && \
    git init && \
    git add * && \
    git commit -m "Dummy."
