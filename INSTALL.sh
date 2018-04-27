# Boot something like ubuntu-wily-15.10-amd64-server-20160222 (ami-05384865).

# install required Debian packages
sudo apt update && \
     sudo apt -y install python-dev python-pip python-virtualenv

# activate virtualenv
python -m virtualenv env
. env/bin/activate

# get spacegraphcats
git clone https://github.com/spacegraphcats/spacegraphcats.git

# install requirements
pip install -r spacegraphcats/requirements.txt

# install
python setup.py develop
or if you are not a developer:
python setup.py install

