
set -e

wget -O v3.0.0.tar.gz https://github.com/yukiteruono/pbsim3/archive/refs/tags/v3.0.0.tar.gz
tar -xvzf v3.0.0.tar.gz
cd pbsim3-3.0.0 && sudo ./configure && sudo make && sudo make install