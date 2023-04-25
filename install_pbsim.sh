
set -e

wget -O v3.0.0.zip https://github.com/yukiteruono/pbsim3/archive/refs/tags/v3.0.0.zip
unzip v3.0.0.zip
cd pbsim3-3.0.0 && sudo ./configure && sudo make && sudo make install