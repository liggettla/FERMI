wget https://repo.continuum.io/archive/Anaconda2-4.3.0-Linux-x86_64.sh
bash Anaconda2-4.3.0-Linux-x86_64.sh
conda env create -f fermi.yml
cd ../main
x=pwd
PATH=$PATH:$x
echo 'export PATH='$PATH':'$x >> ~/.profile