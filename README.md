AuToGraFS
=========

Automatic Topological Generator for Framework Structures

Manual coming soon...

AuToGraFs is now accessible via pip!

https://pypi.python.org/pypi/AuToGraFS/2.0.4b0

If you'd stil like to install manually: Here's how to install AuToGraFS on a clean Ubuntu (tested with 16.04) install:
```
sudo apt-get install python-dev
sudo apt-get install python-pip
sudo pip install numpy
```
then unpack the attached tarball in your home directory:
```
tar -zxf autografs_0.1.1.tar.gz
```

add the following two lines to your .bashrc
```
export PYTHONPATH=$HOME/autografs/:$PYTHONPATH
export PATH=$HOME/bin/:$HOME/autografs/ase/tools:$PATH
```
in the autografs directory that was just created in your homedir, you should see a subdir called autografs_examples

To run one of the examples, copy it to a new directory and (for sra example):
```
~/autografs/mofgen.py -c control_sra.txt -o sra -e cif >sra.out
```

