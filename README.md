# ss_structure_db
Beta Version for a Tool to generate RNA secondary structures across the transcriptome, 
  should work fine if one wishes to generate secondary structures for the complete human transcriptome

To begin utilizing this tool, one must have either RNAStructure or ViennnaRNA downloaded. By default, this CLI will use
  ViennaRNA.
  
ViennaRNA installation:
 
 wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_5_x/ViennaRNA-2.5.0.tar.gz
 tar -zxvf ViennaRNA-2.5.0.tar.gz
 cd ViennaRNA-2.5.0
 ./configure --with-python
 
At this point, you should get an output indicating the install directories, make a not of these because ViennaRNA
  can be a pain to get running in python. You should get an output which looks like this:
  
Install Directories
-------------------
  * Executables               : /usr/local/bin
  * Libraries                 : /usr/local/lib
  * Header files              : /usr/local/include
  * Extra Data                : /usr/local/share
  * Man pages                 : /usr/local/share/man
  * Documentation             : /usr/local/share/doc/ViennaRNA
      (HTML)                  : /usr/local/share/doc/ViennaRNA/html
      (PDF)                   : /usr/local/share/doc/ViennaRNA
  * Perl5 Interface           :
      (binaries)              : /usr/local/lib/x86_64-linux-gnu/perl/5.30.0
      (scripts)               : /usr/local/share/perl/5.30.0
  * Python 3.x Interface      :
      (binaries)              : /usr/local/lib/python3.8/site-packages
      (scripts)               : /usr/local/lib/python3.8/site-packages
  * Python 2.x Interface      : Not to be installed python2 executable missing
      (binaries)              :
      (scripts)               :
 
 
The python 3.x Interface binaries and scripts are both to be installed in the /usr/local/lib/python3.8/site-packages directory
 make
 sudo make install
 
Depending on your version of libc, this should execute seamlessly. 
  If it doesn't and you have root privileges, it shouldn't be anything you shouldn't be able to resolve, 
  but if not, you may want to work something out with a sysadmin.
  
To append the scripts and binary to the appropriate python path. Change line 2 of the 'functions.py' to match the
  binaries and scripts directory.
  
I know this is kinda crude, I intend to change this in the future to 
  execute the command from bash to bypass ViennaRNA's idiosyncratic and difficult python installation procedure.

Speaking of executing from Bash, RNAStructure is utilized via the os.system() command, so RNAStructure should be installed.
  
  wget http://rna.urmc.rochester.edu/Releases/current/RNAstructureLinuxTextInterfaces64bit.tgz
  tar xvf RNAstructureLinuxTextInterfaces64bit.tgz
  cd RNAstructure
  make all
  
Note the present working directory for RNAStructure's installation. 
  You will need to add the executables to the PATH variable so,
  
  cd 
  
