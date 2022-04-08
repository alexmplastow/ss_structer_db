# ss_structure_db
Beta Version for a Tool to generate RNA secondary structures across the transcriptome, 
  should work fine if one wishes to generate secondary structures for the complete human transcriptome

To begin utilizing this tool, one must have either RNAStructure or ViennnaRNA downloaded. By default, this CLI will use
  ViennaRNA.
  
ViennaRNA installation:

` wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_5_x/ViennaRNA-2.5.0.tar.gz
 tar -zxvf ViennaRNA-2.5.0.tar.gz
 cd ViennaRNA-2.5.0
 ./configure --with-python`

 
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
 `make
 sudo make install`
 
Depending on your version of libc, this should execute seamlessly. 
  If it doesn't and you have root privileges, it shouldn't be anything you shouldn't be able to resolve, 
  but if not, you may want to work something out with a sysadmin.
  
To append the scripts and binary to the appropriate python path. Change line 2 of the 'functions.py' to match the
  binaries and scripts directory.
  
I know this is kinda crude, I intend to change this in the future to 
  execute the command from bash to bypass ViennaRNA's idiosyncratic and difficult python installation procedure.

Speaking of executing from Bash, RNAStructure is utilized via the os.system() command, so RNAStructure should be installed.
  
  
 `wget http://rna.urmc.rochester.edu/Releases/current/RNAstructureLinuxTextInterfaces64bit.tgz
  tar xvf RNAstructureLinuxTextInterfaces64bit.tgz
  cd RNAstructure
  make all`
 
Note the present working directory for RNAStructure's installation. 
  You will need to add the executables to the PATH variable so,
  
  `cd ; vim .bashrc`
  
At the beginning of the file or anywhere else, really, write:

  `export PATH="/home/spectre/tools/RNAstructure/exe:$PATH"
  export DATAPATH="/home/spectre/tools/RNAstructure/data_tables/"`
  
 After that, one should install the python libraries: 
 
 `pip install {matplotlib,numpy}`
or 
  `pip3 install {matplotlib, numpy}`
  
That is the end of the installation.

  `python ss_structure_db.py -h`
  
To begin the calibration to estimate the amount of time needed to fold the transcriptome

usage: ss_structure_db.py [-h] [-i INPUT] [-c] [-s [SAMPLE_SIZE]] [-ff [FAST_FOLD_TIME]] [-cof [CALIBRATION_OUTPUT]] [-ft] [-l [NUCLEIC_ACID_LENGTH]] [-rs]

A CLI for generating or improving on transcriptome-wide secondarystructure datasets

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input transcriptome_file
  -c, --calibrate       Determine how long folding will take depending on folding time of a set of nucleotides
  -s [SAMPLE_SIZE], --sample_size [SAMPLE_SIZE]
                        nucleic acid count for folding (default 50)
  -ff [FAST_FOLD_TIME], --fast_fold_time [FAST_FOLD_TIME]
                        nucleic acid count for folding (default 50)
  -cof [CALIBRATION_OUTPUT], --calibration_output [CALIBRATION_OUTPUT]
                        output file for calibration
  -ft, --fold_transcriptome
                        indicate if whole transcriptome should begin folding
  -l [NUCLEIC_ACID_LENGTH], --nucleic_acid_length [NUCLEIC_ACID_LENGTH]
                        maximum length of nucleic acids, anything above will be folded in iterations of 120 nucleotides
  -rs, --RNAStructure   switches the folding algorithm to RNAStructure
 
Calibrating indicates how long transcriptome folding will take depending on folding algorithm and preferred length.

Calibrate like so:
  
 `python ss_structure_db.py -i  hg38* -c -s 30 -ff 100 -cof calibration_output.csv`

This script is configured to use multiprocessing to switch folding modes based on the amount of time the script has taken to fold the nucleic acid sequence.

The user must identify the nucleotide length at which the folding algorithm performs optimally. Execute the following:

  `python folding_time.py -it hg38* -ico calibration_output.csv`
  
  
One should get a file output like this
  
![folding_time_plot_test](https://user-images.githubusercontent.com/79552389/162453802-cae32cbd-4d93-43a1-9efa-d4b578d0d1ed.png)

It looks like the local minima is around 1800 nucleotides, so the cutoff for the slow folding mode is 1800 nucleotides. This is the default value in any case.

To begin folding the transcriptome, one must execute the following command.

   `python ss_structure_db.py -i *0221* -s 10 -ft`
   
By default, this script performs folding via ViennaRNA, but the folding mode may be changed to RNAStructure. To execute the folding with RNAStructure:

  `python ss_structure_db.py -i *0221* -s 10 -ft -rs` 

If RNAStructure has been configured correctly, the following should output to the terminal:

Initializing nucleic acids...done.
 95% [================================================  ] |                     done.
Writing output ct file...done.
Single strand folding complete.
Converting CT file...CT file conversion complete.

RNAStructure seems a bit slower than ViennaRNA. In the future, RNAStructure will be configured to run alongside ViennaRNA so a database of secondary structure which compares multiple folding algorithm outputs.

The input for the -s variable determines how many transcripts will be folded. This will be reconfigured to allow the user to fold the transcriptome indefinately.

The output for the folding algorithms will be placed in a csv file in the following format:


![Screenshot 2022-04-08 103733](https://user-images.githubusercontent.com/79552389/162459404-fbfde93d-2283-405c-87de-d6e4a2e2aae0.png)

The script is configured to check the CSV each time it executes to ensure nothing is being refolded, so take care to ensure you are not indiscriminately folding nucleic acids without regard for algorithm types.


