#!/bin/bash 	

Fold RNAStructure_nucleic_acid.txt RNAStructure_nucleic_acid_output.txt

ct2dot RNAStructure_nucleic_acid_output.txt -1 RNAStructure_bracket_output.txt
