# primer2comment
## create executable
pyinstaller primer2comment.py

## required subfolder / files
.
+-- p2c-async.py
+-- data
|    +-- primer
|    |   +-- <primer_collection>.csv
|    |   +--...
|    +-- test_sequence.txt

## primer-files
probe;sequence
primer_name;ATGC - Sequence 

Note: Squence can contain wildcards as used in regular expressions

## sequence file (fasta format)
>sequence-name
ATGC-sequence
