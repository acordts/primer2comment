# primer2comment
## create executable
pyinstaller --onefiles --console primer2comment.py

## required subfolder / files
```
primer2comment
│   README.md
│   p2c-async.py
│
└───data
│   │   test_sequence.txt
│   │
│   └───primer
│       │   primer_collection1.csv
│       │   primer_collection2.csv
│       │   ...
│   
└───final (generated)
    │   hits_primer_collection1.csv
    │   hits_primer_collection2.csv
    │   ...
    
```

## primer-files
```
probe;sequence
primer_name;ATGC - Sequence 
```

Note: Sequence can contain regular expressions

## sequence file (fasta format)
```
>sequence-name
ATGC-sequence
```
