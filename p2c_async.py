import csv
import glob
import multiprocessing
import os
import re
import time

SEQUENCE_FILE = './data/test_sequences.txt'
PRIMER_FILES = glob.glob('./data/primer/*.csv')
RES_COLUMNS = ['contig name', 'primer', 'start', 'end', 'length', 'requested', 'located']
RES_FOLDER = os.path.abspath('final')
RES_PREFIX = 'hits_'

MAX_PARALLEL_PROC = 0

# progress Bar variables
DISPLAY_PROGRESS = True

PRIMER_AMOUNT_TOTAL = 0
PRIMER_AMOUNT_PROCESSED = 0
PRIMER_CALC_STARTING_TIMER = None

class contig(object):
    """contig / primer structure including option for a position parameter
    
    """
    
    def __init__(self, name, sequence, counter = None):
        if name.startswith('>'):
            name = name[1:]
        self.name       = name
        self.sequence   = sequence.upper()

class primer_reader(object):
    """read primer_file as csv
    
    """
    def __init__(self, primer_file):
        """ prepare primer file and create ordered primer list
        
        """
        self._primer_file = primer_file
        self._primer_list = []
        self._read_csv()
        
    def _read_csv(self, delimiter = ';'):
        """ try to create primer object incl. probe, sequence [,position]
        
        """
        stream = csv.DictReader(open(self._primer_file, 'r'), delimiter = delimiter)
        for line in stream:
            primer_name     = line.get('probe')
            primer_sequence = line.get('sequence')
            if primer_name and primer_sequence: 
                self._primer_list.append(contig(primer_name, primer_sequence))
                
    def get_primer(self):
        """

        :return: primer collection of primer / contig objects
        :rtype: list
        """
        return self._primer_list
    

class contig_reader(object):
    """ get contigs from txt file. 
    contig file format: >name\nsequence\n>name ... aso.
    """
    def __init__(self, sequencefile):
        """

        """
        self.sequencefile   = sequencefile
        self.contig_list    = []
        self._read_contigs_list()
        
    def get_contigs(self):
        """ get contigs collection
        
        :return: contig object collection
        :rtype: list
        """
        return self.contig_list
    
    def _read_contigs_list(self):
        """ create contig list out of fasta formated textfile
        
        """
        lines = [line.strip() for line in open(self.sequencefile, 'r')]
        contig_name = ''
        contig_sequence = ''
        for line in lines:
            if line.startswith('>'):
                if contig_sequence:
                    self.contig_list.append(contig(contig_name, contig_sequence))
                contig_name     = line.strip()
                contig_sequence = ''
            else:
                contig_sequence += line.strip()
        self.contig_list.append(contig(contig_name, contig_sequence))


def progress_bar(result):
    """ display progress of current sequence / primer matching calculations 
    using global variables

    :param result: flag represents calc of primer is done
    :type result: bool
    
    :returns: None
    :rtype: NoneType
    """
    if not DISPLAY_PROGRESS:
        return None

    globals()['PRIMER_AMOUNT_PROCESSED'] += 1

    runtime = time.time() - PRIMER_CALC_STARTING_TIMER
    appr_runtime = (runtime / PRIMER_AMOUNT_PROCESSED) * (PRIMER_AMOUNT_TOTAL - PRIMER_AMOUNT_PROCESSED)

    print('primer calculated {0} / {1} - est. runtime: {2:.0f}s'
        .format(
            PRIMER_AMOUNT_PROCESSED, 
            PRIMER_AMOUNT_TOTAL, 
            appr_runtime
        ).ljust(80), 
        end = '\r'
    )

    # final status line
    if PRIMER_AMOUNT_TOTAL == PRIMER_AMOUNT_PROCESSED:
        print('primer calculated {0} / {1} - total runtime: {2:.2f}s'
            .format(
                PRIMER_AMOUNT_PROCESSED, 
                PRIMER_AMOUNT_TOTAL, 
                runtime
            ).ljust(80), 
            end = '\n\r'
        )        

def write_matches(contig_entry, primer_list, result_file):
    """ find primer matches in contig entry and write them in result_file

    """
    proc_count = multiprocessing.cpu_count()
    if MAX_PARALLEL_PROC > 0:
        proc_count = min(proc_count, MAX_PARALLEL_PROC)

    pool = multiprocessing.Pool(proc_count)

    if DISPLAY_PROGRESS:
        globals()['PRIMER_AMOUNT_TOTAL'] = len(primer_list)
        globals()['PRIMER_AMOUNT_PROCESSED'] = 0
        globals()['PRIMER_CALC_STARTING_TIMER'] = time.time()

    results= list()
    for primer in primer_list:
        results.append(
            pool.apply_async(
                find_match, 
                args = (contig_entry, primer), 
                callback = progress_bar
            )
        )
    pool.close()
    pool.join()

    lines = [p.get() for p in results]
    lines = sum(lines, [])

    with open(result_file, 'w', newline = '') as csvfile:
        writer = csv.DictWriter(csvfile, delimiter = ';', fieldnames=RES_COLUMNS)
        writer.writeheader()
        writer.writerows(lines)

def find_match(contig_entry, primer):
    """

    """
    text = contig_entry.sequence
    compiled_pattern = re.compile(r'(%s)'%(primer.sequence))
    
    hits = []
    for m in compiled_pattern.finditer(text):
        if m.group():
            start, end = m.span()
            hits.append((m.group(), start, end))

    if hits:
        return prepare_hit_lines(contig_entry, primer, hits)
    return list()

def prepare_hit_lines(contig, primer, hit_list):
    """ create entry for each hit
    
    """
    lines = list()
    for hit in hit_list:
        start = hit[1]
        end = hit[2]
        length = end-start
        contig_name = contig.name 
        primer_name = primer.name
        requested = primer.sequence
        located = hit[0]
        line = {
            'contig name': contig_name, 
            'primer': primer_name, 
            'start': start, 
            'end': end, 
            'length': length, 
            'requested': requested, 
            'located': located
        }

        lines.append(line)
    return lines

def main():
    if not os.path.exists(RES_FOLDER):
        os.makedirs(RES_FOLDER)
    
    seq = contig_reader(SEQUENCE_FILE)
    
    for i, primer_file in enumerate(PRIMER_FILES):
        primer = primer_reader(primer_file)
        res_file = os.path.join(
            RES_FOLDER, 
            RES_PREFIX + os.path.basename(primer_file)
        )

        for contig in seq.get_contigs():
            print('sequence: {0} - primer collection: {1} ({2} / {3})'.
                format(
                    contig.name, 
                    os.path.basename(primer_file), 
                    i + 1, 
                    len(PRIMER_FILES)
                )
            )
            write_matches(contig, primer.get_primer(), res_file)

if __name__ == '__main__':
    multiprocessing.freeze_support()

    print('primer 2 comment (v2.3)\n-----------------------', )

    start = time.time()
    try:
        main()
    except Exception as e:
        print(e)
        exit()

    end = time.time()
    print('absolute runtime: {:.2f}s'.format(end - start))
