import csv
from concurrent.futures import ProcessPoolExecutor

PRIMER_FILE = './data/primer/primer1collection.csv'
TEST_SEQUENCE_FILE = './data/test_sequences.txt'
RESULT_FILE = './final/result.csv'

def KMP_search(primer_seq, test_sequence):
    """
    Function for performing KMP search for a single primer sequence
    in a single test sequence
    """
    matches = []
    n, m = len(test_sequence), len(primer_seq)
    lps = [0] * m
    j = 0  # index for primer_seq
    computeLPSArray(primer_seq, m, lps)
    i = 0  # index for test_sequence
    while i < n:
        if primer_seq[j] == test_sequence[i]:
            i += 1
            j += 1
        if j == m:
            matches.append((i-j, j))
            j = lps[j-1]
        elif i < n and primer_seq[j] != test_sequence[i]:
            if j != 0:
                j = lps[j-1]
            else:
                i += 1
    return matches

def computeLPSArray(pattern, m, lps):
    len = 0
    lps[0] = 0
    i = 1
    while i < m:
        if pattern[i] == pattern[len]:
            len += 1
            lps[i] = len
            i += 1
        else:
            if len != 0:
                len = lps[len-1]
            else:
                lps[i] = 0
                i += 1

def search_sequences(primer_name, primer_seq, sequences):
    """
    Function for searching a primer sequence in a list of test sequences
    """
    matches = []
    for name, sequence in sequences:
        match = KMP_search(primer_seq, sequence)
        if match:
            for start, length in match:
                matches.append((primer_name, name, sequence[start:start+length], primer_seq, start+1, length))
    return matches

def main():
    import os
    for filename in os.listdir('./data/primer/'):
        if filename.endswith(".csv"):
            primers = {}
            # Öffnen der Primer-Datei
            with open(os.path.join('./data/primer/', filename)) as primer_file:
                for line in primer_file:
                    # Überspringen der Kopfzeile
                    if "probe" in line:
                        continue
                    # Aufteilen der Zeile in Probe-Name und Sequenz
                    parts = line.strip().split(";")
                    primers[parts[0]] = parts[1]

            test_sequences = []
            current_sequence = ""
            current_name = ""
            with open(TEST_SEQUENCE_FILE) as test_file:
                for line in test_file:
                    if line.startswith(">"):
                        test_sequences.append((current_name, current_sequence))
                        current_name = line[1:].strip()
                        current_sequence = ""
                    else:
                        current_sequence += line.strip()
            test_sequences.append((current_name, current_sequence))
            
            # Erstellen einer eigenen Ergebnisdatei für jede Primer-Datei
            result_file = os.path.join('./final/', os.path.splitext(filename)[0] + '_result.csv')
            with open(result_file, 'w', newline='') as result_file:
                result_writer = csv.writer(result_file, delimiter=';')
                result_writer.writerow(['probe','sequence','matched_sequence', 'regular_expression', 'position', 'match_length'])

                with ProcessPoolExecutor() as executor:
                    results = [executor.submit(search_sequences, primer_name, primer_seq, test_sequences) for primer_name, primer_seq in primers.items()]

                    for f in results:
                        matches = f.result()
                        for match in matches:
                            result_writer.writerow(match)

if __name__ == "__main__":
    main()
