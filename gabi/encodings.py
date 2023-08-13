# Encodings
GC = {
   'A': '0',
   'C': '1',
   'G': '1', 
   'T': '0',
   'bits': 1,
   'id': 'gc',
}
GCME = {
   'A': '00',
   'C': '01',
   'G': '10', 
   'T': '11',
   'bits': 2,
   'id': 'gcme',
}
GCMM = {
   'A': '001',
   'C': '011',
   'G': '110', 
   'T': '100',
   'N': '000',
   'bits': 3,
   'id': 'gcmm',
}
CATEGORICAL = {
   'A': '0001',
   'C': '0010',
   'G': '0100', 
   'T': '1000',
   'N': '0000',
   'bits': 4,
   'id': 'categorical',
}
ALL = [GC, GCME, GCMM, CATEGORICAL]
