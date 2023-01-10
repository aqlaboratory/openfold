# from openfold.data.parsers import parse_a3m
#
# def recover_sequences_from_a3m_entries(a3m_text):
#     a3m_entries = parse_a3m(a3m_text)
#     sanitized = []
#     print(a3m_entries)
#     for entry in a3m_entries:
#         # 1. Convert insertions ('e') to sequence ('E')
#         # 2. Remove deletions (delete '-'s)
#         # 3. No need to delete gaps ('.') since A3M can omit them.
#         sanitized_entry = entry.upper().replace('-', '')
#         sanitized.append(sanitized_entry)
#
#     return sanitized
#
# if __name__ == '__main__':
#     a3m = '>1LLD_1_A\n'\
# 'AETTVKPTKLAVIGAGAVGSTLAFAAAQRGIAREIVLEDIAKERVEAEVLDMQHGSSFYPTVSIDGSDDPEICRDADMVVITAGPRQKPGQSRLELVGATVNILKAIMPNLVKVAPNAIYMLITNPVDIATHVAQKLTGLPENQIFGSGTNLDSARLRFLIAQQTGVNVKNVHAYIAGEHGDSEVPLWESATIGGVPMSDWTPLPGHDPLDADKREEIHQEVKNAAYKIINGKGATNYAIGMSGVDIIEAVLHDTNRILPVSSMLKDFHGISDICMSVPTLLNRQGVNNTINTPVSDKELAALKRSAETLKETAAQFGF\n'\
# '>UniRef100_A0A0 L-lactate dehydrogenase n=1 Tax=Microbacterium sp. MRS-1 TaxID=1451261 RepID=A0A011U3P1_9MICO\n'\
# '---PRRNSKLAIVGAGSVGSSLAYAALIRGSAREVALYDINAAKVEAEVLDLAHGTQFTPASSVTGGDDIEVCAGSDVVVITAGAKQKPGQSRLDLAGANVRILRSLMPQLVEVAPDAVYVLVTNPVDVLTYAAQKISGLPPNRVFGSGTVLDSSRLRWLLAQRAGVAVQSVHAYIVGEHGDSEFPLWSSATIGGVPLLDWRPLDGRPPFTEEERDEIAHEVVNAAYEIIEGKGATNYAIGLSGARIVEAILRDEHRVLPVSTVLDGYHGISGVALSVPSVVGRGGVERVLEVPLSDEELAGLRASAEALREVARSLGF\n'
#     recovered_sequences = parse_a3m(a3m)
#     print(recovered_sequences)