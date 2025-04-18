from Bio.Affy import CelFile

with open("sample.CEL", "r") as handle:
    c = CelFile.read(handle)
    for i, row in enumerate(c):
        print(row)  # ogni riga è un intensità spot
        if i > 10: break