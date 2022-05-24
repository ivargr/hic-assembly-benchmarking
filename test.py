from bionumpy.files import bnp_open

data = bnp_open("test.gfa", chunk_size=10000000)

for line in data:
    print(line)
