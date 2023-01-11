import sys


def parse_mummer():
    mummer = sys.argv[2]
    all_contig = sys.argv[1]
    contig_len = dict()
    with open(all_contig, "r") as f:
        lines = f.readlines()
        contig = ""
        for i in range(0, len(lines), 2):
            line = lines[i]
            contig = (line.split(" ")[0])[1:]
            current_line = lines[i+1]
            contig_len[contig] = len(current_line)
    f.close()
    # print(contig_len)
    mummer_dict = dict()
    with open(mummer, "r") as f:
        lines = f.readlines()
        contig = ""
        # build a dictionary using the mummer output
        for line in lines:
            if line[0] == ">":
                contig = line[1:].strip()
                mummer_dict[contig] = 0
            else:
                line = line.strip()
                arr = line.split()
                # print(arr)
                if len(arr) == 3:
                    og_pos, contig_pos, length = arr[0],arr[1],arr[2]
                else:
                    ref, og_pos, contig_pos, length = arr[0],arr[1],arr[2],arr[3]
                mummer_dict[contig] += int(length)
                # somehow save the next few lines into dict until the next contig
        f.close()
        
    if len(sys.argv) > 4 and sys.argv[3].isnumeric():
        thresh = float(sys.argv[3])
    else: thresh = 0.5
    for contig in mummer_dict:
        if mummer_dict[contig] <= thresh * contig_len[contig]:
            print(contig)


if __name__ == "__main__":
    # print("usage: python filter_contigs.py <contigs> <mummer output>")
    parse_mummer()
