import sys

print("central_sample_id,cluster")
for fname in sys.argv[1:]:
    cluster = fname.split("_")[-1]
    with open(fname) as f:
        for l in f:
            if "/" in l:
                print(l.split("/")[1], cluster, sep=",")
            else:
                print(l.rstrip(), cluster, sep=",")
