import sys

intervals = []
for line in sys.stdin:
    line = line.strip().split()
    intervals.append((line[0], line[1], int(line[2]), int(line[3]), int(line[4])))
intervals.sort(key=lambda x: (x[1], x[2]))

selected = []
for interval in intervals:
    if len(selected) == 0:
        selected.append(interval)
    else:
        if interval[1] != selected[-1][1] or interval[2] >= selected[-1][3]:
            selected.append(interval)

for interval in selected:
    print(*interval)


