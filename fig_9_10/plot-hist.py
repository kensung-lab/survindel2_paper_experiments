# read numbers from stdin, one per line, and plot them as a histogram between 0 and N
# numbers greated than N are excluded

import sys, math
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description='Plot numbers from stdin as a barchart.')
parser.add_argument('output', help='Output file.')
parser.add_argument('--x-label', help='X axis label.', default='')
parser.add_argument('--y-label', help='Y axis label.', default='')
parser.add_argument('--y-log', help='Y axis log.', action='store_true')
parser.add_argument('--bin-size', help='Bin size.', type=float, default=1)
parser.add_argument('--cum', help='Plot a line showing the cumulative frequency.', action='store_true')
parser.add_argument('--edge-color', help='Edge color.', default='none')
parser.add_argument('--legend', help='Legend.')
args = parser.parse_args()

numbers = []
for line in sys.stdin:
    numbers.append(float(line))
max_number = max(numbers)

plt.hist(numbers, bins=math.ceil(max_number/args.bin_size), range=(0, max_number), label=args.legend, edgecolor=args.edge_color)
if args.cum:
    plt.hist(numbers, bins=max_number, range=(0, max_number), cumulative=True, histtype='step', color='r', alpha=0.5)
plt.xlabel(args.x_label)
plt.ylabel(args.y_label)
if args.y_log:
    plt.yscale('log')
plt.savefig(args.output, bbox_inches='tight')
