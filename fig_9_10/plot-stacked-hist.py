import sys
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser(description='Plot numbers from stdin as a barchart.')
parser.add_argument('output', help='Output file.')
parser.add_argument('--x-label', help='X axis label.')
parser.add_argument('--y-label', help='Y axis label.')
parser.add_argument('--x-line', help='X=arg line.', type=float, default=0)
parser.add_argument('--y-line', help='Y=arg line.', type=float, default=0)
parser.add_argument('--legend', help='Legend.', nargs='*', default=[])
args = parser.parse_args()

series = defaultdict(list)
for line in sys.stdin:
    category, value = int(line.split()[0]), float(line.split()[1])
    series[category].append(value)
    
if args.y_line > 0:
    plt.axhline(y=args.y_line, color='r', alpha=0.5)
if args.x_line > 0:
    plt.axvline(x=args.x_line, color='r', alpha=0.5)
plt.xlabel(args.x_label)
plt.ylabel(args.y_label)
plt.hist([series[i] for i in range(len(series))], bins=50, range=(0, 50), stacked=True, label=args.legend)
plt.legend()
plt.savefig(args.output, bbox_inches='tight')
