import argparse
import sys
import pylab as pl
import numpy as np

class Line:
	def __init__(self, x1, y1, x2, y2, col):
		self._x = [x1, x2]
		self._y = [y1, y2]
		self._col = tuple(c / 255.0 for c in col)
		self._a = (y2 - y1) / (x2 - x1) if x2 != x1 else float('inf') 
		self._theta = np.math.atan(self._a)
		self._c = np.cos(self._theta)
		self._s = np.sin(self._theta)
		self._r = self._c * x1 + self._s * y1
	
	def dist(self, x, y):
		return self._c * x + self._s * y
		
	def plot(self):
		pl.plot(self._x, self._y, color=self._col)

def main():
	parser = argparse.ArgumentParser(description='Redraw the perspective lines of a SVG correctly')
	parser.add_argument('filename', type=str,
						help='The path to the SVG file to redraw')
	args = parser.parse_args()
	try:
		f = open(args.filename, "r")
	except FileNotFoundError:
		print("File not found:", args.filename, file=sys.stderr)
		exit(1)
	l = f.readline().replace(' ', '')
	size = {}
	for a in ['width', 'height']:
		start = l.find(a) + len(a) + 2
		end = l.find('"', start)
		v = a[0].upper()
		size[v] = float(l[start:end])
	W, H = size['W'], size['H']
	pl.xlim(0, W)
	pl.ylim(H, 0)
	print(W, H)
	f.readline()
	groups = {}
	while True:
		l = f.readline().replace(' ', '').replace('\t', '')
		if not l.startswith("<line"): break
		cs = {}
		for a in ['x1', 'x2', 'y1', 'y2']:
			start = l.find(a) + len(a) + 2
			end = l.find('"', start)
			cs[a] = float(l[start:end])
		start = l.find('rgb(') + 4
		end = l.find(',', start)
		r = int(l[start:end])
		start = end+1
		end = l.find(',', start)
		g = int(l[start:end])
		start = end+1
		end = l.find(')', start)
		b = int(l[start:end])
		li = Line(cs['x1'], cs['y1'], cs['x2'], cs['y2'], (r, g, b))
		g = r + 256*(g + 256*b)
		if g in groups: groups[g].append(li)
		else: groups[g] = [li]
		li.plot()
	horizon = groups[g]
	dots = {}
	while True:
		if not l.startswith("<circle"): break
		cs = {}
		for a in ['cx', 'cy']:
			start = l.find(a) + len(a) + 2
			end = l.find('"', start)
			cs[a] = float(l[start:end])
		best = float('inf')
		bg = 0
		for g in groups:
			li = groups[g][0]
			dis = li.dist(cs['cx'], cs['cy'])
			if dis < best:
				best = dis
				bg = g
		dots[g] = (cs['cx'], cs['cy'])
		pl.scatter(cs['cx'], cs['cy'])
		l = f.readline().replace(' ', '').replace('\t', '')
	pl.show()

if __name__ == "__main__":
	main()