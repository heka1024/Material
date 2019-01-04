from pylab import figure, axes, pie, title, savefig
import matplotlib.pyplot as plt

#그래프 그리기
def draw_graph(x, y):
    plt.axhline(0, color='black', linewidth = 0.2)
    plt.plot(x, y, color = 'red')
    plt.xlabel('x')
    plt.ylabel('Shear Force')
    plt.title('SFD')

def generate():
	x = []
	sfd = []
	bmd = []
	f = open('data.txt', "rt")
	lines = f.read().splitlines()
	for line in lines:
		k = list(map(float, line.split(', ')))
		x.append(k[0])
		sfd.append(k[1])
	f.close()
	draw_graph(x, sfd)
	savefig('output/sfd.pdf')

generate()
