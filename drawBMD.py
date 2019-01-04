from pylab import figure, axes, pie, title, savefig
import matplotlib.pyplot as plt

#그래프 그리기
def draw_graph(x, y):
    plt.axhline(0, color='black', linewidth = 0.2)
    plt.plot(x, y, color = 'blue')
    plt.xlabel('x')
    plt.ylabel('Moment')
    plt.title('BMD')

def generate():
	xs = []
	bmds = []
	with open("data.txt", "r") as file:
		for line in file :
			(x, sfd, bmd) = map(float, line.strip().split(", "))
			xs.append(x)
			bmds.append(bmd)
	draw_graph(xs, bmds)
	savefig('output/bmd.pdf')

generate()
