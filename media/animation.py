from utils_plot import *

from IPython.display import HTML
from matplotlib import animation

start_points = [([-0.4205748172100388, -0.19077787500541843, 0.8869727309997525], np.array([-0.57843339+0.7803212j ,  0.20781569-0.11543948j]))]

step = 0.1
t0 = 0
t = 2*np.pi
dXdt = lambda x, t: np.array([-x[1], (x[0] - x[2]), x[1]])

# # Hadamard Animation
X0, Z0 = start_points[0]
dXdt = lambda x, t: np.array([-x[1], (x[0] - x[2]), x[1]])
dataSet, ts = rk4(X0, t0, t, step, dXdt)
dataSet = dataSet.T
numDataPoints = len(ts)

def animate_func(num, axis_off=True):
    ax.clear()
    
    # Plots sphere
    plot_H_bloch_regions(ax)

    # Plots trajectory line
    ax.plot(dataSet[0, :num+1], dataSet[1, :num+1], 
              dataSet[2, :num+1], color='blue')

    # Updating Point Location
    ax.scatter(dataSet[0, num], dataSet[1, num], dataSet[2, num], 
               color='blue', marker='o', s=100)
    
    # Add start point
    origin = str((dataSet[0, 0], dataSet[1, 0], dataSet[2, 0]))
    ax.plot3D(dataSet[0, 0], dataSet[1, 0], dataSet[2, 0],
              c='black', marker='o', label=origin)

    # Setting Axes Limits
    ax.set_xlim3d([-1, 1])
    ax.set_ylim3d([-1, 1])
    ax.set_zlim3d([-1, 1])

    # Adding Figure
    if axis_off:
        ax.axis('off')
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        ax.get_zaxis().set_visible(False)
    else:
        ax.set_title('Rotated Cartesian Evolution\nTime = ' + str(np.round(ts[num], decimals=2)) + ' sec')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')

fig = plt.figure()
ax = fig.add_subplot(111,projection='3d',box_aspect=(1,1,1))
ax.set_facecolor((1,1,1,0.1))
plt.tight_layout()

# Plotting the Animation
if 1:
    line_ani = animation.FuncAnimation(fig,
                                       animate_func,
                                       interval=100,   
                                       frames=numDataPoints)
    anim = HTML(line_ani.to_jshtml())
    plt.close(fig)
    anim
    print("Animation made")
    writergif = animation.PillowWriter(fps=numDataPoints/(2*np.pi))
    line_ani.save(path + "/Animations/had.gif",writer=writergif)
    print("Animation stored")
