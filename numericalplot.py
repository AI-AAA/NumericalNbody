#%matplotlib qt5
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import axes3d

#G= 2.95912208286 * (10**-4)
colors = ["yellow", "red", "green", "blue", "cyan", "purple",
          "pink","navy","black","brown","turquoise","magenta","lavender"]
def plot_orbit(planet_name, all_positions, method,type,title):
    # Define colors for each planet
    
    fig = plt.figure()
    if type=="3d":
        ax = fig.gca(projection='3d')
        
    for i, name in enumerate(planet_name):
        # Extract x, y, and z positions for each planet
        x_data = np.array([pos[i][0] for pos in all_positions])
        y_data = np.array([pos[i][1] for pos in all_positions])
        
        if type=="2d":
            if len(planet_name)>20:
                first_point_size = 100
                other_point_size = 1
                sizes = [first_point_size] + [other_point_size] * (len(x_data) - 1)
                ax.scatter(x_data, y_data, s=sizes, color="black")
            else:
                # Plot the orbit and the planet position
                plt.plot(x_data, y_data, linestyle='-', color=colors[i])
                plt.plot(x_data[-1], y_data[-1], marker='o', markersize=8, 
                         linestyle='-', label=name, color=colors[i])

                # Add legend and planet name text
                plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
                plt.text(x_data[-1], y_data[-1], name[0], color='black')

            # Set axis labels
            plt.xlabel("x")
            plt.ylabel("y")
        else:
            z_data = np.array([pos[i][2] for pos in all_positions])
            if len(planet_name)>20:
                first_point_size = 100
                other_point_size = 1
                sizes = [first_point_size] + [other_point_size] * (len(x_data) - 1)
                #plt.scatter3D(x_data, y_data,z_data,s=sizes, color="black")
                ax.scatter3D(x_data, y_data,z_data, s=sizes, color='red', edgecolor='black')
            else:
                # Plot the orbit and the planet position
                ax.plot(x_data, y_data, z_data, linestyle='-', color=colors[i])
                ax.plot(x_data[-2:-1], y_data[-2:-1], z_data[-2:-1], marker='o', 
                        markersize=8, linestyle='-', label=name, color=colors[i])

                # Add legend and planet name text
                ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
                ax.text(all_positions[-1][i][0], all_positions[-1][i][1], 
                        all_positions[-1][i][2], name[0], color='black')
            # Set axis labels
            ax.set_xlabel("x")
            ax.set_ylabel("y")
            ax.set_zlabel("z")
            
    # Set plot title and display the plot
    plt.title(f"{title} ({method} Method)")
    plt.savefig(f"image/{title}({method} Method)-{type}{np.random.randint(0,1000000)}img.png")
    plt.show()
    

#q:position p:momentum
from numerical import SimulateOribt
from data import loadOuterSolarData   
def plot_energy(methods,Ts,step):
    # Calculate the energy loss
    #methods=["Explicity Euler","Symplectic Euler","Stormer–Verlet","4th order Runge-Kutta"]
    fig, ax = plt.subplots(figsize=(10, 6))
    #ax.set_yscale('log')
    for i in range(len(methods)):
        allenergy=SimulateOribt(loadOuterSolarData,Ts,step, methods[i])[3]
        
        #Energyloss = [np.abs((allenergy[0]-element)/allenergy[0]) for element in allenergy]
        # Create a time array with a step of 10
        time = np.arange(0,step+1,1)
        # Select the corresponding elements of Energyloss based on the time array
        #Energyloss_selected = [Energyloss[i] for i in time]
        # Plot the energy loss over time
        ax.plot(time, allenergy, label=methods[i], marker='.', linestyle='-', markersize=3)
        #plt.plot(time, Energyloss_selected,color=colors[i])
    ax.set_title('Energy over time')
    ax.set_xlabel('time [Days]')
    ax.set_ylabel('Energy')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    plt.savefig(f"image/Energy over time.png")
    plt.show()

def plot_energyloss(methods,Ts,step,log):
    # Calculate the energy loss
    #method=["Explicity Euler","Symplectic Euler","Stormer–Verlet","4th order Runge-Kutta"]
    fig, ax = plt.subplots(figsize=(10, 6))
    if log==1:
        ax.set_yscale('log')
    for i in range(len(methods)):
        allenergy=SimulateOribt(loadOuterSolarData,Ts,step, methods[i])[3]
        
        Energyloss = [np.abs((allenergy[0]-element)/allenergy[0]) for element in allenergy]
        # Create a time array with a step of 10
        #time = np.arange(0, step+1,1)
        time = np.arange(0, len(Energyloss),1)
        # Select the corresponding elements of Energyloss based on the time array
        #Energyloss_selected = [Energyloss[i] for i in time]
        # Plot the energy loss over time
        ax.plot(time, Energyloss, label=methods[i], marker='.', linestyle='-', markersize=6)
        #plt.plot(time, Energyloss_selected,color=colors[i])
    ax.set_title('Energy loss over time')
    ax.set_xlabel('time [Days]')
    ax.set_ylabel('% loss from initial energy')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    plt.savefig(f"image/Energy loss over time.png")
    plt.show()


from matplotlib.animation import FuncAnimation
def update(frame, planet_name, all_positions, method, type, ax,title):
    ax.clear()
    plt.title(f"{title} ({method} Method)")
    ax.set_xlim((np.array(all_positions)[:,:,0].min()*2, 
                 np.array(all_positions)[:,:,0].max()*2))
    ax.set_ylim((np.array(all_positions)[:,:,1].min()*2, 
                 np.array(all_positions)[:,:,1].max()*2))
    ax.grid(False)
    for i, name in enumerate(planet_name):
        x_data = np.array([pos[i][0] for pos in all_positions[:frame+1]])
        y_data = np.array([pos[i][1] for pos in all_positions[:frame+1]])
        
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        if type == "2d":
            ax.plot(x_data, y_data, linestyle='-', color=colors[i])
            if frame > 0:
                ax.plot(x_data[-1], y_data[-1], marker='o', markersize=8, 
                        linestyle='-', label=name, color=colors[i])
                
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
        else:
            z_data = np.array([pos[i][2] for pos in all_positions[:frame+1]])
            ax.set_zlim((np.array(all_positions)[:,:,2].min()*2, np.array(all_positions)[:,:,2].max()*2))
            ax.plot(x_data, y_data, z_data, linestyle='-', color=colors[i])
            if frame > 0:
                ax.plot(x_data[-1], y_data[-1], z_data[-1], marker='o', markersize=8,
                        linestyle='-', label=name, color=colors[i])
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
            ax.set_zlabel("z")

def animate_orbit(planet_name, all_positions, method, type,title):
    fig = plt.figure()
    if type == "2d":
        ax = fig.add_subplot()
    else:
        ax = fig.add_subplot( projection='3d')
    all_positions = all_positions[::20]
    ani = FuncAnimation(fig, update, frames=range(len(all_positions)), 
                        fargs=(planet_name, all_positions, method, type, ax,title), interval=10, blit=False)
    ani.save(f'{title}({method}){type}animation.gif', writer='imagemagick')
    #writer = animation.FFMpegWriter(fps=30)
    #ani.save('animation.mp4', writer=writer)
    return ani