import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib import animation
from itertools import combinations
from mpl_toolkits.axes_grid1 import make_axes_locatable


STATUS = ['susceptible', 'infected', 'recovered', 'death', 'hospitalization']
COLORS = ['lightskyblue', 'orangered', 'green', 'gray']


class Particle:
    """A class representing a two-dimensional particle."""

    def __init__(self, x, y, vx, vy, radius, status, color):
        """Initialize the particle's position, velocity, and radius.

        Any key-value pairs passed in the styles dictionary will be passed
        as arguments to Matplotlib's Circle patch constructor.

        """

        self.r = np.array((x, y))
        self.v = np.array((vx, vy))
        self.radius = radius
        self.status = status
        self.color = color
        self.infect_time = 0
        self.hospital_time = 0
        self.detect_time = 0

    # For convenience, map the components of the particle's position and
    # velocity vector onto the attributes x, y, vx and vy.
    @property
    def x(self):
        return self.r[0]
    @x.setter
    def x(self, value):
        self.r[0] = value
    @property
    def y(self):
        return self.r[1]
    @y.setter
    def y(self, value):
        self.r[1] = value
    @property
    def vx(self):
        return self.v[0]
    @vx.setter
    def vx(self, value):
        self.v[0] = value
    @property
    def vy(self):
        return self.v[1]
    @vy.setter
    def vy(self, value):
        self.v[1] = value

    def overlaps(self, other):
        """Does the circle of this Particle overlap that of other?"""

        return np.hypot(*(self.r - other.r)) < self.radius + other.radius

    def draw(self, ax):
        """Add this Particle's Circle patch to the Matplotlib Axes ax."""

        circle = Circle(xy=self.r, radius=self.radius, color=self.color)
        ax.add_patch(circle)
        return circle

    def advance(self, dt):#move function
        """Advance the Particle's position forward in time by dt."""
        if self.status == STATUS[3] or self.status == STATUS[4]:
            return

        # ix = self.vx * np.random.uniform(-2,2)
        # iy = self.vy * np.random.uniform(-2,2)
        # if abs(ix) < 0.01:
        #     self.x += 0 * dt
        # else:
        #     self.x += ix * dt
        # if abs(iy) < 0.01:
        #     self.y += 0 * dt
        # else:
        #     self.y += iy * dt

        self.r += self.v * dt

class Simulation:
    """A class for a simple hard-circle molecular dynamics simulation.

    The simulation is carried out on a square domain: 0 <= x < 1, 0 <= y < 1.

    """

    ParticleClass = Particle

    def __init__(self, n, infectn, radius, iterations, hospital=False):
        """Initialize the simulation with n Particles with radii radius.

        radius can be a single value or a sequence with n values.

        Any key-value pairs passed in the styles dictionary will be passed
        as arguments to Matplotlib's Circle patch constructor when drawing
        the Particles.

        """

        self.dt = 0.5
        self.death_ratio = 0.002
        self.infect_ratio = 0.9
        self.contact_distance = 0.03
        self.iterations = iterations
        self.current_infections = 0
        self.hospital = hospital
        self.hospital_capacity = 15
        self.hospital_num = 0
        self.init_particles(n, infectn, radius)

    def place_particle(self, rad, status, color):
        # # Choose x, y so that the Particle is entirely inside the
        # # domain of the simulation.
        if self.hospital:
            x = np.random.uniform(0.1+rad, 1-rad)
            y = np.random.uniform(rad, 1-rad)
        # x, y = rad + (1 - 2*rad) * np.random.random(2)
        else:
            x = np.random.uniform(rad, 1-rad)
            y = np.random.uniform(rad, 1-rad)
        # x, y = self.random_position()

        # Choose a random velocity (within some reasonable range of
        # values) for the Particle.
        vr = 0.1 * np.sqrt(np.random.random()) + 0.05
        vphi = 2*np.pi * np.random.random()
        vx, vy = vr * np.cos(vphi), vr * np.sin(vphi) * 0.5
        # vx = np.random.random() * 0.05
        # vy = np.random.random() * 0.05

        particle = self.ParticleClass(x, y, vx, vy, rad, status, color)
        # Check that the Particle doesn't overlap one that's already
        # been placed.
        # for p2 in self.particles:
        #     if p2.overlaps(particle):
        #         break
        # else:
        self.particles.append(particle)
        #     return True
        # return False

        

    def init_particles(self, n, infectn, radius):
        """Initialize the n Particles of the simulation.

        Positions and velocities are chosen randomly; radius can be a single
        value or a sequence with n values.

        """


        self.n = n
        self.infectn = infectn
        self.particles = []

        # for j in range(self.infectn):
        #     self.place_particle(radius, STATUS[1], COLORS[1])
        for i in np.arange(self.infectn):
            self.place_particle(radius, STATUS[1], COLORS[1])
        for i in np.arange(self.n - len(self.particles)):
            self.place_particle(radius, STATUS[0], COLORS[0])
        
        
    def infectionDetect(self, a1, a2):
        
        
        if a1.status == STATUS[0] and a2.status == STATUS[1]:
            infect_test = np.random.random()
            if infect_test < self.infect_ratio:
                a1.status = STATUS[1]
                a1.color = COLORS[1]
                a1.draw(self.ax)
                self.current_infections += 1
        
    def handle_collisions(self):
        """Detect and handle any collisions between the Particles.

        When two Particles collide, they do so elastically: their velocities
        change such that both energy and momentum are conserved.

        """ 

        # We're going to need a sequence of all of the pairs of particles when
        # we are detecting collisions. combinations generates pairs of indexes
        # into the self.particles list of Particles on the fly.
        pairs = combinations(range(self.n), 2)
        for i,j in pairs:
            # if self.particles[i].status == STATUS[4] or self.particles[j].status == STATUS[4]:
            #     return
            if self.particles[i].overlaps(self.particles[j]):
                self.infectionDetect(self.particles[i], self.particles[j])
                self.infectionDetect(self.particles[j], self.particles[i])


    def handle_boundary_collisions(self, p):
        """Bounce the particles off the walls elastically."""

        if p.x - p.radius < 0.1:
            #p.x = p.radius
            p.vx = -p.vx
        if p.x + p.radius > 1:
            #p.x = 1-p.radius
            p.vx = -p.vx
        if p.y - p.radius < 0:
            #p.y = p.radius
            p.vy = -p.vy
        if p.y + p.radius > 1:
            #p.y = 1-p.radius
            p.vy = -p.vy
        

    def update(self, p):
        # status change after being infected (except susceptible status)
        if p.status == STATUS[3]:
            return

        if self.hospital:
            if p.status == STATUS[1]:
                p.infect_time += 1
                p.detect_time += 1
                death_test = np.random.random()
                if death_test < self.death_ratio:
                    p.status = STATUS[3]
                    p.color = COLORS[3]
                    p.draw(self.ax)
                    p.infect_time = 0
                    p.detect_time = 0
                    self.current_infections -= 1
                    return

                if p.infect_time > 24:
                    p.status = STATUS[2]
                    p.color = COLORS[2]
                    p.draw(self.ax)
                    p.infect_time = 0
                    self.current_infections -= 1  

                hospital_test = np.random.random()
                if p.detect_time > 14 and hospital_test < 0.8 and self.hospital_num < self.hospital_capacity:
                    p.status = STATUS[4]
                    p.x = np.random.uniform(p.radius, 0.1 - p.radius)
                    p.y = np.random.uniform(p.radius, 1 - p.radius)
                    self.hospital_num += 1
                    


            if p.status == STATUS[4]:
                p.hospital_time += 1
                p.infect_time += 1
                death_testh = np.random.random()
                if death_testh < self.death_ratio / 4:
                    p.status = STATUS[3]
                    p.color = COLORS[3]
                    p.draw(self.ax)
                    p.infect_time = 0
                    p.hospital_time = 0
                    self.current_infections -= 1
                    self.hospital_num -= 1
                    return

                if p.infect_time > 24 or p.hospital_time > 10:
                    p.status = STATUS[2]
                    p.color = COLORS[2]
                    p.x = np.random.uniform(0.1 + p.radius, 1 - p.radius)
                    p.y = np.random.uniform(p.radius, 1 - p.radius)
                    p.draw(self.ax)
                    p.infect_time = 0
                    p.hospital_time = 0
                    self.current_infections -= 1
                    self.hospital_num -= 1
            

        else:
            if p.status == STATUS[1]:
                p.infect_time += 1
                death_test = np.random.random()
                if death_test < self.death_ratio:
                    p.status = STATUS[3]
                    p.color = COLORS[3]
                    p.draw(self.ax)
                    p.infect_time = 0
                    self.current_infections -= 1
                    return
                
                if p.infect_time > 24:
                    p.status = STATUS[2]
                    p.color = COLORS[2]
                    p.draw(self.ax)
                    p.infect_time = 0
                    self.current_infections -= 1

        # if p.status == STATUS[1]:
        #     p.infect_time += 1
        #     death_test = np.random.random()
        #     if death_test < self.death_ratio:
        #         p.status = STATUS[3]
        #         p.color = COLORS[3]
        #         p.draw(self.ax)
        #         p.infect_time = 0
        #         self.current_infections -= 1
        #         return
            
        #     if p.infect_time > 24:
        #         p.status = STATUS[2]
        #         p.color = COLORS[2]
        #         p.draw(self.ax)
        #         p.infect_time = 0
        #         self.current_infections -= 1


    def advance_animation(self):
        """Advance the animation by dt, returning the updated Circles list."""

        for p in self.particles:
            p.advance(self.dt)
            self.handle_boundary_collisions(p)
            self.update(p)
              
        self.handle_collisions()
        # print(self.current_infections, self.hospital_num)
        return self.circles


    def init(self):
        """Initialize the Matplotlib animation."""

        self.circles = []
        for particle in self.particles:
            self.circles.append(particle.draw(self.ax))
        

        return self.circles

    def animate(self, i):
        """The function passed to Matplotlib's FuncAnimation routine."""

        self.advance_animation()
        return self.circles

    def setup_animation(self):
        self.fig, self.ax = plt.subplots()
        for s in ['top','bottom','left','right']:
            self.ax.spines[s].set_linewidth(2)
        self.ax.set_aspect('equal', 'box')
        self.ax.set_xlim(0, 1)
        self.ax.set_ylim(0, 1)
        self.ax.set_title('ABM')

        self.ax.axvline(x=0.1,ls="-",c="black")
        # self.ax.xaxis.set_ticks([])
        # self.ax.yaxis.set_ticks([])

    def save_or_show_animation(self, anim, save, filename='hospital.mp4'):
        if save:
            Writer = animation.writers['ffmpeg']
            writer = Writer(fps=4, bitrate=1800)
            anim.save(filename, writer=writer)
        else:
            plt.show()
            
    def do_animation(self, save=False, filename='hospital.mp4'):
        """Set up and carry out the animation of the molecular dynamics.

        To save the animation as a MP4 movie, set save=True.
        """

        self.setup_animation()
        anim = animation.FuncAnimation(self.fig, self.animate,
                init_func=self.init, frames=self.iterations, interval=250, blit=False)
        self.save_or_show_animation(anim, save, filename)
        

if __name__ == '__main__':
    nparticles = 100
    radius = 0.012
    infectn = 2
    iterations = 120
    sim = Simulation(nparticles, infectn, radius, iterations, hospital=True)
    sim.do_animation(save=False)
