import sys

import gsd.hoomd  # Export for Ovito
import numpy as np


# Function to dump simulation frame that is readable in Ovito
# Also stores radii and velocities in a compressed format which is nice


def create_frame(radii, velocities, sigma, box, n_particle, frame, filename):
    # Set up GSD Simulation
    t = gsd.hoomd.open(name='GSD/Simulation_' + str(filename) + '.gsd', mode='wb')

    # Particle positions, velocities, diameter
    radii = box * radii
    partpos = radii.tolist()
    velocities = velocities.tolist()
    diameter = sigma * np.ones((n_particle,))
    diameter = diameter.tolist()

    # Now make gsd file
    s = gsd.hoomd.Snapshot()
    s.configuration.step = frame
    s.particles.N = n_particle
    s.particles.position = partpos
    s.particles.velocity = velocities
    s.particles.diameter = diameter
    s.configuration.box = [box, box, box, 0, 0, 0]

    # Append results
    t.append(s)

    return 0  # Finished execution of GSD Output


# create_frame(radii, velocities,sigma, box, n_particle, step / n_dump))


total = len(sys.argv)
cmdargs = str(sys.argv)

print("The total numbers of args passed to the script: %d " % total)
print("Args list: %s " % cmdargs)


# Parsing args one by one
# "Script name: %s" % str(sys.argv[0]))

print("First argument: %s" % str(sys.argv[1]))
print("Second argument: %s" % str(sys.argv[2]))
