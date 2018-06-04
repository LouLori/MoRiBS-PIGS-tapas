from MMTK import *
from MMTK.ForceFields import LennardJonesForceField
from MMTK.Environment import NoseThermostat, AndersenBarostat
import string
from Scientific.IO import ArrayIO
from Scientific.Statistics import *
from Scientific.Statistics.Histogram import *
from numpy import *
import math

pos1=Vector(0.0,1.0,0.0)
pos2=Vector(0.0,0.0,0.0)
pos3=Vector(1.0,0.0,0.0)
pos4=Vector(1.0,-1.0,1.0)

universe =InfiniteUniverse(LennardJonesForceField())
atom1=Atom('Ar', position=pos1)
atom2=Atom('Ar', position=pos2)
atom3=Atom('Ar', position=pos3)
atom4=Atom('Ar', position=pos4)
universe.addObject(atom1)
universe.addObject(atom2)
universe.addObject(atom3)
universe.addObject(atom4)
print universe.dihedral(atom1,atom2,atom3,atom4)*180/math.pi
