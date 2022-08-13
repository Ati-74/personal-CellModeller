import random
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
from CellModeller.GUI import Renderers
import numpy
import math

delta = 1.9
delta_sig = 0.45


def setup(sim):
    # Set biophysics module
    # controls the physical properties and interactions of the bacterial cells.
    # This module is connected to the simulator object, and is initialized with the default parameters
    # jitter_z=False keeps the simulation constrained to 2 dimensions by preventing movement in the z-axis.
    # another parameters:
    # max_cells=100000
    # compNeighbours=False
    biophys = CLBacterium(sim, jitter_z=False, compNeighbours=True)

    # Set up regulation module
    # This module connects this model file to the simulator object.
    regul = ModuleRegulator(sim, sim.moduleName)
    # Only biophysics and regulation
    sim.init(biophys, regul, None, None)

    # Specify the initial cell and its location in the simulation
    # Susceptible
    sim.addCell(cellType=0, pos=(-8, 0, 0), dir=(random.gauss(1, -1), random.gauss(1, -1), 0),
                length=random.gauss(1, 0.5))
    # Attacker
    sim.addCell(cellType=1, pos=(8, 0, 0), dir=(random.gauss(1, -1), random.gauss(1, -1), 0),
                length=random.gauss(1, 0.5))

    # Add some objects to draw the models
    therenderer = Renderers.GLBacteriumRenderer(sim)
    sim.addRenderer(therenderer)

    # Specify how often data is saved
    sim.pickleSteps = 50


def init(cell):
    #  state of the cell when it is created – either at the beginning of the simulation
    #  as an initial condition or when a cell is born through division.
    # Specify mean and distribution of initial cell size
    # Specify mean and distribution of initial cell size
    cell.targetVol = cell.length + random.gauss(delta, delta_sig)
    # Specify growth rate of cells
    cell.growthRate = 1.0
    if cell.cellType == 0:
        cell.color = (1.0, 1.0, 0.0)  # Susceptible
    elif cell.cellType == 1:
        cell.color = (1.0, 0.0, 0.0)  # Attacker
    elif cell.cellType == 2:
        cell.color = (0.0, 0.0, 0.0)


def update(cells):
    # Iterate through each cell and flag cells that reach target size for division
    for (id, cell) in cells.items():
        colliding = False
        if cell.cellType == 0:  # if a cell is a Susceptible,
            if cell.volume > cell.targetVol:
                cell.divideFlag = True
            for index in cell.neighbours:  # loop through all contacts
                if cells[index].cellType == 1:  # if Attacker is in contact
                    colliding = True
            if colliding:
                cell.cellType = 2  # become transconjugant
                cell.growthRate = 0.0
                cell.divideFlag = False
                cell.cellAge = 0
        elif cell.cellType == 1:  # if a cell is an Attacker,
            if cell.volume > cell.targetVol:
                cell.divideFlag = True
        elif cell.cellType == 2:
            life_threshold = 10.0 + random.uniform(0.0, 3)
            if cell.cellAge > life_threshold:
                cell.deathFlag = True

        if cell.cellType == 0:
            cell.color = (1.0, 1.0, 0.0)  # Susceptible
        elif cell.cellType == 1:
            cell.color = (1.0, 0.0, 0.0)  # Attacker
        elif cell.cellType == 2:
            cell.color = (0.0, 0.0, 0.0)


def divide(parent, d1, d2):
    # Specify target cell size that triggers cell division
    d1.targetVol = d1.length + random.gauss(delta, delta_sig)
    d2.targetVol = d2.length + random.gauss(delta, delta_sig)
