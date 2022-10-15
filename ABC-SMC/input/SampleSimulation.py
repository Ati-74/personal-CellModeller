import random
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium_reg_param import CLBacterium_reg_param
from CellModeller.GUI import Renderers
import numpy
from scipy.stats import lognorm


def setup(sim):
    sim.dt = 0.025
    # Set biophysics module
    biophys = CLBacterium_reg_param(sim, jitter_z=False, gamma=10 , reg_param=0.1)

    # Set up regulation module
    regul = ModuleRegulator(sim, sim.moduleName)
    # Only biophysics and regulation
    sim.init(biophys, regul, None, None)

    # Specify the initial cell and its location in the simulation
    sim.addCell(cellType=0,
                pos=(37.684918033, 39.40642623, 0),
                length=5.061029829,
                dir=(-0.658, 0.753, 0),
                rad=0.657250262)

    # Add some objects to draw the models
    therenderer = Renderers.GLBacteriumRenderer(sim)
    sim.addRenderer(therenderer)

    # Specify how often data is saved
    sim.pickleSteps = 50
    sim.saveOutput = True


def init(cell):
    # Specify mean and distribution of initial cell size
    cell.targetVol = lognorm.rvs(0.266, -0.681, 6.22, size=1)[0]
    # Specify growth rate of cells
    cell.growthRate = lognorm.rvs(0.355, -0.015, 0.053, size=1)[0]
    cell.color = (1.0, 1.0, 0.0)


def update(cells):
    # Iterate through each cell and flag cells that reach target size for division
    for (id, cell) in cells.items():
        cell.growthRate = lognorm.rvs(0.355, -0.015, 0.053, size=1)[0]
        if cell.volume > cell.targetVol:
            cell.divideFlag = True


def divide(parent, d1, d2):
    # Specify target cell size that triggers cell division
    d1.targetVol = lognorm.rvs(0.266, -0.681, 6.22, size=1)[0]
    d2.targetVol = lognorm.rvs(0.266, -0.681, 6.22, size=1)[0]
    # radius
    d1.rad = lognorm.rvs(0.248, 0.184, 0.252, size=1)[0]
    d2.rad = lognorm.rvs(0.248, 0.184, 0.252, size=1)[0]
    # growth rate
    d1.growthRate = lognorm.rvs(0.355, -0.015, 0.053, size=1)[0]
    d2.growthRate = lognorm.rvs(0.355, -0.015, 0.053, size=1)[0]

