# This simulation code is compatible with CompuCell3D ver. 3.7.9
# Vedang Narain, 29 April 2019

from PySteppables import *
import CompuCell
import sys
import numpy as np

from PySteppablesExamples import MitosisSteppableBase
            
# set starting target volume and lambda volume
class ConstraintInitializerSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        
    def start(self):
        for cell in self.cellList:
            cell.targetVolume=25
            cell.lambdaVolume=2.0
 
class GrowthSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    
    # define plots
    def start(self):
        self.pW = self.addNewPlotWindow(_title='Cell Populations',_xAxisTitle='MonteCarlo Step (MCS)',_yAxisTitle='Number', _xScaleType='linear',_yScaleType='linear')
        self.pW.addPlot('population_PROLIFERATING',_style='Lines',_color='green',_size=3)        
        self.pW.addPlot('population_QUIESCENT',_style='Lines',_color='yellow',_size=3)        
        self.pW.addPlot('population_NECROTIC',_style='Lines',_color='red',_size=3)        
        self.pW.addPlot('population_total',_style='Lines',_color='grey',_size=1)
        self.pW.addPlot('population_CTL',_style='Lines',_color='blue',_size=3)  
    
    # initialize field variables
    def step(self,mcs):
        
        # Make sure Secretion plugin is loaded
        # make sure this field is defined in one of the PDE solvers
        field_survival = self.getConcentrationField('Survival')      
        secretor_survival = self.getFieldSecretor('Survival')
        field_growth = self.getConcentrationField('Growth')      
        secretor_growth = self.getFieldSecretor('Growth')
        field_attractant = self.getConcentrationField('Attractant')      
        secretor_attractant = self.getFieldSecretor('Attractant')   
                   
        # track cell populations
        if not mcs%10: # every 10 mcs steps       
            total_population_PROLIFERATING = 0; total_population_QUIESCENT = 0; total_population_NECROTIC = 0; total_population_CTL = 0; total_population_tumor = 0
            for cell in self.cellListByType(self.PROLIFERATING):
                total_population_PROLIFERATING+=1
            self.pW.addDataPoint('population_PROLIFERATING', mcs, total_population_PROLIFERATING)
            for cell in self.cellListByType(self.QUIESCENT):
                total_population_QUIESCENT+=1
            self.pW.addDataPoint('population_QUIESCENT', mcs, total_population_QUIESCENT)
            for cell in self.cellListByType(self.NECROTIC):
                total_population_NECROTIC+=1
            self.pW.addDataPoint('population_NECROTIC', mcs, total_population_NECROTIC)
            for cell in self.cellListByType(self.CTL):
                total_population_CTL+=1
            self.pW.addDataPoint('population_CTL', mcs, total_population_CTL)
            total_population_tumor = total_population_PROLIFERATING + total_population_QUIESCENT + total_population_NECROTIC
            self.pW.addDataPoint('population_total', mcs, total_population_tumor)
        
        # outline nutrient-dependent behavior
        for cell in self.cellListByType(self.PROLIFERATING,self.QUIESCENT): # iterate over PROLIFERATING and QUIESCENT cells
            res_growth = secretor_growth.uptakeInsideCellTotalCount(cell, 1000.0, 0.002) # arguments are: cell, max uptake, relative uptake
            res_survival = secretor_survival.uptakeInsideCellTotalCount(cell, 1000.0, 0.002) # arguments are: cell, max uptake, relative uptake
            growth_concentration_within = -(res_growth.tot_amount) # since amount is negative by default
            survival_concentration_within = -(res_survival.tot_amount) # since amount is negative by default
            if cell.type == 1 and growth_concentration_within >=0.2: # if cell is PROLIFERATING with sufficient growth factor
                cell.targetVolume+=(1*growth_concentration_within)/(2+growth_concentration_within) # increase the target volume
            if cell.type == 2 and growth_concentration_within >= 0.2: # if cell is QUIESCENT with sufficient growth factor
                cell.type = 1 # QUIESCENT cell becomes PROLIFERATING
                cell.targetVolume+=(1*growth_concentration_within)/(2+growth_concentration_within) # increase the target volume
            if cell.type == 1 and growth_concentration_within < 0.2: # if cell is PROLIFERATING with insufficient growth factor
                cell.type = 2 # PROLIFERATING cell becomes QUIESCENT
                cell.targetVolume = cell.volume # freeze the current volume
            if survival_concentration_within < 0.05: # if there is insufficient survival factor
                cell.type = 3 # cell becomes NECROTIC
                cell.targetVolume *= 0.98 # make the cell shrink

        # iterate over NECROTIC cells and make them shrink   
        for cell in self.cellListByType(self.NECROTIC):
            cell.targetVolume*=0.98 # update the shrinking volume

        # seed CTL cells
        immune_strength = 5.0                       # enter desired strength of immune response, i.e., number of CTLs arriving per wave
        CTL_size = 4                                # enter desired CTL size
        lattice_size = 250                          # enter square lattice size here
        edge = lattice_size - CTL_size              # used in cell seeding
        if not mcs%50:                              # if there is no remainder, i.e., every x mcs steps
            border_attractant_raw = field_attractant[0, 0, 0] + field_attractant[125, 0, 0] + field_attractant[250, 0, 0] + field_attractant[250, 125, 0] \
                                    + field_attractant[250, 250, 0] + field_attractant[125, 250, 0] + field_attractant[0, 250, 0] + field_attractant[0, 125, 0]
            border_attractant_float = (immune_strength*border_attractant_raw)/(5.0+border_attractant_raw)
            border_attractant = int(round(border_attractant_float)) # get rough estimate of border values
            if border_attractant >= 1:                                           # once there's a sufficient amount of attractant at border
                for count in range(border_attractant):                           # seed a proportional number of CTLs
                    while True:                                                  # infinite loop
                        x = np.random.random() * lattice_size                    # random_factor * lattice_size
                        y = np.random.random() * lattice_size                    # random_factor * lattice_size
                        cell_0 = self.cellField[x, y, 0]                         # attempt to choose random cell field
                        if not cell_0 and not ((CTL_size < x < edge) and (CTL_size < y < edge)): # if there is no cell there and the coordinates are not away from the edge
                            break                                                # exit loop
                    cell = self.newCell(self.CTL)                                # create a new cell of type PROLIFERATING
                    self.cellField[x:x + CTL_size, y:y + CTL_size, 0] = cell     # dimensions of cell will be size x size x 1
                    cell.targetVolume = 25.0                                     # set target volume
                    cell.lambdaVolume = 40.0                                     # set lambda volume
                    cell.dict['kills'] = 0                                       # give CTL new attribute 'kills' to track number of tumor cells killed
         
class MitosisSteppable(MitosisSteppableBase):
    def __init__(self,_simulator,_frequency=1):
        MitosisSteppableBase.__init__(self,_simulator, _frequency)    
    def step(self,mcs):
        cells_to_divide = []       
        
        # divide cell once volume > 50
        for cell in self.cellList:
            if cell.volume > 50:
                cells_to_divide.append(cell)                
        for cell in cells_to_divide:
            self.divideCellRandomOrientation(cell)

    def updateAttributes(self):
        self.parentCell.targetVolume /= 2.0 # reduce parent target volume                 
        self.cloneParent2Child()            
               
        # when PROLIFERATING cells divide, they form two QUIESCENT cells
        if self.parentCell.type==1:
            self.childCell.type=2
            self.parentCell.type=2
        
class DeathSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def step(self,mcs):
        
        chemo_mcs = 2000 # enter point of chemotherapy administration
       
       # CTL cells make tumor cells necrotic    
        for cell in self.cellListByType(self.CTL):
            kills = cell.dict['kills']
            prob_death = kills/(10.0+kills)     # probability of death is directly proportional to kill count
            if mcs < chemo_mcs:                 # if chemotherapy not yet administered
                prob_kill = 1.0 - prob_death    # probability of killing is inversely proportional to kill count
            elif mcs >= chemo_mcs:              # if chemotherapy administered
                prob_kill = 0.5 * (1.0 - prob_death) # probability of killing is halved
            for neighbor, commonSurfaceArea in self.getCellNeighborDataList(cell): # iterate over cell's neighbors
                if (neighbor and neighbor.type < 3) and (np.random.random() < prob_kill): # if neighbor is present and type PROLIFERATING or QUIESCENT, roll the killing dice
                    neighbor.type = 3 # make cell NECROTIC
                    cell.dict['kills'] += 1 # increment CTL kill count by 1
                    if np.random.random() < prob_death: # roll the death dice
                        cell.targetVolume=0           
                        cell.lambdaVolume=100
                    break # move onto the next cell
                
        # chemotherapy
        if mcs==chemo_mcs: # simulate chemotherapeutic perturbation
            for cell in self.cellListByType(self.CTL):
                if np.random.random() < 0.9: # kill all cells with 90% probability
                    cell.targetVolume=0
                    cell.lambdaVolume=100
            for cell in self.cellListByType(self.PROLIFERATING, self.QUIESCENT):  
                if np.random.random() < 0.9: # kill all cells with 90% probability
                    cell.type = 3 # cell becomes NECROTIC
                    cell.targetVolume *= 0.98
