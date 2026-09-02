"""Script to build the composant of the protein"""
from loguru import logger
import math 

class protein:
    "Class to build the protein composant"
    def __init__(self, numres):
        self.numres = numres

class residues:
    "Class to build the residue component"
    def __init__(self, name, number):
        self.name = name
        self.number = number
        self.atomlist = []

class sphere:
    "Class to build the sphere component"
    # We set for every atom the water radius to 1.4 Angstroms
    water_vdw = 1.4
    def __init__(self, radius, nbpoints = 400):
        self.radius = radius + self.water_vdw
        self.nbpoints = nbpoints
        self.pointlits = []
    
    def compute_points_coordinate(self):
        # We define the first angle
        first_heigth = -1
        logger.info(f"Point n° 1 | heigth: {first_heigth}")
        first_teta_angle = math.pi
        logger.info(f"Point n° 1 | heigth: {first_teta_angle}")
        first_phi_angle = 0
        logger.info(f"Point n° 1 | heigth: {first_phi_angle}")
        
        # We compute the coordinate of the first point 
        x = self.radius *  math.sin(first_teta_angle) * math.cos(first_phi_angle)
        y = self.radius *  math.sin(first_teta_angle) * math.sin(first_phi_angle)
        z = self.radius *  math.cos(first_teta_angle) 
        
        self.pointlits.append((x, y, z))
        logger.info(f"Added point n° 1 | coordinates : {(x, y, z)}")
        
        previous_phi_angle = first_phi_angle
        # We treat the other points
        for k in range(2, self.nbpoints):
            logger.info(f"Working on the point n° {k}")
            # First we compute the heigth of the sub spheres
            curent_heigth = -1 + 2 * (k - 1)/(self.nbpoints - 1)
            logger.info(f"Point n° {k} | heigth: {curent_heigth}")
            # Second the compute the teta metric
            current_teta_angle = math.acos(curent_heigth)
            logger.info(f"Point n° {k} | téta angle : {current_teta_angle}")
            # Then we calculate the rotation of the angle
            current_phi_angle = (previous_phi_angle + 3.6 / math.sqrt(self.nbpoints) * 
                                 1 / math.sqrt(1-curent_heigth**2)) % (2*math.pi)
            previous_phi_angle = current_phi_angle
            logger.info(f"Point n° {k} | phi angle : {current_phi_angle}")
            
            x = self.radius *  math.sin(current_teta_angle) * math.cos(current_phi_angle)
            y = self.radius *  math.sin(current_teta_angle) * math.sin(current_phi_angle)
            z = self.radius *  math.cos(current_teta_angle) 
        
            self.pointlits.append((x, y, z))
            logger.info(f"Added point n° {k} | coordinates : {(x, y, z)}")
            
        # We take care of the last point
        last_heigth = -1 + 2 * (self.nbpoints - 1) / (self.nbpoints - 1)
        logger.info(f"Point n° {self.nbpoints} | heigth: {last_heigth}")
        last_teta_angle = math.acos(last_heigth)
        logger.info(f"Point n° {self.nbpoints} | heigth: {last_teta_angle}")
        last_phi_angle = 0
        logger.info(f"Point n° {self.nbpoints} | heigth: {last_phi_angle}")
        
        x = self.radius *  math.sin(last_teta_angle) * math.cos(last_phi_angle)
        y = self.radius *  math.sin(last_teta_angle) * math.sin(last_phi_angle)
        z = self.radius *  math.cos(last_teta_angle) 
        
        self.pointlits.append((x, y, z))        
        
        logger.info(f"Added point n° {self.nbpoints} | coordinates : {(x, y, z)}")
        

        
        
        
class atom:
    "Class to build the atom component"
    def __init__(self, type, atonum, atomres, x, y , z, sphere):
        self.type = type
        self.atonum = atonum
        self.atomres = atomres
        self.x = x
        self.y = y
        self.z = z
        self.sphere = sphere
