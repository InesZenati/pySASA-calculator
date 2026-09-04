"""Script to build the composant of the protein"""
from loguru import logger
import math 

class Protein:
    "Class to build the protein composant"
    def __init__(self, pdbname):
        self.pdbname = pdbname
        self.residuelist = []

class Residues:
    "Class to build the residue component"
    def __init__(self, name, number):
        self.name = name
        self.number = number
        self.atomlist = []
        
    

class Sphere:
    "Class to build the sphere component"
    # We set for every atom the water radius to 1.4 Angstroms
    water_vdw = 1.4
    def __init__(self, radius, nbpoints = 400):
        self.radius = radius + self.water_vdw
        self.nbpoints = nbpoints
        self.pointlits = []
        
    def compute_points_coordinate(self):
        """Compute the point coordinate on the sphere surface."""
        # We define the first angle
        first_heigth = -1
        logger.info(f"Point n° 1 | heigth: {first_heigth}")
        first_teta_angle = math.pi
        logger.info(f"Point n° 1 | heigth: {first_teta_angle}")
        first_phi_angle = 0
        logger.info(f"Point n° 1 | heigth: {first_phi_angle}")
        
        # We compute the coordinate of the first point 
        x_point = self.radius *  math.sin(first_teta_angle) * math.cos(first_phi_angle)
        y_point = self.radius *  math.sin(first_teta_angle) * math.sin(first_phi_angle)
        z_point = self.radius *  math.cos(first_teta_angle) 

        self.pointlits.append((x_point, y_point, z_point))
        logger.info(f"Added point n° 1 | coordinates : {(x_point, y_point, z_point)}")
        
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
            
            x_point = self.radius *  math.sin(current_teta_angle) * math.cos(current_phi_angle)
            y_point = self.radius *  math.sin(current_teta_angle) * math.sin(current_phi_angle)
            z_point = self.radius *  math.cos(current_teta_angle) 
        
            self.pointlits.append((x_point, y_point, z_point))
            logger.info(f"Added point n° {k} | coordinates : {(x_point, y_point, z_point)}")
            
        # We take care of the last point
        last_heigth = -1 + 2 * (self.nbpoints - 1) / (self.nbpoints - 1)
        logger.info(f"Point n° {self.nbpoints} | heigth: {last_heigth}")
        last_teta_angle = math.acos(last_heigth)
        logger.info(f"Point n° {self.nbpoints} | heigth: {last_teta_angle}")
        last_phi_angle = 0
        logger.info(f"Point n° {self.nbpoints} | heigth: {last_phi_angle}")
        
        # We translate the sphere points on the atom
        x_point = self.radius *  math.sin(last_teta_angle) * math.cos(last_phi_angle)
        y_point = self.radius *  math.sin(last_teta_angle) * math.sin(last_phi_angle)
        z_point = self.radius *  math.cos(last_teta_angle) 
        
        self.pointlits.append((x_point, y_point, z_point))
        
        logger.info(f"Added point n° {self.nbpoints} | coordinates : {(x_point, y_point, z_point)}")
        

class Atom:
    "Class to build the atom component"
    def __init__(self, type, atomres, x, y , z, sphere):
        self.type = type
        self.atomres = atomres
        self.x = x
        self.y = y
        self.z = z
        self.sphere = sphere
        
    def translate_points_on_atom(self):
        """We translate the sphere point on the atom."""
        for sphere_point in self.sphere.pointlits:
            # For each point we add the atom coordinate 
            sphere_point[0] += self.x
            sphere_point[1] += self.y
            sphere_point[2] += self.z
                  

    def compute_distance_from_point_to_atom(self, point, atom):
        """Compute the distance from a point to the atom center."""
        distance =  math.sqrt((atom.x - point[0])**2 
                   +(atom.y - point[1])**2
                   +(atom.z - point[2])**2)
        logger.info(f"Distance found : {distance}")
        return distance
        
        
    def is_occluded_atom(self, point, atom):
        """Find if an atom is occluded or not by another atom."""
        # We assume that if a point is closer to an atom center than its own point then 
        # the point of our atom is occluded
        distance = self.calculate_point_distance_to_atom(point, atom)
        if distance < atom.sphere.radius:
            return True
        return False
    
    def count_occluded_points(self, atom):    
        """Calculate the number of occluded points on the sphere by another atom."""
        occluded_points = 0
        for point in (self.sphere.nbpoints):
            if self.is_point_occluded(point, atom):
                logger.success(f"Point {point} is occluded by {atom.type}")
                occluded_points +=1
        return occluded_points
                
            

            
