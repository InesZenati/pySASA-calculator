"""Script to build the composant of the protein"""
from loguru import logger
import numpy as np
import math 

class Protein:
    """
    Class used to represent a protein.
    
    Instance Attributes
    -------------------
    pdbname : str
        The code of the protein extracted from the PDB database
    residuelist : list
        A list of Residues objects that make up the protein
    
    Methods
    -------
    get_all_atom(self)
        Returns a list of all Atom objects in the protein by iterating 
        through each Residues object in residuelist and collecting their atomlist.
    """
    def __init__(self, pdbname):
        """
        Construct a protein.
        
        Parameters
        ----------
        pdbname : str
            The code of the protein.
        residuelist : list
            A list of Residues objects that make up the protein.
        """
        self.pdbname = pdbname
        self.residuelist = []
    
    def get_all_atom(self):
        """
        Create a list of all Atom objects in the protein.
        
        Parameters
        ----------
        None
        
        Returns
        -------
        list:
            A list of all Atom objects in the protein.
        """
        all_atoms = []
        for residue in self.residuelist:
            for atom in residue.atomlist:
                all_atoms.append(atom)
        return all_atoms              

class Residues:
    """
    Class used to represent a residue.
    
    Instance Attributes
    -------------------
    name : str
        The name of the residue.
    number : int
        The number of the residue.
    atomlist : list
        A list of Atom objects that make up the residue.
    """
    def __init__(self, name, number, chain_id):
        """
        Construct a residue.
        
        Parameters
        ----------
        name : str
            The name of the residue.
        number : int
            The number of the residue.
        chain_id : str
            The chain identifier of the residue.
        atomlist : list
            A list of Atom objects that make up the residue.
        """
        self.name = name
        self.number = number
        self.chain_id = chain_id
        self.atomlist = []


class Sphere:
    """
    Class used to represent a sphere around an atom.
    
    Class Attributes
    ----------------
    water_vdw : float
        The van der Waals radius of water in Angstroms.
    
    Instance Attributes
    -------------------
    radius : float
        The radius of the sphere in Angstroms.
    nbpoints : int
        The number of points to be generated on the sphere surface.
    pointlist : list
        A list of tuples representing the coordinates of points on the sphere surface.
    occluded_points : list
        A boolean list indicating whether each point in pointlist is occluded or not.
        
    Methods
    -------
    compute_points_coordinate(self)
        Computes the coordinates of points on the sphere surface using spherical coordinates.
    """
    # We set for every atom the water radius to 1.4 Angstroms
    water_vdw = 1.4
    def __init__(self, radius, nbpoints = 92):
        """
        Construct a sphere around an atom.
        
        Parameters
        ----------
        radius : float
            The radius of the sphere in Angstroms.
        nbpoints : int
            The number of points to be generated on the sphere surface. Default is 92.
        pointlist : list
            A list of tuples representing the coordinates of points on the sphere surface.
        occluded_points : list
            A boolean list indicating whether each point in pointlist is occluded or not.
        """
        self.radius = radius + self.water_vdw
        self.nbpoints = nbpoints
        self.pointlist = []
        self.occluded_points = []
        self.compute_points_coordinate()
        
    def compute_points_coordinate(self):
        """
        Compute the point coordinate on the sphere surface.
        
        Generates evenly distributed points on a sphere's surface using the Fibonacci 
        sphere algorithm (golden angle spiral). For each point, it computes a height 
        value, derives the polar angle (theta) and azimuthal angle (phi), then converts 
        these spherical coordinates into (x, y, z) Cartesian coordinates centered at the origin.
                
        """
        # We define the first angle
        first_heigth = -1
        logger.debug(f"Point n° 1 | heigth: {first_heigth}")
        first_teta_angle = math.pi
        logger.debug(f"Point n° 1 | heigth: {first_teta_angle}")
        first_phi_angle = 0
        logger.debug(f"Point n° 1 | heigth: {first_phi_angle}")
        
        # We compute the coordinate of the first point 
        x_point = self.radius *  math.sin(first_teta_angle) * math.cos(first_phi_angle)
        y_point = self.radius *  math.sin(first_teta_angle) * math.sin(first_phi_angle)
        z_point = self.radius *  math.cos(first_teta_angle) 

        self.pointlist.append((x_point, y_point, z_point))
        logger.debug(f"Added point n° 1 | coordinates : {(x_point, y_point, z_point)}")
        
        previous_phi_angle = first_phi_angle
        # We treat the other points
        for k in range(2, self.nbpoints):
            logger.debug(f"Working on the point n° {k}")
            # First we compute the heigth of the sub spheres
            curent_heigth = -1 + 2 * (k - 1)/(self.nbpoints - 1)
            logger.debug(f"Point n° {k} | heigth: {curent_heigth}")
            # Second the compute the teta metric
            current_teta_angle = math.acos(curent_heigth)
            logger.debug(f"Point n° {k} | téta angle : {current_teta_angle}")
            # Then we calculate the rotation of the angle
            current_phi_angle = (previous_phi_angle + 3.6 / math.sqrt(self.nbpoints) * 
                                 1 / math.sqrt(1-curent_heigth**2)) % (2*math.pi)
            previous_phi_angle = current_phi_angle
            logger.debug(f"Point n° {k} | phi angle : {current_phi_angle}")
            
            x_point = self.radius *  math.sin(current_teta_angle) * math.cos(current_phi_angle)
            y_point = self.radius *  math.sin(current_teta_angle) * math.sin(current_phi_angle)
            z_point = self.radius *  math.cos(current_teta_angle) 
        
            self.pointlist.append((x_point, y_point, z_point))
            logger.debug(f"Added point n° {k} | coordinates : {(x_point, y_point, z_point)}")
            
        # We take care of the last point
        last_heigth = -1 + 2 * (self.nbpoints - 1) / (self.nbpoints - 1)
        logger.debug(f"Point n° {self.nbpoints} | heigth: {last_heigth}")
        last_teta_angle = math.acos(last_heigth)
        logger.debug(f"Point n° {self.nbpoints} | heigth: {last_teta_angle}")
        last_phi_angle = 0
        logger.debug(f"Point n° {self.nbpoints} | heigth: {last_phi_angle}")
        
        # We translate the sphere points on the atom
        x_point = self.radius *  math.sin(last_teta_angle) * math.cos(last_phi_angle)
        y_point = self.radius *  math.sin(last_teta_angle) * math.sin(last_phi_angle)
        z_point = self.radius *  math.cos(last_teta_angle) 
        
        self.pointlist.append((x_point, y_point, z_point))
        
        logger.debug(f"Added point n° {self.nbpoints} | "
                     f"coordinates : {(x_point, y_point, z_point)}")
        

class Atom:
    """
    Class used to represent an atom.
    
    Instance attributes
    -------------------
    type : str
        The type of the atom.
    atomres : Residues
        The residue object associated with the atom.
    x : float
        X position.
    y : float
        Y position.
    z : float
        Z position.
    sphere : Sphere
        An instance of the Sphere class representing the atom's surrounding sphere.
    
    Methods
    -------
    translate_points_on_atom(self)
        Translates the points on the sphere to be centered on the atom's coordinates.
    compute_distance_from_point_to_atom(self, point, atom)
        Computes the Euclidean distance from a given point of our atom to the another
        atom's position.
    is_occluded_atom(self, atom)
        Determines if the atom is occluded by another atom based on their positions and 
        radiuss.
    count_occluded_points(self, atom)
        Counts the number of points on the atom's sphere that are occluded by another 
        atom.
    detect_occluded_point(self, atomlist)
        Detects which points on the atom's sphere are occluded by a list of other atoms
        and fill the occluded_pointlist attribute with the corresponding boolean 
        (True / False).        
    """
    def __init__(self, type, atomres, x, y , z, sphere):
        """
        Constructor for the Atom class.
        
        Parameters
        ----------
        type : str
            The type of the atom (e.g., 'C', 'O', 'N').
        atomres : Residues
            The residue name associated with the atom (e.g., 'ALA', 'GLY').
        x : float
            X position.
        y : float
            Y position.
        z : float
            Z position.
        sphere : Sphere
            An instance of the Sphere class representing the atom's surrounding sphere.
        """
        self.type = type
        self.atomres = atomres
        self.x = x
        self.y = y
        self.z = z
        self.sphere = sphere
        self.translate_points_on_atom()
        
    def translate_points_on_atom(self):
        """
        Translates the points on the sphere to be centered on the atom's coordinates.
        
        Translates all sphere points from the origin to the atom's actual position 
        in 3D space, by adding the atom's (x, y, z) coordinates to each point. This 
        is called after the sphere's points have been computed, so that each 
        atom ends up with a sphere correctly centered on its own coordinates.
        """
        points_on_atom = []
        for sphere_point in self.sphere.pointlist:
            # For each point we add the atom coordinate 
            x_point_on_atom = sphere_point[0] + self.x
            y_point_on_atom = sphere_point[1] + self.y
            z_point_on_atom = sphere_point[2] + self.z
            points_on_atom.append((x_point_on_atom, y_point_on_atom, z_point_on_atom)) 
        self.sphere.pointlist = points_on_atom
                  

    def compute_distance_from_point_to_atom(self, point, atom):
        """
        Compute the distance from a point to the atom center.
        
        Computes the Euclidean distance between a 3D point (typically a sphere surface 
        point) and the center of a given atom. Used to determine whether that point 
        lies inside another atom's van der Waals sphere (i.e. whether it's occluded).
        """
        distance =  math.sqrt((atom.x - point[0])**2 
                   +(atom.y - point[1])**2
                   +(atom.z - point[2])**2)
        return distance
        
        
    def is_occluded_atom(self, point, otheratom):
        """
        Detect if an atom is occluded or not by another atom.
        
        Checks whether a single sphere point is occluded by a specific neighboring 
        atom, by comparing the point-to-atom distance against that neighbor's sphere 
        radius. Returns True if the point falls inside the neighbor's sphere (meaning 
        it's buried), False otherwise.
        """
        # We assume that if a point is closer to an atom center than its own point then 
        # the point of our atom is occluded
        distance = self.compute_distance_from_point_to_atom(point, otheratom)
        if distance < otheratom.sphere.radius:
            return True
        return False
    
    def count_occluded_points(self, otheratom):    
        """
        Calculate the number of occluded points on the current atom sphere by another 
        atom. [WIP]
        
        Counts how many points on the atom's own sphere are occluded by a single given 
        neighboring atom, by testing each point with is_occluded_atom. Meant to be 
        called once per neighbor and summed.
        """
        occluded_points = 0
        for point in self.sphere.pointlist:
            if self.is_occluded_atom(point, otheratom):
                occluded_points +=1
        return occluded_points
    
    def detect_occluded_point(self, atomlist):
        """
        Detect which points on the atom's sphere are occluded by a list of other atoms.
        
        Iterates over every point on the atom's own sphere and checks each one against 
        all atoms in atomlist, marking it as occluded (True) if any of them covers it, 
        free (False) otherwise. The result is stored point-by-point in 
        sphere.occluded_points, aligned by index with sphere.pointlits.
        """
        self.sphere.occluded_points = []  
        
        for point in self.sphere.pointlist:
            occluded = False
            for atom in atomlist:
                if atom is self: 
                    continue
                    
                if self.is_occluded_atom(point, atom):
                    occluded = True
                    break

            self.sphere.occluded_points.append(occluded)
            
    def get_coord(self):
        """
        Return the atom's coordinates as a numpy array.

        Used so that Bio.PDB.NeighborSearch can work directly with our own
        Atom objects, since it expects a get_coord method on each object.

        Return
        ------
        numpy.ndarray :
            The (x, y, z) coordinates of the atom.
        """
        return np.array([self.x, self.y, self.z])
            
        
                    
            

            
