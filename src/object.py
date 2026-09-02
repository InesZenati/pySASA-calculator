"""Script to build the composant of the protein"""

class protein:
    "Class to build the protein composant"
    def __init__(self, numres):
        self.numres = numres

class residues:
    "Class to build the residue component"
    def __init__(self, name, number, atomlist):
        self.name = name
        self.number = number
        self.atomlist = atomlist

class sphere:
    "Class to build the sphere component"
    # We set for every atom the water radius to 1.4 Angstroms
    water_vdw = 1.4
    def __init__(self, radius, pointlits):
        self.radius = radius + self.water_vdw
        self.pointlits = pointlits
    
    # def compute_points_coordinate(radius):
        
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
