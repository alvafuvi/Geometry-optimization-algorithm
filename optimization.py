import numpy as np
import math



class Atom:
    def __init__(self,id,element,x,y,z):
        """Atom class contains information for each atom

        Args:
            id (integer): Number lable of each atom in the molecule
            element (string): Either C or H in this case
            x (array(1,3)): x coordinate of each atom
            y (array(1,3)): y coordinate of each atom
            z (array(1,3)): z coordinate of each atom
        """
        self.id = id
        self.element = element
        self.x = x
        self.y = y
        self.z = z



def read_file(filename):
    """Reads .mol2 file and extracts coordinates, atoms and bonds

    Returns:
        bnd (array(bonds,2)): bond list. Format: [at1.id,at2.id]
        coord (array(nat,3)): cartesian coordinates.
        element (list(nat)): string list with atomic labels
        nat (integer): number of atoms  
    """
    

    mol = open(filename,'r')
    line = mol.readline()


    nat = int(line[1:3])
    bonds = int(line[4:7])
    print(f'This molecule has {nat} atoms and {bonds} bonds')
    mol.close()

    # Load the coordinates (X,Y,Z) and the list of atoms for each position

    coord = np.loadtxt(filename,skiprows=1,max_rows=nat,usecols=(0,1,2))
    
    element = np.loadtxt(filename,skiprows=1,max_rows=nat,usecols=3,dtype=str)
    





    # Load the bonds (made for a file with no exit line).
    #   First column = atom nº1
    #   Second column = atom nº2

    bnd = np.loadtxt(filename,skiprows=nat+1,usecols=(0,1),dtype=int)

    return bnd,coord,element,nat
    






#####################################
# Important parameters for geom opt #
#####################################


# Bond force constants (kcal/molA-2)

kcc = 300.
kch = 350.

# Angle force constants (kcal/moldegree-2)

khch = 35.
kcch = 35.
kccc = 60.

# Dihedral force constants (kcal/mol)

kxccx = 0.3 

# Equilibrium bond distances (A)

r0cc = 1.53
r0ch = 1.11

# Equilibrium angles (Degrees) 

a0hch = 109.5 * np.pi/180
a0cch = 109.5 * np.pi/180
a0ccc = 109.5 * np.pi/180

#  Oscillation constant

n = 3

# Van-Der-Waals parameters (kcal/mol and A)

epsh = 0.03
epsc = 0.07

epscc = epsc
epshh = epsh
epsch = np.sqrt(epsh*epsc)

sigmah = 1.20
sigmac = 1.75

sigmacc = 2*sigmac
sigmahh = 2*sigmah
sigmach = 2*np.sqrt(sigmah*sigmac)

Ahh =  4*epshh*sigmahh**12 
Ach = 4*epsch*sigmach**12
Acc = 4*epscc*sigmacc**12

Bhh = 4*epshh*sigmahh**6 
Bch = 4*epsch*sigmach**6
Bcc = 4*epscc*sigmacc**6


#########################################
# Internal coordinate system definition #
#########################################


# Cross product

def cross(r1,r2):
    """Cross product between two vectors

    Args:
        r1 (array(1,3)): x,y and z components of vector 1
        r2 (array(1,3)): x,y and z components of vector 2

    Returns:
        array(1,3): cross product vector
    """

    cp = np.array([r1[1]*r2[2] -r1[2]*r2[1],r1[2]*r2[0] -r1[0]*r2[2],r1[0]*r2[1] -r1[1]*r2[0]])
    return cp



def bond_length(i,j):
    """Calculates bond lengths between two atoms i and j

    Args:
        i (object atom): atom1, introduced from atoms (list of objects)
        j (object atom): atom2, introduced from atoms (list of objects)

    Returns:
        Float: array of scalar norms
        array(1,3): vector between atoms
    """
    rij = np.array([j.x-i.x, j.y-i.y, j.z-i.z])
    r = np.sqrt((rij[0])**2 + (rij[1])**2 + (rij[2])**2 ) # Scalar distance
    return r,rij # r=norm rij=vector


def bond_read(bnd):
    """Recognises bonds and calculates bond length and direction

    Args:
        bnd (array(bonds,2)): bond list

    Returns:
        rab_vect (array(3)): Vectorial direction of the bond
        bond_list (list of objects): list with information of the atoms involved in the bond and scalar distance
    """
    rab_vect = []
    bond_list = [] # Format: [at1,at2,distance] 

    for i in bnd:
        ra = atoms[i[0]-1]
        rb = atoms[i[1]-1]

        rab_vect.append(bond_length(ra,rb)[1])
        bond_list.append([ra,rb,bond_length(ra,rb)[0]])



    # for i in bond_list:
    #     print(f'{i[0].element,i[0].id} - {i[1].element,i[1].id} \t Distance in Angstroms = {i[2]}')

    return rab_vect,bond_list




def angle(at1,at2,at3):
    """Calculate angle (in radians) between three atoms and two bonds

    Args:
        at1 (object atom): atom1
        at2 (object atom): atom2. Central atom
        at3 (object atom): atom3

    Returns:
        theta (float): calculated angle in radians
    """
    dist21 = bond_length(at2,at1)[0]
    dist23= bond_length(at2,at3)[0]

    vect21 = bond_length(at2,at1)[1]
    vect23 = bond_length(at2,at3)[1]

    theta = np.arccos(np.dot(vect21,vect23) / (dist21*dist23))
    return theta


def angle_read(bnd):
    """Obtain angles 

    Args:
        bnd (array(bonds,2)): bond list

    Returns:
        angle_list (list of objects): list with information of the atoms involved in the angle and its value in radians
    """
    angle_list = [] # Format: [at1,at2,at3,angle] with at2 being the vertex

    for i in range(len(bnd)):
        vari = bnd[i]
        for j in range(i+1,len(bnd)):
            varj = bnd[j]
            
            if vari[0] == varj[0]:  # Central atom in the first column
                angle_list.append([atoms[vari[1]-1],atoms[vari[0]-1], atoms[varj[1]-1], angle(atoms[vari[1]-1],atoms[vari[0]-1], atoms[varj[1]-1])])
            
            elif vari[1] == varj[0]:  # Central atom in the first and second column
                angle_list.append([atoms[vari[0]-1],atoms[vari[1]-1], atoms[varj[1]-1], angle(atoms[vari[0]-1],atoms[vari[1]-1], atoms[varj[1]-1])])

            elif vari[1] == varj[1]:  # Central atom in the second column
                angle_list.append([atoms[vari[0]-1],atoms[vari[1]-1], atoms[varj[0]-1], angle(atoms[vari[0]-1],atoms[vari[1]-1], atoms[varj[0]-1])])


    # for i in angle_list:
    #     print(f'{i[0].element,i[0].id} - {i[1].element,i[1].id} - {i[2].element,i[2].id} \t Angle in radians = {i[3]}, angle in degrees = {i[3]*180/np.pi}')

    return angle_list



# Dihedral angle between four atoms. 

def dihedral(at1,at2,at3,at4):
    """Calculate dihedral angles

    Args:
        at1 (object atom): atom1
        at2 (object atom): atom2, central atom 
        at3 (object atom): atom3, central atom 
        at4 (object atom): atom4

    Returns:
        Float: dihedral angle (positive or negative)
    """

    dist23= bond_length(at2,at3)[0]


    vect12 = bond_length(at1,at2)[1]
    vect23 = bond_length(at2,at3)[1]
    vect34 = bond_length(at3,at4)[1]

    t = cross(vect12,vect23)
    u = cross(vect23,vect34)
    v = cross(t,u)

    tnorm = np.sqrt((t[0])**2 + (t[1])**2 + (t[2])**2 )
    unorm = np.sqrt((u[0])**2 + (u[1])**2 + (u[2])**2 )

    sinphi = np.dot(vect23,v) / (dist23*tnorm*unorm)
    cosphi = np.dot(t,u) / (tnorm*unorm)
    
    phi = np.arctan2(sinphi,cosphi)

    return phi



def dihedral_read(angle_list):
    """Obtain dihedral angles

    Args:
        angle_list (list of objects): previously calculated list of angles in angle_read

    Returns:
        dihedral_list (list of objects): list with information of the atoms involved in the dihedral and its value in radians
    """
    dihedral_list = [] # Format: [at1,at2,at3,at4,angle] 

    for i in range(len(angle_list)):
        vari = angle_list[i]
        for j in range(i+1,len(angle_list)):
            varj = angle_list[j]

            if vari[0].id == varj[1].id and vari[1].id == varj[0].id: # Central atoms in the two first columns
                dihedral_list.append([atoms[vari[2].id-1],atoms[vari[1].id-1],atoms[vari[0].id-1],atoms[varj[2].id-1],
                dihedral(atoms[vari[2].id-1],atoms[vari[1].id-1],atoms[vari[0].id-1],atoms[varj[2].id-1])])

            if vari[1].id == varj[0].id and vari[2].id == varj[1].id: # Central atoms in the two last columns in the first iterator and in the two first columns in the second
                dihedral_list.append([atoms[vari[0].id-1],atoms[vari[1].id-1],atoms[vari[2].id-1],atoms[varj[2].id-1],
                dihedral(atoms[vari[0].id-1],atoms[vari[1].id-1],atoms[vari[2].id-1],atoms[varj[2].id-1])])
            
            if vari[1].id == varj[2].id and vari[2].id == varj[1].id: # Central atoms in the two last columns
                dihedral_list.append([atoms[vari[0].id-1],atoms[vari[1].id-1],atoms[vari[2].id-1],atoms[varj[0].id-1],
                dihedral(atoms[vari[0].id-1],atoms[vari[1].id-1],atoms[vari[2].id-1],atoms[varj[0].id-1])])
            

    # for i in dihedral_list:
    #     print(f'{i[0].element,i[0].id} - {i[1].element,i[1].id} - {i[2].element,i[2].id} - {i[3].element,i[3].id} \t Angle in radians = {i[4]}, angle in degrees = {i[4]*180/np.pi}')
    
    return dihedral_list



def vdw_read(angle_list,bnd):
    """Obtain pairs of atoms that differ in more than two atoms

    Args:
        angle_list (list of objects): previously obtained list in angle_read
        bnd (array(bonds,2)): bond list 

    Returns:
        vdw_list (list of objects): list with information about the paired atoms and the Lennard-Jones interaction energy
    """
    
    connect_matrix = np.arange(0.5,nat*nat).reshape(nat,nat)
    for i in range(len(connect_matrix)):
        for j in angle_list:
            for k in bnd:
                connect_matrix[i][i] = 0  # Atom with itself
                connect_matrix[k[0]-1][k[1]-1] = connect_matrix[k[1]-1][k[0]-1] = 0 # Bonded atoms
                connect_matrix[j[0].id - 1][j[2].id - 1] = connect_matrix[j[2].id - 1][j[0].id - 1] = 0 # Angle between atoms

                
    vdw_list = []  # Format [at1,at2,VdWenergy]

    for i in range(len(connect_matrix)):
        for j in range(len(connect_matrix)):
            if connect_matrix[i][j] != 0:
                if atoms[i].element == 'C' and atoms[j].element == 'C':
                    energy = float(Acc/(bond_length(atoms[i],atoms[j])[0])**12 - Bcc/(bond_length(atoms[i],atoms[j])[0])**6)
                    connect_matrix[i][j] = energy
                    vdw_list.append([atoms[i],atoms[j],energy])

                if atoms[i].element == 'C' and atoms[j].element == 'H':
                    energy = float(Ach/(bond_length(atoms[i],atoms[j])[0])**12 - Bch/(bond_length(atoms[i],atoms[j])[0])**6)
                    connect_matrix[i][j] = energy
                    vdw_list.append([atoms[i],atoms[j],energy])

                if atoms[i].element == 'H' and atoms[j].element == 'C':
                    energy = float(Ach/(bond_length(atoms[i],atoms[j])[0])**12 - Bch/(bond_length(atoms[i],atoms[j])[0])**6)
                    connect_matrix[i][j] = energy
                    vdw_list.append([atoms[i],atoms[j],energy])

                if atoms[i].element == 'H' and atoms[j].element == 'H':
                    energy = float(Ahh/(bond_length(atoms[i],atoms[j])[0])**12 - Bhh/(bond_length(atoms[i],atoms[j])[0])**6)
                    connect_matrix[i][j] = energy
                    vdw_list.append([atoms[i],atoms[j],energy])

    # for i in vdw_list:
    #     print(f'{i[0].element,i[0].id} - {i[1].element,i[1].id} \t Van Der Waals energy (kcal/mol) = {i[2]}')

    return vdw_list


def internal_read(bnd):
    bond_list = bond_read(bnd)[1]
    angle_list = angle_read(bnd)
    dihedral_list = dihedral_read(angle_list) 

    internal_list = []
    for i in bond_list:
        internal_list.append(i[2])

    for i in angle_list:
        internal_list.append(i[3])

    for i in dihedral_list:
        internal_list.append(i[4])

    return bond_list, angle_list, dihedral_list, internal_list 


################################
# Potential energy calculation #
################################


def calc_vstretch(bond_list):
    """Calculate bonding contribution to potential energy

    Args:
        bond_list (list of objects): previously obtained list of bonds in bond_read

    Returns:
        vstretch(list(bonds)): list with bonding energy contribution of every bond
    """
    vstretch = []
    for i in bond_list:
        if i[0].element and i[1].element == 'C':
            vbond = kcc*(i[2] - r0cc)**2
            vstretch.append(vbond)
        else:
            vbond = kch*(i[2] - r0ch)**2
            vstretch.append(vbond)

    return vstretch



def calc_vbend(angle_list):
    """Calculate bending contribution to potential energy

    Args:
        angle_list (list of objects): previously obtained list of angles in angle_read

    Returns:
        vbend(list(len(angle_list))): list with bending energy contribution of every angle
    """
    vbend = []
    for i in angle_list:
        if i[0].element == i[1].element == i[2].element == 'C':
            vangle = kccc*(i[3]-a0ccc)**2
            vbend.append(vangle)
        elif i[0].element == i[1].element == 'C' and i[2].element == 'H':
            vangle = kcch*(i[3]-a0cch)**2
            vbend.append(vangle)
        elif i[2].element == i[1].element == 'C' and i[0].element == 'H':
            vangle = kcch*(i[3]-a0cch)**2
            vbend.append(vangle)
        elif i[0].element == i[2].element == 'H' and i[1].element == 'C':
            vangle = khch*(i[3]-a0hch)**2
            vbend.append(vangle)  

    return vbend


def calc_vtors(dihedral_list):
    """Calculate torsional contribution to potential energy

    Args:
        dihedral_list (list of objects): previously obtained list of dihedral angles in dihedral_read

    Returns:
        vtors(list(len(dihedral_list))): list with torsional energy contribution of every dihedral angle
    """
    vtors = []
    for i in dihedral_list:
        vdih = kxccx*(1+np.cos(n*i[4]))
        vtors.append(vdih)

    return vtors





def calc_vvdw(vdw_list):
    """Calculate dispersive interactions contribution to potential energy

    Args:
        vdw_list (list of objects): previously obtained list of paired atoms in vdw_read

    Returns:
        vvdw(list(len(vdw_list))): list with dispersive interactions energy contribution of every paired atom
    """    
    vvdw = []

    for i in vdw_list:
        vvdw.append(i[2])

    return vvdw



def getV(bnd):
    """Calculate total potential energy 

    Args:
        bnd (array(nat,2)): list of bonds

    Returns:
        Float: Potential energy
    """
    bond_list = bond_read(bnd)[1]
    angle_list = angle_read(bnd)
    dihedral_list = dihedral_read(angle_list)
    vdw_list = vdw_read(angle_list,bnd)

    vstretch = calc_vstretch(bond_list)
    vbend = calc_vbend(angle_list)
    vtors = calc_vtors(dihedral_list)
    vvdw = calc_vvdw(vdw_list)


    V = np.sum(vstretch) + np.sum(vbend) + np.sum(vtors) + np.sum(vvdw)/2

    return V



#########################
# Derivative definition #
#########################


def dr(at1,at2):
    """Derivative of a stretch coordinate

    Args:
        at1 (object atom): atom1 
        at2 (object atom): atom2

    Returns:
        array(1,3): derivative of the bond distance wrt cartesian coord 
    """
    r21 = bond_length(at2,at1)[1]
    r_norm = bond_length(at2,at1)[0]
    dd1 = r21/r_norm
    dd2 = -dd1
    return dd1,dd2


def dtheta(at1,at2,at3):
    """Derivative of a bend coordinate

    Args:
        at1 (object atom): atom1
        at2 (object atom): atom2 (central)
        at3 (object atom): atom3
    Returns:
        array(1,3): derivative wrt to cartesian coordinates of each atom 
    """
    dist21 = bond_length(at2,at1)[0]
    dist23= bond_length(at2,at3)[0]

    vect21 = bond_length(at2,at1)[1]
    vect23 = bond_length(at2,at3)[1]
    
    p = cross(vect21,vect23)
    p_norm = np.sqrt(p[0]**2 + p[1]**2 + p[2]**2)

    dd2 = -cross(vect21,p)/((dist21**2)*p_norm) + cross(vect23,p)/((dist23**2)*p_norm) # Central atom
    dd1 = cross(vect21,p)/((dist21**2)*p_norm) 
    dd3 = -cross(vect23,p)/((dist23**2)*p_norm)

    return dd1,dd2,dd3


def dphi(at1,at2,at3,at4):
    """Derivative of a torsion angle

    Args:
        at1 (oject atom): at1
        at2 (oject atom): at2
        at3 (oject atom): at3
        at4 (oject atom): at4

    Returns:
        array(1,3): derivative wrt to cartesian coordinates of each atom 
    """
    dist23= bond_length(at2,at3)[0]
    

    vect12 = bond_length(at1,at2)[1]
    vect23 = bond_length(at2,at3)[1]
    vect34 = bond_length(at3,at4)[1]
    vect13 = bond_length(at1,at3)[1]
    vect24 = bond_length(at2,at4)[1]

    t = cross(vect12,vect23)
    u = cross(vect23,vect34)

    tnorm = np.sqrt((t[0])**2 + (t[1])**2 + (t[2])**2 )
    unorm = np.sqrt((u[0])**2 + (u[1])**2 + (u[2])**2 )

    dd1 = cross(cross(t,vect23)/((tnorm**2)*dist23), vect23)
    dd2 = cross(vect13, cross(t,vect23)/((tnorm**2)*dist23)) + cross(cross(-u,vect23)/((unorm**2)*dist23), vect34)
    dd3 = cross(cross(t,vect23)/((tnorm**2)*dist23),vect12) + cross(vect24,cross(-u,vect23)/((unorm**2)*dist23))
    dd4 = cross(cross(-u,vect23)/((unorm**2)*dist23),vect23)

    return dd1,dd2,dd3,dd4




#######################################
# Energy gradient wrt cartesian coord #
#######################################


def grad_stretch(bond_list):
    """Calculates the gradient matrix of bending coordinates wrt cartesian coordinates

    Args:
        bond_list (list of objects): previously obtained list of bonds in bond_read

    Returns:
        array(nat,3): gradient matrix
    """
    gstretch = np.zeros((nat,3),dtype=float)

    for i in bond_list:
        if i[0].element and i[1].element == 'C':
            gstretch[i[0].id-1] += 2*kcc*(i[2] - r0cc) * dr(i[0],i[1])[0] # Stretch gradient of atom A wrt cart. coord.
            gstretch[i[1].id-1] += 2*kcc*(i[2] - r0cc) * dr(i[0],i[1])[1] # Stretch gradient of atom B wrt cart. coord.
            
        else:
            gstretch[i[0].id-1] += 2*kch*(i[2] - r0ch) * dr(i[0],i[1])[0] # Stretch gradient of atom A wrt cart. coord.
            gstretch[i[1].id-1] += 2*kch*(i[2] - r0ch) * dr(i[0],i[1])[1] # Stretch gradient of atom B wrt cart. coord.
    
    return gstretch



def grad_bend(angle_list):
    """Calculates the gradient matrix of bending coordinates wrt cartesian coordinates

    Args:
        angle_list (list of objects): previously obtained list of angles in angle_read

    Returns:
        array(nat,3): gradient matrix
    """
    gbend = np.zeros((nat,3),dtype=float)

    for i in angle_list:
        if i[0].element == i[1].element == i[2].element == 'C':
            gbend[i[0].id-1] += 2*kccc*(i[3]-a0ccc)*dtheta(i[0],i[1],i[2])[0]  # Gradient of atom1
            gbend[i[1].id-1] += 2*kccc*(i[3]-a0ccc)*dtheta(i[0],i[1],i[2])[1]  # Gradient of central atom
            gbend[i[2].id-1] += 2*kccc*(i[3]-a0ccc)*dtheta(i[0],i[1],i[2])[2]  # Gradient of atom2

        elif i[0].element == i[1].element == 'C' and i[2].element == 'H':
            gbend[i[0].id-1] += 2*kcch*(i[3]-a0cch)*dtheta(i[0],i[1],i[2])[0]  # Gradient of atom1
            gbend[i[1].id-1] += 2*kcch*(i[3]-a0cch)*dtheta(i[0],i[1],i[2])[1]  # Gradient of central atom
            gbend[i[2].id-1] += 2*kcch*(i[3]-a0cch)*dtheta(i[0],i[1],i[2])[2]  # Gradient of atom2

        elif i[2].element == i[1].element == 'C' and i[0].element == 'H':
            gbend[i[0].id-1] += 2*kcch*(i[3]-a0cch)*dtheta(i[0],i[1],i[2])[0]  # Gradient of atom1
            gbend[i[1].id-1] += 2*kcch*(i[3]-a0cch)*dtheta(i[0],i[1],i[2])[1]  # Gradient of central atom
            gbend[i[2].id-1] += 2*kcch*(i[3]-a0cch)*dtheta(i[0],i[1],i[2])[2]  # Gradient of atom2

        elif i[0].element == i[2].element == 'H' and i[1].element == 'C':
            gbend[i[0].id-1] += 2*khch*(i[3]-a0hch)*dtheta(i[0],i[1],i[2])[0]  # Gradient of atom1
            gbend[i[1].id-1] += 2*khch*(i[3]-a0hch)*dtheta(i[0],i[1],i[2])[1]  # Gradient of central atom
            gbend[i[2].id-1] += 2*khch*(i[3]-a0hch)*dtheta(i[0],i[1],i[2])[2]  # Gradient of atom2
    
    return gbend



def grad_tors(dihedral_list):
    """Calculates the gradient matrix of torsional coordinates wrt cartesian coordinates

    Args:
        dihedral_list (list of objects): previously obtained list of dihedral angles in dihedral_read

    Returns:
        array(nat,3): gradient matrix
    """
    gtors = np.zeros((nat,3),dtype=float)

    for i in dihedral_list:
        gtors[i[0].id-1] += (-n)*kxccx*np.sin(n*i[4])*dphi(i[0],i[1],i[2],i[3])[0]  # Gradient of atom 1
        gtors[i[1].id-1] += (-n)*kxccx*np.sin(n*i[4])*dphi(i[0],i[1],i[2],i[3])[1]  # Gradient of atom 2 (central carbon)
        gtors[i[2].id-1] += (-n)*kxccx*np.sin(n*i[4])*dphi(i[0],i[1],i[2],i[3])[2]  # Gradient of atom 3 (central carbon)
        gtors[i[3].id-1] += (-n)*kxccx*np.sin(n*i[4])*dphi(i[0],i[1],i[2],i[3])[3]  # Gradient of atom 4

    return gtors



def grad_vdw(vdw_list):
    """Calculates the gradient matrix of Lennard-Jones potential wrt cartesian coordinates

    Args:
        vdw_list (list of objects): previously obtained list of paired atoms in vdw_read

    Returns:
        array(nat,3): gradient matrix
    """
    gvdw = np.zeros((nat,3),dtype=float)

    for i in vdw_list:
        if i[0].element == i[1].element == 'C':
            gvdw[i[0].id-1] += np.array([i[0].x-i[1].x, i[0].y-i[1].y, i[0].z-i[1].z]) * (-12*Acc/(bond_length(i[0],i[1])[0]**14) + 6*Bcc/(bond_length(i[0],i[1])[0]**8))

        if i[0].element == 'H' and i[1].element == 'C':
            gvdw[i[0].id-1] += np.array([i[0].x-i[1].x, i[0].y-i[1].y, i[0].z-i[1].z]) * (-12*Ach/(bond_length(i[0],i[1])[0]**14) + 6*Bch/(bond_length(i[0],i[1])[0]**8))

        if i[1].element == 'H' and i[0].element == 'C':
            gvdw[i[0].id-1] += np.array([i[0].x-i[1].x, i[0].y-i[1].y, i[0].z-i[1].z]) * (-12*Ach/(bond_length(i[0],i[1])[0]**14) + 6*Bch/(bond_length(i[0],i[1])[0]**8))

        if i[0].element == i[1].element == 'H':
            gvdw[i[0].id-1] += np.array([i[0].x-i[1].x, i[0].y-i[1].y, i[0].z-i[1].z]) * (-12*Ahh/(bond_length(i[0],i[1])[0]**14) + 6*Bhh/(bond_length(i[0],i[1])[0]**8))

    return gvdw


def getgrad(bond_list,angle_list,dihedral_list,vdw_list):
    """Obtain gradient matrix

    Args:
        bond_list (list of objects): previously obtained list of bonds in bond_read
        angle_list (list of objects): previously obtained list of angles in angle_read
        dihedral_list (list of objects): previously obtained list of dihedral angles in dihedral_read
        vdw_list (list of objects): previously obtained list of paired atoms in vdw_read

    Returns:
        Array(nat,3): gradient matrix 
    """
    gstretch = grad_stretch(bond_list)
    gbend = grad_bend(angle_list)
    gtors = grad_tors(dihedral_list)
    gvdw = grad_vdw(vdw_list)

    gradV = gstretch + gbend + gtors + gvdw
    

    return gradV




#######################################
# Important functions for geom. opt. #
#######################################



def linear_search(coord,p_k,atoms,g_initial,V0):
    """Line search for an optimal value of alpha using Wolve's first condition

    Args:
        coord (array(nat,3)): set of coordinates to optimize
        p_k (array): predicted structure change
        atoms (list of objects): list of objects (atoms)
        g_initial (array): initial gradient matrix
        V0 (float): potential energy of the initial set of coordinates

    Returns:
        float: optimized alpha
    """
    alpha = 1
        
        
    while True:    # Linear search
        alpha = 0.8*alpha
        Vnew = step(coord,alpha,p_k,atoms,bnd)[0]
            

        if Vnew <= V0 + 0.1*alpha*np.dot(p_k.reshape((1,3*len(atoms))), g_initial.reshape((3*len(atoms),1))):
            
            return alpha
            


def step(coord,alpha,p_k,atoms,bnd):
    """_summary_

    Args:
        coord (array(nat,3)): set of coordinates to optimize
        alpha (float): alpha value for step scaling
        p_k (array): predicted structure change
        atoms (list of objects): list of objects (atoms)
        bnd (array(nat,2)): bond list

    Returns:
        Float: potential energy of new structure
        Array(nat,3): updated set of coordinates
    """
    rnew = coord + alpha*p_k.reshape((len(atoms),3))
    

    for i in range(len(atoms)):  # Update objects
        atoms[i].x = rnew[i,0]
        atoms[i].y = rnew[i,1]
        atoms[i].z = rnew[i,2]

    Vnew = getV(bnd) # Get potential energy in the new set

    return Vnew,rnew



def Bmat(bond_list,angle_list,dihedral_list ):
    B = np.zeros((internal_coordinates,3*nat),dtype=float)

    row = 0 # Counts which row to substitute data
    for i in bond_list:
        B[row,3*i[0].id-3:3*i[0].id] = dr(i[0],i[1])[0]
        B[row,3*i[1].id-3:3*i[1].id] = dr(i[0],i[1])[1]
        row += 1

    for i in angle_list:
        B[row,3*i[0].id-3:3*i[0].id] = dtheta(i[0],i[1],i[2])[0]
        B[row,3*i[1].id-3:3*i[1].id] = dtheta(i[0],i[1],i[2])[1]
        B[row,3*i[2].id-3:3*i[2].id] = dtheta(i[0],i[1],i[2])[2]
        row += 1

    for i in dihedral_list:
        B[row,3*i[0].id-3:3*i[0].id] = dphi(i[0],i[1],i[2],i[3])[0]
        B[row,3*i[1].id-3:3*i[1].id] = dphi(i[0],i[1],i[2],i[3])[1]
        B[row,3*i[2].id-3:3*i[2].id] = dphi(i[0],i[1],i[2],i[3])[2]
        B[row,3*i[3].id-3:3*i[3].id] = dphi(i[0],i[1],i[2],i[3])[3]
        row += 1



    return B



def inverse_G(B):
    G = np.dot(B,B.T)

        # Diagonalize G and get the inverse

    eigenvalues = sorted(np.linalg.eig(G)[0]) 
    V = np.linalg.eig(G)[1] #10x10

    l_mat = np.dot(V.T,np.dot(G,V)) # Lambda matrix with the diagonal elements being eigenvalues

        # Set small values to true 0 and invert the non-zero eigenvalues
    inv_l_mat = np.zeros((internal_coordinates,internal_coordinates),dtype=float)

    for i in range(internal_coordinates):
        for j in range(internal_coordinates):
            if abs(l_mat[i][j].real) <= 1.0E-6:               
                inv_l_mat[i][j] = 0.0
            else:
                inv_l_mat[i][j] = 1/l_mat[i][j].real

    
    inv_G = np.dot(V.real,np.dot(inv_l_mat, V.T.real))
        
    return inv_G




def update_Hessian(Mq,sqk,yqk,vqk):

    sqk = np.array(sqk).reshape((internal_coordinates,1))

    Mq1 = Mq + (np.dot(sqk.T,yqk) + np.dot(yqk.T,vqk)) / ((np.dot(sqk.T,yqk))**2)  * np.outer(sqk,sqk)   - (np.outer(vqk,sqk) + np.outer(sqk,vqk)) / np.dot(sqk.T,yqk)

    print(np.dot(sqk.T,yqk))

    # dot = 0
    # for i in range(len(sqk)):
    #     dot += sqk[i]*yqk[i]
    # print(f'Dot product sk*yk = {np.dot(sqk,yqk)}')
    # print(f'Try function: {dot}')

    return Mq1



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#



if __name__ == '__main__':
    
    filename = input('Choose a .mol2 file: \n')  # Load molecule


    # Next read  bonded atoms list, coordinates, element of each atom and number of atoms

    bnd = read_file(filename)[0]
    coord = read_file(filename)[1]
    element = read_file(filename)[2]
    nat = read_file(filename)[3]

    # Assign values for each object in the class and append in a list

    atoms = []
    for i in range(len(coord)):
        atoms.append(Atom(i+1,element[i],coord[i,0],coord[i,1],coord[i,2]))


    # Get bonds,angles, dihedral angles and dispersive interactions

    bond_list = bond_read(bnd)[1]
    angle_list = angle_read(bnd)
    dihedral_list = dihedral_read(angle_list)
    vdw_list = vdw_read(angle_list,bnd)


    internal_coordinates = len(bond_list)+len(angle_list)+len(dihedral_list)  #Number of total internal coordinates

    print(f'Internal coordinates')
    print('--------------------')  
    print(f' Stretching : {len(bond_list)} \n Bending : {len(angle_list)} \n Torsion : {len(dihedral_list)} \n Total : {internal_coordinates}' )
    print(' ')
    print(f'Cartesian coordinates')
    print('--------------------')
    print(f' Total: {3*nat}')



    
    
    ##################################################
    # Geometry optimization in cartesian coordinates #
    ##################################################


    # Start by calculating potential energy contributions and total potential energy

    vstretch = calc_vstretch(bond_list)
    vbend = calc_vbend(angle_list)
    vtors = calc_vtors(dihedral_list)
    vvdw = calc_vvdw(vdw_list)
    
    V0 = getV(bnd)

    print(' ')
    print(f'Potential energy of input geometry = {V0}')


    # Get initial gradient matrix

    g_initial = getgrad(bond_list,angle_list,dihedral_list,vdw_list)
 

    # print(' ')
    # print('Initial gradient matrix')
    # print('------------------------')
    # [print(atoms[i].element, g_initial[i]) for i in range(len(atoms))]


    # # Define the inverse Hessian as a diagonal matrix with values 1/300

    # M_k = np.zeros((3*nat,3*nat),dtype=float) # Approximate inverse Hessian
    # for i in range(len(M_k)):
    #     M_k[i][i] += 1/300


    # # Set initial direction of max displacement 

    # p_k = -np.dot(M_k,g_initial.reshape(3*nat,1)) 
    

    

    
        
    ################# Optimization algorithm #######################
 
    # Set previously calculated initial variables so that they enter the loop and can be updated

    # rold = coord
    # Vold = V0
    # gold = g_initial
    # Mold = M_k
    # count = 1 # Index

    # while True:
    #     print(' ')
    #     print('#############################################')
    #     print(f'Geometry optimization: cycle number {count}')
    #     print('#############################################')

    #     print(' ')
    #     print('Initial geometry')
    #     print('------------------------')
    #     print(rold)

    #     print(' ')
    #     print('Cartesian displacement')
    #     print('------------------------')
    #     print(p_k.reshape(nat,3))

    #     alpha = linear_search(rold,p_k,atoms,gold,Vold) #Calculate alpha at each cycle
    #     print(f'Linear search: alpha = {alpha}')

    #     Vnew = step(rold,alpha,p_k,atoms,bnd)[0]
    #     rnew = step(rold,alpha,p_k,atoms,bnd)[1]

    #     # After updating the coordinates for each object, define new lists

    #     bond = bond_read(bnd)[1]
    #     ang = angle_read(bnd)
    #     dih = dihedral_read(angle_list)
    #     vdw = vdw_read(angle_list,bnd)

    #     gnew = getgrad(bond,ang,dih,vdw)
        

    #     grms = np.sqrt(np.dot(gnew.reshape((1,3*nat)),gnew.reshape((3*nat,1))) / (3*nat)) # Gradient root mean square. 
        

    #     print(' ')
    #     print(f'New set of coordinates')
    #     print('------------------------')
    #     print(rnew)

        
    #     print(' ')
    #     print('Updated gradient')
    #     print('------------------------')
    #     print(gnew)

    #     print(' ')
    #     print(f'Energy difference = {Vnew-Vold} kcal/mol')       
    #     print(f'Gradient root mean square deviation (GRMS) = {grms}')





    #     if grms > 0.001:  
            
    #         s_k = alpha*p_k  # Scaled step (3*nat,1)
            
    #         y_k = gnew - gold  # Gradient difference (nat,3)
            
    #         v_k = np.dot(Mold,y_k.reshape((3*nat,1))) # (3*nat,1)
            
    #         # Calculate the new Hessian matrix

    #         M_new = Mold + (np.outer((np.dot(s_k.T,y_k.reshape((3*nat,1))) + np.dot(y_k.reshape((1,3*nat)),v_k) )* s_k, s_k)) / (np.dot(s_k.T,y_k.reshape((3*nat,1)))**2) - (np.outer(v_k,s_k.T) + np.outer(s_k,v_k.T)) / (np.dot(s_k.T,y_k.reshape((3*nat,1))))

    #         # Update variables so that the calculated become the old ones

    #         rold = rnew
    #         Mold = M_new 
    #         gold = gnew
    #         Vold = Vnew

    #         p_k = -np.dot(Mold,gold.reshape(3*nat,1))

    #         count += 1




    #     else: #Condition is fulfilled
    #         print(' ')
    #         print('Optimization converged :)')

    #         print(' ')
    #         print(' ')
    #         print(f'Final energy = {Vnew} kcal/mol')
    #         print(' ')
    #         print('Optimized geometry in Å: ')
    #         print('------------------------')
    #         print(rnew)
    #         break

        
        

    #################################################
    # Geometry optimization in internal coordinates #
    #################################################
   
   # First, build the B matrix

    B = Bmat(internal_read(bnd)[0],internal_read(bnd)[1],internal_read(bnd)[2])

    print(' ')
    print('Wilson B matrix: \n') 
    print(B)

    # Now, obtain the G matrix as G=B*B^T

    inv_G = inverse_G(B)

    print(' ')
    print('Inverse G matrix')
    print(inv_G)

    # With the Wilson B matrix and the inverse of G, we can obtain the gradient in internal coordinates as gq = G^-*B*gx

    gx = g_initial.reshape((3*nat,1))

    gq0 = np.dot(inv_G, np.dot(B,gx))

    print(' ')
    print('Gradient wrt internal coordinates')
    print(gq0.reshape((1,internal_coordinates)))


    # Internal coordinates system


    q0 = internal_read(bnd)[3]
    print(' ')
    print('Internal coordinates: \n')
    print(q0)



    # Initial guess for Inverse Hessian matrix

    M_q = np.zeros((internal_coordinates,internal_coordinates),dtype=float)

    for i in range(len(bond_list)):
        M_q[i][i] = 1/600

    for i in range(len(bond_list),len(bond_list)+len(angle_list)):
        M_q[i][i] = 1/150 

    for i in range(len(bond_list) + len(angle_list),len(bond_list)+len(angle_list)+len(dihedral_list)):
        M_q[i][i] = 1/80

    x0 = coord

    cycle = 1 # Iterator



    flag = True

    while flag: 
        print(' ')
        print('#############################################')
        print(f'Geometry optimization: cycle number {cycle}')
        print('#############################################')


        # Initial direction of max displacement

        p_q = -np.dot(M_q,gq0)

    ################# Optimization algorithm #######################

        rmax = 0.02

        rms = np.sqrt(np.dot(p_q.T,p_q)/internal_coordinates)
        print(' ')
        print('Predicted displacement: \n')
        print(p_q.reshape((1,internal_coordinates)))
        print(' ')


        if rms >= rmax:
            print('RMS is too high. Rescaling: \n')
            scaling = np.sqrt((0.02**2)*internal_coordinates / np.dot(p_q.T,p_q)) # Scaling factor so that rms = 0.02
            p_q *= scaling 
        
        print('Scaled step: \n')
        print(p_q.reshape((1,internal_coordinates)))

        # Update internal coordinates

        qk = q0 + p_q.reshape((1,internal_coordinates))
        
        print(' ')
        print('Updated set of internal coordinates: \n')
        print(qk)


        # Search for updated cartesian coordinates

        print(' ')
        print('Search for optimally updated cartesian coordinates')
        s_qk = qk - q0 #desired change
        print(s_qk)

        for i in range(len(bond_list) + len(angle_list), len(bond_list)+len(angle_list)+len(dihedral_list)):
            if np.abs(s_qk[0][i]) >= math.pi:
                if s_qk[0][i] < 0:
                    s_qk[0][i] += 2*math.pi
                elif s_qk[0][i] > 0:
                    s_qk[0][i] -= 2*math.pi # Ensure that the change is not more than 360º when changing from negative to positive

        
        xj1 = x0.reshape((nat*3,1)) + np.dot(B.T, np.dot(inv_G, s_qk.T)) 

    

        V0 = getV(bnd)

        diff = 0.00001
        displacement = 1
        iter = 0

        while displacement > diff:

            iter+=1
            print(' \n')
            print(f'Starting iteration number {iter}')
            print('----------------------------')
            
            

            # Obtain new set of internals

            for i in range(len(atoms)):  # Update objects
                atoms[i].x = xj1.reshape((nat,3))[i,0]
                atoms[i].y = xj1.reshape((nat,3))[i,1]
                atoms[i].z = xj1.reshape((nat,3))[i,2]

            qk1 = internal_read(bnd)[3]

            print(' ')
            print('New set of internals: \n')
            print(qk1)

            print(' ')
            print('Difference between internals: \n')
            s_qk1 = qk-qk1
            for i in range(len(bond_list) + len(angle_list),len(bond_list)+len(angle_list)+len(dihedral_list)):  # Check dihedral errors
                if np.abs(s_qk1[0][i]) >= math.pi:
                    if s_qk1[0][i] < 0:
                        s_qk1[0][i] += 2*math.pi
                    elif s_qk1[0][i] > 0:
                        s_qk1[0][i] -= 2*math.pi
            print(s_qk1)

            # Calculate new set of cartesians

            xj2 = xj1 + np.dot(B.T, np.dot(inv_G, s_qk1.T))

            print(' ')
            print('New set of cartesians: \n')       
            print(xj2.T)

            print(' ')
            print('Maximum cartesian displacement: \n')
            displacement = max(np.abs(xj2-xj1))
            print(displacement)

            xj1 = xj2

        print(' ')
        print('Optimal cartesians found! \n')

        print(' ')
        print('Updated set of cartesian coordinates: \n') 
        
        for i in range(nat):    
            print(atoms[i].element, xj2.reshape((nat,3))[i])

        for i in range(len(atoms)):  # Update objects
            atoms[i].x = xj2.reshape((nat,3))[i,0]
            atoms[i].y = xj2.reshape((nat,3))[i,1]
            atoms[i].z = xj2.reshape((nat,3))[i,2]
            

        qk1 = internal_read(bnd)[3]

        print(' ')
        print('Updated set of internal coordinates: \n') 
        print(qk1)


        # Update Wilson B matrix

        Bnew = Bmat(internal_read(bnd)[0],internal_read(bnd)[1],internal_read(bnd)[2])
        print(' ')
        print('Updated Wilson B matrix: \n') 
        print(Bnew)

        inv_Gnew = inverse_G(Bnew)
        print(' ')
        print('Updated inverse G matrix: \n')
        print(inv_Gnew)

        # Obtain gradient wrt internal coords

        gxk1 = getgrad(internal_read(bnd)[0],internal_read(bnd)[1],internal_read(bnd)[2],vdw_read(internal_read(bnd)[1],bnd))
        gqk1 = np.dot(inv_Gnew,np.dot(Bnew,gxk1.reshape((3*nat,1))))

        
        print(' ')
        print('Gradient of internal coordinates: \n')
        print(gqk1.T)

        yq = gqk1 - gq0
        vq = np.dot(M_q,yq)

        print(' ')
        print('Vector yq \n')
        print(yq)
        

        print(' ')
        print('Vector vq:')
        print(vq)

        s_qk_new = np.array(qk1) - np.array(q0)
        for i in range(len(bond_list) + len(angle_list),len(bond_list)+len(angle_list)+len(dihedral_list)):  # Check dihedral errors
                if np.abs(s_qk_new[i]) >= math.pi:
                    if s_qk_new[i] < 0:
                        s_qk_new[i] += 2*math.pi
                    elif s_qk_new[i] > 0:
                        s_qk_new[i] -= 2*math.pi
        print(' ')
        print('Vector sq:')
        print(s_qk_new)


        new_Mq = update_Hessian(M_q,s_qk_new,yq,vq)

        print(' ')
        print('Updated inverse Hessian Mq: \n')
        print(new_Mq)


        V = getV(bnd)

        print(' ')
        print(f'Old energy = {V0} \t New energy = {V}')
        
        rmse = np.sqrt(np.dot(gxk1.reshape((1,3*nat)),gxk1.reshape((3*nat,1))) / (3*nat))
        print(f'RMSE = {rmse}')

        M_q = new_Mq
        gq0 = gqk1

        B = Bnew
        q0 = qk1
        inv_G = inv_Gnew
        x0 = xj2

        if rmse >= 0.00001:
            cycle += 1

        else:
            flag = False

    print(' ')
    print('Optimization converged :)')
    print('------------------------')
    print(' ')

    print(f'Final energy = {V} kcal/mol \n')
    print('Optimized geometry: ')
    for i in range(nat):    
        print(atoms[i].element, xj2.reshape((nat,3))[i])



