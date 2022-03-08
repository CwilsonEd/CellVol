##############################################################################################
#                                    Cell Vol Script                                         #
#                            Created @ The University of Edinburgh                           #
#                         with functionality designed for the CCDC API                       #
#                                                                                            #
#                      _____ _____  __    __  _   _  ______  __                              #
#                     / ___// ___/ / /   / /  \\  \\ \  __ \ \ \                             #
#                    / /   / /__  / /   / /    \\  \\ \ \ \ \ \ \                            #
#	                / /   / ___/ / /   / /      \\  \\ \ \ \ \ \ \                           #
#	               / /_  / /__  / /_  / /_       \\ ˯\\ \ \_\ \ \ \__                        #
#                 /__ / /____/ /___/ /___/        \ __ \ \_____\ \___\                       #
##############################################################################################

"""This code features functionality to run both through Mercury and directly from the command line on single or multiple cif files.
   The necessary order in which the code must progress has restricted many global definitions in many instances and multiple local definitions are required.
   Command line inputs should be formulated >>python CellVol cifname pointsperrun no.runs precision
   Comments are included to allow the user to follow the code and for backend modification."""   

# Global imports of modules used for calculations #
   
from __future__ import division, print_function, absolute_import
import six
import ccdc.utilities
from ccdc.diagram import DiagramGenerator
import numpy as np
from numpy import random
import math
import sys
import os


# This import is required for scripts running from Mercury but will throw an error if attempts are made to import it from the command line #
# The try function is used to identify where the script is being run from for the calculation #  

try:
    from mercury_interface import MercuryInterface
    running_from_mercury = True
except ModuleNotFoundError:
    running_from_mercury = False

# List of van der Waals radii ordered by atomic number, Alvarez 2013, DOI: 10.1039/c3dt50599e #
# These are used to define the extent of the network volume around discrete atom coordinates #
	
vdw_radii = [1.2, 1.43, 2.12, 1.98, 1.91, 1.77, 1.66, 1.5, 1.46, 1.58, 2.5, 2.51, 2.25, 2.19, 1.9, 1.89, 1.82, 1.83, 2.73, 2.62, 2.58, 2.46,
                 2.42, 2.45, 2.45, 2.44, 2.40, 2.40, 2.38, 2.39, 2.32, 2.29, 1.88, 1.82, 1.86, 2.25, 3.21, 2.84, 2.75, 2.52, 2.56, 2.45, 2.44,
    			 2.46, 2.44, 2.15, 2.53, 2.49, 2.43, 2.42, 2.47, 1.99, 2.04, 2.06, 3.48, 3.03, 2.98, 2.88, 2.92, 2.95, -1, 2.9, 2.87, 2.83, 2.79, 2.87,
    			 2.81, 2.83, 2.79, 2.8, 2.74, 2.63, 2.53, 2.57, 2.49, 2.48, 2.41, 2.29, 2.32, 2.45, 2.47, 2.6, 2.54, -1, -1, -1, -1, -1, 2.8, 2.93, 2.88, 2.71, 2.82,
    			 2.81, 2.83, 3.05, 3.4, 3.05, 2.7]
				 

# This function calculates the magnitude of the interatomic vector between a random point and all atoms within the calculation #
# This is used to calculate whether it is within the van der Waals radius of any atom # 
  
def length_interatomic_matrix(atom_fract_coord, random_point, abc):
    delta = atom_fract_coord[:,:-1]-random_point
    magnitude = np.sum(np.square(abc[:3]*delta), axis=1) + \
                2*(abc[1]*abc[2]*np.cos(abc[3])*delta[:, 1]*delta[:, 2]+ \
                abc[0]*abc[2]*np.cos(abc[4])*delta[:, 0]*delta[:, 2]+ \
                abc[0]*abc[1]*np.cos(abc[5])*delta[:, 0]*delta[:, 1])
    magnitude = np.sqrt(magnitude)
    diff = atom_fract_coord[:, 3] - magnitude[:]
    num_vdw_radius_in = np.sum(np.array(diff) >= 0)
    return int(num_vdw_radius_in >= 1)

# Function to create a structure diagram in Mercury, taken from the CCDC #	
	
def make_diagram(mol, identifier, destdir, add_hydrogens=True):
    """Generate a diagram for the given molecule.

    :param mol: (:obj:`ccdc.molecule.Molecule`) The molecule to generate the diagram for.
    :param identifier: (:obj:`str`) The identifier for the molecule to use as a file name.
    :param destdir: (:obj:`str`) The target directory to save the diagram to.
    :param add_hydrogens: (:obj:`bool`) Whether to add hydrogens to the molecule for the diagram.

    :returns: (:obj:`str`) The absolute path to the saved diagram.
    """
    molecule_diagram_generator = DiagramGenerator()
    molecule_diagram_generator.settings.line_width = 1.6
    molecule_diagram_generator.settings.font_size = 12
    molecule_diagram_generator.settings.image_height = 300
    molecule = mol.copy()

    if add_hydrogens:
        molecule.add_hydrogens()

    img = molecule_diagram_generator.image(molecule)
    fname = os.path.join(destdir, '%s_diagram.png' % identifier)
    if img:
        img.save(fname)
    return os.path.abspath(fname)	
	
	
# Empty dictionaries and lists defined globally for use in the calculation # 
	
atom_dictionary = {}	
list_network_volume = []
list_void_volume = []

# The section of code contained within the 'True' return, runs the calculation from Mercury # 

if running_from_mercury == True:
   
   # Assign modules for the Mercury calculation #
    
    helper = MercuryInterface()
    entry = helper.current_entry
    crystal = entry.crystal
    molecule = entry.molecule
    label = entry.identifier
    html_file = helper.output_html_file
    
    # Import Tk classes for GUI dialogs. This are for the user to input the number of runs, points per run and precision required #
	
    from tkinter import *
	# Check if the user is running python 2 #
    if six.PY2:
        import Tkinter
        import tkSimpleDialog as SimpleDialog
    else:
        # import Python 3 variants to match Python 2 names for compatibility #
        import tkinter as Tkinter
        import tkinter.simpledialog
    
    # The following tells enures OS compatibility #
	# It tells the CSD Python API not to create a QApplication object on macOS. This script doesn't require any QApplication-dependent functionality #
	# On macOS, an Apple bug will cause Python crashes when using Tk if a QApplication object exists in the same Python script. #
	
    if sys.platform.startswith("darwin"):
        os.environ['CCDC_PYTHON_API_NO_QAPPLICATION'] = '1'
    
    # GUI dialog formatting #
    
    class MyDialog(tkinter.simpledialog.Dialog):
    
        def body(self, master):
    
            Label(master, text="Maximum Number of Runs (e.g. 16):").grid(row=0)
            Label(master, text="Points per Run (millions e.g. 1 = 1000000):").grid(row=1)
            Label(master, text="Precision % (e.g. 0.1 = 0.1% error):").grid(row=2)
    
            self.e1 = Entry(master)
            self.e2 = Entry(master)
            self.e3 = Entry(master)
    
            self.e1.grid(row=0, column=1)
            self.e2.grid(row=1, column=1)
            self.e3.grid(row=2, column=1)
            return self.e1 # initial focus
    
        def apply(self):
            first = int(self.e1.get())
            second = float(self.e2.get())
            third = float(self.e3.get())
            self.result = first, second, third
    
	# Withdraw inputs from the GUI #
	
    root = Tk()
    root.withdraw()
    d = MyDialog(root)
    
    # Assign the returned inputs and scale the number of points by a million #
	# The number of runs assigned is a maximum, a low enough standard deviation will allow the program to kick out any time after 3 runs #
    
    runs = d.result[0]
    n = int(1000000*d.result[1])
    precision_percent = d.result[2]
    precision = precision_percent*0.01
	
    # Check radiation source and normalise hydrogen containing bonds for X-ray studies #
    # This accounts for the artificial shortening of hydrogen containing bonds in X-ray studies # 
    
    if entry.radiation_source == 'Neutron':
            print('\n')
            print('Neutron radiation source detected')
    else:
        print('\n')
        print('X-ray radiation source detected, normalising hydrogens...')
        molecule.normalise_hydrogens()
    
	# Pack molecules to account for all contributing atoms to the unit cell by expanding to the largest VdW radii #
	# Convert to fractional coordinates #
    
    vdw_radii_in_unit = []
	
    for atom in molecule.atoms:
        Z = atom.atomic_number
        vdw_radii_in_unit.append(vdw_radii[(Z-1)])

    max_radii = max(vdw_radii_in_unit)
    a_length = crystal.cell_lengths[0]
    b_length = crystal.cell_lengths[1]
    c_length = crystal.cell_lengths[2]
    a_extra = (max_radii)/a_length
    b_extra = (max_radii)/b_length
    c_extra = (max_radii)/c_length
	
    a_extra = a_extra*100
    a_extra = math.ceil(a_extra)
    a_extra = a_extra/100
    a_minus = -a_extra
    a_plus = a_extra +1
	
    b_extra = b_extra*100
    b_extra = math.ceil(b_extra)
    b_extra = b_extra/100
    b_minus = -b_extra
    b_plus = b_extra +1
	
    c_extra = c_extra*100
    c_extra = math.ceil(c_extra)
    c_extra = c_extra/100
    c_minus = -c_extra
    c_plus = c_extra +1
	
    crystal.molecule = molecule
    packed_crystal = crystal.packing(box_dimensions=((a_minus, b_minus, c_minus), (a_plus, b_plus, c_plus)), inclusion='OnlyAtomsIncluded')
    
    # Assign an atom dictionary of the packed (and normalised) fractional coordinates # 
		
    for atom in packed_crystal.atoms:
        x,y,z = atom.fractional_coordinates
        Z = atom.atomic_number
        if vdw_radii[(Z-1)] != -1:
            atom_dictionary[atom] = [x,y,z,vdw_radii[(Z-1)]]
        else:
            f = open(html_file, "w")
            f.write('</i> <b>There is no defined van der Waals radius for the element {}</b><br>'.format(atom.atomic_symbol))
            f.write('</i> <b>Please provide one in the van der Waals radii list in the CellVol script at index {}'.format(Z-1))
            exit()
    
    # Mercury's void volume probe method and a number of other useful fields are calculated or extracted here for comparison #
	# The Mercury void method maps space available for a probe of defined radius, sampling at intervals defined by a grid spacing #
    # Defaults here map the smallest limits allowed by the calculation, these are set by the CCDC #	
    
    def packing_characteristics(molecule, crystal, probe_radius=0.2, grid_spacing=0.1):			
        
		# 18 Å^3 rule calcualtion #
		
        heavy = 0
        for atom in molecule.atoms:
            if atom.atomic_symbol != 'H' and atom.atomic_symbol != 'D':
                heavy += 1
        theoretical_volume = heavy * 18 * crystal.z_value if crystal.z_value else 0.0
        
		# Experimentally measured unit cell volume based on axes and angles #
		
        experimental_volume = round(crystal.cell_volume, 3)
		
		# Mercury calculated packing coefficient compared to typical values for the structure type #
		
        packing_coefficient = round(crystal.packing_coefficient, 3)
        if molecule.is_organic:
            moltype = "organic"
            av_pack_coeff = "0.68(4)"
        elif molecule.is_organometallic:
            moltype = "organometallic"
            av_pack_coeff = "0.67(5)"
        else:
            moltype = ""
            av_pack_coeff = ""
    
	    # Mercury void calculation #
	
        void_perc = crystal.void_volume(probe_radius, grid_spacing, mode='contact')
        void_vol = (void_perc * crystal.cell_volume) / 100
    
        return (theoretical_volume,
                experimental_volume,
                packing_coefficient,
                moltype,
                av_pack_coeff,
                void_perc,
                void_vol)
    
    # Function to run the network and void volume calculation and print out HTML tables # 
    			
    def run_report():
	
        # Generate a HTML report on the currently selected entry's crystal structure #
		
        interface = ccdc.utilities.ApplicationInterface(parse_commandline=False)
        interface.parse_commandline()
        entry_id = interface.identifier
    
        # Run packing charecteristics function, an alternative probe radius and grid spacing may be entered here when calling the function #
        # Adding extra variables to the function call will overite defaults. Minimum values are: probe radius = 0.2, grid spacing = 0.1 #
		
        theoretical_volume, experimental_volume, packing_coefficient, moltype,\
            av_pack_coeff, void_perc, void_vol =\
            packing_characteristics(molecule, crystal)
        
        # Open a txt file labelled by ref code to write your results to. This will write to a 'temp' folder on your C-drive. #
		# If this folder does not exist, it shall create it #
    	# With higher numbers of runs, a time-out error will prevent viewing HTML tables in Mercury but results will write to here regardless #
		
        filepath = r'c:/temp'
    	
        if not os.path.exists(filepath):
            os.makedirs(filepath)
    	
        f = open('c:/temp/%s.txt' % entry_id, 'w')
    
        # Assign parameters for calculation #
    	
        cell_volume = crystal.cell_volume
        lengths_abc = crystal.cell_lengths
        angles = crystal.cell_angles
        sum_percentage_points_in_envelope = 0
        sum_network_volume = 0
        sum_void_volume = 0
        mol_in_cell = crystal.z_value
    
    	# Convert angles to radians for numpy module #
    		
        deg_to_rad = 45. / math.atan(1.)
        abc = [lengths_abc[0], lengths_abc[1], lengths_abc[2], float(angles[0])/deg_to_rad, float(angles[1])/deg_to_rad, float(angles[2])/deg_to_rad]
    		
        """ Calculation loop, loops through the required number of random points and checks if they are within the vdw radius of any atom by calling the length_interatomic_matrix function. This is repeated for the required number of runs. This is then used to calculate network and void volume, averages are then taken over all runs and a standard deviation error calculated."""
		
		# Write numpy arrays for the calculation and run it#
		
        atom_fract_coord = atom_dictionary.values()
        list_atom_fract_coord = list(atom_fract_coord)
        atom_fract_matrix = np.array(list_atom_fract_coord)
        NUM_CALCULATIONS = int(n)
        abc_matrix = np.array(abc) 

        for i in range(0,runs):
            count_in_vdw_radius = 0
            random_matrix = np.random.rand(NUM_CALCULATIONS, 3)	
            count_in_vdw_radius = np.array([length_interatomic_matrix(atom_fract_matrix, random_point, abc_matrix) for random_point in random_matrix])
            count_in_vdw_radius = np.sum(count_in_vdw_radius)
            percentage_points_in_envelope = (count_in_vdw_radius/n)*100
            sum_percentage_points_in_envelope += percentage_points_in_envelope
            network_volume = (count_in_vdw_radius/n)*cell_volume
            list_network_volume.append(network_volume)
            sum_network_volume += network_volume
            void_volume = cell_volume - network_volume
            list_void_volume.append(void_volume)
            sum_void_volume += void_volume
            print('% Points in envelope:', '%.2f' % percentage_points_in_envelope, 'Network Volume:', '%.3f' % network_volume, 'Void Volume:', '%.3f' % void_volume, sep=' ',file=f)
            average_network_volume = sum_network_volume/(i+1)
            average_points_in_envelope = sum_percentage_points_in_envelope/(i+1)	
            average_void_volume = sum_void_volume/(i+1)
            if i != int(0):
                standard_deviation_network = 0
                standard_deviation_void = 0
                for vol in list_network_volume:
                    standard_deviation_network += (vol - average_network_volume)**2
                standard_deviation_network = math.sqrt((standard_deviation_network / (i+1 -1)))
                for vol in list_void_volume:  
                    standard_deviation_void += (vol - average_void_volume)**2
                standard_deviation_void = math.sqrt((standard_deviation_void / (i+1 -1)))
            else:
                standard_deviation_network = float(network_volume)
                standard_deviation_void = float(void_volume)
                number_of_runs_used = 1
            if standard_deviation_network >= average_network_volume*precision or i < 2:
                continue
            else:
                number_of_runs_used = i+1
                break
				
    
	    # Calculation to see if the standard deviation falls below 0.1% of the average value after the max number of runs, if not it will print an warning to concider increasing the number of points per run # 
	
        if standard_deviation_network >= average_network_volume*precision:
            print('\n', file=f)
            print('Standard deviation is >', precision_percent,'% of the volume. Consider increasing the number of points.', sep=' ',file=f)
            warning_to_screen = True
        else:
            warning_to_screen = False
		
		# This section writes averages and standard deviations to the output file, it also scales results to the number of molecules in a unit cell for analysis across phase transitions #   
		
        if standard_deviation_network == float(network_volume):
            standard_deviation_network = 'N/A'
            print('\n', file=f)
            print('CellVol packing coefficient:', '%.3f' % (average_points_in_envelope/100), file=f)
            print('Average Network Volume:', '%.3f' % average_network_volume,'Å^3', file=f)
            print('Average Void Volume:', '%.3f' % average_void_volume,'Å^3', file=f)
            print('Average Network Volume Scaled to Z:', '%.3f' % (average_network_volume/mol_in_cell),'Å^3', file=f)
            print('Average Void Volume Scaled to Z:', '%.3f' % (average_void_volume/mol_in_cell),'Å^3', file=f)	
        else:
            print('\n', file=f)
            print('CellVol packing coefficient:', '%.3f' % (average_points_in_envelope/100), file=f)
            print('Average Network Volume:', '%.3f' % average_network_volume, '+/-', '%.3f' % standard_deviation_network,'Å^3', file=f)
            print('Average Void Volume:', '%.3f' % average_void_volume, '+/-', '%.3f' % standard_deviation_void,'Å^3', file=f)
            print('Average Network Volume Scaled to Z:', '%.3f' % (average_network_volume/mol_in_cell),'+/-', '%.3f' % (standard_deviation_network/mol_in_cell),'Å^3', file=f)
            print('Average Void Volume Scaled to Z:', '%.3f' % (average_void_volume/mol_in_cell), '+/-', '%.3f' % (standard_deviation_void/mol_in_cell),'Å^3', file=f)		
        f.close()	
    		
        # Write output report to HTML tables in Mercury window #
    	
        tableno = 1
        with interface.html_report(title='CellVol Report for %s' % entry_id) as report:
    
            # Write structure overview section #
				        
            report.write_figure(file_name=make_diagram(molecule, crystal.identifier,
                                                   interface.output_directory_path),
                            alt_text='Figure 1. Diagram for %s' % entry_id,
                            caption='Figure 1. Diagram for %s' % entry_id)
            report.write_section_header('Crystal Structure Analysis')
            report.write(ccdc.utilities.html_table(
                data=[
                    ['Identifier', entry_id],
                    ['Formula', crystal.formula],
                    ['Space Group', crystal.spacegroup_symbol],
                    ['Cell Lengths (&#197;)', '<b>a</b> %.3f <b>b</b> %.3f <b>c</b> %.3f' % crystal.cell_lengths],
                    ['Cell Angles (&#176;)', '<b>&#945;</b> %.2f <b>&#946;</b> %.2f <b>&#947;</b> %.2f' % crystal.cell_angles],
                    ['Cell Volume (&#197;)', '%.2f' % crystal.cell_volume],
                    ['R-Factor', '%.2f' % entry.r_factor if entry.r_factor else 'unknown']
                 ],
                caption='Table %d. Selected Crystal Structure Information' % tableno
            ))
            tableno += 1
    
            # Write packing data table #
			
            report.write_section_header('Volume and Packing Analysis')
            report.write(ccdc.utilities.html_table(
                data=[
                    ['Estimated Volume from 18 &#197;&sup3; Rule', theoretical_volume if theoretical_volume else 'unknown'],
                    ['Experimental Volume', experimental_volume],
                    ['Classically Calculated Packing Coefficient', packing_coefficient],
                    ['CSD Average Packing Coefficient for %s Molecules' % moltype, av_pack_coeff],
                    ['Rolling Probe Void Percentage %', '%.2f' % void_perc],
                    ['Rolling Probe Void Volume /&#197;&sup3;', '%.3f' % void_vol],
                ],
                caption='Table %d. Crystal Packing Information' % tableno
            ))
            tableno += 1
    		
    		# Write Cell Vol Analysis table #
			
            if standard_deviation_network == 'N/A':
                report.write_section_header('Cell Vol Analysis')
                report.write(ccdc.utilities.html_table(
                    data=[
                        ['Z', mol_in_cell],
                        ['Number of Runs Used', number_of_runs_used],
                        ['Number of Points Per Run', n],
                        ['CellVol Packing Coefficient', '%.3f' % (average_points_in_envelope/100)],
                        ['Network Volume /&#197;&sup3;','%.3f' % average_network_volume],
                        ['Network Volume Scaled to Z /&#197;&sup3;','%.3f' % (average_network_volume/mol_in_cell)],
                        ['Void Volume /&#197;&sup3;', '%.3f' % average_void_volume],
                        ['Void Volume Scaled to Z /&#197;&sup3;', '%.3f' % (average_void_volume/mol_in_cell)],
                        ['Standard Deviation /&#197;&sup3', standard_deviation_network],
                    ],
                    caption='Table %d. Cell Vol Calculations' % tableno
                ))
                tableno += 1
            else:
                report.write_section_header('Cell Vol Analysis')
                report.write(ccdc.utilities.html_table(
                    data=[
                        ['Z', mol_in_cell],
                        ['Number of Runs Used', number_of_runs_used],
                        ['Number of Points Per Run', n],
                        ['CellVol Packing Coefficient', '%.3f' % (average_points_in_envelope/100)],
                        ['Network Volume /&#197;&sup3;','%.3f' % average_network_volume],
                        ['Network Volume Scaled to Z /&#197;&sup3;','%.3f' % (average_network_volume/mol_in_cell)],
                        ['Void Volume /&#197;&sup3;', '%.3f' % average_void_volume],
                        ['Void Volume Scaled to Z /&#197;&sup3;', '%.3f' % (average_void_volume/mol_in_cell)],
                        ['Standard Deviation /&#197;&sup3', '%.3f' % standard_deviation_network],
                    ],
                    caption='Table %d. Cell Vol Calculations' % tableno
                ))
                tableno += 1
            if warning_to_screen == True:
                report.write('Standard deviation is larger than the defined precision. Consider increasing the number of points.')
				
# The Mercury code ends here, the following section of code is contained within the 'False' return for the 'if' statement #
# It is only called if the program is run from the command line #

else:
    # Import and assign relevent modules for the command line calculation #
    
    from ccdc import io
    from ccdc.io import EntryReader
    csd_reader = EntryReader('CSD')
 
    # Assign the returned inputs and scale the number of points by a million #
	# The number of runs assigned is a maximum, a low enough standard deviation will allow the program to kick out any time after 3 runs #
    
    ref = sys.argv[1]	
    n = float(sys.argv[2])
    n = int(n*1000000)
    runs = int(sys.argv[3])
    try:
        precision_percent = float(sys.argv[4])
        precision = precision_percent*0.01
    except IndexError:
        precision_percent = 0.1
        precision = 0.001
	
    # Extract the directory path the command line is being run from, this is the location results will be written to #
    
    dir_path = os.path.dirname(os.path.realpath(__file__))
    
    """ Attempt first to search for an unpublished cif within the directory of the active command line with the defined name. This is not a case sensitive. 
     The program can run on a single cif file or a file containing multiple cifs. If running on a file with multiple cifs, the calculation is entirely
     completed within the 'for loop' below within the 'try' function. However a single unpublished cif or a CSD cif will be calculated by the code outside of 
     the 'try' function. If unsuccessful in directory searches, the program then attempts to read-in crystal and molecule paramaters via a ref code search of
     the CSD. This is case sensitive. Note, an unpublished cif in the directory will prevent a CSD ref code search of the same name. If still unsuccesful we
	 then check for any capitilisation errors in the typed ref code and search again. """
    
    cif_counter = 0
    counter = 0
	
    try: 
        filepath = open('{}\{}.cif'.format(dir_path, str(ref)), 'r')
    except OSError:
        print('\n')
        print('Filepath {}\{}.cif'.format(dir_path, str(ref)), 'does not exist')
        print('Searching CSD for', ref, '...')
        try:
            entry = csd_reader.entry(ref)
        except RuntimeError:
            print('\n')
            print(ref, 'not found within CSD')
            print('Checking for capitilisation error...')
            try:
                ref_up = ref.upper()
                entry = csd_reader.entry(ref_up)
            except RuntimeError:
                print('\n')
                print('Cif could not be located in directory of command line or in the CSD')
                print('Check you have named the CIF correctly and placed it in the correct folder')
                sys.exit()
                # Unsuccessful searches will terminate the script here and print a handled error #
                # If you see this error, enusre the unpublished cif is in the same directory as the command line or check naming in the call. #
            else:
                print('\n')
                print('CIF taken from the CSD')
                crystal = entry.crystal
                molecule = entry.molecule
                label = entry.identifier
                ref = ref_up
                published_cif = True
        else:
            print('\n')
            print('CIF taken from the CSD')
            crystal = entry.crystal
            molecule = entry.molecule
            label = entry.identifier
            published_cif = True
    else:
        published_cif = False
        print('\n')
        print('CIF taken from the command line directory')
        filepath = r"{}\{}.cif".format(dir_path, str(ref))
        f = open('{}\{}.cif'.format(dir_path, str(ref)), 'r')
        data = f.read().splitlines()
        f.close()
        for line in data:
            if line.strip()[0:5] == 'data_':
                cif_counter += 1
            else:
                continue
        if cif_counter == 1:
            crystal_reader = io.EntryReader(filepath)
            unpublished_data = crystal_reader[0]
            crystal = unpublished_data.crystal
            molecule = unpublished_data.molecule
        else:
            for line in data:
                if line.strip()[0:5] == 'data_':
                    counter +=1
                if counter != 0:
                    f = open('{}\{}.cif'.format(dir_path, str(counter)), 'a')
                    print(line, file=f)
                    f.close()
            for cif in range(1, counter+1):
                filepath = r"{}\{}.cif".format(dir_path, str(cif))
                crystal_reader = io.EntryReader(filepath)
                unpublished_data = crystal_reader[0]
                crystal = unpublished_data.crystal
                molecule = unpublished_data.molecule
                ref = cif
                
				# This section of the code runs on multi-CIF files #
				# Despite being defined globally, for multi-CIF files these paramaters must be reset for each calculation #
				
                atom_dictionary = {}	
                list_network_volume = []
                list_void_volume = []
				
                # Check radiation source and normalise hydrogen containing bonds for X-ray studies #
                # This accounts for the artificial shortening of hydrogen containing bonds in X-ray studies #
				
                if unpublished_data.radiation_source == 'Neutron':
                    print('\n')
                    print('Neutron radiation source detected')
                else:
                    print('\n')
                    print('X-ray radiation source detected, normalising hydrogens...')
                    molecule.normalise_hydrogens()
				
              	# Pack molecules to account for all contributing atoms to the unit cell by expanding to the largest VdW radii #
              	# Convert to fractional coordinates #
				
                vdw_radii_in_unit = []
	            
                for atom in molecule.atoms:
                    Z = atom.atomic_number
                    vdw_radii_in_unit.append(vdw_radii[(Z-1)])
                
                max_radii = max(vdw_radii_in_unit)
                a_length = crystal.cell_lengths[0]
                b_length = crystal.cell_lengths[1]
                c_length = crystal.cell_lengths[2]
                a_extra = (max_radii)/a_length
                b_extra = (max_radii)/b_length
                c_extra = (max_radii)/c_length
	            
                a_extra = a_extra*100
                a_extra = math.ceil(a_extra)
                a_extra = a_extra/100
                a_minus = -a_extra
                a_plus = a_extra +1
	            
                b_extra = b_extra*100
                b_extra = math.ceil(b_extra)
                b_extra = b_extra/100
                b_minus = -b_extra
                b_plus = b_extra +1
	            
                c_extra = c_extra*100
                c_extra = math.ceil(c_extra)
                c_extra = c_extra/100
                c_minus = -c_extra
                c_plus = c_extra +1
	            
                crystal.molecule = molecule
                packed_crystal = crystal.packing(box_dimensions=((a_minus, b_minus, c_minus), (a_plus, b_plus, c_plus)), inclusion='OnlyAtomsIncluded') 
				
				# Assign an atom dictionary of the packed (and normalised) fractional coordinates #
					
                for atom in packed_crystal.atoms:
                    x,y,z = atom.fractional_coordinates
                    Z = atom.atomic_number
                    if vdw_radii[(Z-1)] != -1:
                        atom_dictionary[atom] = [x,y,z,vdw_radii[(Z-1)]]
                    else:
                        print('There is no defined van der Waals radius for the element {}'.format(atom.atomic_symbol))
                        print('please provide one in the van der Waals radii list in the CellVol script at index {}'.format(Z-1))
                        exit()
				
			    # Mercury's void volume probe method function for comparison #
	            # The Mercury void method maps space available for a probe of defined radius, sampling at intervals defined by a grid spacing #	
	            # 18 Å^3 heavy atom calculation has been omittted from the CMD line calculation #
					
                def packing_characteristics(molecule, crystal, probe_radius, grid_spacing):
                    void_perc = crystal.void_volume(probe_radius, grid_spacing, mode='contact')
                    void_vol = (void_perc * crystal.cell_volume) / 100
                    return (void_vol)	
					
				# Parent function to run network and void volume calculation #
					
                def run_report():
                
                    # Set probe size and grid spacing for the rolling probe void volume calculation and run it # 
					# Defaults here map the smallest limits allowed by the calculation, these are set by the CCDC #
                
                    probe_radius = 0.2
                    grid_spacing = 0.1
                	
                    rolling_probe_vol = packing_characteristics(molecule, crystal, probe_radius, grid_spacing)
                    print('Rolling probe void volume =', '%.3f' % rolling_probe_vol, 'Å^3') 
                
                	# Open a file to write results to, we open in apend mode to allow extended studies to print to the same file #
                	
                    f = open('{}/cell_vol_results.txt'.format(dir_path), 'a')
                	
                	# Print calculation set-up variables to command line #
                	
                    print('\n')
                    print(ref)
                    print('Number of runs (maximum):', runs)
                    print('Number of points per run:', n)
                    print('Precision:', precision_percent, '%')
                	
                    # Assign parameters for calculation #
                	
                    cell_volume = crystal.cell_volume
                    lengths_abc = crystal.cell_lengths
                    angles = crystal.cell_angles
                    sum_percentage_points_in_envelope = 0
                    sum_network_volume = 0
                    sum_void_volume = 0
                
	            	# Read in the number of molcules in the unit cell to allow scaling to Z. If undefined, set to 1. #
	            	
                    if crystal.z_value != 0:
                        mol_in_cell = crystal.z_value
                    else:
                        mol_in_cell = 1                				
				
                 	# Convert angles to radians for numpy module #
                	
                    deg_to_rad = 45. / math.atan(1.)
                    abc = [lengths_abc[0], lengths_abc[1], lengths_abc[2], float(angles[0])/deg_to_rad, float(angles[1])/deg_to_rad, float(angles[2])/deg_to_rad]
                		
                    """ Calculation loop, loops through the required number of random points and checks if they are within the vdw radius of any atom by calling the length_interatomic_matrix function. This is repeated for the required number of runs. This is then used to calculate network and void volume, averages are then taken over all runs and a standard deviation error calculated."""
		
		            # Write numpy arrays for the calculation and run it#
					
                    atom_fract_coord = atom_dictionary.values()
                    list_atom_fract_coord = list(atom_fract_coord)
                    atom_fract_matrix = np.array(list_atom_fract_coord)
                    NUM_CALCULATIONS = int(n)
                    abc_matrix = np.array(abc)                

                    for i in range(0,runs):
                        count_in_vdw_radius = 0
                        random_matrix = np.random.rand(NUM_CALCULATIONS, 3)	
                        count_in_vdw_radius = np.array([length_interatomic_matrix(atom_fract_matrix, random_point, abc_matrix) for random_point in random_matrix])
                        count_in_vdw_radius = np.sum(count_in_vdw_radius)
                        percentage_points_in_envelope = (count_in_vdw_radius/n)*100
                        sum_percentage_points_in_envelope += percentage_points_in_envelope
                        network_volume = (count_in_vdw_radius/n)*cell_volume
                        list_network_volume.append(network_volume)
                        sum_network_volume += network_volume
                        void_volume = cell_volume - network_volume
                        list_void_volume.append(void_volume)
                        sum_void_volume += void_volume
                        average_network_volume = sum_network_volume/(i+1)
                        average_points_in_envelope = sum_percentage_points_in_envelope/(i+1)	
                        average_void_volume = sum_void_volume/(i+1)
                        if i < 9:
                            print('Run  ', i+1,'% Points in envelope:', '%.2f' % percentage_points_in_envelope, 'Network Volume:', '%.3f' % network_volume, 'Void Volume:', '%.3f' % void_volume)
                        else:
                            print('Run ', i+1,'% Points in envelope:', '%.2f' % percentage_points_in_envelope, 'Network Volume:', '%.3f' % network_volume, 'Void Volume:', '%.3f' % void_volume)
                        if i != int(0):
                            standard_deviation_network = 0
                            standard_deviation_void = 0
                            for vol in list_network_volume:
                                standard_deviation_network += (vol - average_network_volume)**2
                            standard_deviation_network = math.sqrt((standard_deviation_network / (i+1 -1)))
                            for vol in list_void_volume:  
                                standard_deviation_void += (vol - average_void_volume)**2
                            standard_deviation_void = math.sqrt((standard_deviation_void / (i+1 -1)))
                        else:
                            standard_deviation_network = float(network_volume)
                            standard_deviation_void = float(void_volume)
                        if standard_deviation_network >= average_network_volume*precision or i < 2:
                            continue
                        else:
                            break
                        
                
                    
                    print('\n')
                
	                # Calculation to see if the final standard deviation falls below 0.1% of the average value, if not it will print an warning to concider increasing the number of points per run #
	            
                    if standard_deviation_network >= average_network_volume*precision:
                        print('Standard deviation is >', precision_percent,'% of the volume. Consider increasing the number of points.', sep=' ',file=f)
                        print('Standard deviation is >', precision_percent,'% of the volume. Consider increasing the number of points.')
                        print('\n')
                	
                	# Print results to the command line and to the file (file=f) #
					
                    if standard_deviation_network == float(network_volume):
                        standard_deviation_network = 'N/A'
                        standard_deviation_void = 'N/A'
                        print(ref,'CellVol packing coefficient:', '%.3f' % (average_points_in_envelope/100),'Average Network Volume:', '%.3f' % average_network_volume,
                        '+/-', standard_deviation_network, 'Å^3', 'Average Void Volume:', '%.3f' % average_void_volume, '+/-', standard_deviation_void, 
                        'Å^3', 'Rolling probe void volume:', '%.3f' % rolling_probe_vol, 'Å^3', 'Network Volume Scaled to Z:', '%.3f' % (average_network_volume/mol_in_cell), 
                        '+/-', standard_deviation_network, 'Å^3', 'Void Volume Scaled to Z:', '%.3f' % (average_void_volume/mol_in_cell),
                        '+/-', standard_deviation_void, 'Å^3', 'Rolling probe void volume Scaled to Z:', '%.3f' % (rolling_probe_vol/mol_in_cell), 'Å^3', file=f)
                        
                        print('CellVol packing coefficient:', '%.3f' % (average_points_in_envelope/100))
                        print('Average Network Volume:', '%.3f' % average_network_volume, '+/-', standard_deviation_network, 'Å^3')
                        print('Average Void Volume:', '%.3f' % average_void_volume, '+/-', standard_deviation_void, 'Å^3')
                        print('\n')
                        if crystal.z_value == 0:
                            print('Z value not defined within CIF')
                        else:
                            print('Network Volume Scaled to Z:', '%.3f' % (average_network_volume/mol_in_cell), '+/-', standard_deviation_network, 'Å^3')
                            print('Void Volume Scaled to Z:', '%.3f' % (average_void_volume/mol_in_cell), '+/-', standard_deviation_void, 'Å^3')
                        print('\n')
                        print('Results were written to {}\cell_vol_results.txt'.format(dir_path))
                        print('\n')
                    else:
                        print(ref,'CellVol packing coefficient:', '%.3f' % (average_points_in_envelope/100),'Average Network Volume:', '%.3f' % average_network_volume,
                        '+/-', '%.3f' % standard_deviation_network, 'Å^3', 'Average Void Volume:', '%.3f' % average_void_volume, '+/-', '%.3f' % standard_deviation_void, 
                        'Å^3', 'Rolling probe void volume:', '%.3f' % rolling_probe_vol, 'Å^3', 'Network Volume Scaled to Z:', '%.3f' % (average_network_volume/mol_in_cell), 
                        '+/-', '%.3f' % (standard_deviation_network/mol_in_cell), 'Å^3', 'Void Volume Scaled to Z:', '%.3f' % (average_void_volume/mol_in_cell),
                        '+/-', '%.3f' % (standard_deviation_void/mol_in_cell), 'Å^3', 'Rolling probe void volume Scaled to Z:', '%.3f' % (rolling_probe_vol/mol_in_cell), 'Å^3', file=f)
                        
                        print('CellVol packing coefficient:', '%.3f' % (average_points_in_envelope/100))
                        print('Average Network Volume:', '%.3f' % average_network_volume, '+/-', '%.3f' % standard_deviation_network, 'Å^3')
                        print('Average Void Volume:', '%.3f' % average_void_volume, '+/-', '%.3f' % standard_deviation_void, 'Å^3')
                        print('\n')
                        if crystal.z_value == 0:
                            print('Z value not defined within CIF')
                        else:
                            print('Network Volume Scaled to Z:', '%.3f' % (average_network_volume/mol_in_cell), '+/-', '%.3f' % (standard_deviation_network/mol_in_cell), 'Å^3')
                            print('Void Volume Scaled to Z:', '%.3f' % (average_void_volume/mol_in_cell), '+/-', '%.3f' % (standard_deviation_void/mol_in_cell), 'Å^3')
                        print('\n')
                        print('Results were written to {}\cell_vol_results.txt'.format(dir_path))
                        print('\n')
                    f.close()
				
                if __name__ == '__main__':
                    run_report()

                
                if os.path.isfile(filepath):
                    os.remove(filepath)	
            sys.exit()				
	
	# Successful calculations on files containing multiple CIFs will terminate here. The following code is used to treat single CIF files and CSD searches only # 
    
    # Check radiation source and normalise hydrogen containing bonds for X-ray studies #
    # This accounts for the artificial shortening of hydrogen containing bonds in X-ray studies # 
    
    if published_cif == True:
        if entry.radiation_source == 'Neutron':
            print('\n')
            print('Neutron radiation source detected')
        else:
            print('\n')
            print('X-ray radiation source detected, normalising hydrogens...')
            molecule.normalise_hydrogens()
    else:
        if unpublished_data.radiation_source == 'Neutron':
            print('\n')
            print('Neutron radiation source detected')
        else:
            print('\n')
            print('X-ray radiation source detected, normalising hydrogens...')
            molecule.normalise_hydrogens()
    		
    
    # Pack molecules to account for all contributing atoms to the unit cell by expanding to the largest VdW radii #
    # Convert to fractional coordinates #
    
    vdw_radii_in_unit = []
	
    for atom in molecule.atoms:
        Z = atom.atomic_number
        vdw_radii_in_unit.append(vdw_radii[(Z-1)])

    max_radii = max(vdw_radii_in_unit)
    a_length = crystal.cell_lengths[0]
    b_length = crystal.cell_lengths[1]
    c_length = crystal.cell_lengths[2]
    a_extra = (max_radii)/a_length
    b_extra = (max_radii)/b_length
    c_extra = (max_radii)/c_length
	
    a_extra = a_extra*100
    a_extra = math.ceil(a_extra)
    a_extra = a_extra/100
    a_minus = -a_extra
    a_plus = a_extra +1
	
    b_extra = b_extra*100
    b_extra = math.ceil(b_extra)
    b_extra = b_extra/100
    b_minus = -b_extra
    b_plus = b_extra +1
	
    c_extra = c_extra*100
    c_extra = math.ceil(c_extra)
    c_extra = c_extra/100
    c_minus = -c_extra
    c_plus = c_extra +1
	

	
    crystal.molecule = molecule
    packed_crystal = crystal.packing(box_dimensions=((a_minus, b_minus, c_minus), (a_plus, b_plus, c_plus)), inclusion='OnlyAtomsIncluded') 
    
    # Assign an atom dictionary of the packed (and normalised) fractional coordinates # 	
    
	
    for atom in packed_crystal.atoms:
        x,y,z = atom.fractional_coordinates
        Z = atom.atomic_number
        if vdw_radii[(Z-1)] != -1:
            atom_dictionary[atom] = [x,y,z,vdw_radii[(Z-1)]]
        else:
            print('There is no defined van der Waals radius for the element {}'.format(atom.atomic_symbol))
            print('please provide one in the van der Waals radii list in the CellVol script at index {}'.format(Z-1))
            exit()
		
    
    # Mercury's void volume probe method function for comparison #
    # The Mercury void method maps space available for a probe of defined radius, sampling at intervals defined by a grid spacing #	
    # 18 Å^3 heavy atom calculation has been omittted from the CMD line calculation #
    	
    def packing_characteristics(molecule, crystal, probe_radius, grid_spacing):
        void_perc = crystal.void_volume(probe_radius, grid_spacing, mode='contact')
        void_vol = (void_perc * crystal.cell_volume) / 100
        return (void_vol)
    
    # Parent function to run network and void volume calculation # 
    			
    def run_report():
    
        # Set probe size and grid spacing for the rolling probe void volume calculation and run it # 
	    # Defaults here map the smallest limits allowed by the calculation, these are set by the CCDC #
    
        probe_radius = 0.2
        grid_spacing = 0.1
    	
        rolling_probe_vol = packing_characteristics(molecule, crystal, probe_radius, grid_spacing)
        print('Rolling probe void volume =', '%.3f' % rolling_probe_vol, 'Å^3') 
    
    	# Open a file to write results to, we open in apend mode to allow extended studies to print to the same file #
    	
        f = open('{}/cell_vol_results.txt'.format(dir_path), 'a')
    	
    	# Print calculation set-up variables to command line #
    	
        print('\n')
        print(ref)
        print('Number of runs (maximum):', runs)
        print('Number of points per run:', n)
        print('Precision:', precision_percent,'%')
		
        # Assign parameters for calculation #
    	
        cell_volume = crystal.cell_volume
        lengths_abc = crystal.cell_lengths
        angles = crystal.cell_angles
        sum_percentage_points_in_envelope = 0
        sum_network_volume = 0
        sum_void_volume = 0

		# Read in the number of molcules in the unit cell to allow scaling to Z. If undefined, set to 1. #
		
        if crystal.z_value != 0:
            mol_in_cell = crystal.z_value
        else:
            mol_in_cell = 1
    
        # Convert angles to radians for numpy module #
    	
        deg_to_rad = 45. / math.atan(1.)
        abc = (lengths_abc[0], lengths_abc[1], lengths_abc[2], float(angles[0])/deg_to_rad, float(angles[1])/deg_to_rad, float(angles[2])/deg_to_rad)
    		
        """ Calculation loop, loops through the required number of random points and checks if they are within the vdw radius of any atom by calling the length_interatomic_matrix function. This is repeated for the required number of runs. This is then used to calculate network and void volume, averages are then taken over all runs and a standard deviation error calculated."""
		
		# Write numpy arrays for the calculation and run it#
    	
        atom_fract_coord = atom_dictionary.values()
        list_atom_fract_coord = list(atom_fract_coord)
        atom_fract_matrix = np.array(list_atom_fract_coord)
        NUM_CALCULATIONS = int(n)
        abc_matrix = np.array(abc)

        for i in range(0,runs):
            count_in_vdw_radius = 0
            random_matrix = np.random.rand(NUM_CALCULATIONS, 3)	
            count_in_vdw_radius = np.array([length_interatomic_matrix(atom_fract_matrix, random_point, abc_matrix) for random_point in random_matrix])
            count_in_vdw_radius = np.sum(count_in_vdw_radius)
            percentage_points_in_envelope = (count_in_vdw_radius/n)*100
            sum_percentage_points_in_envelope += percentage_points_in_envelope
            network_volume = (count_in_vdw_radius/n)*cell_volume
            list_network_volume.append(network_volume)
            sum_network_volume += network_volume
            void_volume = cell_volume - network_volume
            list_void_volume.append(void_volume)
            sum_void_volume += void_volume
            average_network_volume = sum_network_volume/(i+1)
            average_points_in_envelope = sum_percentage_points_in_envelope/(i+1)	
            average_void_volume = sum_void_volume/(i+1)
            if i < 9:
                print('Run  ', i+1,'% Points in envelope:', '%.2f' % percentage_points_in_envelope, 'Network Volume:', '%.3f' % network_volume, 'Void Volume:', '%.3f' % void_volume)
            else:
                print('Run ', i+1,'% Points in envelope:', '%.2f' % percentage_points_in_envelope, 'Network Volume:', '%.3f' % network_volume, 'Void Volume:', '%.3f' % void_volume)
            if i != int(0):
                standard_deviation_network = 0
                standard_deviation_void = 0
                for vol in list_network_volume:
                    standard_deviation_network += (vol - average_network_volume)**2
                standard_deviation_network = math.sqrt((standard_deviation_network / (i+1 -1)))
                for vol in list_void_volume:  
                    standard_deviation_void += (vol - average_void_volume)**2
                standard_deviation_void = math.sqrt((standard_deviation_void / (i+1 -1)))
            else:
                standard_deviation_network = float(network_volume)
                standard_deviation_void = float(void_volume)
            if standard_deviation_network >= average_network_volume*precision or i < 2:
                continue
            else:
                break
            
    
        
        print('\n')
    
	    # Calculation to see if the final standard deviation falls below 0.1% of the average value, if not it will print an warning to concider increasing the number of points per run #
	
        if standard_deviation_network >= average_network_volume*precision:
            print('Standard deviation is >', precision_percent,'% of the volume. Consider increasing the number of points.', sep=' ',file=f)
            print('Standard deviation is >', precision_percent,'% of the volume. Consider increasing the number of points.')
            print('\n')
    	
    	# Print results to the command line and to the file (file=f) #
    	
        if standard_deviation_network == float(network_volume):
            standard_deviation_network = 'N/A'
            standard_deviation_void = 'N/A'
            if len(ref) > 6:
                print(ref,'CellVol packing coefficient:', '%.3f' % (average_points_in_envelope/100),'Average Network Volume:', '%.3f' % average_network_volume,
                '+/-', standard_deviation_network, 'Å^3', 'Average Void Volume:', '%.3f' % average_void_volume, '+/-', standard_deviation_void, 
                'Å^3', 'Rolling probe void volume:', '%.3f' % rolling_probe_vol, 'Å^3', 'Network Volume Scaled to Z:', '%.3f' % (average_network_volume/mol_in_cell), 
                '+/-', standard_deviation_network, 'Å^3', 'Void Volume Scaled to Z:', '%.3f' % (average_void_volume/mol_in_cell),
                '+/-', standard_deviation_void, 'Å^3', 'Rolling probe void volume Scaled to Z:', '%.3f' % (rolling_probe_vol/mol_in_cell), 'Å^3', file=f)
            else:
                print(ref,' ','CellVol packing coefficient:', '%.3f' % (average_points_in_envelope/100),'Average Network Volume:', '%.3f' % average_network_volume,
                '+/-', standard_deviation_network, 'Å^3', 'Average Void Volume:', '%.3f' % average_void_volume, '+/-', standard_deviation_void, 
                'Å^3', 'Rolling probe void volume:', '%.3f' % rolling_probe_vol, 'Å^3', 'Network Volume Scaled to Z:', '%.3f' % (average_network_volume/mol_in_cell), 
                '+/-', standard_deviation_network, 'Å^3', 'Void Volume Scaled to Z:', '%.3f' % (average_void_volume/mol_in_cell),
                '+/-', standard_deviation_void, 'Å^3', 'Rolling probe void volume Scaled to Z:', '%.3f' % (rolling_probe_vol/mol_in_cell), 'Å^3', file=f)                
            print('CellVol packing coefficient:', '%.3f' % (average_points_in_envelope/100))
            print('Average Network Volume:', '%.3f' % average_network_volume, '+/-', standard_deviation_network, 'Å^3')
            print('Average Void Volume:', '%.3f' % average_void_volume, '+/-', standard_deviation_void, 'Å^3')
            print('\n')
            if crystal.z_value == 0:
                print('Z value not defined within CIF')
            else:
                print('Network Volume Scaled to Z:', '%.3f' % (average_network_volume/mol_in_cell), '+/-', standard_deviation_network, 'Å^3')
                print('Void Volume Scaled to Z:', '%.3f' % (average_void_volume/mol_in_cell), '+/-', standard_deviation_void, 'Å^3')
            print('\n')
            print('Results were written to {}\cell_vol_results.txt'.format(dir_path))
            print('\n')
        else:
            if len(ref) > 6:
                print(ref,'CellVol packing coefficient:', '%.3f' % (average_points_in_envelope/100),'Average Network Volume:', '%.3f' % average_network_volume,
                '+/-', '%.3f' % standard_deviation_network, 'Å^3', 'Average Void Volume:', '%.3f' % average_void_volume, '+/-', '%.3f' % standard_deviation_void, 
                'Å^3', 'Rolling probe void volume:', '%.3f' % rolling_probe_vol, 'Å^3', 'Network Volume Scaled to Z:', '%.3f' % (average_network_volume/mol_in_cell), 
                '+/-', '%.3f' % (standard_deviation_network/mol_in_cell), 'Å^3', 'Void Volume Scaled to Z:', '%.3f' % (average_void_volume/mol_in_cell),
                '+/-', '%.3f' % (standard_deviation_void/mol_in_cell), 'Å^3', 'Rolling probe void volume Scaled to Z:', '%.3f' % (rolling_probe_vol/mol_in_cell), 'Å^3', file=f)
            else:
                print(ref,' ','CellVol packing coefficient:', '%.3f' % (average_points_in_envelope/100),'Average Network Volume:', '%.3f' % average_network_volume,
                '+/-', '%.3f' % standard_deviation_network, 'Å^3', 'Average Void Volume:', '%.3f' % average_void_volume, '+/-', '%.3f' % standard_deviation_void, 
                'Å^3', 'Rolling probe void volume:', '%.3f' % rolling_probe_vol, 'Å^3', 'Network Volume Scaled to Z:', '%.3f' % (average_network_volume/mol_in_cell), 
                '+/-', '%.3f' % (standard_deviation_network/mol_in_cell), 'Å^3', 'Void Volume Scaled to Z:', '%.3f' % (average_void_volume/mol_in_cell),
                '+/-', '%.3f' % (standard_deviation_void/mol_in_cell), 'Å^3', 'Rolling probe void volume Scaled to Z:', '%.3f' % (rolling_probe_vol/mol_in_cell), 'Å^3', file=f)			
            print('CellVol packing coefficient:', '%.3f' % (average_points_in_envelope/100))
            print('Average Network Volume:', '%.3f' % average_network_volume, '+/-', '%.3f' % standard_deviation_network, 'Å^3')
            print('Average Void Volume:', '%.3f' % average_void_volume, '+/-', '%.3f' % standard_deviation_void, 'Å^3')
            print('\n')
            if crystal.z_value == 0:
                print('Z value not defined within CIF')
            else:
                print('Network Volume Scaled to Z:', '%.3f' % (average_network_volume/mol_in_cell), '+/-', '%.3f' % (standard_deviation_network/mol_in_cell), 'Å^3')
                print('Void Volume Scaled to Z:', '%.3f' % (average_void_volume/mol_in_cell), '+/-', '%.3f' % (standard_deviation_void/mol_in_cell), 'Å^3')
            print('\n')
            print('Results were written to {}\cell_vol_results.txt'.format(dir_path))
        f.close()
    
# Call the cell vol report and run it # 		
    	    
if __name__ == '__main__':
    run_report()	
