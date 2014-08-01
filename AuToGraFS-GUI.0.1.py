####################################################################################################################################################################        WHAT IS THIS ?
###############################################################################################################################################################
##      This is the GUI for the novel AuToGraFs software (Automated Topological GeneRAtor for Framework Structures)
##      It enables the users t0 experience some capabilities of the base software while simplifying the use through graphical tools and
##      3D visualization.
##
####################################################################################################################################################################        IMPORTS
#############################################################################################################################################################
import Tkinter as tk
import tkMessageBox, tkFileDialog
import os
import ase
from AuToGraFS-GUI-dependencies import topodict, elements_dict
from ase import Atom, Atoms
from subprocess import check_call
from os.path import exists
from vtk import * 
from vtk.tk.vtkTkRenderWindowInteractor import vtkTkRenderWindowInteractor
from math import cos, sin


#############################################################################################################################################################
########        DECLARED STUFF
#############################################################################################################################################################

pi = 3.141592
cdic={}
ldic={}

#############################################################################################################################################################
########        FUNCTIONS
#############################################################################################################################################################

def makeSphere(radius, position, color):
    # create the object to be drawn
    sphere = vtkSphereSource()
    sphere.SetRadius( radius )
    sphere.SetThetaResolution(30)
    sphere.SetPhiResolution(30)
    sphere.SetCenter(position[0], position[1], position[2])
    # Convert the sphere into polygons
    sphereMapper = vtkPolyDataMapper()
    sphereMapper.SetInput(sphere.GetOutput())
    #Create an actor for the sphere 
    sphereActor = vtkActor()
    sphereActor.SetMapper(sphereMapper)
    (sphereActor.GetProperty()).SetColor(color[0], color[1], color[2])
    Radius = sphereActor.MakeProperty()
    return sphereActor

    
def functionalize(obj, event, *args):                   # when pointer used to select functionalizable atom, returns a functionalize dialog
                                                        # comment: the pointer is prone to error. multiple trials are often necessary
    C = var_center.get().split('.')
    L = var_linker.get().split('.')
    T= var_topo.get().split('.')
    supercell = var_supercell.get()
    if supercell != '1,1,1': 
        tkMessageBox.showerror('Error', 'Discrete functionalization of supercells is not yet implemented. Try again on unit cell')
        return
    picker = vtk.vtkPicker()
    picker.SetTolerance(0.001)
    renWin_disp.MakeCurrent()
    X,Y = disp_widget.GetEventPosition() 
    picker.Pick([X, Y,0], ren_disp)
    sphereActor = picker.GetActor()
    sphere_index_raw = sphereActor.GetProperty().GetEdgeColor()         #the atom number and name are stored into Edge color and actor color
    sphere_index = int(str(sphere_index_raw[0]).split('.')[1])
    print sphere_index
    with open('frags.out','r') as fragfile:                             #get corresponding info from fragsfile
        line_frag = fragfile.readlines()[sphere_index]
        atomtype_picked = line_frag.split()[2].strip()
        print atomtype_picked
        if atomtype_picked != 'H' and  atomtype_picked != 'X':
            tkMessageBox.showerror('Error', 'The picking process has encountered an error (discordant atom types). Please change you picking angle and try again. ')
            return
    color = sphereActor.GetProperty().GetAmbientColor()
    selected_atom = []
    for element, prop in elements_dict.iteritems():
        rgb = list(color)
        if prop[1] == rgb :
            selected_atom.append(element)
        else:
            pass
    if (selected_atom[0] == 'X' or selected_atom[0] == 'H') and len(choices_func) > 0:      #pops open the functionalize dialog if the atom picked
        top = tk.Toplevel(root)                                                             #is functionalizable
        top.title('{0}-{1}-{2} functionalization'.format(T[0],C[0],L[0]))
        frame_top = tk.Frame(top)
        frame_top.grid(row=0, column=1, sticky='n'+'w'+'e'+'s')
        myLabel = tk.Label(frame_top, text='To functionalize selected atom\n please choose a funcional group')
        myLabel.grid(row=0, column = 0, columnspan = 3)
        choose_func = tk.OptionMenu(frame_top, var_func, *choices_func)
        choose_func.grid(row=1, column=0)
        cancelbutton = tk.Button(frame_top, text='CANCEL/FINISH', command=lambda:top.destroy())
        cancelbutton.grid(row=2, column=0, columnspan=2)
        applybutton = tk.Button(frame_top, text='APPLY', command=lambda:newmake(sphere_index))
        applybutton.grid(row=1, column=1)
        tk.Grid.columnconfigure(top,1,weight=1)
        tk.Grid.rowconfigure(top,0,weight=1)



def newmake(sphere_index, *args):                                                         # makes a mof using AuToGraFs, taking functionalization into
    path = var_SBUpath.get()                                                              # account
    F = var_func.get().split('.')
    C = var_center.get().split('.')
    L = var_linker.get().split('.')
    T= var_topo.get().split('.')
    with open('control.txt','r') as checkfile:
        linecheck = checkfile.readlines()
        if 'model' in linecheck[0]:
            controlfile = 'control.txt'
        else:
            controlfile = 'control-mofgen.txt'
    with open('frags.out', 'r') as fragsfile:
        line_frag = fragsfile.readlines()[sphere_index]
        id_link = int(line_frag.split()[6])
        id_atom = int(line_frag.split()[1])
    with open(controlfile,'r') as tempfile:
        orig_lines = tempfile.readlines()
        if 'func' in orig_lines[id_link+1]:
            orig_lines[id_link+1] = orig_lines[id_link+1].strip()+':{0},{1}\n'.format(str(id_atom),F[0])
        else: 
            orig_lines[id_link+1] = orig_lines[id_link+1].strip()+' func={0},{1}\n'.format(str(id_atom),F[0])
    with open('control.txt','w') as newfile:
        for line in orig_lines:
            newfile.write(line)
    try:
       check_call(['~/source_codes/mofgen -p {3} -o {0}-{1}-{2}-func-{4} &> {0}-{1}-{2}.log'.format(T[0], C[0], L[0],path, F[0] )], shell= True)
    except CalledProcessError:
        tkMessageBox.showerror('Error', 'an error occured during functionalization. Check log file for details')
    top.destroy()


def askdirSBU(*args):                                               #self explaining
    var_SBUpath.set(tkFileDialog.askdirectory())


def updateSBU(*args):                                               #when changing directory or topology, the SBUs are updated
    global choices_func
    choices_func = [os.path.splitext(f)[0] for f in os.listdir('{0}/functional_groups/'.format(var_SBUpath.get())) if os.path.isfile(os.path.join('{0}/functional_groups/'.format(var_SBUpath.get()), f)) and os.path.splitext(f)[1]=='.inp']
    var_func.set(choices_func[0])

    centers = [f for f in os.listdir('{0}/centers/'.format(var_SBUpath.get())) if os.path.isfile(os.path.join('{0}/centers/'.format(var_SBUpath.get()), f)) and os.path.splitext(f)[1]=='.inp']
    for topo in list(topodict.keys()):
        geom = topodict[topo]
        c = [] 
        for center in centers:
            centername = os.path.splitext(center)[0]
            with open('{0}/centers/{1}'.format(var_SBUpath.get(),center), 'r') as centerfile:
                lines=centerfile.readlines()
                for line in lines:
                    if 'Data: shape =' in line:
                        v = line.split()
                        if v[3] == geom[0]:
                            c.append(center)
                        else:
                            pass
                        cdic.update({topo:c})
                    else:
                        pass

    linkers = [f for f in os.listdir('{0}/linkers/'.format(var_SBUpath.get())) if os.path.isfile(os.path.join('{0}/linkers/'.format(var_SBUpath.get()), f)) and os.path.splitext(f)[1]=='.inp']
    for topo in list(topodict.keys()):
        geom = topodict[topo]
        l = []
        for linker in linkers:
            linkername = os.path.splitext(linker)[0]
            with open('{0}/linkers/{1}'.format(var_SBUpath.get(),linker), 'r') as linkerfile:
                lines=linkerfile.readlines()
                for line in lines:
                    if 'Data: shape =' in line:
                        v = line.split()
                        if v[3] == geom[1]:
                            l.append(linker)
                        else:
                            pass
                        ldic.update({topo:l})
                    else:
                        pass

def make(*args):                                                        #will create a MOF using the selected parameters using AuToGraFs
    T= var_topo.get().split('.')
    C = var_center.get().split('.')
    L = var_linker.get().split('.')
    path = var_SBUpath.get()
    supercell_tmp = var_supercell.get()
    if supercell_tmp != '1,1,1':
        supercell = '\nsupercell = {0}'.format(supercell_tmp)
    else:
        supercell = ''
    with open('control.txt','w') as controlfile:
        controlfile.write('topology = {0}\ncenter = {1}\nlinker = {2}{3}'.format(T[0], C[0], L[0],supercell))
    try:
        check_call(['~/source_codes/mofgen -p {3} -o {0}-{1}-{2} &> {0}-{1}-{2}.log'.format(T[0], C[0], L[0],path )], shell= True)
        tkMessageBox.showinfo('Info', '{0}-{1}-{2} is created'.format(T[0], C[0], L[0]))
    except CalledProcessError:
        tkMessageBox.showerror('Make Error', '{0}-{1}-{2} was not created, check log file for precisions'.format(T[0], C[0], L[0]))
    
def optimize(*args):                        #uses Gulp optimization software to optimize geometry of the MOF using MM and UFF4MOF library
    T= var_topo.get().split('.')
    C = var_center.get().split('.')
    L = var_linker.get().split('.')
    F = var_func.get().split('.')
    if exists('{0}-{1}-{2}-func-{3}.gin'.format(T[0], C[0], L[0], F[0])):
        MOF_to_optimize = '{0}-{1}-{2}-func-{3}'.format(T[0], C[0], L[0], F[0])
    else:
        MOF_to_optimize = '{0}-{1}-{2}'.format(T[0], C[0], L[0])
    with open('{0}.gin'.format(MOF_to_optimize), 'r') as infile:
        line=infile.readlines()
        with open('{0}-opti.gin'.format(MOF_to_optimize), 'w') as outfile:
            outfile.write('opti conp molq noautobond cartesian fix orthorhombic\n')
            for i in line:
                outfile.write(i)
            outfile.write('\n\nlibrary uff4mof.lib\noutput cif {0}-opti.cif\noutput xyz {0}-opti.xyz'.format(MOF_to_optimize))
    try:
        check_call(['gulp {0}-opti'.format(MOF_to_optimize)], shell=True)
        tkMessageBox.showinfo('Info', 'Optimisation of {0} is finished'.format(MOF_to_optimize))
    except CalledProcessError:
        tkMessageBox.showerror('Make Error', 'error while optimizing {0}, check log file for precisions'.format(MOF_to_optimize))

def visualize(renWin_disp,ren_disp, *args):     #renders the result of MOF creation and optimization
    ren_disp.RemoveAllViewProps()
    ren_disp.ResetCamera()
    T= var_topo.get().split('.')
    C = var_center.get().split('.')
    L = var_linker.get().split('.')
    F = var_func.get().split('.')
    MOFname = '{0}-{1}-{2}'.format(T[0], C[0], L[0])
    CIFfile = '{0}-opti.cif'.format(MOFname)
    CIFfile_func ='{0}-func-{1}-opti.cif'.format(MOFname, F[0] ) 
    XYZfile = '{0}.xyz'.format(MOFname)
    XYZfile_func = '{0}-func-{1}.xyz'.format(MOFname, F[0] )
    if exists(CIFfile_func) and exists(XYZfile_func):
        MOF = ase.io.read(CIFfile_func)
        print 'CIF FUNC'
    elif exists(XYZfile_func) and not exists(CIFfile_func):
        MOF = ase.io.read(XYZfile_func)
        print 'XYZ FUNC'
    elif exists(CIFfile) and exists(XYZfile):
        MOF = ase.io.read(CIFfile)
        print 'CIF'
    elif exists(XYZfile) and not exists(CIFfile):
        MOF = ase.io.read(XYZfile)
        print 'XYZ'
    else:    
        tkMessageBox.showerror('Error', 'No files found for {0}. Please create it.'.format(MOFname))   
        MOF = Atoms('H', positions=[(0,0,0)], cell=(1,1,1))     
    i = 0
    elements = MOF.get_chemical_symbols()
    for coordinates in MOF.get_positions():
        scaled_coordinates = (5 * coordinates)
        element_data = elements_dict[elements[i]]
        radius = (element_data[0]/25) 
        color= element_data[1]
        sphereActor = makeSphere(radius, scaled_coordinates,color)
        sphereActor.GetProperty().SetEdgeColor(float('0.'+str(i)),0,0)
        ren_disp.AddActor(sphereActor)
        i += 1
    ren_disp.ResetCamera()
    renWin_disp.Render()

def visualize_linp(*args):          #automatically renders the selected linker
    
    inpfile = '{0}/linkers/{1}'.format(var_SBUpath.get(),var_linker.get())
    linker = ase.io.read(inpfile)
    i = 0
    ren_link.RemoveAllViewProps()
    ren_disp.ResetCamera()
    elements = linker.get_chemical_symbols()
    for coordinates in linker.get_positions():
        scaled_coordinates = (5 * coordinates)
        element_data = elements_dict[elements[i]]
        radius = (element_data[0]/25)
        color= element_data[1]
        linkactor = makeSphere(radius, scaled_coordinates,color)
        linkactor.GetProperty().SetEdgeColor(float('0.'+str(i)),0,0)
        ren_link.AddActor(linkactor)
        i += 1
    ren_link.ResetCamera()
    renWin_link.Render()


def visualize_cinp(*args):          #automatically renders the selected center

    inpfile = '{0}/centers/{1}'.format(var_SBUpath.get(),var_center.get())
    center = ase.io.read(inpfile)
    i = 0
    ren_cent.ResetCamera()
    ren_cent.RemoveAllViewProps()
    elements = center.get_chemical_symbols()
    for coordinates in center.get_positions():
        scaled_coordinates = (5 * coordinates)
        element_data = elements_dict[elements[i]]
        radius = (element_data[0]/25)
        color= element_data[1]
        centactor = makeSphere(radius, scaled_coordinates,color)
        centactor.GetProperty().SetEdgeColor(float('0.'+str(i)),0,0)
        ren_cent.AddActor(centactor)
        i += 1
    ren_cent.ResetCamera()
    renWin_cent.Render()

def updateoptions(*args):       #updates the visible options (available linkers and centers) when changing topology


    actualised_centers = sorted(cdic[var_topo.get()])
    if len(actualised_centers) != 0:
        var_center.set(actualised_centers[0])
    else:
        var_center.set('No center of compatible geometry')

    menu2 = optionCenter.children['menu']
    menu2.delete(0,'end')
    for center in actualised_centers:
        menu2.add_command(label=center, command=lambda v=var_center,l=center: v.set(l))
    
    actualised_linkers = sorted(ldic[var_topo.get()])
    if len(actualised_linkers) !=0 :
        var_linker.set(actualised_linkers[0])
    else:
        var_linker.set('No linker of compatible geometry')

    menu3 = optionLinker.children['menu']
    menu3.delete(0,'end')
    for linker in actualised_linkers:
        menu3.add_command(label=linker, command=lambda v=var_linker,l=linker: v.set(l))


    
    
#####################################################################################################################################################################       ROOTLOOP                 
#############################################################################################################################################################


root = tk.Tk()
root.title("AuToGraFS GUI")

#declaration of tk variables
var_topo = tk.StringVar(root)
var_center = tk.StringVar(root)
var_linker = tk.StringVar(root)
var_SBUpath = tk.StringVar(root)
var_SBUpath.set(os.getcwd())
var_supercell = tk.StringVar(root)
var_supercell.set('1,1,1')
var_updatedisp = tk.DoubleVar(root)
var_func = tk.StringVar(root)



var_SBUpath.trace('w', updateSBU)

#declaration of frames for a beautiful, well organized GUI
frame_buttons = tk.Frame(root)
frame_buttons.grid(row=0,rowspan=2, column=0)
frame_dropdown = tk.Frame(frame_buttons, relief='raised')
frame_dropdown.grid(row=6,rowspan =3, column=0,padx=20, pady=10)
frame_make = tk.Frame(frame_buttons)
frame_make.grid(row=9, column=0)

#creates option menus (drop down) for topologies, and their related centers and linkers
var_topo.set('')
choices_topology = sorted(list(topodict.keys())) 
optionTopology = tk.OptionMenu(frame_dropdown, var_topo, *choices_topology)
optionTopology.grid(column=0, row=0, padx=20, pady=10)
if var_topo.get() != '':
    GEOMS = topodict[var_topo.get()]
else:
    GEOMS = ''

if var_topo.get() != '':
    choices_centers = cdic[var_topo.get()]
else:
    choices_centers = ['No center available']
optionCenter = tk.OptionMenu(frame_dropdown, var_center, *choices_centers)
optionCenter.grid(column=0, row=1, padx=25, pady=10)


if var_topo.get() != '':
    choices_linkers = ldic[var_topo.get()]
else:
    choices_linkers = ['No linker available']
optionLinker = tk.OptionMenu(frame_dropdown, var_linker, *choices_linkers)
optionLinker.grid(column=0, row=2, padx=25, pady=10)

var_topo.trace('w', updateoptions)

#create the buttons
#declares the linked functions
button_make = tk.Button(frame_make, text="Create the selected MOF", command=make)
button_make.grid(column=0, row=0, padx=25, pady=10)
choices_supercell = ['1,1,1','1,1,2','2,2,2','3,3,3','4,4,4']
option_supercell = tk.OptionMenu(frame_make, var_supercell, *choices_supercell)
option_supercell.grid(column=1,row=0,padx=25, pady=10)
button_opti = tk.Button(frame_buttons, text = "Optimize current working MOF", command=optimize)
button_opti.grid(column=0, row=10, padx=25, pady=10)


button_vis = tk.Button(frame_buttons, text = "Refresh visualization", command=lambda:visualize(renWin_disp, ren_disp))
button_vis.grid(column=0, row=11, padx=25, pady=10)
button_askdirSBU = tk.Button(frame_buttons, text='Location of sub-building units (SBU = linkers, centers,function groups)', command=askdirSBU) 
button_askdirSBU.grid(column=0, row = 2,padx=25, pady=10)
label_SBU = tk.Label(frame_buttons, textvariable = var_SBUpath, relief='sunken')
label_SBU.grid(column=0, row=3, padx=15, pady=10)

#declaration of rendering frames for 3D visual effects. 

frame_sbu = tk.Frame(root)
frame_sbu.grid(row=0, column=1, sticky='n'+'w'+'e'+'s')
frame_cent = tk.Frame(frame_sbu, bg='white', relief='sunken')
frame_cent.grid(column=0,row=0,rowspan=2, sticky='n'+'w'+'e'+'s')
frame_link = tk.Frame(frame_sbu, bg='white', relief='sunken')
frame_link.grid(column=1,row=0,rowspan=2, sticky='n'+'w'+'e'+'s')
frame_disp = tk.Frame(root, bg='white', relief='sunken')
frame_disp.grid(column=1, row=1,rowspan=2, sticky='n'+'w'+'e'+'s')

#configuration of rows and column, their behaviour when subjected to resizing
tk.Grid.columnconfigure(root,1,weight=1)
tk.Grid.rowconfigure(root,0,weight=1)
tk.Grid.rowconfigure(root,1,weight=1)
tk.Grid.columnconfigure(frame_sbu,0,weight=1)
tk.Grid.columnconfigure(frame_sbu,1,weight=1)
tk.Grid.rowconfigure(frame_sbu,0,weight=1)
tk.Grid.columnconfigure(frame_disp,0,weight=1)
tk.Grid.rowconfigure(frame_disp,0,weight=1)
tk.Grid.columnconfigure(frame_link,0,weight=1)
tk.Grid.rowconfigure(frame_link,0,weight=1)
tk.Grid.columnconfigure(frame_cent,0,weight=1)
tk.Grid.rowconfigure(frame_cent,0,weight=1)
tk.Grid.columnconfigure(root,0,weight=1)

#creation of renderers for 3D visualization
ren_disp = vtkRenderer()                                #for the created mofs
renWin_disp = vtkRenderWindow()       
renWin_disp.AddRenderer(ren_disp)
disp_widget = vtkTkRenderWindowInteractor(frame_disp,rw=renWin_disp)
disp_widget.grid(row=0,column=0, sticky='n'+'s'+'e'+'w')
disp_widget.Initialize()
disp_widget.AddObserver("StartPickEvent", functionalize)
disp_widget.Start()

ren_cent = vtkRenderer()                                #for the selected center
renWin_cent = vtkRenderWindow()
renWin_cent.AddRenderer(ren_cent)
center_widget = vtkTkRenderWindowInteractor(frame_cent,rw=renWin_cent)
center_widget.grid(row=0,column=0, sticky='n'+'s'+'e'+'w')
center_widget.Initialize()
center_widget.Start()

ren_link = vtkRenderer()                                 #for the selected linker
renWin_link = vtkRenderWindow()
renWin_link.AddRenderer(ren_link)
linker_widget = vtkTkRenderWindowInteractor(frame_link,rw=renWin_link)
linker_widget.grid(row=0,column=0, sticky='n'+'s'+'e'+'w')
linker_widget.Initialize()
linker_widget.Start()

#updates the visualizatyion for linkers and centers
var_linker.trace('w', visualize_linp)
var_center.trace('w', visualize_cinp)


root.mainloop()
