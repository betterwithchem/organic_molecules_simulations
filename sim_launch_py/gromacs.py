import shutil


class _Simulation():
       
    def __init__(self, name: str, mdrun_options: str, topology: str,
                 path: str):

        self._name=name
        self._mdrun_options=mdrun_options
        self._topology=topology
        #self._run_options=run_options
        self._path=path

        @property
        def name(self):
            return self._name

        @property
        def mdrun_options(self):
            return self._mdrun_options

        @property
        def topology(self):
            return self._topology

        #@property
        #def run_options(self):
        #    return self._run_options

        @property
        def path(self):
            return self._path

        @name.setter
        def name(self,n):
            self._name=n


class MD(_Simulation):

    def __init__(self,  name: str, mdrun_options='', coord='', topology='',
                 path_mdp='', path_input='', path_output='', temperature=300, thermostat='vv',
                 pressure=1, barostat='Parrinello-Rahman', nsteps=100, plumed='',
                 print_bash=False, maxwarn=0):

        super().__init__(name, mdrun_options, topology, path_output)

        self.name=name
        # copy mdp files

        shutil.copy(path_mdp,path_input)

        if print_bash:
            mdp=path_mdp
            self.bash_command="gmx grompp -f {0} -p {1} -c {2} -o {3}.tpr -maxwarn {4}\n\n".format(mdp,topology,coord,name,maxwarn)+\
                "gmx mdrun -deffnm {} {} {}\n".format(name,mdrun_options,plumed)
            

class EnergyMinimization(_Simulation):

    def __init__(self,name: str,
                 mdrun_options='', coord='start.pdb', topology='topol.top',
                 path_mdp='', path_input='',path_output='',
                 print_bash=False,maxwarn=0):

        super().__init__(name, mdrun_options, topology, path_output)

        self.name=name
        # copy mdp files

        shutil.copy(path_mdp,path_input)

        if print_bash:            
            mdp=path_mdp
            self.bash_command="gmx grompp -f {0} -p {1} -c {2} -o {3}.tpr -maxwarn {4}\n\n".format(mdp,topology,coord,name,maxwarn)+\
                "gmx mdrun -deffnm {} {}\n".format(name,mdrun_options)
                
            

        



"""
class Metad(MD, _Simulation):

    def __init__(self,plumedfile='plumed.dat'):

        super().__init__()

        self._plumedfile=plumedfile

    @property
    def plumedfile(self):
        return self._plumedfile

    @plumedfile.setter
    def plumedfile(self,p):
        self._plumedfile=p
        
    def add_cv(self):

        newcv=CV()
        self.append(newcv)


    def add_metad(self):

        newmetad=Metad()
        self.append(newmetad)

    def writePlumedFile(self):

        with open(self.plumedfile,'w') as f:
            f.write("# do nothing")

    def writeRunFile(self,runfile):

        with open(runfile,'w') as f:
            f.write("gmx mdrun -deffnm md -v -plumed {0} ".format())
"""

        
        
        


