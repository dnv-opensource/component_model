from collections import OrderedDict
from collections.abc import Iterable
from abc import ABC, abstractmethod
from pint import UnitRegistry
import tempfile
from pathlib import Path
import zipfile
import datetime
import numpy as np
from .variable import Variable, Variable_NP
from math import log
import uuid
import xml.etree.ElementTree as ET
from pythonfmu import Fmi2Slave, FmuBuilder, DefaultExperiment, __version__ as PythonFMU_version
from pythonfmu.fmi2slave import FMI2_MODEL_OPTIONS
from pythonfmu.builder import get_model_description
from pythonfmu.enums import Fmi2Causality as Causality, Fmi2Initial as Initial, Fmi2Variability as Variability


class ModelInitError(Exception):
    '''Special error indicating that something is wrong with the boom definition'''
    pass
class ModelOperationError(Exception):
    '''Special error indicating that something went wrong during crane operation (rotations, translations,calculation of CoM,...)'''
    pass
class ModelAnimationError(Exception):
    '''Special error indicating that something went wrong during crane animation'''
    pass


class Model(Fmi2Slave):
    """ Defines a model including some common model concepts, like variables and units.
    The model interface and the inner working of the model is missing here and must be defined in a given application.
    For a fully defined model instance shall
    
    * define a full set of interface variables
    * set the current variable values
    * run the model in isolation for a time interval
    * retrieve updated variable values
    
    The following FMI concepts are (so far) not implemented
    
    * TypeDefinitions. Instead of defining SimpleType variables, ScalarVariable variables are always based on the pre-defined types and details provided there
    * DisplayUnit. Variable units contain a Unit (the unit as used for inputs and outputs) and BaseUnit (the unit as used in internal model calculations, i.e. based on SI units).
      Additional DisplayUnit(s) are so far not defined/used. Unit is used for that purpose.
    * Explicit <Derivatives> (see <ModelStructure>) are so far not implemented. This could be done as an additional (optional) property in Variable.
      Need to check how this could be done with Variable_NP
    * <InitialUnknowns> (see <ModelStructure>) are so far not implemented.
      Does overriding of setup_experiment(), enter_initialization_mode(), exit_initialization_mode() provide the necessary functionality?
    
    Args:
        name (str): unique model of the instantiated model
        author (str) = 'anonymous': The author of the model
        version (str) = '0.1': The version number of the model
        unitSystem (str)='SI': The unit system to be used. self.uReg.default_system contains this information for all variables
        license (str)=None: License text or file name (relative to source files).
          If None, the BSD-3-Clause license as also used in the component_model package is used, with a modified copyright line
        copyright (str)=None: Copyright line for use in full license text. If None, an automatic copyright line is constructed from the author name and the file date.
       defaultExperiment (dict) = None: key/value dictionary for the default experiment setup
       guid (str)=None: Unique identifier of the model (supplied or automatically generated)
       non_default_flags (dict)={}: Any of the defined FMI flags with a non-default value (see FMI 2.0.4, Section 4.3.1)
    """
    instances = []
    def __init__(self, name, description:str='A component model', author:str='anonymous', version:str='0.1', unitSystem='SI',
                 license:str|None=None, copyright:str|None=None, defaultExperiment:dict=None, nonDefaultFlags:dict=None, guid=None, **kwargs):
        kwargs.update( { 'modelName':name, 'description':description,
                   'author':author, 'version':version, 'copyright':copyright, 'license':license, 
                   'guid':guid, 'default_experiment':defaultExperiment})
        if 'instance_name' not in kwargs: kwargs['instance_name'] = self.make_instance_name(__name__)
        self.check_and_register_instance_name( kwargs['instance_name'])
        if 'resources' not in kwargs: kwargs['resources'] = None 
        super().__init__( **kwargs) # in addition, OrderedDict vars is initialized
        self.name = name
        self.description = description
        self.author = author
        self.version = version
        self.uReg = UnitRegistry( system=unitSystem, autoconvert_offset_to_baseunit=True) # use a common UnitRegistry for all variables
        self.copyright, self.license = self.make_copyright_license( copyright, license)
        self.default_experiment = DefaultExperiment( None, None, None, None) if defaultExperiment is None else DefaultExperiment( **defaultExperiment)
        self.guid = guid if guid is not None else uuid.uuid4().hex
#        print("FLAGS", nonDefaultFlags)
        self.nonDefaultFlags = self.check_flags( nonDefaultFlags)
        self._units = {} # dict of units and displayUnits (unitName : conversionFacto) used in the model (for usage in UnitDefinitions element)
        self.currentTime = 0 # keeping track of time when dynamic calculations are performed
        self.changedVariables = [] # list of input variables which are changed at any time, i.e. the variables which are taken into account during do_step()
                                   # A changed variable is kept in the list until its value is not changed any more
                                   # Adding/removing of variables happens through change_variable()
        self._eventList = [] # possibility for a list of events that will be activated on time during a simulation
                             # Events consist of tuples of (time, changedVariable)

    def do_step(self, currentTime, stepSize):
        '''Do a simulation step of size 'stepSize at time 'currentTime'''
        pass # this need to be re-defined

    def _ensure_unit_registered(self, candidate:Variable|tuple):
        '''Ensure that the displayUnit of a variable is registered. To register the units of a compound variable, the whole variable is entered and a recursive call to the underlying displayUnits is made'''
        if isinstance( candidate, Variable_NP):
            for i in range( len( candidate)): # recursive call to the components
                if candidate.displayUnit is None: self._ensure_unit_registered( (candidate.unit[i], None))
                else: self._ensure_unit_registered( (candidate.unit[i], candidate.displayUnit[i]))
        elif isinstance( candidate, Variable):
            self._ensure_unit_registered( (candidate.unit, candidate.displayUnit))
        else: # here the actual work is done
            if candidate[0] not in self._units: # the main unit is not yet registered
                self._units[ candidate[0]] = [] # main unit has no factor
            if candidate[1] is not None: # displayUnits are defined
                if candidate[1] not in self._units[ candidate[0]]:
                    self._units[ candidate[0]].append( candidate[1])

    def register_variable(self, var):
        '''Register the variable 'var' as model variable. Add the unit if not yet used. Perform some checks and return the basic index (useable as valueReference)
        Note that only the first element of compound variables includes the variable reference, while the following sub-elements contain None, so that an index is reserved.
        Note that the variable name, _initialVal and _unit must be set before calling this function 
        '''
        for idx, v in self.vars.items():
            if v is not None and v.name == var.name:
                raise ModelInitError("Variable name " +var.name +" is not unique in model " +self.name +" already used as reference " +str(idx))
        idx = len( self.vars)
        self.vars[ idx] = var
        if isinstance( var, Variable_NP):
            for i in range(1, len( var)):
                self.vars[idx+i] = None # marking that this is a sub-element
        self._ensure_unit_registered( var)
        return( idx)
    
    @property
    def units(self): return(self._units)
    
    def add_variable(self, *args, **kwargs):
        '''Convenience method, where the model reference is automatically added to the variable initialization'''
        return( Variable( self, *args, **kwargs))
    
    def add_event(self, time:float|None=None, event:tuple=None):
        '''Register a new event to the event list. Ensure that the list is sorted.
        Note that the event mechanism is mainly used for model testing, since normally events are initiated by input variable changes.

        Args:
            time (float): the time at which the event shall be issued. If None, the event shall happen immediatelly
            event (tuple): tuple of the variable (by name or object) and its changed value
        '''
        if event is None: return # no action
        var = event[0] if isinstance( event[0], Variable) else self.variable_by_name( event[0])
        if var is None:
            raise ModelOperationError("Trying to add event related to unknown variable " +str( event[0]) +". Ignored.")
            return        
        if time is None:
            self._eventList.append(-1, (var, event[1])) # append (the list is sorted wrt. decending time)            
        else:
            if not len( self._eventList):
                self._eventList.append( ( time, (var, event[1])))
            else:
                for i, (t,_) in enumerate( self._eventList):
                    if t<time:
                        self._eventList.insert( i, (time,( var, event[1])))
                        break
                    
    def variable_by_name(self, name:str, errorMsg:str|None=None):
        '''Return Variable object related to name, or None, if not found.
        If errorMsg is not None, an error is raised and the message provided
        Note: So far, this does not handle components of compound variables!'''
        for ref, var in self.vars.items():
            if var is not None and var.name == name:
                return( var)
        if errorMsg is not None:
            raise ModelInitError(errorMsg)
        return( None)
    
    def xml_unit_definitions(self):
        '''Make the xml element for the unit definitions used in the model. See FMI 2.0.4 specification 2.2.2'''
        unitDefinitions = ET.Element('UnitDefinitions')
        for u in self._units:
            uBase = self.uReg(u).to_base_units()
            dim = uBase.dimensionality
            exponents = {}
            for key,value in { 'mass':'kg', 'length':'m', 'time':'s', 'current':'A', 'temperature':'K', 'substance':'mol', 'luminosity':'cd'}.items():
                if '['+key+']' in dim:
                    exponents.update({ value : str(int(dim['['+key+']']))})
            if 'radian' in str(uBase.units): # radians are formally a dimensionless quantity. To include 'rad' as specified in FMI standard this dirty trick is used
                uDeg = str(uBase.units).replace('radian','degree')
#                print("EXPONENT", uBase.units, uDeg, log(uBase.magnitude), log(self.uReg('degree').to_base_units().magnitude))
                exponents.update( {'rad':str(int( log(uBase.magnitude) / log(self.uReg('degree').to_base_units().magnitude)))})

            unit = ET.Element("Unit", {'name':u})
            base = ET.Element("BaseUnit", exponents)
            base.attrib.update( {'factor':str(self.uReg(u).to_base_units().magnitude)})
            unit.append( base)
            for dU in self._units[ u]: # list also the displayUnits (if defined)
                unit.append( ET.Element("DisplayUnit", {'name': dU[0], 'factor': str(dU[1])}))
            unitDefinitions.append( unit)
        return( unitDefinitions)

    def make_instanceName(self,  base):
        '''Make a new (unique) instance name, using 'base_#'''
        ext = []
        for name in Model.instances:
            if name.startswith( base+'_') and name[ len(base)+1:].isnumeric():
                ext.append( int( name[ len(base)+1:]))
        return( base+'_'+str(sorted( ext)[-1]+1))
    
    def check_and_register_instance_name(self, iName):
        if any( name==iName for name in Model.instances):
            raise ModelInitError(f"The instance name {iName} is not unique")
        Model.instances.append( iName)

    def make_copyright_license(self, copyright=None, license=None):
        '''Prepares a copyright notice (one line) and a license text (without copyright line)
        If license is None, the license text of the component_model package is used (BSD-3-Clause)
        If copyright is None, a copyright text is construced from self.author and the file date
        '''
        import os, inspect
        import datetime
        pkgPath = inspect.getfile( Model).split(os.path.sep+'component_model'+os.path.sep+'model.py')[0]
        if license is None:
            with open( pkgPath+os.path.sep+'LICENSE', 'r') as f:
                license = f.read()
            license = license.split('\n',2)[2] #Note: the copyright line of component_model cannot be used
        elif license.startswith('Copyright'):
            [c,license] = license.split('\n',1)
            license = license.strip()
            if copyright is None: copyright = c
            
        if copyright is None:
            copyright = "Copyright (c) " +str(datetime.datetime.fromtimestamp(os.path.getatime(__file__)).year) +" "+self.author
            
        return(copyright, license)

# =====================
# FMU-related functions
# =====================
    def build(scriptFile:str|None=None, project_files:list=[], dest:Path = ".", documentation_folder:Path|None = None):
        if scriptFile is None: scriptFile = self.instance_name+'.py'
        print("BUILD", scriptFile, project_files)
        project_files.append( Path(__file__).parents[0])
        print("PROJECT_FILES", project_files)
        with tempfile.TemporaryDirectory() as documentation_dir:
            doc_dir = Path(documentation_dir)
            license_file = doc_dir / "licenses" / "license.txt"
            license_file.parent.mkdir()
            license_file.write_text("Dummy license")
            index_file = doc_dir / "index.html"
            index_file.write_text("dummy index")
            asBuilt = FmuBuilder.build_FMU( scriptFile, project_files=project_files, dest='.', documentation_folder=doc_dir)#, xFunc=None)
            return( asBuilt)
        
    def to_xml(self, model_options:dict={}) -> ET.Element:
        """Build the XML FMI2 modelDescription.xml tree. (adapted from Fmi2Slave.to_xml())
        
        Args:
            model_options ({[str, str]}) : FMU model options
        
        Returns:
            (xml.etree.TreeElement.Element) XML description of the FMU
        """

        t = datetime.datetime.now(datetime.timezone.utc)
        date_str = t.isoformat(timespec="seconds")

        attrib = dict(
            fmiVersion="2.0",
            modelName=self.modelName,
            guid=f"{self.guid!s}",
            generationTool=f"PythonFMU {PythonFMU_version}",
            generationDateAndTime=date_str,
            variableNamingConvention="structured"
        )
        if self.description is not None:
            attrib["description"] = self.description
        if self.author is not None:
            attrib["author"] = self.author
        if self.license is not None:
            attrib["license"] = self.license
        if self.version is not None:
            attrib["version"] = self.version
        if self.copyright is not None:
            attrib["copyright"] = self.copyright

        root = ET.Element("fmiModelDescription", attrib)

        options = dict()
        for option in FMI2_MODEL_OPTIONS:
            value = model_options.get(option.name, option.value)
            options[option.name] = str(value).lower()
        options["modelIdentifier"] = self.modelName
        options["canNotUseMemoryManagementFunctions"] = "true"

        ET.SubElement(root, "CoSimulation", attrib=options)

        root.append( self.xml_unit_definitions())
        
        if len(self.log_categories) > 0:
            categories = ET.SubElement(root, "LogCategories")
            for category, description in self.log_categories.items():
                categories.append( ET.Element( "Category", attrib={"name": category, "description": description},))

        if self.default_experiment is not None:
            attrib = dict()
            if self.default_experiment.start_time is not None:
                attrib["startTime"] = str(self.default_experiment.start_time)
            if self.default_experiment.stop_time is not None:
                attrib["stopTime"] = str(self.default_experiment.stop_time)
            if self.default_experiment.step_size is not None:
                attrib["stepSize"] = str(self.default_experiment.step_size)
            if self.default_experiment.tolerance is not None:
                attrib["tolerance"] = str(self.default_experiment.tolerance)
            ET.SubElement(root, "DefaultExperiment", attrib)

        variables = self.xml_variables()
        root.append( variables) # append <ModelVariables>

        structure = ET.SubElement( root, "ModelStructure")
        outputs = list(
            filter( lambda v: v is not None and v.causality == Causality.output, self.vars.values())
        )
        if outputs:
            outputs_node = ET.SubElement(structure, "Outputs")
            for v in outputs:
                if len(v) == 1:
                    ET.SubElement(outputs_node, "Unknown", attrib=dict( index=str(v.valueReference)))
                else:
                    for i in range( len(v)):
                        ET.SubElement(outputs_node, "Unknown", attrib=dict( index=str(v.valueReference+i)))
        return root

    def xml_variables(self):
        '''Generate the FMI2 modelDescription.xml sub-tree <ModelVariables>'''
        mv = ET.Element( 'ModelVariables')
        for var in self.vars_iter():
            et = var.to_xml()
            if isinstance( et, list): mv.extend( et)
            else:                     mv.append( et)
        return( mv)        

    @staticmethod
    def check_flags( flags):
        '''Check and collect provided flags and return the non-default flags
        Any of the defined FMI flags with a non-default value (see FMI 2.0.4, Section 4.3.1)
        .. todo:: Check also whether the model actually provides these features
        '''
        def check_flag( fl, typ):
            if flags is not None and fl in flags:
                if isinstance( flags[fl], bool) and flags[fl]: # a nondefault value
                    _flags.update( { fl: flags[fl]})
                elif isinstance( flags[fl], int) and flags[fl]!=0: # nondefault integer
                    _flags.update( { fl: flags[fl]})
        _flags = {}
        check_flag( 'needsExecutionTool', bool)
        check_flag( 'canHandleVariableCommunicationStepSize', bool)
        check_flag( 'canInterpolateInputs', bool)
        check_flag( 'maxOutputDerivativeOrder', int)
        check_flag( 'canRunAsynchchronously', bool)
        check_flag( 'canBeInstantiatedOnlyOncePerProcess', bool)
        check_flag( 'canNotUstMemoryManagementFunctions', bool)
        check_flag( 'canGetAndSetFMUstate', bool)
        check_flag( 'canSerializeFMUstate', bool)
        check_flag( 'providesDirectionalDerivative', bool)
        return( _flags)

    def vars_iter(self, key=None):
        '''Iterator over model variables ('vars'). The returned variables depend on 'key' (see below)

        Args:
            key: filter for returned variables. The following possibilities exist:
            
            * None: All variables are returned
            * type: A type designator (int, float, bool, Enum, str), returning only variables matching on this type
            * causality: A Causality value, returning only one causality type, e.g. Causality.input
            * variability: A Variability value, returning only one variability type, e.g. Variability.fixed
            * callable: Any bool function of model variable object
            
        If typ and causality are None, all variables are included, otherwise only the specified type/causality
        The returned list can be indexed to retrieve given valueReference variables.
        '''
        if key is None: # all variables
            for idx, v in self.vars.items():
                if v is not None:
                    yield v
        elif isinstance( key, type): # variable type iterator
            for idx, v in self.vars.items():
                if v is not None and v.type==key:
                    yield v
        elif isinstance( key, Causality):
            for idx, v in self.vars.items():
                if v is not None and v.causality==key:
                    yield v
        elif isinstance( key, Variability):
            for idx, v in self.vars.items():
                if v is not None and v.variability==key:
                    yield v
        elif callable( key):
            for idx, v in self.vars.items():
                if v is not None and key(v):
                    yield v
            
        else:
            raise KeyError(f"Unknown iteration key {key} in 'vars_iter'")

    # ================
    # Need to over-write the get_ and set_ variable access functions, since we need to deal with compound variables
    def ref_to_var(self, vr):
        '''Find Variable and sub-index (for compound variable), based on a valueReference value'''
        _vr = vr
        while True:
            var = self.vars[ _vr]
            if var is None:
                _vr -= 1
            else: # found the base of the variable
                return( var, vr-_vr)
    
    def _get(self, vrs:list, typ:type) -> list:
        '''Generic get function covering all types. This method is called by get_xxx and translates to fmi2GetXxx'''
        refs = list()
        for vr in vrs:
            if vr >= len(self.vars):
                raise KeyError(f"Variable with valueReference={vr} does not exist in model {self.name}")
            var = self.vars[vr]
            if var is not None and var.type != typ:
                raise TypeError( f"Variable with valueReference={vr} is not of type " +typ.__name__)
            if var is None: # sub-element of a compound variable
                var,sub = self.ref_to_var( vr)
            else:
                sub = 0
            if isinstance( var, Variable_NP):
                refs.append( var.unit_convert( var.value[sub].tolist(), toBase=False))
            else: # non-compound variable. Just append to refs
                if var.displayUnit is None:
                    refs.append( var.value)
                else:
                    refs.append( var.unit_convert( var.value, toBase=False))
        return( refs)
    def get_integer(self, vrs):   return( self._get( vrs, int)) 
    def get_real(self, vrs):      return( self._get( vrs, float))
    def get_boolean(self, vrs):   return( self._get( vrs, bool))
    def get_string(self, vrs):    return( self._get( vrs, str))

    def _set(self, vrs:list, values:list, typ:type):
        '''Generic set function covering all types. This method is called by set_xxx and translates to fmi2SetXxx.
        Variable range check, unit check and unit conversion is performed here.
        Compound variables do not really exist in fmi2 and therefore complicate this function.
        With respect to non-scalar variables we only set the component value,
        but final 'triggers' need to wait until all values are set (which may come in any order)
        '''
        nonScalar = {} 
        for vr, value in zip(vrs, values):
            if vr >= len(self.vars):
                raise KeyError(f"Variable with valueReference={vr} does not exist in model {self.name}")
            var = self.vars[vr]
            if var is None: # sub-element of a compound variable
                var,sub = self.ref_to_var( vr)
            else: # component 0 of a compound variable, or a scalar
                sub = 0
            if var.type != typ:
                raise TypeError( f"Variable with valueReference={vr} is not of type " +typ.__name__)
            if isinstance( var, Variable_NP):
                if var in nonScalar: #already registered
                    nonScalar[var][sub] = var.unit_convert( value, sub)
                else:
                    val = var.value # the current value of the whole array
                    val[sub] = var.unit_convert( val, sub) # change component 'sub
                    nonScalar.update({var: val}) # .. and register
            elif isinstance( var, Variable):
                var.value = var.unit_convert( value)
        # finally we set all compound variables (in one piece)
        for var in nonScalar:
            if isinstance(var, Variable_NP):
                var.value = np.array( nonScalar[var])

    def set_integer(self, vrs:list, values: list):
        self._set( vrs, values, int)
    def set_real(self, vrs:list, values:list):
        self._set( vrs, values, float)
    def set_boolean(self, vrs:list, values:list):
        self._set( vrs, values, bool)
    def set_string(self, vrs:list, values:list):
        self._set( vrs, values, str)

    def _get_fmu_state(self) -> dict:
        state = dict()
        for var in self.vars.values():
            state[var.name] = var.getter()
        return state

    def _set_fmu_state(self, state: dict):
        vars_by_name = dict([(v.name, v) for v in self.vars.values()])
        for name, value in state.items():
            if name not in vars_by_name:
                setattr(self, name, value)
            else:
                v = vars_by_name[name]
                if v.setter is not None:
                    v.setter(value)
