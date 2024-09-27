PythonFMU
=========

The `PythonFMU`package is used heavily in the `component_model` package. 
In this note additional information is collected which is important for the development work here.

Fmi2Slave
---------

* containts the OrderedDict `vars` containing all ScalarVariable of the model as { varRef : ScalarVariable object }
* the model defines also a property with the same name as the variable.name, which contains the current value.
  It can be accessed directly or through getter/setter functions (see below)


ScalarVariable
--------------

* contains getter()/setter() methods which retrieve/set the variable value