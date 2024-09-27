Introduction
============
The package extends the `PythonFMU package <https://github.com/NTNU-IHB/PythonFMU>`_.
It includes the necessary modules to construct a component model according to the fmi, OSP and DNV-RP-0513 standards 
with focus on the following features:

* seamless translation of a Python model to an FMU package with minimal overhead (definition of FMU interface)
* support of vector variables (numpy)
* support of variable units and display units
* support of range checking of variables

Features which facilitate `Assurance of Simulation Models, DNV-RP-0513 <https://standards.dnv.com/explorer/document/6A4F5922251B496B9216572C23730D33/2>`_
shall have a special focus in this package.


Getting Started
---------------
A new model can consist of any python code. To turn the python code into an FMU the following is necessary

#. The model code is wrapped into a Python class which inherits from `Model`
#. The exposed interface variables (model parameters, input- and output connectors) are defined as `Variable` objects
#. The `(model).do_step( time, dt)` function of the model class is extended with model internal code,
   i.e. model evolves from `time` to `time+dt`.
#. Calling the method `Model.build()` will then compile the FMU and package it into a suitable FMU file.

See the files `example_models/bouncing_ball.py` and `tests/test_make_bouncingBall.py` supplied with this package
as a simple example of this process. The first file defines the model class 
and the second file demonstrates the process of making the FMU and using it within fmpy and OSP.


1.	Install the `component_model` package: ``pip install component_model``
2.	Software dependencies: `PythonFMU`, `numpy`, `pint`, `uuid`, `ElementTree`
3.	Latest releases: Version 0.1, based on PythonFMU 0.64

Usage example
-------------
This is another BouncingBall example, using 3D vectors and units.

.. code-block:: Python

    from math import sqrt

    import numpy as np

    from component_model.model import Model
    from component_model.variable import Variable


    class BouncingBall(Model):
        """Another BouncingBall model, made in Python and using Model and Variable to construct a FMU.

        Special features:

        * The ball has a 3-D vector as position and speed
        * As output variable the model estimates the next bouncing point
        * As input variables, the restitution coefficient `e` and the ground angle at the bouncing point can be changed.
        * Internal units are SI (m,s,rad)

        Args:
            pos (np.array)=(0,0,1): The 3-D position in of the ball at time [m]
            speed (np.array)=(1,0,0): The 3-D speed of the ball at time [m/s]
            g (float)=9.81: The gravitational acceleration [m/s^2]
            e (float)=0.9: The coefficient of restitution (dimensionless): |speed after| / |speed before| collision
            min_speed_z (float)=1e-6: The minimum speed in z-direction when bouncing stops [m/s]
        """

        def __init__(
            self,
            name: str = "BouncingBall_3D",
            description="Another BouncingBall model, made in Python and using Model and Variable to construct a FMU",
            pos: tuple = (0, 0, 10),
            speed: tuple = (1, 0, 0),
            g: float = 9.81,
            e: float = 0.9,
            min_speed_z: float = 1e-6,
            **kwargs,
        ):
            super().__init__(name, description, **kwargs)
            self.pos = np.array(pos, dtype=float)
            self.speed = np.array(speed, dtype=float)
            self.a = np.array((0, 0, -g), float)
            self.g = g
            self.e = e
            self.min_speed_z = min_speed_z
            self.stopped = False
            self.time = 0.0
            self.t_bounce, self.p_bounce = self.next_bounce()
            self._interface_variables()

        def _interface_variables(self):
            """Define the FMU2 interface variables, using the variable interface."""
            self._pos = Variable(
                self,
                name="pos",
                description="The 3D position of the ball [m] (height in inch as displayUnit example.",
                causality="output",
                variability="continuous",
                initial="exact",
                start=(str(self.pos[0]) + "m", str(self.pos[1]) + "m", str(self.pos[2]) + "inch"),
                rng=((0, "100 m"), None, (0, "10 m")),
            )
            self._speed = Variable(
                self,
                name="speed",
                description="The 3D speed of the ball, i.e. d pos / dt [m/s]",
                causality="output",
                variability="continuous",
                initial="exact",
                start=tuple(str(x) + "m/s" for x in self.speed),
                rng=((0, "1 m/s"), None, ("-100 m/s", "100 m/s")),
            )
            self._g = Variable(
                self,
                name="g",
                description="The gravitational acceleration (absolute value).",
                causality="parameter",
                variability="fixed",
                start=str(self.g) + "m/s^2",
                rng=(),
            )
            self._e = Variable(
                self,
                name="e",
                description="The coefficient of restitution, i.e. |speed after| / |speed before| bounce.",
                causality="parameter",
                variability="fixed",
                start=self.e,
                rng=(),
            )
            self._p_bounce = Variable(
                self,
                name="p_bounce",
                description="The expected position of the next bounce as 3D vector",
                causality="output",
                variability="continuous",
                start=tuple(str(x) for x in self.p_bounce),
                rng=(),
            )

        def do_step(self, time, dt):
            """Perform a simulation step from `time` to `time + dt`."""
            if not super().do_step(time, dt):
                return False
            self.t_bounce, self.p_bounce = self.next_bounce()
            while dt > self.t_bounce:  # if the time is this long
                dt -= self.t_bounce
                self.pos = self.p_bounce
                self.speed -= self.a * self.t_bounce  # speed before bouncing
                self.speed[2] = -self.speed[2]  # speed after bouncing if e==1.0
                self.speed *= self.e  # speed reduction due to coefficient of restitution
                if self.speed[2] < self.min_speed_z:
                    self.stopped = True
                    self.a[2] = 0.0
                    self.speed[2] = 0.0
                    self.pos[2] = 0.0
                self.t_bounce, self.p_bounce = self.next_bounce()
            self.speed += self.a * dt
            self.pos += self.speed * dt + 0.5 * self.a * dt**2
            if self.pos[2] < 0:
                self.pos[2] = 0
            # print(f"@{time}. pos {self.pos}, speed {self.speed}, bounce {self.t_bounce}")
            return True

        def next_bounce(self):
            """Calculate time until next bounce and position where the ground will be hit,
            based on current time, pos and speed.
            """
            if self.stopped:  # stopped bouncing
                return (1e300, np.array((1e300, 1e300, 0), float))
                # return ( float('inf'), np.array( (float('inf'), float('inf'), 0), float))
            else:
                t_bounce = (self.speed[2] + sqrt(self.speed[2] ** 2 + 2 * self.g * self.pos[2])) / self.g
                p_bounce = self.pos + self.speed * t_bounce  # linear. not correct for z-direction!
                p_bounce[2] = 0
                return (t_bounce, p_bounce)

        def setup_experiment(self, start: float):
            """Set initial (non-interface) variables."""
            super().setup_experiment(start)
            self.stopped = False

The following might be noted:

* The interface variables are defined in a separate local method `_interface_variables`,
  keeping it separate from the model code.
* The `do_step()` method contains the essential code, describing how the ball moves through the air.
  It calls the `super().do_step()` method, which is essential to link it to `Model`.
  The `return True` statement is also essential for the working of the emerging FMU.
* The `next_bounce()` method is a helper method.
* In addition to the extension of `do_step()`, here also the `setup_experiment()` method is extended.
  Local (non-interface) variables can thus be initialized in a convenient way.

It should be self-evident that thorough testing of any model is necessary **before** translation to a FMU.
The simulation orchestration engine (e.g. OSP) used to run FMUs obfuscates error messages, 
such that first stage assurance of a model should aways done using e.g. `pytest`.

The minimal code to make the FMU file package is 

.. code-block:: Python

   from component_model.model import Model
   from fmpy.util import fmu_info

   asBuilt = Model.build("../component_model/example_models/bouncing_ball.py")
   info = fmu_info(asBuilt.name)  # not necessary, but it lists essential properties of the FMU

The model can then be run using `fmpy <https://pypi.org/project/FMPy/>`_ 

.. code-block:: Python

   from fmpy import plot_result, simulate_fmu

   result = simulate_fmu(
       "BouncingBall.fmu",
       stop_time=3.0,
       step_size=0.1,
       validate=True,
       solver="Euler",
       debug_logging=True,
       logger=print,
       start_values={"pos[2]": 2}, # optional start value settings
   )
   plot_result(result)

Similarly, the model can be run using `OSP <https://opensimulationplatform.com/>`_ 
(or rather `libcosimpy <https://pypi.org/project/libcosimpy/>`_ - OSP wrapped into Python):

.. code-block:: Python

   from libcosimpy.CosimEnums import CosimExecutionState
   from libcosimpy.CosimExecution import CosimExecution
   from libcosimpy.CosimSlave import CosimLocalSlave

   sim = CosimExecution.from_step_size(step_size=1e7)  # empty execution object with fixed time step in nanos
   bb = CosimLocalSlave(fmu_path="./BouncingBall.fmu", instance_name="bb")

   print("SLAVE", bb, sim.status())

   ibb = sim.add_local_slave(bb)
   assert ibb == 0, f"local slave number {ibb}"

   reference_dict = {var_ref.name.decode(): var_ref.reference for var_ref in sim.slave_variables(ibb)}

   # Set initial values
   sim.real_initial_value(ibb, reference_dict["pos[2]"], 2.0)

   sim_status = sim.status()
   assert sim_status.current_time == 0
   assert CosimExecutionState(sim_status.state) == CosimExecutionState.STOPPED
   infos = sim.slave_infos()
   print("INFOS", infos)

   # Simulate for 1 second
   sim.simulate_until(target_time=3e9)

This is admittedly more complex than the `fmpy` example,
but it should be emphasised that fmpy is made for single component model simulation (testing),
while OSP is made for multi-component systems.

Contribute
----------
Anybody in the FMU and OSP community is welcome to contribute to this code, to make it better, 
and especially including other features from model assurance, 
as we firmly believe that trust in our models is needed 
if we want to base critical decisions on the support from these models.