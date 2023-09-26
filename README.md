# Modified Deviation-Based Aiming (MDBA) Strategy

This set of Python scripts (mdbapy) applies MDBA (modified deviation-based aiming) strategy to optimise heliostat aiming points to respect flux limiation of solar thermal receivers. It includes heliostat field design, receiver flow path design, receiver thermal performance simulation, aiming points optimisation, and annual performance simulation. 

The MDBA method is presented by [Wang et al (2021)](https://doi.org/10.1016/j.solener.2021.07.059)

This MDBA package also needs the following packages:

	* Receiver thermal model: https://github.com/casselineau/Open_CSPERB
	* SOLSTICE, installation guideline: https://github.com/anustg/solstice-scripts.git

## TODO
	* there are scripts duplicated with solsticepy
	* the receiver thermal model is from [Open_CSPERB](https://github.com/casselineau/Open_CSPERB), but a duplicated copy is within the package
	* the role of dakota in this package is not clear NOTE by shuangw1012: Dakota is no longer needed.

## Instruction
To be added

To add this package into your local Python library:

	$ cd $HOME
	$ git clone https://github.com/yewang726/MDBA-Aiming-Strategy.git mdba-scripts
	$ cd mdba-scripts
	To add it into the default path:
	$ pip install . 
	To add it into a user defined path:
	$ pip install . --prefix=/user/defined/path


