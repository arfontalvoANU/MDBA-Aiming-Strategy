#!/usr/bin/env python

# Dakota will execute this script as
#  The command line arguments will be extracted by dakota.interfacing automatically.

# necessary python modules
import dakota.interfacing as di

# ----------------------------
# Parse Dakota parameters file
# ----------------------------

params, results = di.read_parameters_file()

# -------------------------------
# Convert and send to application
# -------------------------------

# set up the data structures 
# for this simple example, put all the variables into a single hardwired array
# The asv has to be mapped back into an integer
continuous_vars = [params['s1'],params['a1']]
active_set_vector = 0
if results["obj_fn"].asv.function:
    active_set_vector += 1



# set a dictionary for passing to via Python kwargs
dakota_params = {}
dakota_params['cv'] = continuous_vars
dakota_params['asv'] = [active_set_vector,active_set_vector]
dakota_params['functions'] = 2

# execute the analysis as a separate Python module
print "Running simulation..."
from Dakota_aiming import dakota_opt

dakota_results = dakota_opt(**dakota_params)
print "Simulation complete."

# ----------------------------
# Return the results to Dakota
# ----------------------------

# Insert functions from rosen into results
# Results.responses() yields the Response objects.
for i, r in enumerate(results.responses()):
    if r.asv.function:
        r.function = dakota_results['fns'][i]

results.write()
