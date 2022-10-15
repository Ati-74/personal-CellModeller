
def create_simulation_module_file(sample_simulation_file, params, params_patterns):
    """"
        goal: By providing an example of the desired simulation file that includes the parameters to calibrate,
         the user can generate the module needed to run the simulation in ABC.
        @param script_file directory   simulation script python file directory
        @param  params list  Parameters to be calibrated
        @param params_patterns list  Parameters to be calibrated pattern in script file
        @output module_file  str written in script directory
    """
    # read script file
    script_file = open(sample_simulation_file, "r")
    script_file = script_file.read()

    for index, param in enumerate(params):
        # Replace the target string
        replaced_string = param + '''=""" + str(params[''' + str("'") + param + str("'") + ''']) + """ '''
        script_file = script_file.replace(params_patterns[index], replaced_string)

    header = 'def simulation_script(script_name, params): \n    script = """\n'
    footer = '    """ \n\n    # write script \n    f = open("scripts/" + script_name + ".py", "w") \n    ' \
             'f.write(script) \n    f.close()'

    module_file = header + "\n" + script_file + "\n" + footer
    with open("CreateSimulationScript.py", "w") as python_module_file:
        python_module_file.write(module_file)

