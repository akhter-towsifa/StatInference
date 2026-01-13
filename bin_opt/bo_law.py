import law
import os
import yaml
import json

from FLAF.RunKit.run_tools import ps_call
from FLAF.run_tools.law_customizations import Task, HTCondorWorkflow, copy_param

clean_env = {k: os.environ[k] for k in [
    'HOME', 'USER', 'LOGNAME', 'PATH', 'SHELL', 'ANALYSIS_SOFT_PATH', 'FLAF_CMSSW_BASE'
] if k in os.environ}

class ServerTask(HTCondorWorkflow, law.LocalWorkflow):
    max_runtime = copy_param(HTCondorWorkflow.max_runtime, 48.0)
    def __init__(self, *args, **kwargs):
        super(ServerTask, self).__init__(*args, **kwargs)
        self.bin_opt_params = yaml.safe_load(open(os.path.join(os.environ["ANALYSIS_PATH"], "StatInference", "bin_opt", "bin_optimization.yaml")))

    def create_branch_map(self):
        branches = {}
        for channel_id, channel in enumerate(self.bin_opt_params["input"]["channels"]):
            branches[channel_id] = channel
        return branches

    def output(self):
        
        output_path = os.path.join(
            os.environ["ANALYSIS_PATH"], 
            "StatInference", 
            self.bin_opt_params["output"]["directory"], 
            self.branch_data + '_' + str(self.bin_opt_params["input"]["era"])
        )
        return law.LocalFileTarget(output_path)
    
    def run(self):
        print("Test run servertask")
        work_directory = os.path.join(os.environ["ANALYSIS_PATH"], "StatInference")
        print("work_directory", work_directory)

        cmd = "python3 bin_opt/optimize_channel.py --channel {}".format(self.branch_data)
        
        ps_call(
            "bash -c 'source env.sh && cd {} && ".format(work_directory)
             + cmd + "'", 
             shell=True, env=clean_env)