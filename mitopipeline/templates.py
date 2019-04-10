from jinja2 import Template

#template for "import module" and "from module import function"
import_template = Template("{% for module in imports %}import {{ module }}\n{% endfor %}")

#template for directory and output paths
paths_template = Template(
    'STEPS = "./steps"\nPIPELINE_START = "{{ directory }}"\nPIPELINE_STORE = "{{ output }}"\nSLURM_DIR = "{{ output }}slurm"\n')

#template for a task without any requirements
task_template = Template(
    "class {{ task_name }}(luigi.Task):\n    id=luigi.Parameter()\n    def run(self):\n        subprocess.call([STEPS + \'/{{ job_name }}\', self.id, PIPELINE_START, PIPELINE_STORE + \'/{{ task_name }}\'])\n    def output(self):\n        return luigi.LocalTarget(PIPELINE_STORE + \'/{{ task_name }}/\' + \'{}_{{ file_name }}\'.format(self.id))")

#template for a slurm task without requirements
slurm_task_template = Template(
    "class {{ task_name }}(luigi.Task):\n    id=luigi.Parameter()\n    def run(self):\n        subprocess.call([STEPS + \'/submit_job.sh\', self.id, {{ job_name }}, SLURM_DIR, PIPELINE_START, PIPELINE_STORE + \'/{{ task_name }}\'])\n    def output(self):\n        return luigi.LocalTarget(PIPELINE_STORE + \'/{{ task_name }}/\' + \'{}_{{ file_name }}\'.format(self.id))")

#template for task with a required task
task_with_req_template = Template(
    "class {{ task_name }}(luigi.Task):\n    id=luigi.Parameter()\n    def require(self):\n        return {{ req_name }}(id=self.id)\n    def run(self):\n        subprocess.call([STEPS + \'/{{ job_name }}\', self.id, PIPELINE_START + \'/{{ req_name }}\', PIPELINE_STORE + \'/{{ task_name }}\'])\n    def output(self):\n        return luigi.LocalTarget(PIPELINE_STORE + \'/{{ task_name }}/\' + \'{}_{{ file_name }}\'.format(self.id))")

#slurm template for task with a required task
slurm_task_with_req_template = Template(
    "class {{ task_name }}(luigi.Task):\n    id=luigi.Parameter()\n    def require(self):\n        return {{ req_name }}(id=self.id)\n    def run(self):\n        subprocess.call([STEPS + \'/submit_job.sh\', self.id, {{ job_name }}, SLURM_DIR, PIPELINE_START + \'/{{ req_name }}\', PIPELINE_STORE + \'/{{ task_name }}\'])\n    def output(self):\n        return luigi.LocalTarget(PIPELINE_STORE + \'/{{ task_name }}/\' + \'{}_{{ file_name }}\'.format(self.id))")

#template for the wrapper task
wrapper_task_template = Template("class {{ task_name }}(luigi.WrapperTask):\n    _file = luigi.Parameter(default=\"\")\n    def requires(self):\n        if PIPELINE_START:\n            for f in os.listdir(PIPELINE_START):\n                if os.path.isfile(PIPELINE_START + f):\n                    {% for yield in yields %}yield {{ yield }}(id=util.parse_fid(f))\n                    {% endfor %}\n        else:\n            {% for yield in yields %}yield {{ yield }}(id=util.parse_fid(f))\n            {% endfor %}")
