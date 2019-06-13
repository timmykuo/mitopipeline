import shutil, os, unittest
from unittest.mock import patch
from mitopipeline.templates import task_template, task_with_req_template, import_template, paths_template, wrapper_task_template, slurm_task_template, slurm_task_with_req_template
from mitopipeline.pipeline_builder import PipelineBuilder

class TestPB(unittest.TestCase):

    def setUp(self):
        self.builder = PipelineBuilder()
        self.output_path = "./output"
        self.tools_path = "./tools_test"
        self.dir_path = "./dir_test"
        self.refs_path = "./refs_path"
        if not os.path.isdir(self.tools_path):
            os.makedirs(self.tools_path)
        if not os.path.isdir(self.dir_path):
            os.makedirs(self.dir_path)
        if not os.path.isdir(self.refs_path):
            os.makedirs(self.refs_path)
        if not os.path.isdir(self.output_path):
            os.makedirs(self.output_path)

    def cleanUpDirectory(self, directory):
        for f in os.listdir(directory):
            os.remove(os.path.join(directory, f))

    def tearDown(self):
        shutil.rmtree(self.dir_path)
        shutil.rmtree(self.tools_path)
        shutil.rmtree(self.refs_path)
        shutil.rmtree(self.output_path)
        
    def test_buildtemplate_onestep(self):
        steps = ["extractmito"]
        with patch.object(PipelineBuilder, 'get_template', return_value='None') as mock:
            builder = PipelineBuilder()
            builder.build_from_template(self.dir_path, steps, False, self.output_path, self.tools_path, self.refs_path)
        mock.assert_called_once_with(False, "", "extractmito.sh", "extractmito")
    
    def test_buildtemplate_splitgap(self):
        steps = ["extractmito", "splitgap", "removenumts"]
        with patch.object(PipelineBuilder, 'get_template', return_value=self.output_path + "/pipeline.py") as mock:
            builder = PipelineBuilder()
            builder.build_from_template(self.dir_path, steps, False, self.output_path, self.tools_path, self.refs_path)
        mock.return_value = self.output_path + "/pipeline.py"
        mock.assert_any_call(False, "splitgap", "removenumts.sh", "removenumts")
        mock.assert_any_call(False, "extractmito", "splitgap.sh", "splitgap")
        mock.assert_any_call(False, "", "extractmito.sh", "extractmito")

    def test_buildtemplate_nosplitgap(self):
        steps = ["extractmito", "removenumts"]
        with patch.object(PipelineBuilder, 'get_template', return_value=self.output_path + "/pipeline.py") as mock:
            builder = PipelineBuilder()
            builder.build_from_template(
                self.dir_path, steps, False, self.output_path, self.tools_path, self.refs_path)
        mock.return_value = self.output_path + "/pipeline.py"
        mock.assert_any_call(False, "extractmito", "remove_numts_no_split_gap.sh", "removenumts")
        mock.assert_any_call(False, "", "extractmito.sh", "extractmito")

    def test_buildtemplate_softwares(self):
        steps = ["clipping", "removenumts", "gatk", "annovar", "snpeff"]
        with patch.object(PipelineBuilder, 'get_template', return_value=self.output_path + "/pipeline.py") as mock:
            builder = PipelineBuilder()
            builder.build_from_template(
                self.dir_path, steps, False, self.output_path, self.tools_path, self.refs_path)
        mock.assert_any_call(False, "", "clipping.sh", "clipping")
        mock.assert_any_call(False, "clipping", "remove_numts_no_split_gap.sh", "removenumts")
        mock.assert_any_call(False, "removenumts", "gatk.sh", "gatk")
        mock.assert_any_call(False, "gatk", "annovar.sh", "annovar")
        mock.assert_any_call(False, "gatk", "snpeff.sh", "snpeff")
    
    def test_gettemplate_slurm(self):
        with patch.object(slurm_task_template, 'render') as mock:
                self.builder.get_template(True, "", "gatk.sh", "gatk")
        mock.assert_called_once_with(task_name="GATK", job_name="gatk.sh", file_name="gatk.vcf")
        
    def test_gettemplate_slurmreq(self):
        with patch.object(slurm_task_with_req_template, 'render') as mock:
            self.builder.get_template(True, "removenumts", "gatk.sh", "gatk")
        mock.assert_called_once_with(task_name="GATK", req_name="RemoveNuMTs", job_name="gatk.sh", file_name="gatk.vcf")

    def test_gettemplate_noslurm(self):
        with patch.object(task_template, 'render') as mock:
                self.builder.get_template(False, "", "gatk.sh", "gatk")
        mock.assert_called_once_with(task_name="GATK", job_name="gatk.sh", file_name="gatk.vcf")
        
    def test_gettemplate_noslurmreq(self):
        with patch.object(task_with_req_template, 'render') as mock:
            self.builder.get_template(False, "removenumts", "gatk.sh", "gatk")
        mock.assert_called_once_with(task_name="GATK", req_name="RemoveNuMTs", job_name="gatk.sh", file_name="gatk.vcf")


if __name__ == "__main__":
    unittest.main()
        
        
