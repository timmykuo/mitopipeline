import unittest, os, shutil
from mitopipeline.util import parse_fid, correct_format, is_valid_directories, files_correct_format, tools_exist, is_annovar_downloaded, found_loc, get_dir_name, is_exe, is_downloaded, make_subdirectories, get_wrapper_tasks, which, execute, query_yes_no

class TestUtils(unittest.TestCase):

    def setUp(self):
        self.output_path = "./output"
        self.tools_path = "./tools_test"
        self.dir_path = "./dir_test"
        self.refs_path = "./refs_path"
        #folder output of the step, and output extension added onto the fild id of each step, ex. extractmito step --> filename_extractmito.bam
        self.task_template_info = {'extractmito': ['ExtractMito', 'extractmito.bam'],
                                   'splitgap': ['SplitGap', '1_splitgap.fastq'],
                                   'clipping': ['Clipping', '1_clipping.fastq'],
                                   'removenumts': ['RemoveNuMTs', 'removenumts.bam'],
                                   'downsample': ['Downsample', 'downsample.bam'],
                                   'gatk': ['GATK', 'gatk.vcf'],
                                   'snpeff': ['SNPEFF', 'snpeff.vcf'],
                                   'annovar': ['ANNOVAR', 'avinput.hg38_multianno.txt'],
                                   #for empty prevstep
                                   '': ['', '']}
        #software dependencies for each step
        self.dependencies = {'snpeff': ['snpEff.jar'],
                             'gatk': ['GenomeAnalysisTK.jar'],
                             'annovar': ['table_annovar.pl', 'convert2annovar.pl'],
                             'removenumts': ['samtools', 'bwa'],
                             'clipping': ['samtools', 'bwa', 'seqtk'],
                             'extractmito': ['samtools'],
                             'splitgap': ['samtools']}
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
    
    def test_parsefid(self):
        f = "NAME.bam"
        self.assertEqual("NAME", parse_fid(f))
        f = "NA.M.E.bam"
        self.assertEqual("NA.M.E", parse_fid(f))
        f = "NAME"
        self.assertEqual("", parse_fid(f))
    
    def test_correctformat(self):
        f = "CORRECT.bam"
        self.assertTrue(correct_format(f))
        f2 = "CORRECT.bai"
        self.assertTrue(correct_format(f2))
        f_wrong = "WRONG.b"
        self.assertFalse(correct_format(f_wrong))

    def test_isvaliddirectories(self):
        #GATK and removeNuMts in steps, should return true
        self.assertTrue(is_valid_directories(self.dir_path, self.refs_path, ["gatk", "removenumts"]))
        #Return error if input directory is not a directory
        self.assertRaises(ValueError, is_valid_directories, "", self.refs_path, ["gatk", "removenumts"])
        #Return error if refs is not a directory
        self.assertRaises(ValueError, is_valid_directories, self.dir_path, "", ["gatk", "removenumts"])
        #Don't return error if refs is not a directory but gatk/removenumts not in steps
        self.assertTrue(is_valid_directories(self.dir_path, "", [""]))
    
    def test_filescorrectformat(self):
        f = open(self.dir_path + "/test.bam", "w")
        f.close()
        f = open(self.dir_path + "/test.bai", "w")
        f.close()
        self.assertTrue(files_correct_format(self.dir_path))
        f = open(self.dir_path + "/test.wrong", "w")
        f.close()
        self.assertRaises(ValueError, files_correct_format, self.dir_path)
        self.cleanUpDirectory(self.dir_path)

    def test_foundloc_commandline(self):
        softwares = ["samtools", "bwa"]
        for software in softwares:
            if shutil.which(software):
                self.assertTrue(found_loc(software, self.tools_path))
            else:
                self.assertRaises(ValueError, found_loc(software, self.tools_path))
    
    def test_foundloc_software(self):
        software = "GenomeAnalysisTK.jar"
        os.makedirs(self.tools_path + "/gatk")
        f = open(self.tools_path + "/gatk/GenomeAnalysisTK.jar", "w")
        f.close()
        self.assertTrue(found_loc(software, self.tools_path))
        software = "snpEff.jar"
        self.assertFalse(found_loc(software, self.tools_path))

    def test_getdirname(self):
        software = "gatk"
        self.assertEqual(None, get_dir_name(software, self.dir_path))
        os.makedirs(self.dir_path + "/" + software)
        self.assertEqual(self.dir_path + "/" + software, get_dir_name(software, self.dir_path))
    
    def test_isexe(self):
        self.assertFalse(is_exe(self.tools_path + "/doesnt_exist.txt"))
        f = open(self.tools_path + "/test.txt", "w")
        f.close()
        self.assertTrue(is_exe(self.tools_path + "/test.txt"))
    
    def test_makesubdirectories(self):
        folder_name = 0
        steps = ["extractmito", "clipping", "annovar"]
        make_subdirectories(self.output_path, self.task_template_info, steps, True)
        assert os.path.exists(self.output_path + "/" + self.task_template_info["extractmito"][folder_name])
        assert os.path.exists(self.output_path + "/" + self.task_template_info["clipping"][folder_name])
        assert os.path.exists(self.output_path + "/" + self.task_template_info["annovar"][folder_name])
        assert os.path.exists(self.output_path + "/slurm")
        steps = ["gatk", "removenumts", "clipping", "extractmito"]
        make_subdirectories(self.output_path, self.task_template_info, steps, True)
        assert os.path.exists(self.output_path + "/" + self.task_template_info["gatk"][folder_name] + "/gatk_stor")
        assert os.path.exists(self.output_path + "/" + self.task_template_info["removenumts"][folder_name] + "/fastqs")
        assert os.path.exists(self.output_path + "/" + self.task_template_info["removenumts"][folder_name] + "/pileups")
        assert os.path.exists(self.output_path + "/" + self.task_template_info["removenumts"][folder_name] + "/numt_removal_stor")

    def test_getwrappertasks(self):
        softwares = ["gatk", "annovar", "snpeff"]
        steps = ["removenumts", "extractmito", "clipping"]
        #no softwares
        self.assertEqual(["RemoveNuMTs"], get_wrapper_tasks(self.task_template_info, steps, softwares))
        #multiple softwares including gatk
        steps = ["removenumts", "extractmito", "gatk", "annovar"]
        self.assertEqual(["ANNOVAR"], get_wrapper_tasks(self.task_template_info, steps, softwares))
        #only gatk
        steps = ["removenumts", "extractmito", "gatk"]
        self.assertEqual(["GATK"], get_wrapper_tasks(self.task_template_info, steps, softwares))
        #no gatk, but includes snpeff/annovar
        steps = ["removenumts", "snpeff"]
        self.assertRaises(ValueError, get_wrapper_tasks, self.task_template_info, steps, softwares)  

if __name__ == '__main__':
    unittest.main()



            
