import urllib.request, pkg_resources, shutil, subprocess
from mitopipeline.util import is_downloaded
TOOLS = pkg_resources.resource_filename('mitopipeline', "tools")
class Downloader:

    def __init__(self):
        self.dependencies = {'snpeff': ['snpeff'],
                                'annovar': ['annovar'],
                                'gatk': ['gatk', 'picard', 'rCRS'],
                                'removenumts': ['samtools', 'bwa', 'rCRS', 'hg38-nocrs'],
                                'clipping': ['samtools', 'bwa', 'bedtools', 'seqtk'],
                                'extractmito': ['samtools'],
                                'splitgap': ['samtools']}
        self.downloads = {'snpeff': 'http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip/download',
                        'annovar': 'None',
                        'gatk': 'https://software.broadinstitute.org/gatk/download/auth?package=GATK-archive&version=3.1-1-g07a4bf8',
                        'samtools': 'https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2',
                        'bwa': 'git clone https://github.com/lh3/bwa.git',
                        'bam2fastq': 'git clone --recursive https://github.com/jts/bam2fastq',
                        'seqtk': 'git clone https://github.com/lh3/seqtk.git',
                        'rCRS': 'rCRS',
                        'hg38-nocrs': 'hg38-norcrs'}

    def download_dependencies(self, steps, tools=None):
        def is_link(software):
            return software.startswith("https://")
        def is_gitclone(software):
            return software.startswith("git clone")
        for step in steps:
            for software in self.dependencies[step]:
                if not is_downloaded(software):
                    print('Downloading software dependencies for ' + step)
                    if is_link(self.downloads[software]):
                        self.download_link(software, self.downloads[software])
                    elif is_gitclone(self.downloads[software]):
                        self.git_clone(software, self.downloads[software], tools)
                    else:
                        raise ValueError('The software annovar needs registration to be downloaded. You can register here: http://www.openbioinformatics.org/annovar/annovar_download_form.php')
    
    def git_clone(self, software, command, tools):
        if tools:
            subprocess.call(['cd', tools])
        else:
            subprocess.call(['cd', TOOLS])
        subprocess.call([command])
        subprocess.call(['cd', software])
        subprocess.call(['make'])
        
    def download_link(self, software, url):
        # Download the file from `url` and save it locally under `file_name`:
        with urllib.request.urlopen(url) as response, open(str(software), 'wb') as out_file:
            shutil.copyfileobj(response, out_file)
