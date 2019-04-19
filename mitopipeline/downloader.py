import os, platform, urllib.request, pkg_resources, shutil, subprocess, tarfile
from mitopipeline.util import is_downloaded, is_exe
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
                        'picard': 'https://github.com/broadinstitute/picard/releases/download/2.19.1/picard.jar',
                        'rCRS': 'rCRS',
                        'hg38-nocrs': 'hg38-norcrs'}

    def download_dependencies(self, steps, tools=None):
        def in_root_folder(software):
            return is_downloaded(software, TOOLS) or is_downloaded(software, TOOLS + "/" + software)
        def get_root_folder(software):
            return TOOLS if is_downloaded(software,TOOLS) else TOOLS + "/" + software

        for step in steps:
            print('Downloading software dependencies for ' + step)
            for software in self.dependencies[step]:
                print('Checking if ' + software + ' is already downloaded...' )
                if software == 'samtools' or software == 'bwa' and not is_downloaded(software, tools):
                    print(software + ' is not available on the command line. Starting download...')
                    self.download(software, tools)
                elif tools:
                    #if tools directory specified, copy from TOOLS over to tools if it already exists
                    if in_root_folder(software):
                        print(software + " is already downloaded in mitopipeline\'s tool directory. Copying over to specified tools directory")
                        os.makedirs(tools + "/" + software)
                        shutil.copyfile(get_root_folder(software), tools + "/" + software)
                    elif not is_downloaded(software, tools) and not is_downloaded(software, tools + "/" + software):
                        print("Downloading " + software + " to specified tools directory")
                        self.download(software, tools)
                #tools not specified, download to mitopipeline's directory instead
                else:
                    if not in_root_folder(software):
                        print("Downloading " + software + " to mitopipeline's tool directory")
                        self.download(software, TOOLS + "/" + software)

    def download(self, software, tools):
        print("Downloading " + software + ". This may take awhile...")
        def is_link(software):
            return software.startswith("https://")
        def is_gitclone(software):
            return software.startswith("git clone")

        if is_link(self.downloads[software]):
            self.download_link(software, self.downloads[software])
        elif is_gitclone(self.downloads[software]):
            self.git_clone(software, self.downloads[software], tools)
        else:
            raise ValueError('The software annovar needs registration to be downloaded. You can register here: http://www.openbioinformatics.org/annovar/annovar_download_form.php')
                
    def git_clone(self, software, command, tools):
        subprocess.call(['cd', tools])
        subprocess.call([command])
        subprocess.call(['cd', software])
        make = True
        #if make install command is unavailable (samtools)
        try:
            subprocess.check_output(['make', 'install'])
            make = False
        except subprocess.CalledProcessError:
            print("make install command unavailable. Attempting to make " + software)
            subprocess.check_output(['make'])
            make = True

        #move to /usr/local/bin for command line usage
        if make and is_exe(tools + "/" + software + "/" + software) and platform.system() == 'Darwin':
            print("Downloaded successfully and make successfully. Looks like you\'re using a Mac. Requesting permission to copy " + software + " executable over to /usr/local/bin...")
            try:
                subprocess.check_output(['sudo', 'mv', software, '/usr/local/bin'])
            except subprocess.CalledProcessError as e:
                print(e.output)
                print("Try downloading " + software + " manually.")

          
    def download_link(self, software, url):
        # Download the file from `url` and save it locally under `file_name`:
        with urllib.request.urlopen(url) as response, open(str(software), 'wb') as out_file:
            shutil.copyfileobj(response, out_file)
