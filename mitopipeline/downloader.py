import os, platform, urllib.request, pkg_resources, shutil, subprocess, tarfile, zipfile, requests, io
from mitopipeline.util import is_downloaded, is_exe, cd, which
TOOLS = pkg_resources.resource_filename('mitopipeline', "tools")
class Downloader:

    def __init__(self):
        self.msgs = []
        self.dependencies = {'snpeff': ['snpeff'],
                                'annovar': ['annovar'],
                                'gatk': ['gatk', 'rCRS'],
                                'removenumts': ['samtools', 'bwa', 'rCRS', 'hg38-nocrs'],
                                'clipping': ['samtools', 'bwa', 'seqtk'],
                                'extractmito': ['samtools'],
                                'splitgap': ['samtools']}
        self.downloads = {'snpeff': 'http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip/download',
                        'annovar': 'None',
                        'gatk': 'https://software.broadinstitute.org/gatk/download/auth?package=GATK-archive&version=3.1-1-g07a4bf8',
                        'samtools': 'https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2',
                        #'samtools': 'git clone https://github.com/samtools/samtools.git',
                        'bwa': 'git clone https://github.com/lh3/bwa.git',
                        'seqtk': 'git clone https://github.com/lh3/seqtk.git',
                        'rCRS': 'rCRS',
                        'hg38-nocrs': 'hg38-norcrs'}
    
    #def refs(self, steps):

    def download_dependencies(self, steps):
        def get_root_folder(software):
            return TOOLS if is_downloaded(software,TOOLS) else TOOLS + "/" + software

        for step in steps:
            print('Downloading software dependencies for ' + step + ".....")
            for software in self.dependencies[step]:
                print('Checking if ' + software + ' is already downloaded...' )
                if (software == 'samtools' or software == 'bwa'):
                    if not which(software):
                        print(software + ' is not available on the command line. Starting download...')
                        self.download(software, TOOLS)
                        self.add_to_command_line(software, TOOLS + "/" + software)
                    else:
                        print(software + " is available on the command line. Skipping download.")
                #TODO: copy from TOOLS to TOOLS + "/" + software
                elif is_downloaded(software, TOOLS):
                    shutil.copyfile(get_root_folder(software), TOOLS + "/" + software)
                else:
                    print("Downloading " + software + " to mitopipeline's tool directory")
                    self.download(software, TOOLS)
        if self.msgs:
            self.print_msgs()
        else:
            print("Downloads successful.")
    
    def download(self, software, dir):
        print("Downloading " + software + ". This may take awhile...")
        url = self.downloads[software]
        with cd(str(dir)):
            if software == "annovar":
                self.cache_msg(software, "annovar")
            if 'git clone' in url:
                self.git_clone(software, url, dir)
            if 'tar.bz2' in url or software == 'gatk':
                self.download_tar(url)
            elif 'zip' in url:
                self.download_zip(url)
            # else:
            #     raise ValueError('The software annovar needs registration to be downloaded. You can register here: http://www.openbioinformatics.org/annovar/annovar_download_form.php')
                
    def git_clone(self, software, link, dir):
        try: 
            subprocess.check_output(link.split(" "))
        except subprocess.CalledProcessError as e:
            print("git clone unsuccessful")
            print(e.output)
            raise
            
        #if make install command is unavailable (samtools)
        try:
            subprocess.check_output(['make', 'install'])
        except subprocess.CalledProcessError:
            print("'make install' command unavailable. Attempting to make " + software)
            try:
                subprocess.check_output(['make'])
            except subprocess.CalledProcessError as e:
                print(e.output)
                print("'make' command unavailable. Please attempt this manually through the git cloned directory at " + dir)
                raise

    def add_to_command_line(self, software, dir):
        with cd(str(dir)):
            #move to /usr/local/bin for command line usage
            if is_exe(dir + "/" + software + "/" + software):
                if platform.system() == 'Darwin':
                    print("Downloaded successfully and make successfully. Looks like you\'re using a Mac. Requesting permission to copy " + software + " executable over to /usr/local/bin...")
                    try:
                        subprocess.check_output(['sudo', 'mv', software, '/usr/local/bin'])
                    except subprocess.CalledProcessError as e:
                        print(e.output)
                        print("There was an error moving the executable to your /usr/local/bin folder. Try doing it manually or adding the directory to your $PATH variable. The path to cloned git repository directory is " + TOOLS + "/" + software)
                        raise 
                elif platform.system() == 'Linux':
                    #code for linux
                    print('Linux')
                elif platform.system() == 'Windows':
                    #code for windows
                    print('Windows')
                else:
                    self.cache_msg(software, "Unknown operating system, please download manually.")
            else:
                self.cache_msg(software, "Make unsuccessful")

    #TODO write download code for each kind of file, i.e. .zip, .tar, git clone, etc
    def download_zip(self, url):
        r = requests.get(url)
        z = zipfile.ZipFile(io.BytesIO(r.content))
        z.extractall()

    def download_tar(self, url):
        tar = tarfile.open(fileobj=urllib.request.urlopen(url), mode="r|bz2")
        print(tar.extractall(path=os.getcwd()))
        tar.close()

    def cache_msg(self, software, msg):
        self.msgs.append("While downloading " + software + ", the following error occured: " + msg)
    
    def print_msgs(self):
        for msg in self.msgs:
            print(msg)
