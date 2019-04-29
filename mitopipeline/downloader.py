import os, platform, urllib.request, pkg_resources, subprocess, tarfile, zipfile, requests, io
from mitopipeline.util import is_downloaded, is_exe, cd, which, execute, query_yes_no, get_dir_name
from mitopipeline.logging import Download_Logger
logger = Download_Logger(__name__)

class Downloader():

    def __init__(self):
        self.msgs = []
        #hack to create tools folder in cache
        try:
            self.TOOLS = pkg_resources.resource_filename('mitopipeline', "tools")
        except KeyError:
            self.TOOLS = pkg_resources.resource_filename('mitopipeline', "/") + "/tools"
        if not os.path.isdir(self.TOOLS):
            os.makedirs(self.TOOLS)
        
        self.dependencies = {'snpeff': ['snpEff'],
                                'annovar': ['annovar'],
                                'gatk': ['GenomeAnalysisTK'],
                                'removenumts': ['samtools', 'bwa'],
                                'clipping': ['samtools', 'bwa', 'seqtk'],
                                'extractmito': ['samtools'],
                                'splitgap': ['samtools']}
        self.downloads = {'snpEff': 'http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip/download',
                          'annovar': 'http://annovar.openbioinformatics.org/en/latest/user-guide/download/',
                          'GenomeAnalysisTK': 'https://software.broadinstitute.org/gatk/download/auth?package=GATK-archive&version=3.4-46-gbc02625',
                        'samtools': 'https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2',
                        #'samtools': 'git clone https://github.com/samtools/samtools.git',
                        'bwa': 'git clone https://github.com/lh3/bwa.git',
                        'seqtk': 'git clone https://github.com/lh3/seqtk.git'}
    
    #def refs(self, steps):

    def download_dependencies(self, steps):
        for step in steps:
            logger.info('Checking softwar software dependencies for ' + step + "...")
            for software in self.dependencies[step]:
                if get_dir_name(software, self.TOOLS) or get_dir_name(step, self.TOOLS):
                    logger.info(software + " is already downloaded in " + self.TOOLS + ". Skipping download.")
                elif (software == 'samtools' or software == 'bwa'):
                    if not which(software):
                        logger.info(software + ' is not available on the command line. Starting download...')
                        self.download(software, self.TOOLS)
                        dir_name = get_dir_name(software, self.TOOLS)
                        self.make(software, dir_name)
                        self.add_to_command_line(software, dir_name)
                    else:
                        logger.info(software + " is available on the command line. Skipping download.")
                else:
                    logger.info("Downloading " + software + " to mitopipeline's tool directory")
                    self.download(software, self.TOOLS)
                    if software == "seqtk":
                        dir_name = get_dir_name(software, self.TOOLS)
                        self.make(software, dir_name)
                    if software == "GenomeAnalysisTK":
                        if not os.path.isdir(self.TOOLS + "/gatk"):
                            os.makedirs(self.TOOLS + "/gatk")
                        os.rename(self.TOOLS + "/GenomeAnalysisTK.jar", self.TOOLS + "/gatk/GenomeAnalysisTK.jar")

    def download(self, software, dir):
        logger.info("Downloading " + software + ". This may take awhile...")
        url = self.downloads[software]
        with cd(str(dir)):
            if software == "annovar":
                logger.warning("Annovar must be downloaded through website at " + url + " after registration. Skipping download")
            elif 'git clone' in url:
                self.git_clone(software, url, dir)
            elif 'tar.bz2' in url or software == 'GenomeAnalysisTK':
                self.download_tar(url)
            elif 'zip' in url:
                self.download_zip(url)
            else:
                raise ValueError(software + ' url not found')
                
    def git_clone(self, software, link, dir):
        if not os.path.isdir(dir + "/" + software):
            try:
                for path in execute(link.split(" ")):
                    print(path, end="")
            except subprocess.CalledProcessError as e:
                logger.error(e)
                raise e

    def make(self, software, dir):    
        with cd(str(dir)):
            #if make install command is unavailable (samtools)
            try:
                for path in execute(["make"]):
                    print(path, end="")
            except subprocess.CalledProcessError as e:
                logger.error("'make' command unavailable. Please attempt this manually through the git cloned directory at " + dir)
                raise e

    def add_to_command_line(self, software, dir):
        with cd(str(dir)):
            #move to /usr/local/bin for command line usage
            if platform.system() == 'Darwin' or platform.system() == 'Linux':
                logger.info("Downloaded successfully and make successfully. Looks like you\'re using a Unix machine.")
                if query_yes_no("Permission to copy " + software + " executable over to /usr/local/bin?"):
                    try:
                        subprocess.check_output(['mv', software, '/usr/local/bin'])
                    except subprocess.CalledProcessError as e:
                        logger.error("There was an error moving the executable to your /usr/local/bin folder. Try doing it manually or adding the directory to your $PATH variable. The path to cloned git repository directory is " + self.TOOLS + "/" + software)
                        raise e
                else:
                    logger.info("Exiting... " + software + " was downloaded into " + dir)
            elif platform.system() == 'Windows':
                logger.warning("Looks like this is a Windows machine. Please download software manually.")
            else:
                logger.warning("When trying to add to the command line, this is an unknown operating system. Please download manually.")

    def download_zip(self, url):
        r = requests.get(url)
        z = zipfile.ZipFile(io.BytesIO(r.content))
        z.extractall()

    def download_tar(self, url):
        ftpstream = urllib.request.urlopen(url)

        tmpfile = io.BytesIO()
        while True:
            s = ftpstream.read(16384)
            if not s:
                break

            tmpfile.write(s)
        ftpstream.close()

        #seek back to the beginning of the temporary file.
        tmpfile.seek(0)

        tfile = tarfile.open(fileobj=tmpfile, mode="r:bz2")
        tfile.extractall(path=os.getcwd())
        tfile.close()
        tmpfile.close()
