import os, platform, urllib.request, pkg_resources, subprocess, tarfile, zipfile, requests, io, sys
from mitopipeline.util import is_downloaded, is_exe, cd, which, execute, query_yes_no, get_dir_name
from mitopipeline.logging import Download_Logger
logger = Download_Logger(__name__)
NUM_BARS = 50
class Downloader():

    cache = []

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
                                'removenumts': ['samtools', 'bwa', 'bedtools2'],
                                'clipping': ['samtools', 'bwa', 'seqtk', 'bedtools2'],
                                'extractmito': ['samtools'],
                                'splitgap': ['samtools']}
        self.downloads = {'snpEff': 'http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip/download',
                          'annovar': 'http://annovar.openbioinformatics.org/en/latest/user-guide/download/',
                          'GenomeAnalysisTK': 'https://github.com/broadinstitute/gatk/releases/download/4.1.2.0/gatk-4.1.2.0.zip',
                          'samtools': 'https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2',
                          'bedtools2': 'https://github.com/arq5x/bedtools2/releases/download/v2.28.0/bedtools-2.28.0.tar.gz',
                          'bwa': 'git clone https://github.com/lh3/bwa.git',
                          'seqtk': 'git clone https://github.com/lh3/seqtk.git'}

    #Downloads all dependent softwares of a step if not yet downloaded in the cached tools directory or not available on the command line if necessary
    def download_dependencies(self, steps):
        logger.info('Tools directory is ' + self.TOOLS)
        for step in steps:
            logger.info('Checking software dependencies for ' + step + "...")
            for software in self.dependencies[step]:
                if software not in self.cache:
                    if get_dir_name(software, self.TOOLS) or get_dir_name(step, self.TOOLS):
                        logger.info(software + " is already downloaded in " + self.TOOLS + ". Skipping download.")
                    elif (software == 'samtools' or software == 'bwa' or software == 'bedtools2'):
                        #hack for bedtools2 because it's referencned as just bedtools on the command line, not bedtools2
                        if not which(software) or (software =='bedtools2' and not which('bedtools')):
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
                            os.rename(self.TOOLS + "/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar", self.TOOLS + "/gatk/GenomeAnalysisTK.jar")
                    self.cache.append(software)
            logger.info('All software dependencies for ' + step + ' satisfied.')
        self.cache = []

    def download(self, software, dir):
        logger.info("Downloading " + software + ". This may take awhile...")
        url = self.downloads[software]
        with cd(str(dir)):
            if software == "annovar":
                logger.warning("Annovar must be downloaded through website at " + url + " after registration. Skipping download")
            elif 'git clone' in url:
                self.git_clone(software, url, dir)
            elif 'tar.bz2' in url or 'tar.gz' in url:
                self.download_tar(url)
            elif 'zip' in url:
                self.download_zip(url, software)
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
                        if software == "bedtools2":
                            subprocess.check_output(
                                ['cp', '-a', self.TOOLS + '/' + software + '/bin/', '/usr/local/bin'])
                        else:
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

    def download_zip(self, url, software):
        file_name = software + ".data"
        with open(file_name, "wb") as f:
                response = requests.get(url, stream=True)
                total_length = response.headers.get('content-length')

                if total_length is None:  # no content length header
                    f.write(response.content)
                else:
                    dl = 0
                    total_length = int(total_length)
                    for data in response.iter_content(chunk_size=4096):
                        dl += len(data)
                        f.write(data)
                        done = int(NUM_BARS * dl / total_length)
                        sys.stdout.write("\r|%s%s| %s " % ('|' * done, ' ' * (NUM_BARS-done), "{0:.0f}%".format(done / NUM_BARS * 100)))
                        sys.stdout.flush()
                    sys.stdout.write("\n")
                    sys.stdout.flush()
        logger.info("Extracting zipfile contents for " + software)
        z = zipfile.ZipFile(file_name, mode='r')
        z.extractall()
        os.remove(file_name)

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
        if "tar.bz2" in url:
            tfile = tarfile.open(fileobj=tmpfile, mode="r:bz2")
        else:
            tfile = tarfile.open(fileobj=tmpfile, mode="r:gz")
        tfile.extractall(path=os.getcwd())
        tfile.close()
        tmpfile.close()
