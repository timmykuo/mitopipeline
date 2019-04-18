import urllib.request, pkg_resources, shutil
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
        self.links = {'snpeff': 'snpeff',
                        'annovar': 'annovar',
                        'gatk': 'gatk',
                        'samtools': 'samtools',
                        'bwa': 'bwa',
                        'bedtools': 'bedtools',
                        'seqtk': 'seqtk',
                        'rCRS': 'rCRS',
                        'hg38-nocrs': 'hg38-norcrs'}

    def download_dependencies(self, steps):
        for step in steps:
            for software in self.dependencies[step]:
                if not is_downloaded(software):
                    print('Downloading software dependencies for ' + step)
                    self.download(software, self.links[software])

    def download(self, software, url):
        # Download the file from `url` and save it locally under `file_name`:
        with urllib.request.urlopen(url) as response, open(str(software), 'wb') as out_file:
            shutil.copyfileobj(response, out_file)
