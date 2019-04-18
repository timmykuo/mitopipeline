from mitopipeline.cmdline_parser import CommandLineParser
from mitopipeline.downloader import Downloader

def run():
    cmdline_parser = CommandLineParser()
    opts = cmdline_parser.parse_commands()
    if opts.download:
        downloader = Downloader()
        downloader.download(opts.steps)
    cmdline_parser.build_and_run()