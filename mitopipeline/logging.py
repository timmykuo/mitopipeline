import logging

class Download_Logger():

    def __init__(self, name):
        self.logger = logging.getLogger(name)
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s:  %(message)s')

    def warning(self, msg):
        self.logger.warning(msg)

    def info(self, msg):
        self.logger.info(msg)
    
    def error(self, msg):
        self.logger.error(msg)

    def exception(self, msg):
        self.logger.exception(msg)
