import logging

class Download_Logger():

    def __init__(self, name):
        self.logger = logging.getLogger(name)
        logging.basicConfig(level=logging.INFO)
        stream_handler = logging.StreamHandler()
        stream_format = logging.Formatter('%(asctime)s - %(levelname)s:  %(message)s')
        stream_handler.setFormatter(stream_format)
        self.logger.addHandler(stream_handler)


    def warning(self, msg):
        self.logger.warning(msg)

    def info(self, msg):
        self.logger.info(msg)
    
    def error(self, msg):
        self.logger.error(msg)

    def exception(self, msg):
        self.logger.exception(msg)
