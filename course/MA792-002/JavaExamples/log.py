import logging
import logging.config

logging.config.fileConfig("log.conf")
logger = logging.getLogger() # root by default

if __name__ == "__main__":
    logger.debug("debug message")
    logger.info("info message")
    logger.warn("warn message")
    logger.error("error message")
    logger.critical("critical message")
