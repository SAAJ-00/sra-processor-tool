# config_logging.py
import logging
import logging.config
from pathlib import Path

def setup_logging(output_dir="logs"):
    """
    Configura el sistema de logging global
    """
    log_dir = Path(output_dir)
    log_dir.mkdir(exist_ok=True)
    
    logging_config = {
        'version': 1,
        'disable_existing_loggers': False,
        'formatters': {
            'standard': {
                'format': '%(asctime)s [%(levelname)s] %(name)s: %(message)s',
                'datefmt': '%Y-%m-%d %H:%M:%S'
            },
        },
        'handlers': {
            'console': {
                'class': 'logging.StreamHandler',
                'formatter': 'standard',
                'level': 'INFO',
                'stream': 'ext://sys.stdout'
            },
            'file': {
                'class': 'logging.handlers.RotatingFileHandler',
                'formatter': 'standard',
                'filename': str(log_dir / 'sra_processor.log'),
                'maxBytes': 10485760,  # 10MB
                'backupCount': 5,
                'level': 'DEBUG'
            }
        },
        'loggers': {
            '': {  # root logger
                'handlers': ['console', 'file'],
                'level': 'DEBUG',
                'propagate': True
            },
            'sra_processor': {  # nuestro logger principal
                'level': 'DEBUG',
                'propagate': False
            }
        }
    }
    
    logging.config.dictConfig(logging_config)