# Paquete de procesamiento SRA

from .config_logging import setup_logging
import logging

__version__ = "0.2.0"

# Configura logging básico si no está configurado
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(name)s: %(message)s'
)

logger = logging.getLogger(__name__)
logger.info(f"Inicializando sra_processor v{__version__}")