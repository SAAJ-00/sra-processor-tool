"""
Módulo de excepciones personalizadas para el pipeline SRA Processor

Las excepciones siguen una jerarquía lógica para mejor manejo de errores:
SRAProcessingError (base)
├── DownloadError
├── ConversionError
├── TrimmingError
├── UnsupportedDataTypeError
└── ConfigurationError
"""

class SRAProcessingError(Exception):
    """Clase base para todas las excepciones del pipeline"""
    def __init__(self, message="Error en el procesamiento SRA"):
        self.message = message
        super().__init__(self.message)

class DownloadError(SRAProcessingError):
    """Error durante la descarga de datos SRA"""
    def __init__(self, srr_id=None, message=None):
        self.srr_id = srr_id
        full_msg = f"Error descargando {srr_id}: {message}" if srr_id else message
        super().__init__(full_msg)

class ConversionError(SRAProcessingError):
    """Error durante la conversión SRA a FASTQ"""
    def __init__(self, srr_id=None, message=None):
        self.srr_id = srr_id
        full_msg = f"Error convirtiendo {srr_id}: {message}" if srr_id else message
        super().__init__(full_msg)

class TrimmingError(SRAProcessingError):
    """Error durante el proceso de trimming"""
    def __init__(self, tool=None, message=None):
        self.tool = tool
        full_msg = f"Error en {tool}: {message}" if tool else message
        super().__init__(full_msg)

class UnsupportedDataTypeError(SRAProcessingError):
    """Error cuando se encuentra un tipo de datos no soportado"""
    def __init__(self, data_type=None, message=None):
        self.data_type = data_type
        if not message:
            message = f"Tipo de datos no soportado: {data_type}" if data_type else "Tipo de datos no soportado"
        super().__init__(message)

class ConfigurationError(SRAProcessingError):
    """Error en la configuración del pipeline"""
    def __init__(self, parameter=None, message=None):
        self.parameter = parameter
        if not message:
            message = f"Error de configuración en parámetro: {parameter}" if parameter else "Error de configuración"
        super().__init__(message)

class DiskSpaceError(SRAProcessingError):
    """Error por espacio insuficiente en disco"""
    def __init__(self, required=None, available=None):
        message = "Espacio en disco insuficiente"
        if required and available:
            message += f" (Se necesitan {required} bytes, disponibles {available} bytes)"
        super().__init__(message)