class SRAProcessingError(Exception):
    """Base class for SRA processing exceptions"""
    pass

class DownloadError(SRAProcessingError):
    """Error during SRA download"""
    pass

class ConversionError(SRAProcessingError):
    """Error during FASTQ conversion"""
    pass

class TrimmingError(SRAProcessingError):
    """Error during trimming"""
    pass

class UnsupportedDataType(SRAProcessingError):
    """Unsupported data type (e.g., single-end when expecting paired-end)"""
    pass