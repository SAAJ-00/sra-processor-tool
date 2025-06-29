import subprocess
import shutil
from pathlib import Path
from .exceptions import DownloadError, ConversionError
import logging

logger = logging.getLogger(__name__)

class SRADownloader:
    def __init__(self, config):
        self.config = config
        self._validate_tools()

    def _validate_tools(self):
        """Verifica que las herramientas SRA estén instaladas y disponibles"""
        tools = ['prefetch', 'fasterq-dump']
        missing = [tool for tool in tools if not shutil.which(tool)]
        
        if missing:
            error_msg = f"Herramientas SRA no encontradas: {', '.join(missing)}"
            logger.error(error_msg)
            raise DownloadError(message=error_msg)
        
        logger.debug("Herramientas SRA validadas correctamente")

    def download_and_convert(self, srr_id):
        """
        Descarga y convierte un archivo SRA en un solo paso
        Retorna: (tipo_datos, lista_archivos)
        """
        try:
            output_dir = self._download_sra(srr_id)
            return self._convert_to_fastq(srr_id, output_dir)
        except Exception as e:
            logger.error(f"Error procesando {srr_id}: {str(e)}")
            raise

    def _download_sra(self, srr_id):
        """Descarga un archivo SRA usando prefetch"""
        output_dir = self.config['output_dir'] / srr_id
        output_dir.mkdir(exist_ok=True)
        sra_file = output_dir / srr_id / f"{srr_id}.sra"
        if sra_file.exists():
            logger.info(f"Archivo SRA ya existe para {srr_id}, omitiendo descarga.")
            return output_dir

        logger.info(f"Iniciando descarga de {srr_id}...")

        cmd = [
            'prefetch',
            '--max-size', self.config['max_size'],
            '--output-directory', str(output_dir),
            srr_id
        ]

        try:
            subprocess.run(
                cmd, 
                check=True, 
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True
                )
            logger.debug(f"Descarga completada para {srr_id}")
            return output_dir
        except subprocess.CalledProcessError as e:
            error_msg = f"Error en prefetch para {srr_id}: {e.stderr}"
            logger.error(error_msg)
            raise DownloadError(error_msg)

    def _convert_to_fastq(self, srr_id, sra_dir):
        """Convierte SRA a FASTQ usando fasterq-dump"""
        logger.debug(f"Iniciando conversión a FASTQ para {srr_id}...")
        cmd = [
            'fasterq-dump',
            '--outdir', str(sra_dir),
            '--temp', str(sra_dir / 'tmp'),
            '--format', 'fastq',
            '--threads', str(self.config['threads']),
            '--split-files',
            '--skip-technical',  # Nuevo parámetro recomendado
            str(sra_dir / srr_id / f"{srr_id}.sra")
        ]

        try:
            subprocess.run(
                cmd, 
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True
                )
            return self._validate_output(srr_id, sra_dir)
        except subprocess.CalledProcessError as e:
            error_msg = f"Error en fasterq-dump para {srr_id}: {e.stderr}"
            logger.error(error_msg)
            raise ConversionError(error_msg)

    def _validate_output(self, srr_id, sra_dir):
        """Valida los archivos FASTQ generados y determina el tipo de datos"""
        # Verificar paired-end con extensiones alternativas
        extensions = ['.fastq', '.fq', '.fas']
        for ext in extensions:
            cand1 = sra_dir / f"{srr_id}_1{ext}"
            cand2 = sra_dir / f"{srr_id}_2{ext}"
            if cand1.exists() and cand2.exists():
                logger.info(f"Datos paired-end detectados para {srr_id}")
                return 'paired', [cand1, cand2]
        # Verificar single-end con múltiples extensiones
        for ext in extensions:
            single = sra_dir / f"{srr_id}{ext}"
            if single.exists():
                logger.info(f"Datos single-end detectados para {srr_id}")
                return 'single', [single]
        error_msg = f"No se generaron archivos FASTQ válidos para {srr_id}"
        logger.error(error_msg)
        raise ConversionError(error_msg)
    
    def batch_download(self, srr_list_file):
        """Procesa múltiples SRRs desde un archivo de lista"""
        with open(srr_list_file, encoding="utf-8") as f:
            srr_ids = [line.strip() for line in f if line.strip()]
        
        results = {}
        for srr_id in srr_ids:
            try:
                data_type, files = self.download_and_convert(srr_id)
                results[srr_id] = {
                    'status': 'completed',
                    'type': data_type,
                    'files': files
                }
                logger.info(f"Procesamiento completado para {srr_id}")
            except Exception as e:
                results[srr_id] = {
                    'status': 'failed',
                    'error': str(e)
                }
                logger.error(f"Falló el procesamiento de {srr_id}: {str(e)}")
        
        return results