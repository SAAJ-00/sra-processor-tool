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
            output_dir = self.download_sra(srr_id)
            data_type, files = self.convert_to_fastq(srr_id, output_dir)
            
            # Limpiar archivo .sra si no se debe mantener
            if not self.config.get('keep_sra', False):
                self._clean_sra_file(srr_id, output_dir)
            
            return data_type, files
        except Exception as e:
            logger.error(f"Error procesando {srr_id}: {str(e)}")
            raise

    def download_only(self, srr_id):
        """
        Solo descarga el archivo SRA sin conversión
        Retorna: ruta al directorio de salida
        """
        try:
            output_dir = self.download_sra(srr_id)
            logger.info(f"Descarga completada para {srr_id} (sin conversión)")
            return output_dir
        except Exception as e:
            logger.error(f"Error descargando {srr_id}: {str(e)}")
            raise

    def convert_only(self, srr_id):
        """
        Convierte archivo SRA existente a FASTQ
        Retorna: (tipo_datos, lista_archivos)
        """
        output_dir = self.config['output_dir'] / srr_id
        
        # Intentar ubicaciones comunes para archivo SRA
        sra_locations = [
            output_dir / srr_id / f"{srr_id}.sra",  # Estructura típica de prefetch
            output_dir / f"{srr_id}.sra",  # Alternativa: SRA en directorio principal
        ]
        
        sra_file = None
        for location in sra_locations:
            if location.exists():
                sra_file = location
                break
        
        if sra_file is None:
            error_msg = f"Archivo SRA no encontrado para {srr_id} en ninguna ubicación esperada"
            logger.error(error_msg)
            raise ConversionError(error_msg)
        
        try:
            data_type, files = self.convert_to_fastq(srr_id, output_dir)
            
            # Limpiar archivo .sra si no se debe mantener
            if not self.config.get('keep_sra', False):
                self._clean_sra_file(srr_id, output_dir)
            
            return data_type, files
        except Exception as e:
            logger.error(f"Error convirtiendo {srr_id}: {str(e)}")
            raise

    def download_sra(self, srr_id):
        """Descarga un archivo SRA usando prefetch (método público)"""
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

    def convert_to_fastq(self, srr_id, sra_dir):
        """Convierte SRA a FASTQ usando fasterq-dump (método público)"""
        logger.debug(f"Iniciando conversión a FASTQ para {srr_id}...")
        
        # Buscar archivo SRA en ubicaciones comunes
        sra_locations = [
            sra_dir / srr_id / f"{srr_id}.sra",  # Estructura típica de prefetch
            sra_dir / f"{srr_id}.sra",  # Alternativa: SRA en directorio principal
        ]
        
        sra_file = None
        for location in sra_locations:
            if location.exists():
                sra_file = location
                break
        
        if sra_file is None:
            error_msg = f"Archivo SRA no encontrado para conversión: {srr_id}"
            logger.error(error_msg)
            raise ConversionError(error_msg)
        
        cmd = [
            'fasterq-dump',
            '--outdir', str(sra_dir),
            '--temp', str(sra_dir / 'tmp'),
            '--format', 'fastq',
            '--threads', str(self.config['threads']),
            '--split-files',
            '--skip-technical',  # Nuevo parámetro recomendado
            str(sra_file)
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

    def _clean_sra_file(self, srr_id, output_dir):
        """Elimina el archivo .sra después de conversión"""
        # Buscar SRA en ubicaciones comunes
        sra_locations = [
            output_dir / srr_id / f"{srr_id}.sra",
            output_dir / f"{srr_id}.sra",
        ]
        
        for sra_file in sra_locations:
            if sra_file.exists():
                try:
                    sra_file.unlink()
                    logger.debug(f"Archivo SRA eliminado: {sra_file}")
                    # Intentar eliminar directorio del SRA si está vacío
                    sra_dir = sra_file.parent
                    # Comparar paths resueltos para evitar problemas con rutas relativas/absolutas
                    if sra_dir.resolve() != output_dir.resolve() and sra_dir.exists() and not any(sra_dir.iterdir()):
                        sra_dir.rmdir()
                        logger.debug(f"Directorio SRA eliminado: {sra_dir}")
                except Exception as e:
                    logger.warning(f"No se pudo eliminar archivo SRA: {str(e)}")

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