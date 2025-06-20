import subprocess
import shutil
from pathlib import Path
from .exceptions import DownloadError, ConversionError, UnsupportedDataType

class SRADownloader:
    def __init__(self, config):
        self.config = config
        self.sra_tools = self._check_sra_tools()

    def _check_sra_tools(self):
        """Verifica que las herramientas SRA estén disponibles"""
        tools = ['prefetch', 'fasterq-dump']
        missing = [tool for tool in tools if not shutil.which(tool)]
        if missing:
            raise EnvironmentError(f"Herramientas SRA no encontradas: {', '.join(missing)}")
        return True

    def download_sra(self, srr_id):
        """Descarga un archivo SRA usando prefetch"""
        output_dir = self.config['output_dir'] / srr_id
        output_dir.mkdir(exist_ok=True)

        cmd = [
            'prefetch',
            '--max-size', self.config['max_size'],
            '--output-directory', str(output_dir),
            srr_id
        ]

        try:
            subprocess.run(cmd, check=True, capture_output=True)
            return output_dir
        except subprocess.CalledProcessError as e:
            raise DownloadError(f"Error descargando {srr_id}: {e.stderr.decode()}")

    def convert_to_fastq(self, srr_id, sra_dir):
        """Convierte SRA a FASTQ usando fasterq-dump"""
        cmd = [
            'fasterq-dump',
            '--outdir', str(sra_dir),
            '--temp', str(sra_dir / 'tmp'),
            '--format', 'fastq',
            '--threads', str(self.config['threads']),
            '--split-files',
            '--skip-technical',  # Nuevo parámetro recomendado
            '--format', 'fastq',  # Fuerza formato FASTQ
            str(sra_dir / srr_id / f"{srr_id}.sra")
        ]

        try:
            subprocess.run(cmd, check=True)
            return self._check_fastq_files(srr_id, sra_dir)
        except subprocess.CalledProcessError as e:
            raise ConversionError(f"Error convirtiendo {srr_id}: {e.stderr.decode()}")

    def _check_fastq_files(self, srr_id, sra_dir):
        """Verifica los archivos FASTQ generados con múltiples extensiones posibles"""
        extensions = ['.fastq', '.fq', '.fas']  # Todas las extensiones posibles
        
        # Busca archivos paired-end
        for ext in extensions:
            paired_1 = sra_dir / f"{srr_id}_1{ext}"
            paired_2 = sra_dir / f"{srr_id}_2{ext}"
            if paired_1.exists() and paired_2.exists():
                return 'paired', [paired_1, paired_2]
        
        # Busca archivos single-end
        for ext in extensions:
            single = sra_dir / f"{srr_id}{ext}"
            if single.exists():
                return 'single', [single]
        
        raise ConversionError(f"No se encontraron archivos FASTQ para {srr_id} en {sra_dir}")