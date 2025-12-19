import subprocess
import shutil
from pathlib import Path
from .exceptions import TrimmingError
from .utils import open_fastq_file
import logging

logger = logging.getLogger(__name__)

class FastpProcessor:
    def __init__(self, config):
        self.config = config
        self._validate_tools()
    
    def _validate_tools(self):
        """Verifica que las herramientas necesarias estén instaladas"""
        if not shutil.which('fastp'):
            raise TrimmingError(tool='fastp', message="fastp no encontrado en el PATH")
        if not shutil.which('fastplong'):
            logger.warning("ADVERTENCIA: fastplong no encontrado. No se podrá procesar lecturas largas")

    def _get_output_dir(self, srr_id, output_dir=None):
        """
        Determina el directorio de salida apropiado
        
        Args:
            srr_id: ID del SRR o nombre base
            output_dir: Directorio de salida opcional
        
        Returns:
            Path object al directorio de salida
        """
        if output_dir is None:
            result_dir = self.config['output_dir'] / srr_id
        else:
            output_path = Path(output_dir)
            # Usar solo el nombre del archivo para comparación segura
            srr_name = Path(srr_id).name
            # Si el directorio ya termina con el srr_id, no duplicar
            result_dir = output_path if output_path.name == srr_name else output_path / srr_name
        
        result_dir.mkdir(parents=True, exist_ok=True)
        return result_dir
    
    def _is_long_read(self, input_file):
        """Determina si es una lectura larga basado en la longitud promedio"""
        try:
            with open_fastq_file(input_file) as f:
                lengths = []
                for _ in range(25):  # Muestra de 25 secuencias
                    header = next(f, None)
                    seq = next(f, None)
                    plus = next(f, None)
                    qual = next(f, None)
                    if seq is None:
                        break
                    lengths.append(len(seq.strip()))
                
                if not lengths:
                    return False
                
                avg_length = sum(lengths) / len(lengths)
                return avg_length > 500  # Consideramos largas >500bp
        except Exception as e:
            logger.warning(f"Advertencia: No se pudo determinar tipo de lectura: {str(e)}")
            return False

    def _standardize_fastq_extensions(self, srr_id, output_dir, paired=False):
        """Renombra archivos .fq o .fas a .fastq para compatibilidad con fastp"""
        extensions = ['.fastq', '.fq', '.fas']
        if paired:
            for ext in extensions:
                cand1 = output_dir / f"{srr_id}_1{ext}"
                cand2 = output_dir / f"{srr_id}_2{ext}"
                std1 = output_dir / f"{srr_id}_1.fastq"
                std2 = output_dir / f"{srr_id}_2.fastq"
                if cand1.exists() and cand2.exists():
                    if ext != '.fastq':
                        if std1.exists():
                            std1.unlink()
                        if std2.exists():
                            std2.unlink()
                        cand1.rename(std1)
                        cand2.rename(std2)
                    return std1, std2
        else:
            for ext in extensions:
                cand = output_dir / f"{srr_id}{ext}"
                std = output_dir / f"{srr_id}.fastq"
                if cand.exists():
                    if ext != '.fastq':
                        if std.exists():
                            std.unlink()
                        cand.rename(std)
                    return std,
        return None

    def process(self, input_files, srr_id, output_dir=None):
        """
        Procesa cualquier tipo de datos automáticamente
        
        Args:
            input_files: Lista de archivos FASTQ a procesar
            srr_id: Identificador del SRR o nombre base para archivos de salida
            output_dir: Directorio de salida opcional (usa config si no se especifica)
        """
        if len(input_files) == 1 and self._is_long_read(input_files[0]):
            return self._process_long_read(input_files[0], srr_id, output_dir)
        elif len(input_files) == 2:
            return self._process_paired_end(input_files, srr_id, output_dir)
        else:
            return self._process_short_read(input_files[0], srr_id, output_dir)

    def _process_paired_end(self, input_files, srr_id, output_dir=None):
        """Procesa archivos paired-end con fastp"""
        output_dir = self._get_output_dir(srr_id, output_dir)
        
        out1 = output_dir / f"{srr_id}_1_trimmed.fastq.gz"
        out2 = output_dir / f"{srr_id}_2_trimmed.fastq.gz"
        report_html = output_dir / f"report_{srr_id}.html"
        report_json = output_dir / f"report_{srr_id}.json"

        # Verificar si se debe sobrescribir
        force_overwrite = self.config.get('force_overwrite', False)
        if out1.exists() and out2.exists() and not force_overwrite:
            logger.info(f"Archivos de trimming paired-end ya existen para {srr_id}, omitiendo.")
            return True

        # Para archivos externos (rutas absolutas), usar directamente
        if Path(input_files[0]).is_absolute() and Path(input_files[0]).exists():
            std_files = input_files
        else:
            # Estandarizar extensiones antes de trimming
            std_files = self._standardize_fastq_extensions(srr_id, output_dir, paired=True)
            if not std_files:
                logger.error(f"No se encontraron archivos FASTQ para {srr_id} antes del trimming.")
                raise TrimmingError(tool='fastp', message=f"No se encontraron archivos FASTQ para {srr_id}")

        cmd = [
            'fastp',
            '-i', str(std_files[0]),
            '-I', str(std_files[1]),
            '-o', str(out1),
            '-O', str(out2),
            '-h', str(report_html),
            '-j', str(report_json),
            '-w', str(self.config['threads']),
            '--detect_adapter_for_pe',
            '--qualified_quality_phred', str(self.config['trim_params']['quality_phred']),
            '--length_required', str(self.config['trim_params']['min_length']),
            '--cut_front',
            '--cut_tail',
            '--cut_window_size', str(self.config['trim_params']['cut_window_size']),
            '--cut_mean_quality', str(self.config['trim_params']['cut_mean_quality'])
        ]

        if self.config['trim_params']['disable_adapter_trimming']:
            cmd.append('--disable_adapter_trimming')
        if self.config['trim_params']['disable_quality_filtering']:
            cmd.append('--disable_quality_filtering')
        if self.config['trim_params']['disable_length_filtering']:
            cmd.append('--disable_length_filtering')

        return self._run_trimming(cmd, list(std_files))

    def _process_short_read(self, input_file, srr_id, output_dir=None):
        """Procesamiento para single-end cortas"""
        output_dir = self._get_output_dir(srr_id, output_dir)
        
        out = output_dir / f"{srr_id}_trimmed.fastq.gz"
        
        # Verificar si se debe sobrescribir
        force_overwrite = self.config.get('force_overwrite', False)
        if out.exists() and not force_overwrite:
            logger.info(f"Archivo de trimming single-end ya existe para {srr_id}, omitiendo.")
            return True

        # Para archivos externos (rutas absolutas), usar directamente
        if Path(input_file).is_absolute() and Path(input_file).exists():
            std_files = (input_file,)
        else:
            # Estandarizar extensiones antes de trimming
            std_files = self._standardize_fastq_extensions(srr_id, output_dir, paired=False)
            if not std_files:
                logger.error(f"No se encontró archivo FASTQ para {srr_id} antes del trimming.")
                raise TrimmingError(tool='fastp', message=f"No se encontró archivo FASTQ para {srr_id}")

        cmd = [
            'fastp',
            '-i', str(std_files[0]),
            '-o', str(out),
            '-h', str(output_dir / f"report_{srr_id}.html"),
            '-j', str(output_dir / f"report_{srr_id}.json"),
            '-w', str(self.config['threads']),
            '--qualified_quality_phred', str(self.config['trim_params']['quality_phred']),
            '--length_required', str(self.config['trim_params']['min_length']),
            '--cut_front',
            '--cut_tail',
            '--cut_window_size', str(self.config['trim_params']['cut_window_size']),
            '--cut_mean_quality', str(self.config['trim_params']['cut_mean_quality'])
        ]
        
        if self.config['trim_params']['disable_adapter_trimming']:
            cmd.append('--disable_adapter_trimming')
        if self.config['trim_params']['disable_quality_filtering']:
            cmd.append('--disable_quality_filtering')
        if self.config['trim_params']['disable_length_filtering']:
            cmd.append('--disable_length_filtering')

        return self._run_trimming(cmd, list(std_files))

    def _process_long_read(self, input_file, srr_id, output_dir=None):
        """Procesamiento para lecturas largas (Nanopore/PacBio)"""
        output_dir = self._get_output_dir(srr_id, output_dir)
        
        out = output_dir / f"{srr_id}_trimmed.fastq.gz"
        report_html = output_dir / f"report_{srr_id}.html"
        report_json = output_dir / f"report_{srr_id}.json"
        
        # Verificar si se debe sobrescribir
        force_overwrite = self.config.get('force_overwrite', False)
        if out.exists() and not force_overwrite:
            logger.info(f"Archivo de trimming long-read ya existe para {srr_id}, omitiendo.")
            self._remove_empty_tmp(output_dir)
            return True

        cmd = [
            'fastplong',
            '-i', str(input_file),
            '-o', str(out),
            '-j', str(report_json),
            '-h', str(report_html),
            '--mean_qual', str(self.config['trim_params']['long_read_settings']['min_quality']),
            '--length_required', str(self.config['trim_params']['long_read_settings']['min_length']),
            '--thread', str(self.config['threads'])
        ]

        if self.config['trim_params']['long_read_settings']['disable_adapter_trimming']:
            cmd.append('--disable_adapter_trimming')
        if self.config['trim_params']['long_read_settings']['disable_quality_filtering']:
            cmd.append('--disable_quality_filtering')
        
        result = self._run_trimming(cmd, [input_file])
        self._remove_empty_tmp(output_dir)
        
        return result

    def _run_trimming(self, cmd, input_files):
        """Ejecuta el comando de trimming y maneja errores"""
        try:
            subprocess.run(
                cmd,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True
            )
            if not self.config['keep_temp']:
                self._clean_intermediates(input_files)
            # Eliminar tmp vacío después de trimming para todos los casos
            if len(input_files) > 0:
                output_dir = Path(input_files[0]).parent
                self._remove_empty_tmp(output_dir)
            return True
        except subprocess.CalledProcessError as e:
            error_msg = e.stderr if isinstance(e.stderr, str) else e.stderr.decode()
            raise TrimmingError(f"Error en {' '.join(cmd[:2])}: {error_msg}")

    def _clean_intermediates(self, input_files):
        """Elimina archivos intermedios"""
        for f in input_files:
            try:
                Path(f).unlink(missing_ok=True)
            except Exception as e:
                logger.warning(f"Advertencia: No se pudo eliminar {f}: {str(e)}")

    def _remove_empty_tmp(self, output_dir):
        tmp_dir = output_dir / 'tmp'
        if tmp_dir.exists() and tmp_dir.is_dir() and not any(tmp_dir.iterdir()):
            tmp_dir.rmdir()
            