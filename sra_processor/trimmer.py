import subprocess
import shutil
from pathlib import Path
from .exceptions import TrimmingError

class FastpProcessor:
    def __init__(self, config):
        self.config = config
        self._validate_fastp()

    def _validate_fastp(self):
        """Verifica que fastp esté instalado"""
        if not shutil.which('fastp'):
            raise EnvironmentError("fastp no encontrado en el PATH")

    def process_paired_end(self, input_files, srr_id):
        """Procesa archivos paired-end con fastp"""
        output_dir = self.config['output_dir'] / srr_id
        output_dir.mkdir(exist_ok=True)
        
        out1 = output_dir / f"output_{srr_id}_1_paired.fastq.gz"
        out2 = output_dir / f"output_{srr_id}_2_paired.fastq.gz"
        report_html = output_dir / f"output_{srr_id}_fastp.html"
        report_json = output_dir / f"output_{srr_id}_fastp.json"

        cmd = [
            'fastp',
            '-i', str(input_files[0]),
            '-I', str(input_files[1]),
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

        try:
            subprocess.run(cmd, check=True, capture_output=True)
            
            if not self.config['keep_temp']:
                self._clean_intermediates(input_files)
                
            return True
        except subprocess.CalledProcessError as e:
            raise TrimmingError(f"Error en fastp: {e.stderr.decode()}")

    def _clean_intermediates(self, input_files):
        """Elimina archivos intermedios"""
        for f in input_files:
            Path(f).unlink(missing_ok=True)