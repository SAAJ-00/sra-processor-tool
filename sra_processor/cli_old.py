#!/usr/bin/env python3
import sys
import argparse
import shutil
from pathlib import Path
from .downloader import SRADownloader
from .trimmer import FastpProcessor
from .utils import check_srr_status
from .config import DEFAULT_CONFIG
from .exceptions import SRAProcessingError, UnsupportedDataTypeError
import logging
from .config_logging import setup_logging

def process_srr(srr_id, config):
    """Procesa un solo SRR a través de todo el pipeline"""
    downloader = SRADownloader(config)
    trimmer = FastpProcessor(config)

    status = check_srr_status(srr_id, config['output_dir'])

    try:
        if status in ['new', 'sra_downloaded']:
            print(f"Procesando {srr_id} (estado: {status})...")

            if status == 'new':
                data_type, fastq_files = downloader.download_and_convert(srr_id)
            else:
                # Si ya está descargado, solo convertir
                sra_dir = config['output_dir'] / srr_id
                data_type, fastq_files = downloader._convert_to_fastq(srr_id, sra_dir)

            if data_type == 'single':
                print(f"{srr_id} es single-end. Procesando normalmente...")

            trimmer.process(fastq_files, srr_id)
            return True

        elif status == 'fastq_ready':
            print(f"Reanudando procesamiento para {srr_id} (FASTQs listos)...")
            fastq_files = [
                config['output_dir'] / srr_id / f"{srr_id}_1.fastq",
                config['output_dir'] / srr_id / f"{srr_id}_2.fastq"
            ]
            trimmer.process(fastq_files, srr_id)
            return True

        elif status == 'complete':
            print(f"Saltando {srr_id} (ya procesado completamente)")
            return True

        elif status == 'single_end':
            print(f"Reanudando procesamiento para {srr_id} (FASTQ single-end listo)...")
            # Buscar archivo FASTQ single-end con cualquier extensión
            extensions = ['.fastq', '.fq', '.fas']
            for ext in extensions:
                fastq_file = config['output_dir'] / srr_id / f"{srr_id}{ext}"
                if fastq_file.exists():
                    trimmer.process([fastq_file], srr_id)
                    return True
            print(f"No se encontró archivo FASTQ single-end para {srr_id}")
            return False

        else:  # 'unknown'
            print(f"Estado desconocido para {srr_id}. Intentando procesar como nuevo...")
            # Limpia el directorio y reintenta como nuevo
            import shutil
            srr_dir = config['output_dir'] / srr_id
            if srr_dir.exists():
                shutil.rmtree(srr_dir)
            return process_srr(srr_id, config)

    except SRAProcessingError as e:
        print(f"Error procesando {srr_id}: {str(e)}")
        return False

def main():
    parser = argparse.ArgumentParser(
        description='Herramienta para procesamiento de datos desde SRA IDs (Illumina/Nanopore/PacBio)',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Argumentos principales
    parser.add_argument(
        'srr_list',
        help='Archivo con lista de SRRs (uno por línea) o lista de archivos FASTQ'
    )
    parser.add_argument(
        '-o', '--output',
        default='.',
        help='Directorio de salida'
    )
    parser.add_argument(
        '-t', '--threads',
        type=int,
        default=8,
        help='Número de hilos a usar'
    )
    parser.add_argument(
        '--max-size',
        default='30G',
        help='Tamaño máximo para descarga por cada SRR (ej. 30G, 100G)'
    )
    parser.add_argument(
        '--keep-temp',
        action='store_true',
        help='Mantener archivos temporales/intermedios'
    )
    parser.add_argument(
        '--force-long-reads',
        action='store_true',
        help='Forzar tratamiento como lecturas largas (Nanopore/PacBio)'
    )
    
    # Grupo para lecturas cortas (Illumina)
    short_read_group = parser.add_argument_group('Opciones para lecturas cortas (Illumina)')
    short_read_group.add_argument(
        '--quality-phred',
        type=int,
        default=30,
        help='Calidad Phred mínima para lecturas cortas'
    )
    short_read_group.add_argument(
        '--min-length',
        type=int,
        default=50,
        help='Longitud mínima de reads después de trimming (bp)'
    )
    short_read_group.add_argument(
        '--cut-window-size',
        type=int,
        default=4,
        help='Tamaño de ventana para sliding window trimming'
    )
    short_read_group.add_argument(
        '--cut-mean-quality',
        type=int,
        default=25,
        help='Calidad media para sliding window trimming'
    )
    short_read_group.add_argument(
        '--disable-adapter-trimming',
        action='store_true',
        help='Deshabilitar trimming de adaptadores'
    )
    short_read_group.add_argument(
        '--disable-quality-filtering',
        action='store_true',
        help='Deshabilitar filtrado por calidad'
    )
    short_read_group.add_argument(
        '--disable-length-filtering',
        action='store_true',
        help='Deshabilitar filtrado por longitud'
    )
    
    # Grupo para lecturas largas (Nanopore/PacBio)
    long_read_group = parser.add_argument_group('Opciones para lecturas largas (Nanopore/PacBio)')
    long_read_group.add_argument(
        '--long-min-quality',
        type=int,
        default=10,
        help='Calidad mínima para lecturas largas'
    )
    long_read_group.add_argument(
        '--long-min-length',
        type=int,
        default=1000,
        help='Longitud mínima para lecturas largas después de trimming (bp)'
    )
    long_read_group.add_argument(
        '--disable-long-adapter-trimming',
        action='store_true',
        help='Deshabilitar trimming de adaptadores en lecturas largas'
    )
    long_read_group.add_argument(
        '--disable-long-quality-filtering',
        action='store_true',
        help='Deshabilitar filtrado por calidad en lecturas largas'
    )
    
    args = parser.parse_args()
    
    # Configuración completa
    config = {
        'output_dir': Path(args.output).absolute(),
        'threads': args.threads,
        'max_size': args.max_size,
        'keep_temp': args.keep_temp,
        'force_long_reads': args.force_long_reads,
        'trim_params': {
            'quality_phred': args.quality_phred,
            'min_length': args.min_length,
            'cut_window_size': args.cut_window_size,
            'cut_mean_quality': args.cut_mean_quality,
            'disable_adapter_trimming': args.disable_adapter_trimming,
            'disable_quality_filtering': args.disable_quality_filtering,
            'disable_length_filtering': args.disable_length_filtering,
            'long_read_settings': {
                'min_quality': args.long_min_quality,
                'min_length': args.long_min_length,
                'disable_adapter_trimming': args.disable_long_adapter_trimming,
                'disable_quality_filtering': args.disable_long_quality_filtering
            }
        }
    }
    
    # Procesar lista de SRRs o archivos FASTQ
    try:
        with open(args.srr_list) as f:
            inputs = [line.strip() for line in f if line.strip()]
            
        for item in inputs:
            try:
                process_srr(item, config)
            except Exception as e:
                print(f"Error procesando {item}: {str(e)}")
                continue
                
    except FileNotFoundError:
        print(f"Error: Archivo no encontrado - {args.srr_list}")
        sys.exit(1)

if __name__ == '__main__':
    main()