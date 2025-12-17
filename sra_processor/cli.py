#!/usr/bin/env python3
import sys
import argparse
import shutil
from pathlib import Path
from .downloader import SRADownloader
from .trimmer import FastpProcessor
from .utils import check_srr_status, detect_input_type, detect_fastq_type, find_fastq_files
from .config import DEFAULT_CONFIG
from .exceptions import SRAProcessingError, UnsupportedDataTypeError
import logging
from .config_logging import setup_logging


def create_config(args, subcommand=None):
    """Crea configuración desde argumentos de CLI"""
    config = {
        'output_dir': Path(args.output).absolute(),
        'threads': args.threads,
        'keep_temp': args.keep_temp if hasattr(args, 'keep_temp') else False,
        'force_long_reads': args.force_long_reads if hasattr(args, 'force_long_reads') else False,
    }
    
    # Opciones específicas del subcomando
    if subcommand == 'download':
        config['max_size'] = args.max_size
        config['keep_sra'] = args.keep_sra
    elif subcommand in ['trim', 'full']:
        config['force_overwrite'] = args.force if hasattr(args, 'force') else False
        config['input_type'] = args.input_type if hasattr(args, 'input_type') else 'auto'
        config['trim_params'] = {
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
    
    if subcommand == 'full':
        config['max_size'] = args.max_size
        config['keep_sra'] = False  # Para full, no mantener por defecto
    
    return config


def process_download(srr_id, config):
    """Procesa solo la descarga de un SRR"""
    downloader = SRADownloader(config)
    
    try:
        status = check_srr_status(srr_id, config['output_dir'])
        
        if status == 'new':
            print(f"Descargando y convirtiendo {srr_id}...")
            data_type, fastq_files = downloader.download_and_convert(srr_id)
            print(f"✓ {srr_id} completado - Tipo: {data_type}")
            return True
        elif status == 'sra_downloaded':
            print(f"Convirtiendo {srr_id} (SRA ya descargado)...")
            data_type, fastq_files = downloader.convert_only(srr_id)
            print(f"✓ {srr_id} convertido - Tipo: {data_type}")
            return True
        else:
            print(f"✓ {srr_id} ya descargado (estado: {status})")
            return True
            
    except SRAProcessingError as e:
        print(f"✗ Error procesando {srr_id}: {str(e)}")
        return False


def process_trim_srr(srr_id, config):
    """Procesa solo el trimming de un SRR ya descargado"""
    trimmer = FastpProcessor(config)
    
    try:
        status = check_srr_status(srr_id, config['output_dir'])
        srr_dir = config['output_dir'] / srr_id
        
        if status == 'complete' and not config.get('force_overwrite', False):
            print(f"✓ {srr_id} ya procesado completamente")
            return True
        
        if status in ['fastq_ready', 'complete']:
            print(f"Haciendo trimming de {srr_id}...")
            # Buscar archivos FASTQ
            fastq_files = find_fastq_files(srr_dir, srr_id)
            
            if not fastq_files:
                print(f"✗ No se encontraron archivos FASTQ para {srr_id}")
                return False
            
            # Filtrar archivos ya procesados (trimmed)
            fastq_files = [f for f in fastq_files if 'trimmed' not in f.name]
            
            trimmer.process(fastq_files, srr_id)
            print(f"✓ {srr_id} trimming completado")
            return True
        
        elif status == 'single_end':
            print(f"Haciendo trimming de {srr_id} (single-end)...")
            fastq_files = find_fastq_files(srr_dir, srr_id)
            fastq_files = [f for f in fastq_files if 'trimmed' not in f.name]
            
            if fastq_files:
                trimmer.process(fastq_files, srr_id)
                print(f"✓ {srr_id} trimming completado")
                return True
            else:
                print(f"✗ No se encontraron archivos FASTQ para {srr_id}")
                return False
        else:
            print(f"✗ {srr_id} no tiene archivos FASTQ disponibles (estado: {status})")
            print(f"  Ejecute primero 'sra-processor download' para este SRR")
            return False
            
    except SRAProcessingError as e:
        print(f"✗ Error procesando {srr_id}: {str(e)}")
        return False


def process_trim_fastq(fastq_path, config):
    """Procesa el trimming de un archivo FASTQ externo"""
    trimmer = FastpProcessor(config)
    
    try:
        fastq_file = Path(fastq_path).absolute()
        
        if not fastq_file.exists():
            print(f"✗ Archivo no encontrado: {fastq_path}")
            return False
        
        # Detectar tipo de archivo
        file_type = detect_fastq_type(fastq_file)
        
        # Generar ID base desde nombre de archivo
        base_id = fastq_file.stem.replace('.fastq', '').replace('.fq', '').replace('.fas', '')
        base_id = base_id.replace('_1', '').replace('_2', '')
        
        print(f"Procesando {fastq_file.name} (tipo: {file_type})...")
        
        # Para paired-end, buscar el archivo pareja
        if file_type == 'paired' or '_1' in fastq_file.name:
            # Buscar archivo _2
            partner_name = fastq_file.name.replace('_1.', '_2.').replace('_1_', '_2_')
            partner_file = fastq_file.parent / partner_name
            
            if partner_file.exists():
                fastq_files = [fastq_file, partner_file]
                print(f"  Archivo pareja encontrado: {partner_name}")
            else:
                print(f"  Advertencia: No se encontró archivo _2, procesando como single-end")
                fastq_files = [fastq_file]
        else:
            fastq_files = [fastq_file]
        
        # Procesar con output_dir personalizado
        output_dir = config['output_dir']
        trimmer.process(fastq_files, base_id, output_dir=output_dir)
        print(f"✓ {fastq_file.name} trimming completado")
        return True
        
    except Exception as e:
        print(f"✗ Error procesando {fastq_path}: {str(e)}")
        return False


def process_full_pipeline(srr_id, config):
    """Procesa el pipeline completo (descarga + trimming)"""
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
                data_type, fastq_files = downloader.convert_to_fastq(srr_id, sra_dir)
            
            if data_type == 'single':
                print(f"{srr_id} es single-end. Procesando normalmente...")
            
            trimmer.process(fastq_files, srr_id)
            print(f"✓ {srr_id} completado")
            return True
        
        elif status == 'fastq_ready':
            print(f"Reanudando procesamiento para {srr_id} (FASTQs listos)...")
            srr_dir = config['output_dir'] / srr_id
            fastq_files = find_fastq_files(srr_dir, srr_id)
            fastq_files = [f for f in fastq_files if 'trimmed' not in f.name]
            
            if fastq_files:
                trimmer.process(fastq_files, srr_id)
                print(f"✓ {srr_id} completado")
                return True
            else:
                print(f"✗ No se encontraron archivos FASTQ para {srr_id}")
                return False
        
        elif status == 'complete':
            print(f"✓ {srr_id} ya procesado completamente")
            return True
        
        elif status == 'single_end':
            print(f"Reanudando procesamiento para {srr_id} (FASTQ single-end listo)...")
            srr_dir = config['output_dir'] / srr_id
            fastq_files = find_fastq_files(srr_dir, srr_id)
            fastq_files = [f for f in fastq_files if 'trimmed' not in f.name]
            
            if fastq_files:
                trimmer.process(fastq_files, srr_id)
                print(f"✓ {srr_id} completado")
                return True
            else:
                print(f"✗ No se encontraron archivos FASTQ para {srr_id}")
                return False
        
        else:  # 'unknown'
            print(f"Estado desconocido para {srr_id}. Limpiando y reintentando...")
            srr_dir = config['output_dir'] / srr_id
            if srr_dir.exists():
                shutil.rmtree(srr_dir)
            return process_full_pipeline(srr_id, config)
    
    except SRAProcessingError as e:
        print(f"✗ Error procesando {srr_id}: {str(e)}")
        return False


def cmd_download(args):
    """Comando: download - Solo descarga y conversión"""
    config = create_config(args, 'download')
    
    try:
        with open(args.srr_list) as f:
            srr_ids = [line.strip() for line in f if line.strip()]
        
        print(f"Procesando {len(srr_ids)} SRRs...")
        success_count = 0
        
        for srr_id in srr_ids:
            if process_download(srr_id, config):
                success_count += 1
        
        print(f"\n{'='*50}")
        print(f"Completado: {success_count}/{len(srr_ids)} SRRs procesados exitosamente")
        
    except FileNotFoundError:
        print(f"✗ Error: Archivo no encontrado - {args.srr_list}")
        sys.exit(1)


def cmd_trim(args):
    """Comando: trim - Solo trimming"""
    config = create_config(args, 'trim')
    
    try:
        with open(args.input) as f:
            inputs = [line.strip() for line in f if line.strip()]
        
        if not inputs:
            print("✗ Error: El archivo de entrada está vacío")
            sys.exit(1)
        
        # Detectar tipo de input
        input_type = config['input_type']
        if input_type == 'auto':
            # Auto-detectar desde la primera línea
            detected_type, _ = detect_input_type(inputs[0])
            input_type = detected_type
            print(f"Tipo de input auto-detectado: {input_type}")
        
        print(f"Procesando {len(inputs)} entradas...")
        success_count = 0
        
        if input_type == 'srr':
            for srr_id in inputs:
                if process_trim_srr(srr_id, config):
                    success_count += 1
        else:  # fastq
            for fastq_path in inputs:
                if process_trim_fastq(fastq_path, config):
                    success_count += 1
        
        print(f"\n{'='*50}")
        print(f"Completado: {success_count}/{len(inputs)} archivos procesados exitosamente")
        
    except FileNotFoundError:
        print(f"✗ Error: Archivo no encontrado - {args.input}")
        sys.exit(1)


def cmd_full(args):
    """Comando: full - Pipeline completo"""
    config = create_config(args, 'full')
    
    try:
        with open(args.srr_list) as f:
            srr_ids = [line.strip() for line in f if line.strip()]
        
        print(f"Procesando pipeline completo para {len(srr_ids)} SRRs...")
        success_count = 0
        
        for srr_id in srr_ids:
            if process_full_pipeline(srr_id, config):
                success_count += 1
        
        print(f"\n{'='*50}")
        print(f"Completado: {success_count}/{len(srr_ids)} SRRs procesados exitosamente")
        
    except FileNotFoundError:
        print(f"✗ Error: Archivo no encontrado - {args.srr_list}")
        sys.exit(1)


def add_common_args(parser):
    """Añade argumentos comunes a todos los subcomandos"""
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
        '--keep-temp',
        action='store_true',
        help='Mantener archivos temporales/intermedios'
    )
    parser.add_argument(
        '--force-long-reads',
        action='store_true',
        help='Forzar tratamiento como lecturas largas (Nanopore/PacBio)'
    )


def add_trim_args(parser):
    """Añade argumentos de trimming"""
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


def main():
    parser = argparse.ArgumentParser(
        description='Herramienta para procesamiento de datos desde SRA IDs (Illumina/Nanopore/PacBio)',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Subcomandos disponibles')
    subparsers.required = True
    
    # Subcomando: download
    parser_download = subparsers.add_parser(
        'download',
        help='Solo descarga y convierte archivos SRA a FASTQ',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser_download.add_argument(
        'srr_list',
        help='Archivo con lista de SRRs (uno por línea)'
    )
    add_common_args(parser_download)
    parser_download.add_argument(
        '--max-size',
        default='30G',
        help='Tamaño máximo para descarga por cada SRR (ej. 30G, 100G)'
    )
    parser_download.add_argument(
        '--keep-sra',
        action='store_true',
        help='Mantener archivos .sra después de conversión'
    )
    parser_download.set_defaults(func=cmd_download)
    
    # Subcomando: trim
    parser_trim = subparsers.add_parser(
        'trim',
        help='Solo hace trimming de archivos FASTQ existentes',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser_trim.add_argument(
        'input',
        help='Archivo con lista de SRR IDs o rutas a archivos FASTQ'
    )
    add_common_args(parser_trim)
    parser_trim.add_argument(
        '--input-type',
        choices=['auto', 'srr', 'fastq'],
        default='auto',
        help='Tipo de input: auto-detectar, SRR IDs, o rutas a archivos FASTQ'
    )
    parser_trim.add_argument(
        '--force',
        action='store_true',
        help='Sobrescribir archivos de salida existentes'
    )
    add_trim_args(parser_trim)
    parser_trim.set_defaults(func=cmd_trim)
    
    # Subcomando: full
    parser_full = subparsers.add_parser(
        'full',
        help='Pipeline completo: descarga + conversión + trimming',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser_full.add_argument(
        'srr_list',
        help='Archivo con lista de SRRs (uno por línea)'
    )
    add_common_args(parser_full)
    parser_full.add_argument(
        '--max-size',
        default='30G',
        help='Tamaño máximo para descarga por cada SRR (ej. 30G, 100G)'
    )
    add_trim_args(parser_full)
    parser_full.set_defaults(func=cmd_full)
    
    # Parsear argumentos y ejecutar comando
    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
