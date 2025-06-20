#!/usr/bin/env python3
import argparse
import shutil
from pathlib import Path
from .downloader import SRADownloader
from .trimmer import FastpProcessor
from .utils import check_srr_status
from .config import DEFAULT_CONFIG
from .exceptions import SRAProcessingError, UnsupportedDataType

def process_srr(srr_id, config):
    """Procesa un solo SRR a través de todo el pipeline"""
    downloader = SRADownloader(config)
    trimmer = FastpProcessor(config)
    
    status = check_srr_status(srr_id, config['output_dir'])
    
    try:
        if status in ['new', 'sra_downloaded']:
            print(f"Procesando {srr_id} (estado: {status})...")
            
            if status == 'new':
                sra_dir = downloader.download_sra(srr_id)
            else:
                sra_dir = config['output_dir'] / srr_id
            
            data_type, fastq_files = downloader.convert_to_fastq(srr_id, sra_dir)
            
            if data_type == 'single':
                print(f"AVISO: {srr_id} es single-end. Creando archivo log y continuando...")
                with open(config['output_dir'] / "single_end_samples.log", "a") as f:
                    f.write(f"{srr_id}\n")
                return False  # Continuar con el siguiente SRR sin error
            
            trimmer.process_paired_end(fastq_files, srr_id)
            return True
            
        elif status == 'fastq_ready':
            print(f"Reanudando procesamiento para {srr_id} (FASTQs listos)...")
            fastq_files = [
                config['output_dir'] / srr_id / f"{srr_id}_1.fastq",
                config['output_dir'] / srr_id / f"{srr_id}_2.fastq"
            ]
            trimmer.process_paired_end(fastq_files, srr_id)
            return True
            
        elif status == 'complete':
            print(f"Saltando {srr_id} (ya procesado completamente)")
            return True
            
        elif status == 'single_end':
            print(f"Saltando {srr_id} (single-end no soportado)")
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
    parser = argparse.ArgumentParser(description='Procesador de datos SRA a FASTQ con fastp')
    
    parser.add_argument('srr_list', help='Archivo con lista de IDs SRA (uno por línea)')
    parser.add_argument('-o', '--output', default='.', help='Directorio de salida')
    parser.add_argument('-t', '--threads', type=int, default=2, help='Número de hilos a usar (default: 2)')
    parser.add_argument('--max-size', default='30G', help='Tamaño máximo para descarga (default: 30Gb)')
    parser.add_argument('--keep-temp', action='store_true', help='Mantener archivos temporales')
    
    # Opciones de fastp
    trim_group = parser.add_argument_group('Opciones de fastp')
    trim_group.add_argument('--quality-phred', type=int, default=20, 
                          help='Calidad Phred mínima (default: 20)')
    trim_group.add_argument('--min-length', type=int, default=50, 
                          help='Longitud mínima de reads (default: 50)')
    trim_group.add_argument('--cut-window-size', type=int, default=4, 
                          help='Tamaño de ventana para trimming (default: 4)')
    trim_group.add_argument('--cut-mean-quality', type=int, default=20, 
                          help='Calidad media para trimming (default: 20)')
    trim_group.add_argument('--disable-adapter-trimming', action='store_true', 
                          help='Deshabilitar trimming de adaptadores')
    trim_group.add_argument('--disable-quality-filtering', action='store_true', 
                          help='Deshabilitar filtrado por calidad')
    trim_group.add_argument('--disable-length-filtering', action='store_true', 
                          help='Deshabilitar filtrado por longitud')
    
    args = parser.parse_args()
    
    config = DEFAULT_CONFIG.copy()
    config.update({
        'output_dir': Path(args.output).absolute(),
        'threads': args.threads,
        'max_size': args.max_size,
        'keep_temp': args.keep_temp,
        'trim_params': {
            'quality_phred': args.quality_phred,
            'min_length': args.min_length,
            'cut_window_size': args.cut_window_size,
            'cut_mean_quality': args.cut_mean_quality,
            'disable_adapter_trimming': args.disable_adapter_trimming,
            'disable_quality_filtering': args.disable_quality_filtering,
            'disable_length_filtering': args.disable_length_filtering
        }
    })
    
    with open(args.srr_list) as f:
        srr_ids = [line.strip() for line in f if line.strip()]
    
    for srr_id in srr_ids:
        process_srr(srr_id, config)

if __name__ == '__main__':
    main()