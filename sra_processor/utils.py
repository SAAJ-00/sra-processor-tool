import logging
from pathlib import Path
from typing import Literal, Optional

logger = logging.getLogger(__name__)

def check_srr_status(srr_id: str, output_dir: Path) -> Literal['new', 'complete', 'fastq_ready', 'single_end', 'sra_downloaded', 'unknown']:
    """
    Verifica el estado de procesamiento de un SRR
    
    Args:
        srr_id: Identificador del experimento SRA
        output_dir: Directorio base de salida
    
    Returns:
        Estado del procesamiento:
        - 'new': No iniciado
        - 'sra_downloaded': SRA descargado pero no convertido
        - 'fastq_ready': FASTQs generados (paired-end)
        - 'single_end': FASTQ single-end generado
        - 'complete': Procesamiento finalizado (archivos trimmed)
        - 'unknown': Estado no reconocido
    """
    srr_dir = output_dir / srr_id
    
    if not srr_dir.exists():
        logger.debug(f"Estado para {srr_id}: new (directorio no existe)")
        return 'new'
    
    # 1. Verificar si el procesamiento está completo (archivos finales)
    paired_complete = [
        srr_dir / f"{srr_id}_1_trimmed.fastq.gz",
        srr_dir / f"{srr_id}_2_trimmed.fastq.gz"
    ]
    
    single_complete = srr_dir / f"{srr_id}_trimmed.fastq.gz"
    
    if all(f.exists() for f in paired_complete):
        logger.debug(f"Estado para {srr_id}: complete (paired-end)")
        return 'complete'
    if single_complete.exists():
        logger.debug(f"Estado para {srr_id}: complete (single-end)")
        return 'complete'
    
    # 2. Verificar archivos intermedios FASTQ
    extensions = ['.fastq', '.fq', '.fas']
    
    # Paired-end intermedio
    for ext in extensions:
        paired_files = [
            srr_dir / f"{srr_id}_1{ext}",
            srr_dir / f"{srr_id}_2{ext}"
        ]
        if all(f.exists() for f in paired_files):
            logger.debug(f"Estado para {srr_id}: fastq_ready (paired-end)")
            return 'fastq_ready'
    
    # Single-end intermedio
    for ext in extensions:
        single_file = srr_dir / f"{srr_id}{ext}"
        if single_file.exists():
            logger.debug(f"Estado para {srr_id}: single_end")
            return 'single_end'
    
    # 3. Verificar si solo está descargado el SRA (en carpeta principal o subcarpeta)
    sra_file = srr_dir / f"{srr_id}.sra"
    sra_file_sub = srr_dir / srr_id / f"{srr_id}.sra"
    if sra_file.exists() or sra_file_sub.exists():
        logger.debug(f"Estado para {srr_id}: sra_downloaded")
        return 'sra_downloaded'
    
    logger.warning(f"Estado para {srr_id}: unknown (no se reconoce el estado)")
    return 'unknown'

def get_fastq_files(srr_id: str, output_dir: Path) -> Optional[list[Path]]:
    """
    Obtiene los archivos FASTQ para un SRR dado
    
    Args:
        srr_id: Identificador del experimento SRA
        output_dir: Directorio base de salida
    
    Returns:
        Lista de archivos FASTQ encontrados o None si no hay archivos válidos
    """
    status = check_srr_status(srr_id, output_dir)
    srr_dir = output_dir / srr_id
    
    if status in ['complete', 'fastq_ready']:
        # Buscar archivos finales primero
        paired_files = [
            srr_dir / f"{srr_id}_1_trimmed.fastq.gz",
            srr_dir / f"{srr_id}_2_trimmed.fastq.gz"
        ]
        if all(f.exists() for f in paired_files):
            return paired_files
        
        single_file = srr_dir / f"{srr_id}_trimmed.fastq.gz"
        if single_file.exists():
            return [single_file]
        
        # Buscar archivos intermedios
        extensions = ['.fastq', '.fq', '.fas']
        for ext in extensions:
            paired_intermediate = [
                srr_dir / f"{srr_id}_1{ext}",
                srr_dir / f"{srr_id}_2{ext}"
            ]
            if all(f.exists() for f in paired_intermediate):
                return paired_intermediate
            
            single_intermediate = srr_dir / f"{srr_id}{ext}"
            if single_intermediate.exists():
                return [single_intermediate]
    
    return None