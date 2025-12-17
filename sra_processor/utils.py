import logging
from pathlib import Path
from typing import Literal, Optional, List, Tuple

logger = logging.getLogger(__name__)

# Constants
FASTQ_EXTENSIONS = ['.fastq', '.fq', '.fas']
FASTQ_EXTENSIONS_COMPRESSED = ['.fastq.gz', '.fq.gz']
ALL_FASTQ_EXTENSIONS = FASTQ_EXTENSIONS + FASTQ_EXTENSIONS_COMPRESSED
MIN_SRR_ID_LENGTH = 9  # Minimum length for valid SRR/ERR/DRR IDs

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
    # Paired-end intermedio
    for ext in FASTQ_EXTENSIONS:
        paired_files = [
            srr_dir / f"{srr_id}_1{ext}",
            srr_dir / f"{srr_id}_2{ext}"
        ]
        if all(f.exists() for f in paired_files):
            logger.debug(f"Estado para {srr_id}: fastq_ready (paired-end)")
            return 'fastq_ready'
    
    # Single-end intermedio
    for ext in FASTQ_EXTENSIONS:
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
        for ext in FASTQ_EXTENSIONS:
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


def detect_fastq_type(fastq_file: Path) -> Literal['paired', 'single', 'long']:
    """
    Detecta el tipo de archivo FASTQ basado en su contenido
    
    Args:
        fastq_file: Ruta al archivo FASTQ
    
    Returns:
        Tipo detectado: 'paired', 'single', o 'long'
    """
    try:
        with open(fastq_file, 'r') as f:
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
                logger.warning(f"No se pudieron leer secuencias de {fastq_file}")
                return 'single'
            
            avg_length = sum(lengths) / len(lengths)
            
            # Considerar lecturas largas si promedio > 500bp
            if avg_length > 500:
                logger.debug(f"Detectado long-read: longitud promedio {avg_length:.1f}bp")
                return 'long'
            
            # Para paired/single, verificar si hay archivos _1 y _2
            parent = fastq_file.parent
            name = fastq_file.name
            
            # Buscar patrones _1 y _2 (más precisos usando sufijos)
            import re
            if re.search(r'_1[._]', name):
                # Reemplazar solo la última ocurrencia de _1 antes de la extensión
                partner_name = re.sub(r'_1([._])', r'_2\1', name)
                partner = parent / partner_name
                if partner.exists():
                    logger.debug(f"Detectado paired-end: encontrado archivo _2")
                    return 'paired'
            
            logger.debug(f"Detectado single-end: longitud promedio {avg_length:.1f}bp")
            return 'single'
            
    except Exception as e:
        logger.warning(f"Error detectando tipo FASTQ: {str(e)}")
        return 'single'


def find_fastq_files(directory: Path, srr_id: Optional[str] = None) -> List[Path]:
    """
    Busca archivos FASTQ en un directorio
    
    Args:
        directory: Directorio donde buscar
        srr_id: ID opcional para filtrar archivos específicos
    
    Returns:
        Lista de archivos FASTQ encontrados
    """
    if not directory.exists():
        logger.warning(f"Directorio no existe: {directory}")
        return []
    
    # Build glob patterns from extensions
    extensions = [f'*{ext}' for ext in ALL_FASTQ_EXTENSIONS]
    fastq_files = []
    
    for ext in extensions:
        if srr_id:
            # Buscar archivos específicos del SRR
            pattern = f"{srr_id}*{ext[1:]}"  # Quitar el * del inicio
            fastq_files.extend(directory.glob(pattern))
        else:
            # Buscar todos los archivos con la extensión
            fastq_files.extend(directory.glob(ext))
    
    # Ordenar para tener consistencia (_1 antes de _2)
    fastq_files.sort()
    
    logger.debug(f"Encontrados {len(fastq_files)} archivos FASTQ en {directory}")
    return fastq_files


def detect_input_type(input_line: str) -> Tuple[str, str]:
    """
    Detecta si una línea es un SRR ID o una ruta a archivo FASTQ
    
    Args:
        input_line: Línea de entrada a analizar
    
    Returns:
        Tupla (tipo, valor): tipo es 'srr' o 'fastq', valor es el input normalizado
    """
    input_line = input_line.strip()
    
    # Verificar si es una ruta a archivo (tiene separadores de ruta o extensión FASTQ)
    has_path_sep = '/' in input_line or '\\' in input_line
    has_fastq_ext = any(input_line.endswith(ext) for ext in ALL_FASTQ_EXTENSIONS)
    
    if has_path_sep or has_fastq_ext:
        return 'fastq', input_line
    
    # Verificar si es un SRR ID (formato SRR/ERR/DRR seguido de números)
    if input_line.upper().startswith(('SRR', 'ERR', 'DRR')) and len(input_line) >= MIN_SRR_ID_LENGTH:
        return 'srr', input_line
    
    # Por defecto, asumir que es un SRR ID
    return 'srr', input_line