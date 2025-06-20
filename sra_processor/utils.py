from pathlib import Path

def check_srr_status(srr_id, output_dir):
    """Verifica el estado de procesamiento considerando múltiples extensiones"""
    srr_dir = output_dir / srr_id
    
    if not srr_dir.exists():
        return 'new'
    
    # Verifica archivos finales (comprimidos)
    if (srr_dir / f"output_{srr_id}_1_paired.fastq.gz").exists() and \
       (srr_dir / f"output_{srr_id}_2_paired.fastq.gz").exists():
        return 'complete'
    
    # Verifica archivos intermedios con múltiples extensiones
    extensions = ['.fastq', '.fq', '.fas']
    
    # Paired-end
    for ext in extensions:
        if (srr_dir / f"{srr_id}_1{ext}").exists() and \
           (srr_dir / f"{srr_id}_2{ext}").exists():
            return 'fastq_ready'
    
    # Single-end
    for ext in extensions:
        if (srr_dir / f"{srr_id}{ext}").exists():
            return 'single_end'
    
    # SRA descargado
    if (srr_dir / f"{srr_id}.sra").exists():
        return 'sra_downloaded'
    
    return 'unknown'